import signal
from multiprocessing import Pool
import multiprocessing as mp
import sys

import ssw
from Bio import pairwise2
from Bio.SubsMat import MatrixInfo as matlist

# def ssw_alignment_vector( (x_acc, y_acc, x, y) ):
#     result_vector = []
#     for i in range(100):
#         result_vector.append(ssw_alignment( (x_acc, y_acc, x, y) ))
#     return result_vector

def sw_align_sequences(alignments_to_map, alignment_results_dict, x, y, params):
    if params.single_core:
        for j, y_j in enumerate(alignments_to_map.keys()):
            for i, x_i in enumerate(alignments_to_map[y_j].keys()):
                y_j, x_i, stats = ssw_alignment_helper( (x_i, y_j, x[x_i], y[y_j], i,j) )
                if stats:
                    if y_j in alignment_results_dict:
                        alignment_results_dict[y_j][x_i] = stats
                    else:
                        alignment_results_dict[y_j] = {}
                        alignment_results_dict[y_j][x_i] = stats
                else:
                    del alignment_results_dict[y_j][x_i]
    else:
        ####### parallelize alignment #########
        # pool = Pool(processes=mp.cpu_count())
        original_sigint_handler = signal.signal(signal.SIGINT, signal.SIG_IGN)
        signal.signal(signal.SIGINT, original_sigint_handler)
        pool = Pool(processes=mp.cpu_count())
        try:
            res = pool.map_async(ssw_alignment_helper, [(x_i, y_j, x[x_i], y[y_j], i,j) for j, y_j in enumerate(alignments_to_map) for i, x_i in enumerate(alignments_to_map[y_j]) ] )
            alignment_results =res.get(999999999) # Without the timeout this blocking call ignores all signals.
        except KeyboardInterrupt:
            print("Caught KeyboardInterrupt, terminating workers")
            pool.terminate()
            sys.exit()
        else:
            print("Normal termination")
            pool.close()
        pool.join()
        for y_j, x_i, stats in alignment_results:
            if stats:
                if y_j in alignment_results_dict:
                    alignment_results_dict[y_j][x_i] = stats
                else:
                    alignment_results_dict[y_j] = {}
                    alignment_results_dict[y_j][x_i] = stats
            else: 
                del alignment_results_dict[y_j][x_i]


def ssw_alignment_helper(args):
    return ssw_alignment(*args)

def ssw_alignment( x_acc, y_acc, x, y, i,j, ends_discrepancy_threshold = 25 ):
    """
        Aligns two sequences with SSW
        x: query
        y: reference

    """
    if i % 10000 == 0 and j % 10 == 0:
        print("processing alignments on y_j with j={0}".format(j+1))

    score_matrix = ssw.DNA_ScoreMatrix(match=1, mismatch=-1)
    aligner = ssw.Aligner(gap_open=2, gap_extend=1, matrix=score_matrix)

    # for the ends that SSW leaves behind
    bio_matrix = matlist.blosum62
    g_open = -1
    g_extend = -0.5
    ######################################

    # result = aligner.align("GA", "G", revcomp=False)
    # y_alignment, match_line, x_alignment = result.alignment
    # c = Counter(match_line)
    # matches, mismatches, indels = c["|"], c["*"], c[" "]
    # alignment_length = len(match_line)
    # print("matches:{0}, mismatches:{1}, indels:{2} ".format(matches, mismatches, indels))
    # print(match_line)

    result = aligner.align(x, y, revcomp=False)
    y_alignment, match_line, x_alignment = result.alignment

    matches, mismatches, indels = match_line.count("|"), match_line.count("*"), match_line.count(" ")

    # alignment_length = len(match_line)
    
    start_discrepancy = max(result.query_begin, result.reference_begin)  # 0-indexed # max(result.query_begin, result.reference_begin) - min(result.query_begin, result.reference_begin)
    query_end_discrepancy = len(x) - result.query_end - 1
    ref_end_discrepancy = len(y) - result.reference_end - 1
    end_discrepancy = max(query_end_discrepancy, ref_end_discrepancy)  # max(result.query_end, result.reference_end) - min(result.query_end, result.reference_end)
    # print(start_discrepancy, end_discrepancy)
    tot_discrepancy = start_discrepancy + end_discrepancy

    if 0 < start_discrepancy <= ends_discrepancy_threshold:
        # print("HERE")
        matches_snippet = 0
        mismatches_snippet = 0
        if result.query_begin and result.reference_begin:
            query_start_snippet = x[:result.query_begin]
            ref_start_snippet = y[:result.reference_begin]
            alns = pairwise2.align.globalds(query_start_snippet, ref_start_snippet, bio_matrix, g_open, g_extend)
            top_aln = alns[0]
            # print(alns)
            mismatches_snippet = len(list(filter(lambda x: x[0] != x[1] and x[0] != '-' and x[1] != "-", zip(top_aln[0],top_aln[1]))))
            indels_snippet = top_aln[0].count("-") + top_aln[1].count("-")
            matches_snippet = len(top_aln[0]) - mismatches_snippet - indels_snippet
            # print(matches_snippet, mismatches_snippet, indels_snippet)
            query_start_alignment_snippet = top_aln[0]
            ref_start_alignment_snippet = top_aln[1]
        elif result.query_begin:
            query_start_alignment_snippet = x[:result.query_begin]
            ref_start_alignment_snippet = "-"*len(query_start_alignment_snippet)
            indels_snippet = len(ref_start_alignment_snippet)
        elif result.reference_begin:
            ref_start_alignment_snippet = y[:result.reference_begin]
            query_start_alignment_snippet = "-"*len(ref_start_alignment_snippet)
            indels_snippet = len(query_start_alignment_snippet)
        else:
            print("BUG")
            sys.exit()
        matches, mismatches, indels = matches + matches_snippet, mismatches + mismatches_snippet, indels + indels_snippet

        # print(ref_start_alignment_snippet)
        # print(query_start_alignment_snippet)
        y_alignment = ref_start_alignment_snippet + y_alignment
        x_alignment = query_start_alignment_snippet + x_alignment

    if 0 < end_discrepancy <= ends_discrepancy_threshold:
        # print("HERE2")
        matches_snippet = 0
        mismatches_snippet = 0
        if query_end_discrepancy and ref_end_discrepancy:
            query_end_snippet = x[result.query_end+1:]
            ref_end_snippet = y[result.reference_end+1:]
            alns = pairwise2.align.globalds(query_end_snippet, ref_end_snippet, bio_matrix, g_open, g_extend)
            top_aln = alns[0]
            mismatches_snippet = len(list(filter(lambda x: x[0] != x[1] and x[0] != '-' and x[1] != "-", zip(top_aln[0],top_aln[1]))))
            indels_snippet = top_aln[0].count("-") + top_aln[1].count("-")
            matches_snippet = len(top_aln[0]) - mismatches_snippet - indels_snippet
            query_end_alignment_snippet = top_aln[0]
            ref_end_alignment_snippet = top_aln[1]
        elif query_end_discrepancy:
            query_end_alignment_snippet = x[result.query_end+1:]
            ref_end_alignment_snippet = "-"*len(query_end_alignment_snippet)
            indels_snippet = len(ref_end_alignment_snippet)

        elif ref_end_discrepancy:
            ref_end_alignment_snippet = y[result.reference_end+1:]
            query_end_alignment_snippet = "-"*len(ref_end_alignment_snippet)
            indels_snippet = len(query_end_alignment_snippet)

        else:
            print("BUG")
            sys.exit()
        matches, mismatches, indels = matches + matches_snippet, mismatches + mismatches_snippet, indels + indels_snippet

        y_alignment = y_alignment + ref_end_alignment_snippet
        x_alignment = x_alignment + query_end_alignment_snippet

    if start_discrepancy > ends_discrepancy_threshold or end_discrepancy > ends_discrepancy_threshold:
        # print("REMOVING", start_discrepancy, end_discrepancy)
        return (y_acc, x_acc, None)

    else:
        return (y_acc, x_acc, (y_alignment, x_alignment, (matches, mismatches, indels)) )
    # similarity_SW = matches / float(alignment_length + tot_discrepancy)


