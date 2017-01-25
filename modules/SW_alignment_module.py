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

def sw_align_sequences(matches, single_core = False):
    exact_matches = {}
    if single_core:
        for j, s1 in enumerate(matches):
            for i, s2 in enumerate(matches[s1]):
                if s1 in exact_matches:
                    if s2 in exact_matches[s1]:
                        continue
                # print(s1,s2)
                s1, s2, stats = ssw_alignment_helper( ((s1, s2, i, j), {}) )
                if stats:
                    if s1 in exact_matches:
                        exact_matches[s1][s2] = stats
                    else:
                        exact_matches[s1] = {}
                        exact_matches[s1][s2]  = stats
                else:
                    pass
    else:
        ####### parallelize alignment #########
        # pool = Pool(processes=mp.cpu_count())
        original_sigint_handler = signal.signal(signal.SIGINT, signal.SIG_IGN)
        signal.signal(signal.SIGINT, original_sigint_handler)
        pool = Pool(processes=mp.cpu_count())
        try:
            res = pool.map_async(ssw_alignment_helper, [ ((s1, s2, i,j), {}) for j, s1 in enumerate(matches) for i, s2 in enumerate(matches[s1]) ] )
            alignment_results =res.get(999999999) # Without the timeout this blocking call ignores all signals.
        except KeyboardInterrupt:
            print("Caught KeyboardInterrupt, terminating workers")
            pool.terminate()
            sys.exit()
        else:
            print("Normal termination")
            pool.close()
        pool.join()
        for s1, s2, stats in alignment_results:
            if stats:
                if s1 in exact_matches:
                    exact_matches[s1][s2] = stats
                else:
                    exact_matches[s1] = {}
                    exact_matches[s1][s2] = stats
            else:
                pass

    return exact_matches


def sw_align_sequences_keeping_accession(matches, single_core = False):
    exact_matches = {}

    if single_core:
        for j, s1_acc in enumerate(matches):
            for i, s2_acc in enumerate(matches[s1_acc]):
                s1, s2 = matches[s1_acc][s2_acc]
                if s1_acc in exact_matches:
                    if s2_acc in exact_matches[s1_acc]:
                        continue
                s1_acc, s2_acc, stats = ssw_alignment_helper( ((s1, s2, i, j), {"x_acc" : s1_acc, "y_acc" :s2_acc }) )
                if stats:
                    if s1_acc in exact_matches:
                        exact_matches[s1_acc][s2_acc] = stats
                    else:
                        exact_matches[s1_acc] = {}
                        exact_matches[s1_acc][s2_acc]  = stats
                else:
                    pass
    else:
        ####### parallelize alignment #########
        # pool = Pool(processes=mp.cpu_count())
        original_sigint_handler = signal.signal(signal.SIGINT, signal.SIG_IGN)
        signal.signal(signal.SIGINT, original_sigint_handler)
        pool = Pool(processes=mp.cpu_count())

        # for j, s1_acc in enumerate(matches):
        #     for i, s2_acc in enumerate(matches[s1_acc]):
        #         print("lool", matches[s1_acc][s2_acc][0], matches[s1_acc][s2_acc][1], i,j, {"x_acc": s1_acc, "y_acc" : s2_acc} ) 
        try:
            res = pool.map_async(ssw_alignment_helper, [ ((matches[s1_acc][s2_acc][0], matches[s1_acc][s2_acc][1], i,j), {"x_acc": s1_acc, "y_acc" : s2_acc}) for j, s1_acc in enumerate(matches) for i, s2_acc in enumerate(matches[s1_acc]) ] )
            alignment_results =res.get(999999999) # Without the timeout this blocking call ignores all signals.
        except KeyboardInterrupt:
            print("Caught KeyboardInterrupt, terminating workers")
            pool.terminate()
            sys.exit()
        else:
            print("Normal termination")
            pool.close()
        pool.join()
        for s1_acc, s2_acc, stats in alignment_results:
            if stats:
                if s1_acc in exact_matches:
                    exact_matches[s1_acc][s2_acc] = stats
                else:
                    exact_matches[s1_acc] = {}
                    exact_matches[s1_acc][s2_acc] = stats
            else:
                # print("OMG!")
                # print(len(matches[s1_acc][s2_acc][0]), len(matches[s1_acc][s2_acc][0]) )
                pass
    return exact_matches


def ssw_alignment_helper(arguments):
    args, kwargs = arguments
    return ssw_alignment(*args, **kwargs)

def ssw_alignment(x, y, i,j, ends_discrepancy_threshold = 25 , x_acc = "", y_acc = "" ):
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

    if x_acc == y_acc == "":
        if start_discrepancy > ends_discrepancy_threshold or end_discrepancy > ends_discrepancy_threshold:
            # print("REMOVING", start_discrepancy, end_discrepancy)
            return (x, y, None)

        else:
            return (x, y, (x_alignment, y_alignment, (matches, mismatches, indels)) )
    else:
        if start_discrepancy > ends_discrepancy_threshold or end_discrepancy > ends_discrepancy_threshold:
            # print("REMOVING", start_discrepancy, end_discrepancy)
            return (x_acc, y_acc, None)

        else:
            return (x_acc, y_acc, (x_alignment, y_alignment, (matches, mismatches, indels)) )       



def find_best_matches(approximate_matches):
    """
        input: approximate_matches is a dictionary with a string as key and a list of strings as value
        output: dictionary with a string as key and a dictionary as value. 
                The inner dict contains a tuple (edit_distance, s1_alignment, s2_alignment) as value.
                Each string in the inner dict has the same (lowest) edit distance to the key 
    """
    exact_matches = sw_align_sequences(approximate_matches, single_core = False )

    # process the exact matches here
    best_exact_matches = {}
    for s1 in exact_matches:
        for s2 in exact_matches[s1]:
            s1_alignment, s2_alignment, (matches, mismatches, indels) = exact_matches[s1][s2]
            edit_distance = mismatches + indels
            if s1 in best_exact_matches:
                s1_minimizer = best_exact_matches[s1].keys()[0]
                if edit_distance < best_exact_matches[s1][s1_minimizer][0]:
                    best_exact_matches[s1] = {}
                    best_exact_matches[s1][s2] = (edit_distance, s1_alignment, s2_alignment)
                elif edit_distance == best_exact_matches[s1][s1_minimizer][0]:
                    best_exact_matches[s1][s2] = (edit_distance, s1_alignment, s2_alignment)
            else:
                best_exact_matches[s1] = {}
                best_exact_matches[s1][s2] = (edit_distance, s1_alignment, s2_alignment)

            if s2 in best_exact_matches:
                s2_minimizer = best_exact_matches[s2].keys()[0]
                if edit_distance < best_exact_matches[s2][s2_minimizer][0]:
                    best_exact_matches[s2] = {}
                    best_exact_matches[s2][s1] = (edit_distance, s2_alignment, s1_alignment)
                elif edit_distance == best_exact_matches[s2][s2_minimizer][0]:
                    best_exact_matches[s2][s1] = (edit_distance, s2_alignment, s1_alignment)
            else:
                best_exact_matches[s2] = {}
                best_exact_matches[s2][s1] = (edit_distance, s2_alignment, s1_alignment)
    # for s1 in best_exact_matches:
    #     for s2 in best_exact_matches[s1]:
    #         print(best_exact_matches[s1][s2][0])
    return best_exact_matches

def find_best_matches_2set(highest_paf_scores, X, C):
    """
        input: approximate_matches is a dictionary with a string as key and a list of strings as value
        output: dictionary with a string as key and a dictionary as value. 
                the outer dict contains the reads .
                This inner dict dict contains a tuple (edit_distance, s1_alignment, s2_alignment) as value.
                Each string in the inner dict has the same (lowest) edit distance to the key 
    """

    approximate_matches = {}
    for read_acc, best_hits in  highest_paf_scores.items():
        approximate_matches[read_acc] = {}
        # print("NEW")
        for score, t_acc in best_hits:
            # print(score)
            approximate_matches[read_acc][t_acc] = (X[read_acc], C[t_acc])
        # print(len(approximate_matches[read_acc]))
    exact_matches = sw_align_sequences_keeping_accession(approximate_matches, single_core = False )

    # process the exact matches here
    best_exact_matches = {}
    for x_acc in exact_matches:
        for c_acc in exact_matches[x_acc]:
            x_alignment, c_alignment, (matches, mismatches, indels) = exact_matches[x_acc][c_acc]
            edit_distance = mismatches + indels

            if x_acc in best_exact_matches:
                x_minimizer = best_exact_matches[x_acc].keys()[0]
                if edit_distance < best_exact_matches[x_acc][x_minimizer][0]:
                    best_exact_matches[x_acc] = {}
                    best_exact_matches[x_acc][c_acc] = (edit_distance, x_alignment, c_alignment)
                elif edit_distance == best_exact_matches[x_acc][x_minimizer][0]:
                    best_exact_matches[x_acc][c_acc] = (edit_distance, x_alignment, c_alignment)
                else:
                    pass
                    # print("opt:", best_exact_matches[x_acc][x_minimizer][0], "sub:", edit_distance)
            else:
                best_exact_matches[x_acc] = {}
                best_exact_matches[x_acc][c_acc] = (edit_distance, x_alignment, c_alignment)

            # if c_acc in best_exact_matches:
            #     c_minimizer = best_exact_matches[c_acc].keys()[0]
            #     if edit_distance < best_exact_matches[c_acc][c_minimizer][0]:
            #         best_exact_matches[c_acc] = {}
            #         best_exact_matches[c_acc][x_acc] = (edit_distance, c_alignment, x_alignment)
            #     elif edit_distance == best_exact_matches[c_acc][c_minimizer][0]:
            #         best_exact_matches[c_acc][x_acc] = (edit_distance, c_alignment, x_alignment)
            # else:
            #     best_exact_matches[c_acc] = {}
            #     best_exact_matches[c_acc][x_acc] = (edit_distance, s2_alignment, x_alignment)
    
    return best_exact_matches
