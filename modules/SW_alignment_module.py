import signal
from multiprocessing import Pool
import multiprocessing as mp
import sys

import parasail
import re

# import ssw
# from Bio import pairwise2
# from Bio.SubsMat import MatrixInfo as matlist



def cigar_to_seq(cigar, query, ref):
    cigar_tuples = []
    result = re.split(r'[=DXSMI]+', cigar)
    i = 0
    for length in result[:-1]:
        i += len(length)
        type_ = cigar[i]
        i += 1
        cigar_tuples.append((int(length), type_ ))

    r_index = 0
    q_index = 0
    q_aln = []
    r_aln = []
    for length_ , type_ in cigar_tuples:
        if type_ == "=" or type_ == "X":
            q_aln.append(query[q_index : q_index + length_])
            r_aln.append(ref[r_index : r_index + length_])

            r_index += length_
            q_index += length_
        
        elif  type_ == "I":
            # insertion w.r.t. reference
            r_aln.append('-' * length_)
            q_aln.append(query[q_index: q_index + length_])
            #  only query index change
            q_index += length_

        elif type_ == 'D':
            # deletion w.r.t. reference
            r_aln.append(ref[r_index: r_index + length_])
            q_aln.append('-' * length_)
            #  only ref index change
            r_index += length_
        
        else:
            print("error")
            print(cigar)
            sys.exit()

    return  "".join([s for s in q_aln]), "".join([s for s in r_aln])


def parasail_alignment_helper(arguments):
    args, kwargs = arguments
    return parasail_alignment(*args, **kwargs)


def parasail_alignment(s1, s2, i, j, x_acc = "", y_acc = "", match_score = 2, mismatch_penalty = -3, opening_penalty = 2, gap_ext = 0):
    user_matrix = parasail.matrix_create("ACGT", match_score, mismatch_penalty)
    result = parasail.sg_trace_scan_16(s1, s2, opening_penalty, gap_ext, user_matrix)
    if result.saturated:
        print("SATURATED!")
        result = parasail.sg_trace_scan_32(s1, s2, opening_penalty, gap_ext, user_matrix)
    # print(result.cigar.seq)
    # print(result.cigar.decode )
    # print(str(result.cigar.decode,'utf-8') )
    if sys.version_info[0] < 3:
        cigar_string = str(result.cigar.decode).decode('utf-8')
    else:
        cigar_string = str(result.cigar.decode, 'utf-8')

    s1_alignment, s2_alignment = cigar_to_seq(cigar_string, s1, s2)
    mismatches = len([ 1 for n1, n2 in zip(s1_alignment,s2_alignment) if n1 != n2 and n1 != "-" and n2 != "-" ])
    matches = len([ 1 for n1, n2 in zip(s1_alignment,s2_alignment) if n1 == n2 and n1 != "-"])
    indels = len(s1_alignment) - mismatches - matches

    if x_acc == y_acc == "":
        return (s1, s2, (s1_alignment, s2_alignment, (matches, mismatches, indels)) )
    else:
        return (x_acc, y_acc, (s1_alignment, s2_alignment, (matches, mismatches, indels)) ) 


def sw_align_sequences(matches, nr_cores = 1, mismatch_penalty = -1):
    """
        Matches should be a 2D matrix implemented as a dict of dict, the value should be the edit distance.
    """
    exact_matches = {}

    if nr_cores == 1:
        for j, s1 in enumerate(matches):
            for i, s2 in enumerate(matches[s1]):
                if s1 in exact_matches:
                    if s2 in exact_matches[s1]:
                        continue

                ed = matches[s1][s2]
                error_rate = float(ed)/ min(len(s1), len(s2))
                if error_rate <= 0.01:
                    mismatch_penalty = -1
                elif 0.01 < error_rate <= 0.09:
                    mismatch_penalty = -2
                else:
                    mismatch_penalty = -4

                # print(s1,s2)
                s1, s2, stats = parasail_alignment_helper( ((s1, s2, i, j), {"mismatch_penalty" : mismatch_penalty }) )
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
        pool = Pool(processes=nr_cores)

        matches_with_mismatch = {}
        for j, s1 in enumerate(matches):
            matches_with_mismatch[s1] = {}
            for i, s2 in enumerate(matches[s1]):
                ed = matches[s1][s2]
                error_rate = float(ed)/ min(len(s1), len(s2))
                if error_rate <= 0.01:
                    mismatch_penalty = -1
                elif 0.01 < error_rate <= 0.09:
                    mismatch_penalty = -2
                else:
                    mismatch_penalty = -4

                matches_with_mismatch[s1][s2] = mismatch_penalty

        try:
            res = pool.map_async(parasail_alignment_helper, [ ((s1, s2, i,j), {"mismatch_penalty" : mismatch_penalty}) for j, s1 in enumerate(matches_with_mismatch) for i, (s2, mismatch_penalty) in enumerate(matches_with_mismatch[s1].items()) ] )
            alignment_results =res.get(999999999) # Without the timeout this blocking call ignores all signals.
        except KeyboardInterrupt:
            print("Caught KeyboardInterrupt, terminating workers")
            pool.terminate()
            sys.exit()
        else:
            # print("Normal termination")
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


def sw_align_sequences_keeping_accession(matches, nr_cores = 1):
    """
        Matches should be a 2D matrix implemented as a dict of dict, the value should be a tuple (s1,s2, edit_distance) .
    """
    print(len(matches))
    exact_matches = {}
    if nr_cores == 1:
        for j, s1_acc in enumerate(matches):
            for i, s2_acc in enumerate(matches[s1_acc]):
                s1, s2, ed = matches[s1_acc][s2_acc]
                if s1_acc in exact_matches:
                    if s2_acc in exact_matches[s1_acc]:
                        continue
                error_rate = float(ed)/ min(len(s1), len(s2))
                # print("ERROR_RATE", error_rate)
                if error_rate <= 0.01:
                    mismatch_penalty = -1
                elif 0.01 < error_rate <= 0.09:
                    mismatch_penalty = -2
                else:
                    mismatch_penalty = -4

                s1_acc, s2_acc, stats = parasail_alignment_helper( ((s1, s2, i, j), {"x_acc" : s1_acc, "y_acc" :s2_acc, "mismatch_penalty" : mismatch_penalty}) )

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
        pool = Pool(processes=nr_cores)

        matches_with_mismatch = {}
        for j, s1_acc in enumerate(matches):
            matches_with_mismatch[s1_acc] = {}
            for i, s2_acc in enumerate(matches[s1_acc]):
                s1, s2, ed = matches[s1_acc][s2_acc]
                error_rate = float(ed)/ min(len(s1), len(s2))
                # print("ERROR_RATE",error_rate ) #, s1_acc, ed, len(s1), len(s2) )
                if error_rate <= 0.01:
                    mismatch_penalty = -1
                elif 0.01 < error_rate <= 0.09:
                    mismatch_penalty = -2
                else:
                    mismatch_penalty = -4

                matches_with_mismatch[s1_acc][s2_acc] = mismatch_penalty
        print(len(matches_with_mismatch), "lewl")
        # for j, s1_acc in enumerate(matches):
        #     for i, s2_acc in enumerate(matches[s1_acc]):
        #         print("lool", matches[s1_acc][s2_acc][0], matches[s1_acc][s2_acc][1], i,j, {"x_acc": s1_acc, "y_acc" : s2_acc} ) 
        try:
            res = pool.map_async(parasail_alignment_helper, [ ((matches[s1_acc][s2_acc][0], matches[s1_acc][s2_acc][1], i,j), {"x_acc": s1_acc, "y_acc" : s2_acc, "mismatch_penalty" : mismatch_penalty}) for j, s1_acc in enumerate(matches_with_mismatch) for i, (s2_acc, mismatch_penalty) in enumerate(matches_with_mismatch[s1_acc].items()) ] )
            alignment_results =res.get(999999999) # Without the timeout this blocking call ignores all signals.
        except KeyboardInterrupt:
            print("Caught KeyboardInterrupt, terminating workers")
            pool.terminate()
            sys.exit()
        else:
            # print("Normal termination")
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
                # print(len(matches[s1_acc][s2_acc][0]), len(matches[s1_acc][s2_acc][1]) )
                pass
    print(len(exact_matches))

    return exact_matches


# def ssw_alignment_helper(arguments):
#     args, kwargs = arguments
#     return ssw_alignment(*args, **kwargs)

# def ssw_alignment(x, y, i,j, ends_discrepancy_threshold = 25 , x_acc = "", y_acc = "", mismatch_penalty = -3 ):
#     """
#         Aligns two sequences with SSW
#         x: query
#         y: reference

#     """
#     # if i == 100 and j % 1000 == 0:
#     #     print("SW processed alignments:{0}, mismatch_penalty: {1}".format(j+1, mismatch_penalty))

#     score_matrix = ssw.DNA_ScoreMatrix(match=2, mismatch=mismatch_penalty)
#     aligner = ssw.Aligner(gap_open=2, gap_extend=0, matrix=score_matrix)

#     # for the ends that SSW leaves behind
#     bio_matrix = matlist.blosum62
#     g_open = -1
#     g_extend = -0.5
#     ######################################

#     # result = aligner.align("GA", "G", revcomp=False)
#     # y_alignment, match_line, x_alignment = result.alignment
#     # c = Counter(match_line)
#     # matches, mismatches, indels = c["|"], c["*"], c[" "]
#     # alignment_length = len(match_line)
#     # print("matches:{0}, mismatches:{1}, indels:{2} ".format(matches, mismatches, indels))
#     # print(match_line)

#     result = aligner.align(x, y, revcomp=False)
#     y_alignment, match_line, x_alignment = result.alignment

#     matches, mismatches, indels = match_line.count("|"), match_line.count("*"), match_line.count(" ")

#     # alignment_length = len(match_line)
    
#     start_discrepancy = max(result.query_begin, result.reference_begin)  # 0-indexed # max(result.query_begin, result.reference_begin) - min(result.query_begin, result.reference_begin)
#     query_end_discrepancy = len(x) - result.query_end - 1
#     ref_end_discrepancy = len(y) - result.reference_end - 1
#     end_discrepancy = max(query_end_discrepancy, ref_end_discrepancy)  # max(result.query_end, result.reference_end) - min(result.query_end, result.reference_end)
#     # print(start_discrepancy, end_discrepancy)
#     tot_discrepancy = start_discrepancy + end_discrepancy

#     if 0 < start_discrepancy <= ends_discrepancy_threshold:
#         # print("HERE")
#         matches_snippet = 0
#         mismatches_snippet = 0
#         if result.query_begin and result.reference_begin:
#             query_start_snippet = x[:result.query_begin]
#             ref_start_snippet = y[:result.reference_begin]
#             alns = pairwise2.align.globalds(query_start_snippet, ref_start_snippet, bio_matrix, g_open, g_extend)
#             top_aln = alns[0]
#             # print(alns)
#             mismatches_snippet = len(list(filter(lambda x: x[0] != x[1] and x[0] != '-' and x[1] != "-", zip(top_aln[0],top_aln[1]))))
#             indels_snippet = top_aln[0].count("-") + top_aln[1].count("-")
#             matches_snippet = len(top_aln[0]) - mismatches_snippet - indels_snippet
#             # print(matches_snippet, mismatches_snippet, indels_snippet)
#             query_start_alignment_snippet = top_aln[0]
#             ref_start_alignment_snippet = top_aln[1]
#         elif result.query_begin:
#             query_start_alignment_snippet = x[:result.query_begin]
#             ref_start_alignment_snippet = "-"*len(query_start_alignment_snippet)
#             indels_snippet = len(ref_start_alignment_snippet)
#         elif result.reference_begin:
#             ref_start_alignment_snippet = y[:result.reference_begin]
#             query_start_alignment_snippet = "-"*len(ref_start_alignment_snippet)
#             indels_snippet = len(query_start_alignment_snippet)
#         else:
#             print("BUG")
#             sys.exit()
#         matches, mismatches, indels = matches + matches_snippet, mismatches + mismatches_snippet, indels + indels_snippet

#         # print(ref_start_alignment_snippet)
#         # print(query_start_alignment_snippet)
#         y_alignment = ref_start_alignment_snippet + y_alignment
#         x_alignment = query_start_alignment_snippet + x_alignment

#     if 0 < end_discrepancy <= ends_discrepancy_threshold:
#         # print("HERE2")
#         matches_snippet = 0
#         mismatches_snippet = 0
#         if query_end_discrepancy and ref_end_discrepancy:
#             query_end_snippet = x[result.query_end+1:]
#             ref_end_snippet = y[result.reference_end+1:]
#             alns = pairwise2.align.globalds(query_end_snippet, ref_end_snippet, bio_matrix, g_open, g_extend)
#             top_aln = alns[0]
#             mismatches_snippet = len(list(filter(lambda x: x[0] != x[1] and x[0] != '-' and x[1] != "-", zip(top_aln[0],top_aln[1]))))
#             indels_snippet = top_aln[0].count("-") + top_aln[1].count("-")
#             matches_snippet = len(top_aln[0]) - mismatches_snippet - indels_snippet
#             query_end_alignment_snippet = top_aln[0]
#             ref_end_alignment_snippet = top_aln[1]
#         elif query_end_discrepancy:
#             query_end_alignment_snippet = x[result.query_end+1:]
#             ref_end_alignment_snippet = "-"*len(query_end_alignment_snippet)
#             indels_snippet = len(ref_end_alignment_snippet)

#         elif ref_end_discrepancy:
#             ref_end_alignment_snippet = y[result.reference_end+1:]
#             query_end_alignment_snippet = "-"*len(ref_end_alignment_snippet)
#             indels_snippet = len(query_end_alignment_snippet)

#         else:
#             print("BUG")
#             sys.exit()
#         matches, mismatches, indels = matches + matches_snippet, mismatches + mismatches_snippet, indels + indels_snippet

#         y_alignment = y_alignment + ref_end_alignment_snippet
#         x_alignment = x_alignment + query_end_alignment_snippet


#     if x_acc == y_acc == "":
#         if start_discrepancy > ends_discrepancy_threshold or end_discrepancy > ends_discrepancy_threshold:
#             # print("REMOVING", start_discrepancy, end_discrepancy)
#             return (x, y, None)

#         else:
#             return (x, y, (x_alignment, y_alignment, (matches, mismatches, indels)) )
#     else:
#         if start_discrepancy > ends_discrepancy_threshold or end_discrepancy > ends_discrepancy_threshold:
#             # print("REMOVING", start_discrepancy, end_discrepancy)
#             return (x_acc, y_acc, None)

#         else:

#             return (x_acc, y_acc, (x_alignment, y_alignment, (matches, mismatches, indels)) )       

