# from __future__ import unicode_literals
import gzip
import sys
import heapq

from hitemmodules.io import fasta_parser
from hitemmodules import misc_functions
from hitemmodules.io import write_output
from hitemmodules import align
import operator
from collections import defaultdict

def store_x(read_file):
    """
        Initialize reads.
    """
    x = dict([(acc, seq) for acc, seq in fasta_parser.read_fasta(read_file)])
    return x


def initialize_y(reads):
    """
        Initialize transcript instances y as the original reads at time 0.
        That is, on position 0 in the list for each y
    """
    y = {}
    x_to_y = {}
    y_to_x = {}
    for i, (acc, seq) in enumerate(reads.items()):
        y[i] = seq
        x_to_y[acc] = i
        y_to_x[i] = acc
    return y, x_to_y, y_to_x

def initialize_paf_x_y(paf_file_path, x_to_y, params):
    """
        Store container of paf alignments between reads and initial y,
        i.e., all pairwise self alignments as y=x in first step.
    """
    edges_y_x = {}
    params.expected_support = defaultdict(float)
    # temp_reads_with_nbrs = set()
    # temp_paf_nr_cutoff = {} # x : heap[(y, paf_similarity_score), ... ]
    try:
        file_object = gzip.open(paf_file_path)
        file_object.readline()
        file_object.seek(0)
    except IOError:
        file_object = open(paf_file_path)

    highest_paf_scores = defaultdict(list)
    reads_seen_in_paf = set()
    with file_object as paf_file:
        it = 0
        for line in paf_file:
            it += 1
            if it % 1000000 == 0:
                print("parsing line {0}".format(it))
                # if it == 100000000000:
                #     break
            row_info = line.strip().split()
            q_acc = row_info[0].decode('ascii')
            q_len = int(row_info[1])
            t_acc = row_info[5].decode('ascii')
            t_len = int(row_info[6])
            # print(row_info[0][-10:],row_info[1:5], row_info[5][-10:], row_info[6:], nr_matching/float(max(q_len, t_len)))
            # print(q_acc, row_info[1:5], t_acc, row_info[6:], int(row_info[9])/float(max(q_len, t_len)))
            # self mapping -> bad
            if q_acc == t_acc:
                print("SELF MAPPING DETECTED")
                print(q_acc, q_len, row_info)
                print(t_acc, t_len, row_info)
                continue

            reads_seen_in_paf.add(q_acc)
            reads_seen_in_paf.add(t_acc)

            if q_len < t_len:
                diff_ratio = float(t_len)/float(q_len)
            else: # q_len equal or longer
                diff_ratio = float(q_len)/float(t_len)

            # 10% length cutoff for practical reasons
            # its very unlikely that they will be from the same transcript in that case.
            # and if that extremely unlikely event happen, it shouldnt cause big trouble 
            # since they hopefully converge to the same transcript anyway. 
            if diff_ratio > 1.1: 
                continue

            nr_matching = int(row_info[9])
            paf_similarity_score = nr_matching/float(max(q_len, t_len))

            # print(q_acc, paf_similarity_score, nr_matching, q_len, t_len, t_acc)
            # if paf_similarity_score < params.paf_similarity:
            #     continue

            params.expected_support[q_acc] += paf_similarity_score 
            params.expected_support[t_acc] += paf_similarity_score 

            # temp_reads_with_nbrs.add(q_acc)
            # temp_reads_with_nbrs.add(t_acc)

            y = x_to_y[t_acc]
            if y not in edges_y_x:
                edges_y_x[y] = {q_acc : paf_similarity_score}
                heapq.heappush(highest_paf_scores[y], (paf_similarity_score, q_acc) )
            else:
                if q_acc not in edges_y_x[y]:
                    # if we already have paf_y_limit reads mapped to y
                    if len(highest_paf_scores[y]) >= params.paf_y_limit:
                        # current alignment is better than at least one of previous scores, remove the worst one so far
                        if paf_similarity_score > highest_paf_scores[y][0][0]: 
                            edges_y_x[y][q_acc] = paf_similarity_score
                            paf_score, q_acc_out = heapq.heappushpop(highest_paf_scores[y], (paf_similarity_score, q_acc) )
                            del edges_y_x[y][q_acc_out]
                    else:
                        edges_y_x[y][q_acc] = paf_similarity_score
                        heapq.heappush(highest_paf_scores[y], (paf_similarity_score, q_acc))

                else:
                    if edges_y_x[y][q_acc] < paf_similarity_score:
                        edges_y_x[y][q_acc] = paf_similarity_score
                        print("found a better score of the same alignment..")
                        # print(y, t_acc, q_acc )


            y = x_to_y[q_acc]
            if y not in edges_y_x:
                edges_y_x[y] = {t_acc : paf_similarity_score}
                heapq.heappush(highest_paf_scores[y], (paf_similarity_score, t_acc) )
            else:
                if t_acc not in edges_y_x[y]:
                    # if we already have paf_y_limit reads mapped to y
                    if len(highest_paf_scores[y]) >= params.paf_y_limit:
                        # current alignment is better than at least one of previous scores, remove the worst one so far
                        if paf_similarity_score > highest_paf_scores[y][0][0]: 
                            edges_y_x[y][t_acc] = paf_similarity_score
                            paf_score, t_acc_out = heapq.heappushpop(highest_paf_scores[y], (paf_similarity_score, t_acc) )
                            del edges_y_x[y][t_acc_out]
                    else:
                        edges_y_x[y][t_acc] = paf_similarity_score
                        heapq.heappush(highest_paf_scores[y], (paf_similarity_score, t_acc))
                else:
                    # y is already among the best scores
                    if edges_y_x[y][t_acc] < paf_similarity_score:
                        edges_y_x[y][t_acc] = paf_similarity_score
                        print("found a better score of the same alignment..")


    write_output.logger("Filter reads in due to paf_x_limit set to {0}".format(params.paf_x_limit), params.logfile)
    write_output.logger("Filter reads in due to paf_y_limit set to {0}".format(params.paf_y_limit), params.logfile)

    # prune away edges for paf_y_limit: limit nr of reads x that can contribute to each y, at most k reads can contribute to each y_j
    for y_j in edges_y_x:
        nbrs_sorted_highest_score  = sorted(edges_y_x[y_j], key=lambda nbr: edges_y_x[y_j][nbr], reverse=True)
        nbrs_to_save = nbrs_sorted_highest_score[:params.paf_y_limit]
        for x_i in edges_y_x[y_j].keys():
            if x_i not in nbrs_to_save:
                del edges_y_x[y][x_i]


    edges_x_y = misc_functions.transpose(edges_y_x)
    reads_used = set([x_i for x_i in edges_x_y])
    reads_not_supporting_starting_points = set(x_to_y.keys()).difference(reads_used)
    reads_not_observed_in_paf = set(x_to_y.keys()).difference(reads_seen_in_paf)

    write_output.logger("Total number of reads in fasta:{0}".format(len(x_to_y)), params.logfile)
    write_output.logger("Total number of reads seen in paf:{0}".format(len(reads_seen_in_paf)), params.logfile)

    # print(len(reads_not_observed_in_alignment), reads_not_observed_in_alignment)
    # write_output.logger("Number of reads in PAF file that had at least one neighbor over the similarity threshold (either query or target:{0})".format(len(reads_with_alignments)), params.logfile)
    write_output.logger("Number of reads not used, hence also not contributing to starting points:{0}.".format(len(reads_not_supporting_starting_points)), params.logfile)
    write_output.logger("Number of reads not seen in PAF alignment file:{0}.".format(len(reads_not_observed_in_paf)), params.logfile)
    write_output.logger("Number of reads used in support: {0}.".format(len(edges_x_y)), params.logfile)

    alignment_results_dict = edges_y_x 
    # for e1 in edges_y_x:
    #     for e2 in edges_y_x[e1]:
    #         if edges_y_x[e1][e2]  > 0.6:
    #             print("approx sim:", edges_y_x[e1][e2] )
    # sys.exit()
    return alignment_results_dict, reads_not_observed_in_paf

# def initialize_non_zero_p_x_y(edges_x_y, params):
#     """
#         Only consider probabilities where the mapper has found at least a partial alignment.
#         If the aligner has not found at least a partial alignment, the probability that the sequences are from
#         the same transcript is astronomically low.
#     """
#     return non_zero_p_x_y

def initialize_epsilon(x, alignment_results_dict_transposed, y_to_x, params):
    """
        Initialize base uncertainty based on alignments between reads and initial y,
        i.e., all pairwise self alignments as y=x in first step.
    """

    epsilon_x_min = {}
    errors_x_min = {}
    # assert 'gi|768041113|ref|XM_011531492.1|_PREDICTED:_Homo_sapiens_RNA_binding_motif_protein,_Y-linked,_family_1,_member_E_(RBMY1E),_transcript_variant_X1,_mRNA:read_199' not in x
    epsilon_temp = {acc : 1.001 for acc in x}
    errors_temp = {acc : len(x[acc]) for acc in x}

    x_left = set(x.keys())
    x_aligned = set(alignment_results_dict_transposed.keys())
    x_minus_x_aligned = x_left.difference(x_aligned)
    x_aligned_minus_x = x_aligned.difference(x_left)

    print("x minus aligned:")
    for x_i in x_minus_x_aligned:
        print("x minus aligned", x_i)
    print("aligned minus x:")
    for x_i in x_aligned_minus_x:
        print("aligned minus x", x_i)


    for x_i in alignment_results_dict_transposed:
        # x_i_ref = y_to_x[y_j]
        if x_i not in x:
            print("BUG", x_i, "not in x")
            continue

        for y_j in alignment_results_dict_transposed[x_i]:
            y_j_aligned, x_i_aligned, stats = alignment_results_dict_transposed[x_i][y_j]
            matches, mismatches, indels = stats 
            eps = (mismatches + indels) / float(matches + mismatches + indels)
            err = mismatches + indels
            if eps < epsilon_temp[x_i]:
                epsilon_temp[x_i] = eps
                x_minimizer =  y_to_x[y_j]
                epsilon_x_min[x_i] = x_minimizer

            if err < errors_temp[x_i]:
                errors_temp[x_i] = err
                x_minimizer =  y_to_x[y_j]
                errors_x_min[x_i] = x_minimizer

            if err == 0:
                break
        # if eps == 0:
        #     write_output.logger("Epsilon estimated to 0 for read {0}, read length: {1}".format(x_i, len(x[x_i])), params.logfile)

        # if errors_x_min[x_i] != epsilon_x_min[x_i]:
        #     write_output.logger("Min epsilon (eps:{0}, len read:{1}) and min number of errors ({2}) differ!".format(epsilon_temp[x_i], len(x[x_i]), errors_temp[x_i]), params.logfile)


    epsilon_advanced = defaultdict(float)
    errors_advanced = defaultdict(float)
    # print(sorted(epsilon_temp.items(), key=operator.itemgetter(1)))

    for acc in x:
        if acc in epsilon_temp and acc not in epsilon_x_min:
            print("LOOOOL", acc)
    print("x_minus_x_aligned:", x_minus_x_aligned)
            # print("nr of alignments:", len(alignment_results_dict_transposed[acc]))
            # print("Alignment dict:", alignment_results_dict_transposed[acc])

    for x_i_ref, epsilon in sorted(epsilon_temp.items(), key=operator.itemgetter(1)):
        if x_i_ref in x_minus_x_aligned:
            print("did not estimate error or epsilon for read:", x_i_ref)
            continue
        x_minimizer = epsilon_x_min[x_i_ref]
        if x_minimizer not in epsilon_advanced:
            # need baseline value for smallest epsilon. Last term is because errors can occur ofver the same positions
            # hence epsilon_temp[x_i_ref]/2.0 will underestimate error rate.
            epsilon_advanced[x_i_ref] = epsilon_temp[x_i_ref]/2.0 + (epsilon_temp[x_i_ref]/2.0)**2
        else:
            # epsilon is estimated as max(epsilon_observed - epsilon_of_read_with_already_stored_smallest_epsilon, epsilon_of_read_with_already_stored_smallest_epsilon)
            # (first argument is by definition positive), second one ensures that is is not smaller epsilon than its minimizer. This can happen in the following scenario:
            # Assume all transcripts are 10 bp long for easy example
            # ED(x_0,x_1)=0, so e_0, e_1 is set to 0.0 (they are error free and from same isoform). Now consider x_2 which is also error free but is from an allele with a SNP
            # that is, we have e_1, e_2 = 0.0 are the true values. However, we observe ED(x_1,x_2) = 1. Now, lets also say that we have a 
            # read x_3 with one error but from the same copy as x_2. We have ED(x_2,x_3) = 1, with x_3 containng one error.
            # If epsilon is set in this ordder: e_0, e_1, e_2, e_3, we will have: e_0, e_1 = 0. e_2 = 0.1 - e_1 = 0.1. Consecuently, we have e_3 = 0.1 -0.1 = 0.0
            # which is wrong. We can hedge for these cases but enforcing the lower boundary e_2 when estimating e_3, which is what we do.

            # again, Last term is because errors can occur ofver the same positions
            # hence epsilon_temp[x_i_ref]/2.0 will underestimate error rate.
            epsilon_advanced[x_i_ref] = max(epsilon_temp[x_i_ref] - epsilon_advanced[x_minimizer], epsilon_advanced[x_minimizer]) + max(epsilon_temp[x_i_ref] - epsilon_advanced[x_minimizer], epsilon_advanced[x_minimizer])**2

        x_minimizer = errors_x_min[x_i_ref]
        if x_minimizer not in errors_advanced:
            # need baseline value for smallest epsilon. Last term is because errors can occur over the same positions
            # hence errors_temp[x_i_ref]/2.0 will underestimate error rate.
            errors_advanced[x_i_ref] = round(errors_temp[x_i_ref]/2.0)
        else:
            errors_advanced[x_i_ref] = float(max(errors_temp[x_i_ref] - errors_advanced[x_minimizer], errors_advanced[x_minimizer]))


        # if x_i_ref == "m151209_004839_42146_c100926392550000001823199905121692_s1_p0/3554/1665_60_CCS":
        #     print('epsilon minimizer:', epsilon_advanced[x_minimizer], "acc:", x_minimizer)

    for key, val in sorted(epsilon_temp.items(), key=operator.itemgetter(1)):
        if key in epsilon_advanced:
            print(key, "temp:", val, "advanced:", epsilon_advanced[key], "nr errors:", errors_advanced[key])
        else: 
            print(key, " is not in epsilon advanced!!")

    # if "m151209_004839_42146_c100926392550000001823199905121692_s1_p0/3554/1665_60_CCS" in epsilon_temp:
    #     print("m151209_004839_42146_c100926392550000001823199905121692_s1_p0/3554/1665_60_CCS", "temp", epsilon_temp["m151209_004839_42146_c100926392550000001823199905121692_s1_p0/3554/1665_60_CCS"])
    # if "m151209_004839_42146_c100926392550000001823199905121692_s1_p0/3554/1665_60_CCS" in epsilon_advanced:
    #     print("m151209_004839_42146_c100926392550000001823199905121692_s1_p0/3554/1665_60_CCS", "advanced", epsilon_advanced["m151209_004839_42146_c100926392550000001823199905121692_s1_p0/3554/1665_60_CCS"])

    return epsilon_advanced, errors_advanced

def initialize_epsilon_y(epsilon, y, x_to_y, params):
    epsilon_y = {}
    for x_i, y_j in x_to_y.items():
        if y_j in y:
            epsilon_y[y_j] = epsilon[x_i]
    return epsilon_y

def initialize_errors_y(errors_x, y, x_to_y, params):
    errors_y = {}
    for x_i, y_j in x_to_y.items():
        if y_j in y:
            errors_y[y_j] = errors_x[x_i]
    return errors_y


def estimate_epsilon(x, alignment_results_dict, alignment_results_dict_transposed, y_to_x, x_to_y, epsilon_backup, errors_backup, params):
    """
        Initialize base uncertainty based on alignments between sequences
    """

    pairwise_diff_epsilon = {}
    pairwise_diff_errors = {}

    x_left = set(x.keys())
    x_aligned = set(alignment_results_dict_transposed.keys())
    x_minus_x_aligned = x_left.difference(x_aligned)
    x_aligned_minus_x = x_aligned.difference(x_left)

    print("In estimate_epsilon x minus aligned:")
    for x_i in x_minus_x_aligned:
        print("x minus aligned", x_i)
    print("In estimate_epsilon aligned minus x:")
    for x_i in x_aligned_minus_x:
        print("aligned minus x", x_i)


    epsilons = {acc : [] for acc in alignment_results_dict_transposed}
    errors = {acc : [] for acc in alignment_results_dict_transposed}
    for x_i in alignment_results_dict_transposed:
        if x_i not in x:
            print("BUG", x_i, "not in x")
            continue
        almnts = alignment_results_dict_transposed[x_i].items()
        x_3 = x_i
        best_y_to_x_3 = sorted(almnts, key = lambda z: z[1][2][1] + z[1][2][2])[0][0] 
        matches, mismatches, indels = alignment_results_dict_transposed[x_3][best_y_to_x_3][2]
        err_x_1_x_3 = mismatches + indels
        eps_x_1_x_3 = (mismatches + indels) / float(matches + mismatches + indels)

        x_1 =  y_to_x[ best_y_to_x_3 ]
        # print("xxxxx_1",x_1)
        pairwise_diff_epsilon[best_y_to_x_3] = {}
        pairwise_diff_errors[best_y_to_x_3] = {}
        if len(alignment_results_dict[best_y_to_x_3]) == 1:
            errors[x_3] = err_x_1_x_3
            epsilons[x_3] = eps_x_1_x_3
            print("only one y", "errors:", err_x_1_x_3)
            continue

        for x_2_candidate in alignment_results_dict[best_y_to_x_3]:
            x_2_candidate_aligned, x_1_aligned, stats = alignment_results_dict[best_y_to_x_3][x_2_candidate]
            matches, mismatches, indels = stats 
            eps = (mismatches + indels) / float(matches + mismatches + indels)
            err = mismatches + indels
            # print("x_2_candidate", x_2_candidate)
            pairwise_diff_epsilon[best_y_to_x_3][x_2_candidate] = eps
            pairwise_diff_errors[best_y_to_x_3][x_2_candidate] = err
    
        # pairwise_diff_epsilon[x_1] = sorted(pairwise_diff_epsilon[x_1], key=lambda z: z[1])[:2]
        # pairwise_diff_errors[x_1] =  sorted(pairwise_diff_errors[x_1], key=lambda z: z[1])[:2]
        # print("top 2 eps:", pairwise_diff_epsilon[x_1])
        # print("top2 err:", pairwise_diff_errors[x_1])

        # errors
        best_2_matches = sorted(pairwise_diff_errors[best_y_to_x_3].items(), key=lambda z: z[1])[:2]
        x_2 = best_2_matches[0][0] if best_2_matches[0][0] != x_3 else best_2_matches[1][0]
        err_x_1_x_2 = pairwise_diff_errors[best_y_to_x_3][x_2]

        # if x_2 in alignment_results_dict_transposed and x_to_y[x_3] in alignment_results_dict_transposed[x_2]:
        #     y_j_aligned, x_i_aligned, stats = alignment_results_dict_transposed[x_2][x_to_y[x_3]]
        #     print("here1")
        # else:
        x_3_acc, x_2_acc, alignments  = align.ssw_alignment_helper( (x_2, x_3, x[x_2], x[x_3], 2,2) )
        if alignments:
            stats = alignments[2]
        else:
            print(x_1)
            print(x_2_acc, x[x_2])
            print(x_3_acc,x[x_3])
            print("ERRORS:", err_x_1_x_2, err_x_1_x_3)
            print("Nr nbrs", len(pairwise_diff_errors[best_y_to_x_3]), "Neighbours:", pairwise_diff_errors[best_y_to_x_3])
            print("forcing alignment with 200 in end_discrepancy threshold")
            x_3_acc, x_2_acc, alignments  = align.ssw_alignment_helper( (x_2, x_3, x[x_2], x[x_3], 2,2, 200) )
            if alignments:
                stats = alignments[2]
            else:
                # assign error and epsilon from backup here and continue
                errors[x_3] = errors_backup[x_3]
                epsilons[x_3] =  epsilon_backup[x_3]
                print("Assigning epsilon and error from backup!")
                continue

        # print("here2")


        matches, mismatches, indels = stats
        # print(stats) 
        err_x_2_x_3 = mismatches + indels
        error_x_3 = (err_x_2_x_3 + (err_x_1_x_3 - err_x_1_x_2)) / 2
        if (2*error_x_3) % 2 == 1:
            print(err_x_1_x_2, err_x_1_x_3, err_x_2_x_3, error_x_3, x_1, x_2, x_3)
        # error_x_1_biased = err_x_1_x_3 - error_x_3
        # error_x_2 = err_x_1_x_2 - error_x_1_biased

        # errors[x_2].append(error_x_2) 
        errors[x_3] = error_x_3

        # if error_x_3 > 10:
        #     print("x_1", x_1, x[x_1])
        #     print("x_2",x_2_acc, x[x_2])
        #     print("x_3", x_3_acc,x[x_3])
        #     print("ERRORS:", err_x_1_x_2, err_x_1_x_3, err_x_2_x_3)
        #     print("Nr nbrs", len(pairwise_diff_errors[best_y_to_x_3]), "Neighbours:", pairwise_diff_errors[best_y_to_x_3])
        #     print("Nr nbrs x3", len(almnts), "Neighbours:", almnts)


        # epsilon
        # best_2_matches = sorted(pairwise_diff_epsilon[x_1].items(), key=lambda z: z[1])[:2]
        # x_2 = best_2_matches[0][0] if best_2_matches[0][0] != x_3 else best_2_matches[1][0]
        eps_x_1_x_2 = pairwise_diff_epsilon[best_y_to_x_3][x_2]

        # if x_2 in alignment_results_dict_transposed and x_to_y[x_3] in alignment_results_dict_transposed[x_2]:
        #     y_j_aligned, x_i_aligned, stats = alignment_results_dict_transposed[x_2][x_to_y[x_3]]
        # else:
        #     y_temp = x_to_y[x_3]
        #     y_j, x_i, stats = align.ssw_alignment_helper( (x_2, y_temp, x[x_2], y[y_temp], 2,2) )
        # matches, mismatches, indels = stats 

        eps_x_2_x_3 = (mismatches + indels) / float(matches + mismatches + indels)

        epsilon_x_3 = (eps_x_2_x_3 + (eps_x_1_x_3 - eps_x_1_x_2)) / 2
        # epsilon_x_1_biased = eps_x_1_x_3 - epsilon_x_3
        # epsilon_x_2 = eps_x_1_x_2 - epsilon_x_1_biased

        # epsilons[x_2].append(epsilon_x_2) 
        epsilons[x_3] = epsilon_x_3 

    errors_final = {}
    epsilons_final = {}
    for x_i, eps_x_i in sorted(errors.items(), key=operator.itemgetter(1)):
        print("error estimates", x_i, ":", errors[x_i])
        errors_final[x_i] = errors[x_i]
        # print("epsion estimates", x_i, ":", epsilons[x_i])
        epsilons_final[x_i] = epsilons[x_i] 
    # sys.exit()
    return epsilons_final, errors_final
