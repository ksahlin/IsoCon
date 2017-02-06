import signal
from multiprocessing import Pool
import multiprocessing as mp
import sys

import math
from scipy.stats import poisson, binom, norm

from modules.multinomial_distr import multinomial_
from modules.SW_alignment_module import sw_align_sequences_keeping_accession
from modules.functions import create_position_probability_matrix, get_error_rates_and_lambda, get_difference_coordinates_for_candidates, get_supporting_reads_for_candidates, get_invariant_adjustment


def do_statistical_tests(null_hypothesis_references_t, C_seq_to_acc, partition_of_X, partition_of_C, X, C, single_core = False):
    # statistical_test(t, C_seq_to_acc, partition_of_X, partition_of_C, X, C)
    p_values = {}
    if single_core:
        for t in null_hypothesis_references_t:
            p_vals = statistical_test(t, C_seq_to_acc, partition_of_X, partition_of_C, X, C)
            for c_acc, (t_acc, k, p_val, N_t) in p_vals.items():
                assert c_acc not in p_values
                p_values[c_acc] = (t_acc, k, p_val, N_t)
    else:
        ####### parallelize statistical tests #########
        # pool = Pool(processes=mp.cpu_count())
        original_sigint_handler = signal.signal(signal.SIGINT, signal.SIG_IGN)
        signal.signal(signal.SIGINT, original_sigint_handler)
        pool = Pool(processes=mp.cpu_count())
        try:
            res = pool.map_async(statistical_test_helper, [ ( (t, C_seq_to_acc, partition_of_X, partition_of_C, X, C), {}) for t in null_hypothesis_references_t  ] )
            statistical_test_results =res.get(999999999) # Without the timeout this blocking call ignores all signals.
        except KeyboardInterrupt:
            print("Caught KeyboardInterrupt, terminating workers")
            pool.terminate()
            sys.exit()
        else:
            print("Normal termination")
            pool.close()
        pool.join()
        for candidate_dict in statistical_test_results:
            for c_acc, (t_acc, k, p_val, N_t) in candidate_dict.items(): 
                assert c_acc not in p_values
                p_values[c_acc] = (t_acc, k, p_val, N_t)

    return p_values


def statistical_test_helper(arguments):
    args, kwargs = arguments
    return statistical_test(*args, **kwargs)

def statistical_test(t, C_seq_to_acc, partition_of_X, partition_of_C, X, C):
    t_acc = C_seq_to_acc[t]
    print("t length:", len(t))
    reads_to_t = [x_acc for x_acc in  partition_of_X[t_acc] ]
    reads_to_map = set(reads_to_t)
    reads_to_map.update([t_acc])

    for c in partition_of_C[t]:
        print("c len:", len(c))
        c_acc = C_seq_to_acc[c]
        reads_to_c = [x_acc for x_acc in  partition_of_X[c_acc] ]
        # print(type(reads_to_c))
        reads_to_map.update(reads_to_c) 
        reads_to_map.update([c_acc]) # ad the candidate consensus too!
    
    # p_values[t_acc] = (t_acc, len(reads_to_t), -1, len(reads_to_map)) 
    # print(len(reads_to_map))

    # create alignment matrix A of reads to t
    graph_partition_t = {t_acc : reads_to_map }
    partition_alignments_t = get_partition_alignments_2set(graph_partition_t, C, X)
    alignment_matrix_to_t, PFM_to_t = create_position_probability_matrix(t, partition_alignments_t[t_acc])
    print("done", len(alignment_matrix_to_t), "of which is candidates:", len(partition_of_C[t]))


    candidate_accessions = set([C_seq_to_acc[c] for c in partition_of_C[t]])

    # Get all positions in A where c and t differ, as well as the state and character
    delta_t = get_difference_coordinates_for_candidates(t_acc, candidate_accessions, alignment_matrix_to_t) # format: { c_acc1 : {pos:(state, char), pos2:(state, char) } , c_acc2 : {pos:(state, char), pos2:(state, char) },... }

    # get error rates  # e_s, e_i,e_d = ..... from matrix
    # epsilon format: { x_acc1 : {state : prob}, x_acc2 : {state, prob} ,... }
    # also get base pair errors distribution (approximated as poisson), 
    # estimate this from all positions andvariants except for the variant and position combination of the candidates


    # here we do two separate tests: against cluster center and against closest candidate (they may be the same).
    p_values_on_closest_highest_support = test_against_center(delta_t, alignment_matrix_to_t, t_acc, t, candidate_accessions, partition_of_X, partition_of_C, C)


    # ###########################################################################
    # ###########################################################################
    # ###########################################################################

    # # p_values_on_closest = {}
    # all_c_in_partition = candidate_accessions.union({t_acc})
    
    # min_to_t_acc = min(delta_t, key=lambda x: len(delta_t[x]))


    # # for c in all_c_in_partition:
    # #     # find reference here: this should be the closest candidate in partition
    # #     # print(c)
    # #     # print(all_c_in_partition.difference({c}))
    # #     delta_c = get_difference_coordinates_for_candidates(c, all_c_in_partition.difference({c}), alignment_matrix_to_t)

    # #     min_distance = 2**30
    # #     highest_support = 0
    # #     for c2_acc, delta_to_c in delta_c.items():
    # #         c2_distance = len(delta_to_c)
    # #         if c2_distance < min_distance:
    # #             min_distance = c2_distance
    # #             highest_support = len(partition_of_X[c2_acc])
    # #             min_to_t_acc = c2_acc

    # #         elif c2_distance == min_distance:
    # #             support = len(partition_of_X[c2_acc])
    # #             if support > highest_support:
    # #                 min_distance = c2_distance
    # #                 highest_support = support
    # #                 min_to_t_acc = c2_acc



    #     # p_value_min = 2.0

    #     # test against highest supported minimizer
    #     # for c2_acc, delta_to_c in delta_c.items():
    #     #     if len(delta_to_c) == min_distance:   
    #     # print("test", c, "against", min_to_t_acc)                    
    #     # min_to_t_acc = min(delta_c, key=lambda x: len(delta_c[x]))
    # delta_t_single = get_difference_coordinates_for_candidates(min_to_t_acc, {t_acc}, alignment_matrix_to_t)
    #     # print("minimizer:", len(delta_c[min_to_t_acc]))
    # candidate_accessions_single = set({t_acc})
    #     # delta_t_single = {c : delta_c[c]}

    # relevant_accessions = partition_of_X[min_to_t_acc].union(partition_of_X[t_acc])
    # relevant_accessions.update([t_acc, min_to_t_acc ])
    # alignment_matrix_to_t_single = {acc : alignment_matrix_to_t[acc] for acc in alignment_matrix_to_t if acc in relevant_accessions}
    # p_value_on_closest = test_against_closest(delta_t_single, alignment_matrix_to_t_single, min_to_t_acc, C[min_to_t_acc], candidate_accessions_single, partition_of_X, partition_of_C, C)
    # # if p_vaxlue_on_closest[c][2] < p_value_min:
    # # info_tuple =  p_value_on_closest[c] #(t_acc, k, p_value, N_t)
    # p_values_on_closest[t_acc] = p_value_on_closest[t_acc]
    # # p_values_on_center[t] = p_value_on_closest[t]

    # final_p_values = {}
    # for c_acc in candidate_accessions.union({t_acc}):
    #     # if c_acc == 'read34_support_7':
    #         # print("OKKKKKK", t_acc)
    #         # print(p_values_on_center[c_acc])
    #         # print( p_values_on_closest[c_acc])
    #         # print()
    #     # if c_acc != t_acc:
    #     #     (t_acc_center, k_center, p_value_center, N_t_center) =  p_values_on_center[c_acc]
    #     # else:
    #     #     (t_acc_center, k_center, p_value_center, N_t_center) = ("t", 2**30, -100, 2**30)

    #     (t_acc_closest, k_closest, p_value_closest, N_t_closest) =  p_values_on_closest[c_acc]
    #     if False: #p_value_center > p_value_closest:
    #         final_p_values[c_acc] = (t_acc_center, k_center, p_value_center, N_t_center)
    #     else:
    #         # print("P value for closest used:")
    #         # print(t_acc_center, k_center, p_value_center, N_t_center)
    #         # print(t_acc_closest, k_closest, p_value_closest, N_t_closest)
    #         final_p_values[c_acc] = (t_acc_closest, k_closest, p_value_closest, N_t_closest)

    # ###########################################################################
    # ###########################################################################
    # ###########################################################################

    return p_values_on_closest_highest_support


def choose(n, k):
    """
    A fast way to calculate binomial coefficients by Andrew Dalke (contrib).
    """
    if 0 <= k <= n:
        ntok = 1
        ktok = 1
        for t in xrange(1, min(k, n - k) + 1):
            ntok *= n
            ktok *= t
            n -= 1
        return ntok // ktok
    else:
        return 0

def get_partition_alignments_2set(graph_partition, C, X):
    partition_dict = {}
    for t, x_hits in graph_partition.items():
        partition_dict[t] = {}
        for x in x_hits:
            if x in X: # read assigned to alignment matrix
                partition_dict[t][x] = (C[t], X[x]) 
            else: # a candidate we incluse to test, it needs to be aligned w.r.t. the alignment matrix
                partition_dict[t][x] = (C[t], C[x]) 

    exact_alignments = sw_align_sequences_keeping_accession(partition_dict, single_core = True)
    partition_alignments = {} 

    for c in exact_alignments:
        partition_alignments[c] = {}
        for x in exact_alignments[c]:
            aln_c, aln_x, (matches, mismatches, indels) = exact_alignments[c][x]
            edit_dist = mismatches + indels
            partition_alignments[c][x] = (edit_dist, aln_c, aln_x, 1)

        # BUG keyerror later if an exact alignment is not found here. e.g., best alignment of read had a valid exact alignment to a candidate
        # but not when aligned to t... fix here.. actually, implement barcode detection first, this might clean up things for for targeted data
        # at least...

    # for c in partition_alignments:
    #     if len(partition_alignments[c]) < 1: 
    #         del partition_alignments[c]
    print("NR candidates with at least one hit:", len(partition_alignments))
    return partition_alignments


def test_against_center(delta_t, alignment_matrix_to_t, t_acc, t, candidate_accessions, partition_of_X, partition_of_C, C):
    p_values = {}
    epsilon, lambda_S, lambda_D, lambda_I = get_error_rates_and_lambda(t_acc, len(t), candidate_accessions, alignment_matrix_to_t) 
  
    # get number of reads k supporting the given set of variants, they have to support all the variants within a candidate
    candidate_support = get_supporting_reads_for_candidates(t_acc, candidate_accessions, alignment_matrix_to_t, delta_t, partition_of_X) # format: { c_acc1 : [x_acc1, x_acc2,.....], c_acc2 : [x_acc1, x_acc2,.....] ,... }

    # get the factor to which we divide the candidate support with, this is defined as 
    # u_delta = min_{d in Delta} {u_d}  where u_d is the number of possible insertions/deletions to create an identical alignment 
    # format:  { c_acc1 : u_delta, c_acc2 : u_delta, ... , }
    invariant_factors = get_invariant_adjustment(delta_t, alignment_matrix_to_t, t_acc)

    # print()
    # print("INVARIANT FACTORS:", invariant_factors)
    # print()
    # get p_value
    m = len(t)

    for c_acc in candidate_accessions:
        u_c = invariant_factors[c_acc]
        k = len(candidate_support[c_acc]) #/ float(u_c)
        N_t = len(alignment_matrix_to_t) - len(candidate_accessions) - 1 # all reads minus all candidates and the reference transcript
        # print("reads N_t:", N_t)
        # print("varinats:",delta_t[c_acc].items())

        ############################################################################################
        ############################################################################################


        if k <= 1:
            # print("NO support!")
            # print("lengths:", len(t), len(C[c_acc]))                    
            p_value = 1
        else:
            k_I, k_S, k_D = k, k, k
            p_I = poisson.sf(k_I - 1, u_c*lambda_I)  # SF defined as 1-P(S_I > k), we want 1-P(S_I >= k)
            p_S = poisson.sf(k_S - 1, lambda_S)
            p_D = poisson.sf(k_D - 1, u_c*lambda_D)

            x_S, x_D, x_I = 0, 0, 0
            for pos, (state, char) in delta_t[c_acc].items():
                if state == "S":
                    x_S += 1
                elif state == "D":
                    x_D += 1
                if state == "I":
                    x_I += 1

            # special cases that are obvious or can be accurately approximated
            if (p_I == 0.0 and x_I > 0) or (p_S == 0.0 and x_S > 0) or (p_D == 0.0 and x_D > 0):
                # print("here approx", N_t)
                # print("lambda:", lambda_D, lambda_S, lambda_I)
                # print("k:",k, x_S, x_D, x_I,p_S, p_D, p_I)
                # print("lengths:", len(t), len(C[c_acc]))
                p_value = 0.0
            elif (p_I + p_D + p_S)*m >= 10 :
                # approximate with normal
                p_bin = min(0.99, (p_I + p_D + p_S)) # cannot be larger than one, but may be due to approximation errors when lambda is huge
                mu = p_bin*m
                sigma = math.sqrt( (1 - p_bin)* p_bin*m )
                p_value = norm.sf(x_S + x_D + x_I , loc=mu , scale=sigma)
                # print("LOOOOL NORMAL approx:")
                # print("lambda:", lambda_D, lambda_S, lambda_I)
                # print("k:",k, x_S, x_D, x_I,p_S, p_D, p_I)
                # print("Approx p-val: ", p_value)
                # print("lengths:", len(t), len(C[c_acc]))
            elif x_S + x_D + x_I > 10 and (p_I + p_D + p_S)*m < 10 :
                # approximate with poisson
                lambda_prob = p_I + p_D + p_S
                # print("LOOOOL poisson approx:")
                # print("lambda:", lambda_D, lambda_S, lambda_I)
                # print("k:",k, x_S, x_D, x_I,p_S, p_D, p_I)
                p_value = poisson.sf(x_S + x_D + x_I - 1, lambda_prob)
                # print("Approx p-val: ", p_value)
                # print("lengths:", len(t), len(C[c_acc]))

            else:
                #exact!
                # p_val = 1 - \sum_{(i,j,l) s.t., i < k_I, j < k_S , l < k_D } P(i < k_I, j < k_S , l < k_D)
                # print("lambda:", lambda_D, lambda_S, lambda_I)
                # print("k:",k, x_S, x_D, x_I,p_S, p_D, p_I)
                if p_S == 0.0:
                    p_S = 2.2250738585072014e-308
                if p_D == 0.0:
                    p_D = 2.2250738585072014e-308
                if p_I == 0.0:
                    p_I = 2.2250738585072014e-308

                p_value = 1
                # print("lols:", multinomial_( [x_S, x_D, x_I , m - x_S - x_D - x_I], [3*p_S, p_D, 4*p_I, 1 - 3*p_S - p_D - 4*p_I ]) )
                for i in range(x_S + 1):
                    for j in range(x_D + 1):
                        for l in range(x_I + 1):
                            p_value -= multinomial_( [i, j, l, m - i - j - l], [p_S, p_D, p_I, 1 - 3*p_S - p_D - 4*p_I ]) 
                p_value += multinomial_( [x_S, x_D, x_I, m - x_S - x_D - x_I], [p_S, p_D, p_I, 1 - 3*p_S - p_D - 4*p_I ])

                # p_between_bp = 1

                # print("EXACT:")
                # print("lambda:", lambda_D, lambda_S, lambda_I)
                # print("k:",k, x_S, x_D, x_I,p_S, p_D, p_I, m)
                # print("p-val: ", p_value)
                # print("lengths:", len(t), len(C[c_acc]))

        ############################################################################################
        ############################################################################################

        p_values[c_acc] = (t_acc, k, p_value, N_t)

        if math.isnan(p_value):
            print("LOOOOL math is nan!:")
            print("lambda:", lambda_D, lambda_S, lambda_I)
            print("k:",k, x_S, x_D, x_I,p_S, p_D, p_I)
            print("Approx p-val: ", p_value)
            print("lengths:", len(t), len(C[c_acc]))
            sys.exit()

    return p_values


def test_against_closest(delta_t_single, alignment_matrix_to_t_single, t_acc_single, t, candidate_accessions_single, partition_of_X, partition_of_C, C):
    p_value_closest = test_against_center(delta_t_single, alignment_matrix_to_t_single, t_acc_single, t, candidate_accessions_single, partition_of_X, partition_of_C, C)
    return p_value_closest


