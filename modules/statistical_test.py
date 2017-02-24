import signal
from multiprocessing import Pool
import multiprocessing as mp
import sys

import math
from scipy.stats import poisson, binom, norm

from modules.multinomial_distr import multinomial_
from modules.SW_alignment_module import sw_align_sequences_keeping_accession
from modules.edlib_alignment_module import edlib_align_sequences_keeping_accession
from modules.functions import create_position_probability_matrix, get_error_rates_and_lambda, get_difference_coordinates_for_candidates, get_supporting_reads_for_candidates, get_invariant_adjustment, adjust_probability_of_read_to_alignment_invariant


def do_statistical_tests(null_hypothesis_references_t, C_seq_to_acc, partition_of_X, partition_of_C, X, C, previous_round_tests, previous_candidate_p_values, single_core = False):

    p_values = {}

    ###############################################
    ###############################################
    # if the reference support and all candidate supports are identical to previous round, skip test
    # store the performed tests to avoid identical testing in next iteration.
    current_round_tests = {}
    for t in null_hypothesis_references_t:
        t_acc = C_seq_to_acc[t]
        N_t = len(partition_of_X[t_acc])
        for c in partition_of_C[t]:
            c_acc = C_seq_to_acc[c]
            N_t += len(partition_of_X[c_acc])
        current_round_tests[t] = (set([c for c in partition_of_C[t]]), N_t)
        print("CURRENT TEST: candidates {0}, N_t: {1}".format(len(current_round_tests[t][0]), current_round_tests[t][1]))

    for t in current_round_tests:
        if t in previous_round_tests:
            if current_round_tests[t][0] == previous_round_tests[t][0] and current_round_tests[t][1] == previous_round_tests[t][1]:
                print("SKIPPING TESTS BECAUSE IDENTICAL TO PREVIOUS:", len(previous_round_tests[t][0]), previous_round_tests[t][1]) 
                
                for c in partition_of_C[t]:
                    c_acc = C_seq_to_acc[c]
                    p_values[c_acc] =  previous_candidate_p_values[c_acc] 
                
                null_hypothesis_references_t.remove(t)
    ###############################################
    ###############################################

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

    p_values_on_closest_highest_support = test_against_center(delta_t, alignment_matrix_to_t, t_acc, t, candidate_accessions, partition_of_X, partition_of_C, C)

    return p_values_on_closest_highest_support


def choose(n, k):
    """
    A fast way to calculate binomial coefficients by Andrew Dalke (contrib).
    """
    if 0 <= k <= n:
        ntok = 1
        ktok = 1
        for t in range(1, min(k, n - k) + 1):
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

    exact_edit_distances = edlib_align_sequences_keeping_accession(partition_dict, single_core = True)    
    exact_alignments = sw_align_sequences_keeping_accession(exact_edit_distances, single_core = True)
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

    read_indiv_invariant_factors = adjust_probability_of_read_to_alignment_invariant(delta_t, alignment_matrix_to_t, t_acc)

    # print()
    # print("INVARIANT FACTORS:", invariant_factors)
    # print()
    # get p_value
    m = len(t)

    for c_acc in candidate_accessions:
        u_c_S, u_c_D, u_c_I  = invariant_factors[c_acc]
        # print("UNIQUENESS:", u_c_S, u_c_D, u_c_I)
        k = len(candidate_support[c_acc]) #/ float(u_c)
        N_t = len(alignment_matrix_to_t) - len(candidate_accessions) - 1 # all reads minus all candidates and the reference transcript
        # print("reads N_t:", N_t)
        # print("varinats:",delta_t[c_acc].items())

        ############################################################################################
        ############################################################################################

        ############ NEW TEST ##############

        n_S, n_D, n_I = 0, 0, 0
        for pos, (state, char) in delta_t[c_acc].items():
            if state == "S":
                n_S += 1
            elif state == "D":
                n_D += 1
            if state == "I":
                n_I += 1

        ######### INVARIANTS ##########
        u_v_factor = 1
        epsilon_invariant_adjusted = {}
        for q_acc in epsilon:
            p_i = 1
            for pos, (state, char) in delta_t[c_acc].items():
                u_v = read_indiv_invariant_factors[q_acc][pos][(state, char)]
                p_iv = epsilon[q_acc][state]
                if u_v != 1:
                    p_i *= binom.pmf(1, u_v, p_iv)
                    u_v_factor = u_v
                else:
                    p_i *= p_iv

            epsilon_invariant_adjusted[q_acc] = p_i
            # if u_v_factor == 1 and len(delta_t[c_acc]) == 1:
            #     if epsilon[q_acc][state] != p_i:
            #         print("OMG!!!", epsilon[q_acc][state], p_i, delta_t[c_acc])

        lambda_po_approx_inv = sum([ epsilon_invariant_adjusted[q_acc] for q_acc in epsilon_invariant_adjusted])
        mult_factor_inv = (m**n_I)* choose(m, n_D) * choose(m-n_D, n_S)
        p_val_appprox_inv = mult_factor_inv * poisson.sf(k - 1, lambda_po_approx_inv)
        print("lambda inv adjusted", lambda_po_approx_inv, mult_factor_inv, k, len(delta_t[c_acc]), u_v_factor)
        #############################
        
        lambda_po_approx = sum([ ((epsilon[q_acc]["I"])**n_I) * ((epsilon[q_acc]["S"])**n_S) * ((epsilon[q_acc]["D"])**n_D)  for q_acc in epsilon])

        mult_factor = (m**n_I)* choose(m, n_D) * choose(m-n_D, n_S)
        p_val_appprox = (mult_factor/ (u_c_S * u_c_D * u_c_I) ) * poisson.sf(k - 1, lambda_po_approx)
        print("lambda", lambda_po_approx, mult_factor, k, (u_c_S * u_c_D * u_c_I), n_I, n_D, n_S)
        # print("FFFFFFFF", len(epsilon),len(epsilon_invariant_adjusted))
        #################################### 


        if k <= 1:
            # print("NO support!")
            # print("lengths:", len(t), len(C[c_acc]))                    
            p_value = 1
        else:
            k_I, k_S, k_D = k, k, k
            p_I = poisson.sf(k_I - 1, u_c_I*lambda_I)  # SF defined as 1-P(S_I > k), we want 1-P(S_I >= k)
            p_S = poisson.sf(k_S - 1, u_c_S*lambda_S)
            p_D = poisson.sf(k_D - 1, u_c_D*lambda_D)

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
        print("DIFFERENCE:", p_value, p_val_appprox, p_val_appprox_inv, k)
        # p_values[c_acc] = (t_acc, k, p_value, N_t)
        p_values[c_acc] = (t_acc, k, p_val_appprox_inv, N_t)

        if math.isnan(p_value):
            print("LOOOOL math is nan!:")
            print("lambda:", lambda_D, lambda_S, lambda_I)
            print("k:",k, x_S, x_D, x_I,p_S, p_D, p_I)
            print("Approx p-val: ", p_value)
            print("lengths:", len(t), len(C[c_acc]))
            sys.exit()

        # if m == 1881:
        #     sys.exit()
    return p_values


def test_against_closest(delta_t_single, alignment_matrix_to_t_single, t_acc_single, t, candidate_accessions_single, partition_of_X, partition_of_C, C):
    p_value_closest = test_against_center(delta_t_single, alignment_matrix_to_t_single, t_acc_single, t, candidate_accessions_single, partition_of_X, partition_of_C, C)
    return p_value_closest


