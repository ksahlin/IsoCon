
import signal
from multiprocessing import Pool
import multiprocessing as mp
import sys

import math
from scipy.stats import poisson, binom, norm

from modules import functions
from modules.multinomial_distr import multinomial_
from modules.SW_alignment_module import sw_align_sequences_keeping_accession
from modules.edlib_alignment_module import edlib_align_sequences_keeping_accession
from modules.functions import create_position_probability_matrix, get_error_rates_and_lambda, get_difference_coordinates_for_candidates, get_supporting_reads_for_candidates, get_invariant_adjustment, adjust_probability_of_read_to_alignment_invariant, adjust_probability_of_candidate_to_alignment_invariant


def do_statistical_tests(minimizer_graph, C, X, partition_of_X, single_core = False):
    p_values = {}
    actual_tests = 0

    if single_core:
        for t_acc in minimizer_graph:
            p_vals = statistical_test(t_acc, X, C, partition_of_X, minimizer_graph)
            for c_acc, (corrected_p_value, k, N_t) in p_vals.items():
                if corrected_p_value ==  "not_tested":
                    assert c_acc not in p_values # should only be here once
                    p_values[c_acc] = (corrected_p_value, k, N_t)

                elif c_acc in p_values: # new most insignificant p_value
                    actual_tests += 1
                    if corrected_p_value > p_values[c_acc][0]:
                        p_values[c_acc] = (corrected_p_value, k, N_t)

                else: # not tested before
                    actual_tests += 1
                    p_values[c_acc] = (corrected_p_value, k, N_t)

    else:
        ####### parallelize statistical tests #########
        # pool = Pool(processes=mp.cpu_count())
        original_sigint_handler = signal.signal(signal.SIGINT, signal.SIG_IGN)
        signal.signal(signal.SIGINT, original_sigint_handler)
        pool = Pool(processes=mp.cpu_count())
        try:
            res = pool.map_async(statistical_test_helper, [ ( (t_acc, X, C, partition_of_X, minimizer_graph), {}) for t_acc in minimizer_graph  ] )
            statistical_test_results =res.get(999999999) # Without the timeout this blocking call ignores all signals.
        except KeyboardInterrupt:
            print("Caught KeyboardInterrupt, terminating workers")
            pool.terminate()
            sys.exit()
        else:
            print("Normal termination")
            pool.close()
        pool.join()

        for all_tests_to_a_given_target in statistical_test_results:
            for c_acc, (corrected_p_value, k, N_t) in list(all_tests_to_a_given_target.items()): 
                if corrected_p_value ==  "not_tested":
                    assert c_acc not in p_values # should only be here once
                    p_values[c_acc] = (corrected_p_value, k, N_t)

                elif c_acc in p_values: # new most insignificant p_value
                    actual_tests += 1
                    if corrected_p_value > p_values[c_acc][0]:
                        p_values[c_acc] = (corrected_p_value, k, N_t)
                else: # not tested before
                    actual_tests += 1
                    p_values[c_acc] = (corrected_p_value, k, N_t)

    print("Total number of tests performed this round:", actual_tests)

    return p_values


def statistical_test_helper(arguments):
    args, kwargs = arguments
    return statistical_test(*args, **kwargs)


def stat_test(k, m, epsilon, delta_t, candidate_indiv_invariant_factors, t_acc, c_acc, original_mapped_to_c):
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
            u_v = candidate_indiv_invariant_factors[c_acc][pos][(state, char)]
            p_iv = epsilon[q_acc][state]
            p_i *= u_v*p_iv
        epsilon_invariant_adjusted[q_acc] = p_i

    lambda_po_approx_inv = sum([ epsilon_invariant_adjusted[q_acc] for q_acc in epsilon_invariant_adjusted])
    mult_factor_inv = ( (4*(m+1))**n_I ) * functions.choose(m, n_D) * functions.choose( 3*(m-n_D), n_S)
    p_value = poisson.sf(k - 1, lambda_po_approx_inv)
    corrected_p_value = mult_factor_inv * p_value
    print("Corrected p value:", corrected_p_value, "k:", k, "Nr mapped to:", original_mapped_to_c)
    print("lambda inv adjusted", lambda_po_approx_inv, mult_factor_inv, k, len(delta_t[c_acc]), candidate_indiv_invariant_factors[c_acc])
    #############################
    #################################### 
    return  corrected_p_value

def arrange_alignments(t_acc, reads_and_candidates_and_ref, X, C ):
    partition_dict = {t_acc : {}}
    for seq_acc in reads_and_candidates_and_ref:
        if seq_acc in X:
            partition_dict[t_acc][seq_acc] = (C[t_acc], X[seq_acc])
        else:
            partition_dict[t_acc][seq_acc] = (C[t_acc], C[seq_acc])

    exact_edit_distances = edlib_align_sequences_keeping_accession(partition_dict, single_core = True)    
    exact_alignments = sw_align_sequences_keeping_accession(exact_edit_distances, single_core = True)
    partition_alignments = {} 

    for t_acc in exact_alignments:
        partition_alignments[t_acc] = {}
        for x_acc in exact_alignments[t_acc]:
            aln_t, aln_x, (matches, mismatches, indels) = exact_alignments[t_acc][x_acc]
            edit_dist = mismatches + indels
            partition_alignments[t_acc][x_acc] = (edit_dist, aln_t, aln_x, 1)
    print("NR candidates with at least one hit:", len(partition_alignments))
    alignment_matrix_to_t, PFM_to_t = create_position_probability_matrix(C[t_acc], partition_alignments[t_acc])
    return alignment_matrix_to_t, PFM_to_t



def statistical_test(t_acc, X, C, partition_of_X, minimizer_graph):
    significance_values = {}
    t_seq = C[t_acc]
    if len(minimizer_graph[t_acc]) == 0:
        significance_values[t_acc] = ("not_tested", len(partition_of_X[t_acc]), len(partition_of_X[t_acc]) )
        return significance_values

    reads = set([x_acc for c_acc in minimizer_graph[t_acc] for x_acc in partition_of_X[c_acc]] )
    reads.update(partition_of_X[t_acc])

    #### Bug FIX ############
    total_read_in_partition = len(partition_of_X[t_acc])
    reads_in_partition = set(partition_of_X[t_acc])
    # print(partition_of_X[t_acc])
    for c_acc in minimizer_graph[t_acc]:
        # for x_acc in partition_of_X[c_acc]:
        #     if x_acc in reads_in_partition:
        #         print("Already in partition:", x_acc)
        #         print(x_acc in partition_of_X[t_acc])
        #         print("candidate:", c_acc)
        #     else:
        #         reads_in_partition.add(x_acc)
        total_read_in_partition += len(partition_of_X[c_acc])
        # print(partition_of_X[c_acc])

    ############################

    N_t = len(reads)

    print("N_t:", N_t, "reads in partition:", total_read_in_partition, "ref:", t_acc )
    print("Nr candidates:", len(minimizer_graph[t_acc]), minimizer_graph[t_acc])
    assert total_read_in_partition == N_t # each read should be uiquely assinged to a candidate
    reads_and_candidates = reads.union( [c_acc for c_acc in minimizer_graph[t_acc]]) 
    reads_and_candidates_and_ref = reads_and_candidates.union( [t_acc] ) 

    # get multialignment matrix here
    alignment_matrix_to_t, PFM_to_t =  arrange_alignments(t_acc, reads_and_candidates_and_ref, X, C )

    # get parameter estimates for statistical test
    candidate_accessions = set( [ c_acc for c_acc in minimizer_graph[t_acc]] )
    delta_t = functions.get_difference_coordinates_for_candidates(t_acc, candidate_accessions, alignment_matrix_to_t) # format: { c_acc1 : {pos:(state, char), pos2:(state, char) } , c_acc2 : {pos:(state, char), pos2:(state, char) },... }
    epsilon, lambda_S, lambda_D, lambda_I = functions.get_error_rates_and_lambda(t_acc, len(t_seq), candidate_accessions, alignment_matrix_to_t) 
    # get number of reads k supporting the given set of variants, they have to support all the variants within a candidate
    candidate_support = functions.get_supporting_reads_for_candidates(t_acc, candidate_accessions, alignment_matrix_to_t, delta_t, partition_of_X) # format: { c_acc1 : [x_acc1, x_acc2,.....], c_acc2 : [x_acc1, x_acc2,.....] ,... }
    # read_indiv_invariant_factors = adjust_probability_of_read_to_alignment_invariant(delta_t, alignment_matrix_to_t, t_acc)
    candidate_indiv_invariant_factors = functions.adjust_probability_of_candidate_to_alignment_invariant(delta_t, alignment_matrix_to_t, t_acc)

    for c_acc, c_seq in list(minimizer_graph[t_acc].items()):
        original_mapped_to_c = len(partition_of_X[c_acc])
        k = len(candidate_support[c_acc])
        # print("supprot:", k, "diff:", len(delta_t[c_acc]))
        corrected_p_value = stat_test(k, len(t_seq), epsilon, delta_t, candidate_indiv_invariant_factors, t_acc, c_acc, original_mapped_to_c)
        significance_values[c_acc] = (corrected_p_value, k, N_t)

        # actual_tests += 1
        # if c_acc in significance_values:
        #     if corrected_p_value > significance_values[c_acc][0]:
        #         significance_values[c_acc] = (corrected_p_value, k, N_t)
        # else:
        #     significance_values[c_acc] = (corrected_p_value, k, N_t)

    return significance_values

