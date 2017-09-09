
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
from modules.functions import create_position_probability_matrix, get_error_rates_and_lambda, get_difference_coordinates_for_candidates, get_supporting_reads_for_candidates, adjust_probability_of_candidate_to_alignment_invariant


def do_statistical_tests_per_edge(minimizer_graph_transposed, C, X, partition_of_X, params):
    p_values = {}
    actual_tests = 0
    
    # separate partition_of_X and X, C here into subsets for each t in minimizer graph to speed up parallelization when lots of reads
    partition_of_X_per_candidate = {}
    all_X_in_partition = {}
    C_for_minmizer = {}
    candidates_to = {}
    for t_acc in minimizer_graph_transposed:
        candidates_to[t_acc] = {}
        partition_of_X_per_candidate[t_acc] = {}
        C_for_minmizer[t_acc] = {}
        all_X_in_partition[t_acc] = {}
        for c_acc in minimizer_graph_transposed[t_acc]:
            candidates_to[t_acc][c_acc] = {c_acc : minimizer_graph_transposed[t_acc][c_acc] }
            partition_of_X_per_candidate[t_acc][c_acc] = {acc : set([x_acc for x_acc in partition_of_X[acc]]) for acc in [c_acc, t_acc]}
            C_for_minmizer[t_acc][c_acc] = { acc : C[acc] for acc in [c_acc, t_acc] }
            all_X_in_partition[t_acc][c_acc] = { x_acc : X[x_acc] for acc in [t_acc, c_acc] for x_acc in partition_of_X[acc]}


    if params.single_core:
        for t_acc in minimizer_graph_transposed:
            if len(minimizer_graph_transposed[t_acc]) == 0:
                p_values[t_acc] = ("not_tested", "NA", len(partition_of_X[t_acc]), len(partition_of_X[t_acc]), -1 )
                continue

            for c_acc in minimizer_graph_transposed[t_acc]:
                p_vals = statistical_test_CLT(t_acc, all_X_in_partition[t_acc][c_acc], C_for_minmizer[t_acc][c_acc], partition_of_X_per_candidate[t_acc][c_acc], candidates_to[t_acc][c_acc], params.ignore_ends_len)

                for tested_cand_acc, (p_value, mult_factor_inv, k, N_t, delta_size) in p_vals.items():
                    if p_value == "not_tested":
                        assert tested_cand_acc not in p_values # should only be here once
                        p_values[tested_cand_acc] = (p_value, mult_factor_inv, k, N_t, delta_size)

                    elif tested_cand_acc in p_values: # new most insignificant p_value
                        actual_tests += 1
                        if p_value * mult_factor_inv > p_values[tested_cand_acc][0] * p_values[tested_cand_acc][1]:
                            p_values[tested_cand_acc] = (p_value, mult_factor_inv, k, N_t, delta_size)

                    else: # not tested before
                        actual_tests += 1
                        p_values[tested_cand_acc] = (p_value, mult_factor_inv, k, N_t, delta_size)

    else:
        ####### parallelize statistical tests #########
        # pool = Pool(processes=mp.cpu_count())
        original_sigint_handler = signal.signal(signal.SIGINT, signal.SIG_IGN)
        signal.signal(signal.SIGINT, original_sigint_handler)
        pool = Pool(processes=mp.cpu_count())
        try:
            res = pool.map_async(statistical_test_helper, [ ( (t_acc, all_X_in_partition[t_acc][c_acc], C_for_minmizer[t_acc][c_acc], partition_of_X_per_candidate[t_acc][c_acc], candidates_to[t_acc][c_acc], params.ignore_ends_len), {}) for t_acc in minimizer_graph_transposed for c_acc in minimizer_graph_transposed[t_acc]  ] )
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
            for c_acc, (p_value, mult_factor_inv, k, N_t, delta_size) in list(all_tests_to_a_given_target.items()): 
                if p_value ==  "not_tested":
                    assert c_acc not in p_values # should only be here once
                    p_values[c_acc] = (p_value, mult_factor_inv, k, N_t, delta_size)

                elif c_acc in p_values: # new most insignificant p_value
                    actual_tests += 1
                    if p_value * mult_factor_inv > p_values[c_acc][0] * p_values[c_acc][1]:
                        p_values[c_acc] = (p_value, mult_factor_inv, k, N_t, delta_size)
                else: # not tested before
                    actual_tests += 1
                    p_values[c_acc] = (p_value, mult_factor_inv, k, N_t, delta_size)

        for t_acc in minimizer_graph_transposed:
            if t_acc not in p_values:
                p_values[t_acc] = ("not_tested", "NA", len(partition_of_X[t_acc]), len(partition_of_X[t_acc]), -1 )

    print("Total number of tests performed this round:", actual_tests)

    return p_values


def do_statistical_tests_all_c_to_t(minimizer_graph_transposed, C, X, partition_of_X, params):
    p_values = {}
    actual_tests = 0
    
    # separate partition_of_X and X, C here into subsets for each t in minimizer graph to speed up parallelization when lots of reads
    partition_of_X_for_minmizer = {}
    X_for_minmizer = {}
    C_for_minmizer = {}
    candidates_to = {}
    for t_acc in minimizer_graph_transposed:
        candidates_to[t_acc] = minimizer_graph_transposed[t_acc]
        partition_of_X_for_minmizer[t_acc] = {c_acc : set([x_acc for x_acc in partition_of_X[c_acc]]) for c_acc in list(minimizer_graph_transposed[t_acc].keys()) + [t_acc]}
        C_for_minmizer[t_acc] = { c_acc : C[c_acc] for c_acc in list(minimizer_graph_transposed[t_acc].keys()) + [t_acc] }
        X_for_minmizer[t_acc] = { x_acc : X[x_acc] for c_acc in partition_of_X_for_minmizer[t_acc] for x_acc in partition_of_X_for_minmizer[t_acc][c_acc]}

    if params.single_core:
        for t_acc in minimizer_graph_transposed:
            p_vals = statistical_test(t_acc, X_for_minmizer[t_acc], C_for_minmizer[t_acc], partition_of_X_for_minmizer[t_acc], candidates_to[t_acc], params.ignore_ends_len)
            for c_acc, (p_value, mult_factor_inv, k, N_t, delta_size) in p_vals.items():
                if p_value == "not_tested":
                    assert c_acc not in p_values # should only be here once
                    p_values[c_acc] = (p_value, mult_factor_inv, k, N_t, delta_size)

                elif c_acc in p_values: # new most insignificant p_value
                    actual_tests += 1
                    if p_value * mult_factor_inv > p_values[c_acc][0] * p_values[c_acc][1]:
                        p_values[c_acc] = (p_value, mult_factor_inv, k, N_t, delta_size)

                else: # not tested before
                    actual_tests += 1
                    p_values[c_acc] = (p_value, mult_factor_inv, k, N_t, delta_size)

    else:
        ####### parallelize statistical tests #########
        # pool = Pool(processes=mp.cpu_count())
        original_sigint_handler = signal.signal(signal.SIGINT, signal.SIG_IGN)
        signal.signal(signal.SIGINT, original_sigint_handler)
        pool = Pool(processes=mp.cpu_count())
        try:
            res = pool.map_async(statistical_test_helper, [ ( (t_acc, X_for_minmizer[t_acc], C_for_minmizer[t_acc], partition_of_X_for_minmizer[t_acc], candidates_to[t_acc], params.ignore_ends_len), {}) for t_acc in minimizer_graph_transposed  ] )
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
            for c_acc, (p_value, mult_factor_inv, k, N_t, delta_size) in list(all_tests_to_a_given_target.items()): 
                if p_value ==  "not_tested":
                    assert c_acc not in p_values # should only be here once
                    p_values[c_acc] = (p_value, mult_factor_inv, k, N_t, delta_size)

                elif c_acc in p_values: # new most insignificant p_value
                    actual_tests += 1
                    if p_value * mult_factor_inv > p_values[c_acc][0] * p_values[c_acc][1]:
                        p_values[c_acc] = (p_value, mult_factor_inv, k, N_t, delta_size)
                else: # not tested before
                    actual_tests += 1
                    p_values[c_acc] = (p_value, mult_factor_inv, k, N_t, delta_size)

    print("Total number of tests performed this round:", actual_tests)

    return p_values


def statistical_test_helper(arguments):
    args, kwargs = arguments
    return statistical_test_CLT(*args, **kwargs)


def stat_test(k, t_seq, epsilon, delta_t, candidate_indiv_invariant_factors, t_acc, c_acc, original_mapped_to_c):
    m = len(t_seq)
    n_S, n_D, n_I = 0, 0, 0
    for pos, (state, char) in delta_t[c_acc].items():
        if state == "S":
            n_S += 1
        elif state == "D":
            n_D += 1
        if state == "I":
            n_I += 1

    ######### INVARIANTS FACTORS PER VARIANT ##########
    u_v_factor = 1
    epsilon_invariant_adjusted = {}
    for q_acc in epsilon:
        p_i = 1.0
        for pos, (state, char) in delta_t[c_acc].items():
            u_v = candidate_indiv_invariant_factors[c_acc][pos][(state, char)]
            p_iv = epsilon[q_acc][state]
            p_i *= u_v*p_iv
        
        p_i = min(p_i, 1.0) # we cannot have a probability larger than 1, this can happen due to our heuristic degeneracy multiplier "u_v"
        epsilon_invariant_adjusted[q_acc] = p_i
    #################################################

    pobin_mean = sum([ epsilon_invariant_adjusted[q_acc] for q_acc in epsilon_invariant_adjusted])    
    p_value = poisson.sf(k - 1, pobin_mean)
    
    mult_factor_inv = ( (4*(m+1))**n_I ) * functions.choose(m, n_D) * functions.choose( 3*(m-n_D), n_S)

    pobin_var = sum([ epsilon_invariant_adjusted[q_acc] * ( 1 - epsilon_invariant_adjusted[q_acc] ) for q_acc in epsilon_invariant_adjusted])
    pobin_stddev = math.sqrt(pobin_var)
    if pobin_mean > 5.0: # implies that mean is at least 2.2 times larger than stddev, this should give fairly symmetrical distribution
        p_value_norm = norm.sf(float(k-1), loc= pobin_mean, scale= pobin_stddev )

        print(c_acc, "COMPARE P-VALS:", p_value, "norm:", p_value_norm)


    #############################
    #################################### 
    return  p_value, mult_factor_inv

def arrange_alignments(t_acc, reads_and_candidates_and_ref, X, C, ignore_ends_len):
    partition_dict = {t_acc : {}}
    for seq_acc in reads_and_candidates_and_ref:
        if seq_acc in X:
            partition_dict[t_acc][seq_acc] = (C[t_acc], X[seq_acc])
        else:
            partition_dict[t_acc][seq_acc] = (C[t_acc], C[seq_acc])

    exact_edit_distances = edlib_align_sequences_keeping_accession(partition_dict, single_core = True)    
    exact_alignments = sw_align_sequences_keeping_accession(exact_edit_distances, single_core = True, ignore_ends_len = ignore_ends_len)
    partition_alignments = {} 

    assert len(exact_alignments) == 1
    for t_acc in exact_alignments:
        partition_alignments[t_acc] = {}
        for x_acc in exact_alignments[t_acc]:
            aln_t, aln_x, (matches, mismatches, indels) = exact_alignments[t_acc][x_acc]
            edit_dist = mismatches + indels
            partition_alignments[t_acc][x_acc] = (edit_dist, aln_t, aln_x, 1)

    alignment_matrix_to_t, PFM_to_t = create_position_probability_matrix(C[t_acc], partition_alignments[t_acc])
    return alignment_matrix_to_t, PFM_to_t


def statistical_test(t_acc, X, C, partition_of_X, candidates, ignore_ends_len):
    significance_values = {}
    t_seq = C[t_acc]
    if len(candidates) == 0:
        significance_values[t_acc] = ("not_tested", "NA", len(partition_of_X[t_acc]), len(partition_of_X[t_acc]), -1 )
        return significance_values

    reads = set([x_acc for c_acc in candidates for x_acc in partition_of_X[c_acc]] )
    reads.update(partition_of_X[t_acc])

    #### Bug FIX ############
    total_read_in_partition = len(partition_of_X[t_acc])
    reads_in_partition = set(partition_of_X[t_acc])
    # print(partition_of_X[t_acc])
    for c_acc in candidates:
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

    # print("N_t:", N_t, "reads in partition:", total_read_in_partition, "ref:", t_acc )
    # print("Nr candidates:", len(candidates), candidates)
    assert total_read_in_partition == N_t # each read should be uiquely assinged to a candidate
    reads_and_candidates = reads.union( [c_acc for c_acc in candidates]) 
    reads_and_candidates_and_ref = reads_and_candidates.union( [t_acc] ) 

    # get multialignment matrix here
    alignment_matrix_to_t, PFM_to_t =  arrange_alignments(t_acc, reads_and_candidates_and_ref, X, C, ignore_ends_len)

    # cut multialignment matrix first and last ignore_ends_len bases in ends of reference in the amignment matrix
    # these are bases that we disregard when testing varinats
    # We get individual cut positions depending on which candidate is being tested -- we dont want to include ends spanning over the reference or candidate
    # we cut at the start position in c or t that comes last, and the end position in c or t that comes first
    assert len(list(candidates)) == 1
    for c_acc in list(candidates):
        if ignore_ends_len > 0:
            alignment_matrix_to_t = functions.cut_ends_of_alignment_matrix(alignment_matrix_to_t, t_acc, c_acc, ignore_ends_len)

        # get parameter estimates for statistical test
        candidate_accessions = set( [ c_acc for c_acc in candidates] )
        delta_t = functions.get_difference_coordinates_for_candidates(t_acc, candidate_accessions, alignment_matrix_to_t) # format: { c_acc1 : {pos:(state, char), pos2:(state, char) } , c_acc2 : {pos:(state, char), pos2:(state, char) },... }
        epsilon, lambda_S, lambda_D, lambda_I = functions.get_error_rates_and_lambda(t_acc, len(t_seq), candidate_accessions, alignment_matrix_to_t) 
        # get number of reads k supporting the given set of variants, they have to support all the variants within a candidate
        candidate_support = functions.get_supporting_reads_for_candidates(t_acc, candidate_accessions, alignment_matrix_to_t, delta_t, partition_of_X) # format: { c_acc1 : [x_acc1, x_acc2,.....], c_acc2 : [x_acc1, x_acc2,.....] ,... }
        candidate_indiv_invariant_factors = functions.adjust_probability_of_candidate_to_alignment_invariant(delta_t, alignment_matrix_to_t, t_acc)

    for c_acc in list(candidates):
        original_mapped_to_c = len(partition_of_X[c_acc])
        k = len(candidate_support[c_acc])
        # print("supprot:", k, "diff:", len(delta_t[c_acc]))
        if len(delta_t[c_acc]) == 0:
            print("{0} no difference to ref {1} after ignoring ends!".format(c_acc, t_acc))
        p_value, mult_factor_inv = stat_test(k, t_seq, epsilon, delta_t, candidate_indiv_invariant_factors, t_acc, c_acc, original_mapped_to_c)
        delta_size = len(delta_t[c_acc])
        significance_values[c_acc] = (p_value, mult_factor_inv, k, N_t, delta_size)
        print("Tested", c_acc, "to ref", t_acc, "p_val:{0}, mult_factor:{1}, corrected p_val:{2} k:{3}, N_t:{4}, Delta_size:{5}".format(p_value, mult_factor_inv, p_value * mult_factor_inv,  k, N_t, delta_size) )

    return significance_values



def statistical_test_CLT(t_acc, X, C, partition_of_X, candidates, ignore_ends_len):
    significance_values = {}
    t_seq = C[t_acc]
    if len(candidates) == 0:
        significance_values[t_acc] = ("not_tested", "NA", len(partition_of_X[t_acc]), len(partition_of_X[t_acc]), -1 )
        return significance_values

    reads = set([x_acc for c_acc in candidates for x_acc in partition_of_X[c_acc]] )
    reads.update(partition_of_X[t_acc])

    #### Bug FIX ############
    total_read_in_partition = len(partition_of_X[t_acc])
    reads_in_partition = set(partition_of_X[t_acc])
    # print(partition_of_X[t_acc])
    for c_acc in candidates:
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

    # print("N_t:", N_t, "reads in partition:", total_read_in_partition, "ref:", t_acc )
    # print("Nr candidates:", len(candidates), candidates)
    assert total_read_in_partition == N_t # each read should be uiquely assinged to a candidate
    reads_and_candidates = reads.union( [c_acc for c_acc in candidates]) 
    reads_and_candidates_and_ref = reads_and_candidates.union( [t_acc] ) 

    # get multialignment matrix here
    alignment_matrix_to_t, PFM_to_t =  arrange_alignments(t_acc, reads_and_candidates_and_ref, X, C, ignore_ends_len)

    # cut multialignment matrix first and last ignore_ends_len bases in ends of reference in the amignment matrix
    # these are bases that we disregard when testing varinats
    # We get individual cut positions depending on which candidate is being tested -- we dont want to include ends spanning over the reference or candidate
    # we cut at the start position in c or t that comes last, and the end position in c or t that comes first
    assert len(list(candidates)) == 1
    for c_acc in list(candidates):
        if ignore_ends_len > 0:
            alignment_matrix_to_t = functions.cut_ends_of_alignment_matrix(alignment_matrix_to_t, t_acc, c_acc, ignore_ends_len)

        # get parameter estimates for statistical test
        candidate_accessions = set( [ c_acc for c_acc in candidates] )
        delta_t = functions.get_difference_coordinates_for_candidates(t_acc, candidate_accessions, alignment_matrix_to_t) # format: { c_acc1 : {pos:(state, char), pos2:(state, char) } , c_acc2 : {pos:(state, char), pos2:(state, char) },... }
        errors = functions.get_errors_per_read(t_acc, len(t_seq), candidate_accessions, alignment_matrix_to_t) 
        weight = functions.get_weights_per_read(t_acc, len(t_seq), candidate_accessions, errors) 
        invariant_factors_for_candidate = functions.adjust_probability_of_candidate_to_alignment_invariant(delta_t, alignment_matrix_to_t, t_acc)
        probability = functions.get_prob_of_support_per_read(t_acc, len(t_seq), candidate_accessions, errors, invariant_factors_for_candidate) 

        # get number of reads k supporting the given set of variants, they have to support all the variants within a candidate
        x = functions.reads_supporting_candidate(t_acc, candidate_accessions, alignment_matrix_to_t, delta_t, partition_of_X) # format: { c_acc1 : [x_acc1, x_acc2,.....], c_acc2 : [x_acc1, x_acc2,.....] ,... }
        # p_value = CLT_test(probability, weight, x)
        p_value = poisson_approx_test(probability, weight, x)
        correction_factor = calc_correction_factor(t_seq, c_acc, delta_t)

        delta_size = len(delta_t[c_acc])
        if delta_size == 0:
            print("{0} no difference to ref {1} after ignoring ends!".format(c_acc, t_acc))

        print("Tested", c_acc, "to ref", t_acc, "p_val:{0}, mult_factor:{1}, corrected p_val:{2} k:{3}, N_t:{4}, Delta_size:{5}".format(p_value, correction_factor, p_value * correction_factor,  len(x), N_t, delta_size) )
        significance_values[c_acc] = (p_value, correction_factor, len(x), N_t, delta_size)

    return significance_values


def poisson_approx_test(probability, weight, x):

    print([weight[x_i] for x_i in probability if x_i in x ])
    print([weight[x_i] for x_i in probability ])
    min_w = min(weight.values())
    weight_multiplier = 1.0 / min_w

    weight = { x_i : weight[x_i] * weight_multiplier for x_i in weight}
    print([weight[x_i] for x_i in probability ])

    observed_weighted_x = sum([weight[x_i]*1.0 for x_i in probability if x_i in x ])
    po_lambda = sum([ probability[x_i]* weight[x_i] for x_i in probability ])

    k = observed_weighted_x if observed_weighted_x == 0 or not observed_weighted_x.is_integer() else observed_weighted_x - 1
    print("k:", k)
    p_value = poisson.sf(k - 1, po_lambda)

    return p_value


def CLT_test(probability, weight, x):

    observed_weighted_x = sum([weight[x_i]*1.0 for x_i in probability if x_i in x ])
    mu = sum([ probability[x_i]* weight[x_i] for x_i in probability ])
    var = sum([ probability[x_i]*(1.0 - probability[x_i])* weight[x_i]**2 for x_i in probability ])
    sigma = math.sqrt(var)
    print([weight[x_i] for x_i in probability if x_i in x ])
    print([weight[x_i] for x_i in probability ])
    print(mu, sigma, observed_weighted_x)
    p_value_norm = norm.sf(observed_weighted_x, loc= mu, scale= sigma )

    return p_value_norm

def calc_correction_factor(t_seq, c_acc, delta_t):
    m = len(t_seq)
    n_S, n_D, n_I = 0, 0, 0
    for pos, (state, char) in delta_t[c_acc].items():
        if state == "S":
            n_S += 1
        elif state == "D":
            n_D += 1
        if state == "I":
            n_I += 1
    
    correction_factor = ( (4*(m+1))**n_I ) * functions.choose(m, n_D) * functions.choose( 3*(m-n_D), n_S)
    return correction_factor
