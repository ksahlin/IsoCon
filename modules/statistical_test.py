import signal
from multiprocessing import Pool
import multiprocessing as mp
import sys

import math
from scipy.stats import poisson, binom, norm

from modules.multinomial_distr import multinomial_
from SW_alignment_module import sw_align_sequences_keeping_accession
from functions import create_position_probability_matrix, get_error_rates_and_lambda, get_difference_coordinates_for_candidates, get_supporting_reads_for_candidates


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
    p_values = {}
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
    
    p_values[t_acc] = (t_acc, len(reads_to_t), -1, len(reads_to_map)) 
    print(len(reads_to_map))

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
    forbidden_positions = {} # read_acc: set(positions)
    x_to_c_acc = {}
    for c_acc in candidate_accessions:
        for x_acc in partition_of_X[c_acc]:
            x_to_c_acc[x_acc] = c_acc

    for x_acc in partition_of_X[t_acc]:
        x_to_c_acc[x_acc] = t_acc

    for c_acc in delta_t:
        forbidden_estimation_pos_in_c = delta_t[c_acc].keys() 
        for x_acc in partition_of_X[c_acc]:
            forbidden_positions[x_acc] = set(forbidden_estimation_pos_in_c)
    for x_acc in partition_of_X[t_acc]:
        forbidden_positions[x_acc] = set([])
    forbidden_positions[t_acc] = set([])
    epsilon, lambda_S, lambda_D, lambda_I = get_error_rates_and_lambda(t_acc, len(t), candidate_accessions, alignment_matrix_to_t, forbidden_positions, x_to_c_acc) 
    
  
    # get number of reads k supporting the given set of variants, they have to support all the variants within a candidate
    candidate_support = get_supporting_reads_for_candidates(t_acc, candidate_accessions, alignment_matrix_to_t, delta_t, partition_of_X) # format: { c_acc1 : [x_acc1, x_acc2,.....], c_acc2 : [x_acc1, x_acc2,.....] ,... }


    # get p_value
    m = len(t)
    for c_acc in candidate_accessions:
        k = len(candidate_support[c_acc])
        N_t = len(alignment_matrix_to_t) -  len(partition_of_C[t]) - 1 # all reads minus all candidates and the reference transcript
        # print("reads N_t:", N_t)
        # print("varinats:",delta_t[c_acc].items())

        ############################################################################################
        ############################################################################################
        # if c_acc == "read_81_support_3":
        #     print("HEEEERE!!!")
        #     # print("P-VALUE:", p_value )
        #     print("reads N_t:", N_t)
        #     print("lambda:", lambda_D, lambda_S, lambda_I)
        #     print("k:",k)
        #     print("k:",k, x_S, x_D, x_I,p_S, p_D, p_I)
        #     sys.exit()

        if k <= 1:
            print("NO support!")
            print("lengths:", len(t), len(C[c_acc]))                    
            p_value = 1
        else:
            k_I, k_S, k_D = k, k, k
            p_I = poisson.sf(k_I - 1, lambda_I)  # SF defined as 1-P(S_I > k), we want 1-P(S_I >= k)
            p_S = poisson.sf(k_S - 1, lambda_S)
            p_D = poisson.sf(k_D - 1, lambda_D)

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
                print("here approx")
                p_value = 0.0
            elif (p_I + p_D + p_S)*m >= 10 :
                # approximate with normal
                p_bin = min(0.99, (p_I + p_D + p_S)) # cannot be larger than one, but may be due to approximation errors when lambda is huge
                mu = p_bin*m
                sigma = math.sqrt( (1 - p_bin)* p_bin*m )
                p_value = norm.sf(x_S + x_D + x_I , loc=mu , scale=sigma)
                print("LOOOOL NORMAL approx:")
                print("lambda:", lambda_D, lambda_S, lambda_I)
                print("k:",k, x_S, x_D, x_I,p_S, p_D, p_I)
                print("Approx p-val: ", p_value)
                print("lengths:", len(t), len(C[c_acc]))
            elif x_S + x_D + x_I > 10 and (p_I + p_D + p_S)*m < 10 :
                # approximate with poisson
                lambda_prob = p_I + p_D + p_S
                print("LOOOOL poisson approx:")
                print("lambda:", lambda_D, lambda_S, lambda_I)
                print("k:",k, x_S, x_D, x_I,p_S, p_D, p_I)
                p_value = poisson.sf(x_S + x_D + x_I - 1, lambda_prob)
                print("Approx p-val: ", p_value)
                print("lengths:", len(t), len(C[c_acc]))

            else:
                #exact!
                # p_val = 1 - \sum_{(i,j,l) s.t., i < k_I, j < k_S , l < k_D } P(i < k_I, j < k_S , l < k_D)
                # print("lambda:", lambda_D, lambda_S, lambda_I)
                # print("k:",k, x_S, x_D, x_I,p_S, p_D, p_I)
                p_value = 1

                # print("lols:", multinomial_( [x_S, x_D, x_I , m - x_S - x_D - x_I], [3*p_S, p_D, 4*p_I, 1 - 3*p_S - p_D - 4*p_I ]) )
                for i in range(x_S + 1):
                    for j in range(x_D + 1):
                        for l in range(x_I + 1):
                            p_value -= multinomial_( [i, j, l, m - i - j - l], [p_S, p_D, p_I, 1 - 3*p_S - p_D - 4*p_I ]) 
                p_value += multinomial_( [x_S, x_D, x_I, m - x_S - x_D - x_I], [p_S, p_D, p_I, 1 - 3*p_S - p_D - 4*p_I ])

                # p_between_bp = 1

                # print("lols:", multinomial_( [x_S, x_D, x_I , m - x_S - x_D - x_I], [3*p_S, p_D, 4*p_I, 1 - 3*p_S - p_D - 4*p_I ]) )
                # for l in range(x_I + 1):
                #     p_between_bp -= multinomial_( [l, 0, 0, 0, m+1 - l], [p_I, p_I,p_I,p_I, 1 - 4*p_I ]) 
                # p_between_bp += multinomial_( [x_I, 0, 0, 0, m+1 - x_I], [p_I, p_I,p_I,p_I, 1 - 4*p_I ])
                # p_value = p_between_bp*p_on_bp
                # print("P-VALUE:", p_value )

                # p_value_bin = binom.sf(k_D - 1, m , p_D)
                # print("P-VALUE bin:", p_value_bin )

        ############################################################################################
        ############################################################################################

        p_values[c_acc] = (t_acc, k, p_value, N_t)
        # lambda_poisson = 0
        # for x_acc in epsilon:
        #     x_probabilities = epsilon[x_acc]
        #     p_i = 1
        #     for pos, (state, char) in delta_t[c_acc].items():
        #         p_i *= x_probabilities[state]

        #     lambda_poisson += p_i
        # if len(delta_t[c_acc]) < 11:
        #     m_choose_delta = choose(m, len(delta_t[c_acc])) # number of position combinations where delta could happen
        #     # Nt_choose_k = choose(N_t, k) # possible combinations of k reads where the matches between N_t reads could occur
        #     # print("Stats parameters:", lambda_poisson, lambda_poisson * m_choose_delta, k, N_t)
        #     # Use Bin(n,p) with Bin(len(t), lambda_poisson)??
        #     prob_delta = poisson.sf(k-1, lambda_poisson)
            


        #     multi_corrected_lambda = prob_delta* m_choose_delta   #*Nt_choose_k
        #     # print(prob_delta, m_choose_delta, multi_corrected_lambda)
        #     # if multi_corrected_lambda == 0:
        #     p_value = multi_corrected_lambda
        #     # elif multi_corrected_lambda <= 10:
        #     #     p_value = poisson.sf(0, multi_corrected_lambda)
        #     # else:
        #     #     print("here")
        #     #     p_value = 1 # binom.sf(0, m_choose_delta, prob_delta)
        #     print('length of A:', m, "k:", k, "delta size:", len(delta_t[c_acc]),  "nr support pos:", m_choose_delta, "lambda:", lambda_poisson, "prob_delta:", prob_delta, "p val:", p_value)
        # else:
        #     p_value = 0
        if math.isnan(p_value):
            print("LOOOOL math is nan!:")
            print("lambda:", lambda_D, lambda_S, lambda_I)
            print("k:",k, x_S, x_D, x_I,p_S, p_D, p_I)
            print("Approx p-val: ", p_value)
            print("lengths:", len(t), len(C[c_acc]))
            sys.exit()

    return p_values


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

    # for c in partition_alignments:
    #     if len(partition_alignments[c]) < 1: 
    #         del partition_alignments[c]
    print("NR candidates with at least one hit:", len(partition_alignments))

    return partition_alignments