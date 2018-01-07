
#####################################################################################################################  
#####################################################################################################################  
#################################   hypothesis_test_module.py  ######################################################  
#####################################################################################################################  
#####################################################################################################################  
##################################################################################################################### 


# def do_statistical_tests_all_c_to_t(nearest_neighbor_graph_transposed, C, X, partition_of_X, params):
#     p_values = {}
#     actual_tests = 0
    
#     # separate partition_of_X and X, C here into subsets for each t in nearest_neighbor graph to speed up parallelization when lots of reads
#     partition_of_X_for_minmizer = {}
#     X_for_minmizer = {}
#     C_for_minmizer = {}
#     candidates_to = {}
#     for t_acc in nearest_neighbor_graph_transposed:
#         candidates_to[t_acc] = nearest_neighbor_graph_transposed[t_acc]
#         partition_of_X_for_minmizer[t_acc] = {c_acc : set([x_acc for x_acc in partition_of_X[c_acc]]) for c_acc in list(nearest_neighbor_graph_transposed[t_acc].keys()) + [t_acc]}
#         C_for_minmizer[t_acc] = { c_acc : C[c_acc] for c_acc in list(nearest_neighbor_graph_transposed[t_acc].keys()) + [t_acc] }
#         X_for_minmizer[t_acc] = { x_acc : X[x_acc] for c_acc in partition_of_X_for_minmizer[t_acc] for x_acc in partition_of_X_for_minmizer[t_acc][c_acc]}

#     if params.single_core:
#         for t_acc in nearest_neighbor_graph_transposed:
#             p_vals = statistical_test(t_acc, X_for_minmizer[t_acc], C_for_minmizer[t_acc], partition_of_X_for_minmizer[t_acc], candidates_to[t_acc], params.ignore_ends_len)
#             for c_acc, (p_value, mult_factor_inv, k, N_t, delta_size) in p_vals.items():
#                 if p_value == "not_tested":
#                     assert c_acc not in p_values # should only be here once
#                     p_values[c_acc] = (p_value, mult_factor_inv, k, N_t, delta_size)

#                 elif c_acc in p_values: # new most insignificant p_value
#                     actual_tests += 1
#                     if p_value * mult_factor_inv > p_values[c_acc][0] * p_values[c_acc][1]:
#                         p_values[c_acc] = (p_value, mult_factor_inv, k, N_t, delta_size)

#                 else: # not tested before
#                     actual_tests += 1
#                     p_values[c_acc] = (p_value, mult_factor_inv, k, N_t, delta_size)

#     else:
#         ####### parallelize statistical tests #########
#         # pool = Pool(processes=mp.cpu_count())
#         original_sigint_handler = signal.signal(signal.SIGINT, signal.SIG_IGN)
#         signal.signal(signal.SIGINT, original_sigint_handler)
#         pool = Pool(processes=mp.cpu_count())
#         try:
#             res = pool.map_async(statistical_test_helper, [ ( (t_acc, X_for_minmizer[t_acc], C_for_minmizer[t_acc], partition_of_X_for_minmizer[t_acc], candidates_to[t_acc], params.ignore_ends_len), {}) for t_acc in nearest_neighbor_graph_transposed  ] )
#             statistical_test_results =res.get(999999999) # Without the timeout this blocking call ignores all signals.
#         except KeyboardInterrupt:
#             print("Caught KeyboardInterrupt, terminating workers")
#             pool.terminate()
#             sys.exit()
#         else:
#             print("Normal termination")
#             pool.close()
#         pool.join()

#         for all_tests_to_a_given_target in statistical_test_results:
#             for c_acc, (p_value, mult_factor_inv, k, N_t, delta_size) in list(all_tests_to_a_given_target.items()): 
#                 if p_value ==  "not_tested":
#                     assert c_acc not in p_values # should only be here once
#                     p_values[c_acc] = (p_value, mult_factor_inv, k, N_t, delta_size)

#                 elif c_acc in p_values: # new most insignificant p_value
#                     actual_tests += 1
#                     if p_value * mult_factor_inv > p_values[c_acc][0] * p_values[c_acc][1]:
#                         p_values[c_acc] = (p_value, mult_factor_inv, k, N_t, delta_size)
#                 else: # not tested before
#                     actual_tests += 1
#                     p_values[c_acc] = (p_value, mult_factor_inv, k, N_t, delta_size)

#     print("Total number of tests performed this round:", actual_tests)

#     return p_values



# def stat_test(k, t_seq, epsilon, delta_t, candidate_indiv_invariant_factors, t_acc, c_acc, original_mapped_to_c):
#     m = len(t_seq)
#     n_S, n_D, n_I = 0, 0, 0
#     for pos, (state, char) in delta_t[c_acc].items():
#         if state == "S":
#             n_S += 1
#         elif state == "D":
#             n_D += 1
#         if state == "I":
#             n_I += 1

#     ######### INVARIANTS FACTORS PER VARIANT ##########
#     u_v_factor = 1
#     epsilon_invariant_adjusted = {}
#     for q_acc in epsilon:
#         p_i = 1.0
#         for pos, (state, char) in delta_t[c_acc].items():
#             u_v = candidate_indiv_invariant_factors[c_acc][pos][(state, char)]
#             p_iv = epsilon[q_acc][state]
#             p_i *= u_v*p_iv
        
#         p_i = min(p_i, 1.0) # we cannot have a probability larger than 1, this can happen due to our heuristic degeneracy multiplier "u_v"
#         epsilon_invariant_adjusted[q_acc] = p_i
#     #################################################

#     pobin_mean = sum([ epsilon_invariant_adjusted[q_acc] for q_acc in epsilon_invariant_adjusted])    
#     p_value = poisson.sf(k - 1, pobin_mean)
    
#     mult_factor_inv = ( (4*(m+1))**n_I ) * functions.choose(m, n_D) * functions.choose( 3*(m-n_D), n_S)

#     pobin_var = sum([ epsilon_invariant_adjusted[q_acc] * ( 1 - epsilon_invariant_adjusted[q_acc] ) for q_acc in epsilon_invariant_adjusted])
#     pobin_stddev = math.sqrt(pobin_var)
#     if pobin_mean > 5.0: # implies that mean is at least 2.2 times larger than stddev, this should give fairly symmetrical distribution
#         p_value_norm = norm.sf(float(k-1), loc= pobin_mean, scale= pobin_stddev )

#         print(c_acc, "COMPARE P-VALS:", p_value, "norm:", p_value_norm)


#     #############################
#     #################################### 
#     return  p_value, mult_factor_inv




# def statistical_test(t_acc, X, C, partition_of_X, candidates, ignore_ends_len):
#     significance_values = {}
#     t_seq = C[t_acc]
#     if len(candidates) == 0:
#         significance_values[t_acc] = ("not_tested", "NA", len(partition_of_X[t_acc]), len(partition_of_X[t_acc]), -1 )
#         return significance_values

#     reads = set([x_acc for c_acc in candidates for x_acc in partition_of_X[c_acc]] )
#     reads.update(partition_of_X[t_acc])

#     #### Bug FIX ############
#     total_read_in_partition = len(partition_of_X[t_acc])
#     reads_in_partition = set(partition_of_X[t_acc])
#     # print(partition_of_X[t_acc])
#     for c_acc in candidates:
#         # for x_acc in partition_of_X[c_acc]:
#         #     if x_acc in reads_in_partition:
#         #         print("Already in partition:", x_acc)
#         #         print(x_acc in partition_of_X[t_acc])
#         #         print("candidate:", c_acc)
#         #     else:
#         #         reads_in_partition.add(x_acc)
#         total_read_in_partition += len(partition_of_X[c_acc])
#         # print(partition_of_X[c_acc])

#     ############################

#     N_t = len(reads)

#     # print("N_t:", N_t, "reads in partition:", total_read_in_partition, "ref:", t_acc )
#     # print("Nr candidates:", len(candidates), candidates)
#     assert total_read_in_partition == N_t # each read should be uiquely assinged to a candidate
#     reads_and_candidates = reads.union( [c_acc for c_acc in candidates]) 
#     reads_and_candidates_and_ref = reads_and_candidates.union( [t_acc] ) 

#     # get multialignment matrix here
#     alignment_matrix_to_t, PFM_to_t =  arrange_alignments(t_acc, reads_and_candidates_and_ref, X, C, ignore_ends_len)

#     # cut multialignment matrix first and last ignore_ends_len bases in ends of reference in the amignment matrix
#     # these are bases that we disregard when testing varinats
#     # We get individual cut positions depending on which candidate is being tested -- we dont want to include ends spanning over the reference or candidate
#     # we cut at the start position in c or t that comes last, and the end position in c or t that comes first
#     assert len(list(candidates)) == 1
#     for c_acc in list(candidates):
#         if ignore_ends_len > 0:
#             alignment_matrix_to_t = functions.cut_ends_of_alignment_matrix(alignment_matrix_to_t, t_acc, c_acc, ignore_ends_len)

#         # get parameter estimates for statistical test
#         candidate_accessions = set( [ c_acc for c_acc in candidates] )
#         delta_t = functions.get_difference_coordinates_for_candidates(t_acc, candidate_accessions, alignment_matrix_to_t) # format: { c_acc1 : {pos:(state, char), pos2:(state, char) } , c_acc2 : {pos:(state, char), pos2:(state, char) },... }
#         epsilon, lambda_S, lambda_D, lambda_I = functions.get_error_rates_and_lambda(t_acc, len(t_seq), candidate_accessions, alignment_matrix_to_t) 
#         # get number of reads k supporting the given set of variants, they have to support all the variants within a candidate
#         candidate_support = functions.get_supporting_reads_for_candidates(t_acc, candidate_accessions, alignment_matrix_to_t, delta_t, partition_of_X) # format: { c_acc1 : [x_acc1, x_acc2,.....], c_acc2 : [x_acc1, x_acc2,.....] ,... }
#         candidate_indiv_invariant_factors = functions.adjust_probability_of_candidate_to_alignment_invariant(delta_t, alignment_matrix_to_t, t_acc)

#     for c_acc in list(candidates):
#         original_mapped_to_c = len(partition_of_X[c_acc])
#         k = len(candidate_support[c_acc])
#         # print("supprot:", k, "diff:", len(delta_t[c_acc]))
#         if len(delta_t[c_acc]) == 0:
#             print("{0} no difference to ref {1} after ignoring ends!".format(c_acc, t_acc))
#         p_value, mult_factor_inv = stat_test(k, t_seq, epsilon, delta_t, candidate_indiv_invariant_factors, t_acc, c_acc, original_mapped_to_c)
#         delta_size = len(delta_t[c_acc])
#         significance_values[c_acc] = (p_value, mult_factor_inv, k, N_t, delta_size)
#         print("Tested", c_acc, "to ref", t_acc, "p_val:{0}, mult_factor:{1}, corrected p_val:{2} k:{3}, N_t:{4}, Delta_size:{5}".format(p_value, mult_factor_inv, p_value * mult_factor_inv,  k, N_t, delta_size) )

#     return significance_values


 
#####################################################################################################################  
##### The following three functions are all ways to calculate a p-value from a set of different bernoillis
##### Thr exact test does however not support the -log p_i weights which makes it more or less useless, 
##### The other two are approximations that each have their flaws and cases where they are good approximations
#####################################################################################################################  


# def exact_test(probability, weight, x):
#     """
#         This is for unweighted sum of Bernoullis only!
#         https://stats.stackexchange.com/questions/5347/how-can-i-efficiently-model-the-sum-of-bernoulli-random-variables
#         Could be optimized with this answer:
#             https://stats.stackexchange.com/a/5482
#     """
#     probs = list(probability.values())
#     tmp_distr1 = [Decimal(1.0) - Decimal(probs[0]), Decimal(probs[0])]
#     for i in range(1, len(probs)):
#         tmp_distr2 = [Decimal(1.0) - Decimal(probs[i]), Decimal(probs[i])]
#         distr = [0]*(len(tmp_distr1) + 1)
#         for l in range(len(tmp_distr1)):
#             for m in range(2):
#                 distr[l+m] += tmp_distr1[l] * tmp_distr2[m]

#         tmp_distr1 = distr

#     p_distr = tmp_distr1

#     observed_count = sum([1 for x_i in probability if x_i in x ])    
#     p_value = Decimal(1.0) - sum(p_distr[:observed_count])
#     return float(p_value)



# def poisson_approx_test(probability, weight, x):

#     print([weight[x_i] for x_i in probability if x_i in x ])
#     print([weight[x_i] for x_i in probability ])
#     min_w = min([weight[x_i] for x_i in probability if weight[x_i] > 0 ])
#     weight_multiplier = 1.0 / min_w

#     weight = { x_i : weight[x_i] * weight_multiplier for x_i in weight}
#     print([weight[x_i] for x_i in probability ])

#     observed_weighted_x = sum([weight[x_i]*1.0 for x_i in probability if x_i in x ])
#     po_lambda = sum([ probability[x_i]* weight[x_i] for x_i in probability ])

#     observed_x = sum([1.0 for x_i in probability if x_i in x ])
#     po_lambda = sum([ probability[x_i] for x_i in probability ])
#     if observed_x == 0:
#         k = -1
#     elif observed_x.is_integer():
#         k = observed_x - 1
#     else:
#         k = math.floor(observed_x)
#     # k = observed_weighted_x if observed_weighted_x == 0 or not observed_weighted_x.is_integer() else observed_weighted_x - 1
#     print("k:", k)
#     p_value = poisson.sf(k, po_lambda)

#     return p_value


# def CLT_test(probability, weight, x):

#     observed_weighted_x = sum([weight[x_i]*1.0 for x_i in probability if x_i in x ])
#     mu = sum([ probability[x_i]* weight[x_i] for x_i in probability ])
#     var = sum([ probability[x_i]*(1.0 - probability[x_i])* weight[x_i]**2 for x_i in probability ])
#     sigma = math.sqrt(var)
#     print([weight[x_i] for x_i in probability if x_i in x ])
#     print([weight[x_i] for x_i in probability ])
#     print(mu, sigma, observed_weighted_x)
#     p_value_norm = norm.sf(observed_weighted_x, loc= mu, scale= sigma )

#     return p_value_norm





#####################################################################################################################  
#####################################################################################################################  
##########################################   functions.py  ##########################################################  
#####################################################################################################################  
#####################################################################################################################  
#####################################################################################################################  




# def get_min_uncertainty_per_read(target_accession, segment_length, candidate_accessions, alignment_matrix, invariant_factors_for_candidate):
#     probability = {}
#     assert len(invariant_factors_for_candidate) == 1
#     c_acc = list(invariant_factors_for_candidate.keys())[0]
#     delta_size = float(len(invariant_factors_for_candidate[c_acc]))

#     for q_acc in alignment_matrix:
#         if q_acc == target_accession:
#             continue
#         if q_acc == c_acc:
#             continue 

#     # for q_acc in errors:
#         probability[q_acc] = 1.0
#         p_S =  (delta_size / float(segment_length) ) / 3.0   # p = 0.0 not allowed, min_p is 1/(3*len(seq))
#         p_I =  (delta_size / float(segment_length) ) / 4.0   # p = 0.0 not allowed, min_p is 1/(4*len(seq))
#         p_D =  (delta_size / float(segment_length) )         # p = 0.0 not allowed, min_p is 1/(len(seq))

#         for pos in invariant_factors_for_candidate[c_acc]:
#             for (state, char) in invariant_factors_for_candidate[c_acc][pos]:
#                 u_v = invariant_factors_for_candidate[c_acc][pos][(state, char)]
#                 if state == "S":
#                     probability[q_acc] *= p_S*u_v # *(1.0/u_v)
#                 elif state == "I":
#                     probability[q_acc] *= p_I*u_v #**(1.0/u_v)
#                 elif state == "D":
#                     probability[q_acc] *= p_D*u_v #**(1.0/u_v)
#     return probability

# def get_prob_of_support_per_read(target_accession, segment_length, candidate_accessions, errors, invariant_factors_for_candidate):
#     probability = {}
#     assert len(invariant_factors_for_candidate) == 1
#     c_acc = list(invariant_factors_for_candidate.keys())[0]
#     delta_size = float(len(invariant_factors_for_candidate[c_acc]))

#     for q_acc in errors:
#         probability[q_acc] = 1.0
#         p_S = ( max(errors[q_acc], delta_size) / float(segment_length) ) / 3.0   # p = 0.0 not allowed, min_p is 1/(3*len(seq))
#         p_I = ( max(errors[q_acc], delta_size) / float(segment_length) ) / 4.0   # p = 0.0 not allowed, min_p is 1/(4*len(seq))
#         p_D = ( max(errors[q_acc], delta_size) / float(segment_length) )         # p = 0.0 not allowed, min_p is 1/(len(seq))

#         for pos in invariant_factors_for_candidate[c_acc]:
#             for (state, char) in invariant_factors_for_candidate[c_acc][pos]:
#                 u_v = invariant_factors_for_candidate[c_acc][pos][(state, char)]
#                 if state == "S":
#                     probability[q_acc] *= p_S*u_v # *(1.0/u_v)
#                 elif state == "I":
#                     probability[q_acc] *= p_I*u_v #**(1.0/u_v)
#                 elif state == "D":
#                     probability[q_acc] *= p_D*u_v #**(1.0/u_v)
#     return probability


# def calculate_homopolymenr_lengths(t_seq):
#     homopolymenr_length_numbers = {}

#     h_len = 1
#     for char1, char2 in zip(t_seq[:-1], t_seq[1:]):

#         if char1 != char2:
#             if h_len in homopolymenr_length_numbers:
#                 homopolymenr_length_numbers[h_len] += 1
#             else:
#                 homopolymenr_length_numbers[h_len] = 1
#             h_len = 1

#         else:
#             h_len += 1

#     # end case
#     if h_len in homopolymenr_length_numbers:
#         homopolymenr_length_numbers[h_len] += 1
#     else:
#         homopolymenr_length_numbers[h_len] = 1

#     return homopolymenr_length_numbers


# def get_supporting_reads_for_candidates(target_accession, candidate_accessions, alignment_matrix, Delta_t, partition_of_X):
#     # candidate_support = { c : [] for c in candidate_accessions }
#     # target_alignment = alignment_matrix[target_accession]
#     candidate_support = {}
#     for c in candidate_accessions:
#         candidate_support[c] = []

#         # candidate_alignment = alignment_matrix[c]
#         # for q_acc in alignment_matrix:
#         #     if q_acc == target_accession or q_acc in candidate_accessions:
#         #         continue
#         # print("LEEN", len(partition_of_X[c]), len(partition_of_X[target_accession]) , len(partition_of_X[c].union(partition_of_X[target_accession]) ))
#         # print(len(partition_of_X[c] and  partition_of_X[target_accession]) )
#         # assert len(partition_of_X[c] &  partition_of_X[target_accession]) == 0
#         for q_acc in partition_of_X[c]: #.union(partition_of_X[target_accession]):
#             if q_acc not in  alignment_matrix:
#                 print("READ {0} ALIGNED TO {0} BUT FAILED TO ALIGN TO {1}".format(q_acc, c, target_accession) )
#                 continue
#             query_alignment = alignment_matrix[q_acc]    
#             support = 1
#             for delta in Delta_t[c]:
#                 q_base = query_alignment[delta]
#                 c_state, c_base = Delta_t[c][delta]
#                 if q_base != c_base:
#                     support = 0

#             if support:
#                 candidate_support[c].append(q_acc)
#         # print(candidate_support[c])
#     return candidate_support


# def get_weights_per_read(target_accession, segment_length, candidate_accessions, errors):
#     sum_of_inverse_errors = {}
#     sum_of_inverse_errors["I"] = sum([ 1.0/errors[q_acc]["I"] for q_acc in errors if errors[q_acc]["I"] > 0])
#     sum_of_inverse_errors["D"] = sum([ 1.0/errors[q_acc]["D"] for q_acc in errors if errors[q_acc]["D"] > 0])
#     sum_of_inverse_errors["S"] = sum([ 1.0/errors[q_acc]["S"] for q_acc in errors if errors[q_acc]["S"] > 0])

#     sum_of_inverse_errors = sum([sum_of_inverse_errors["I"], sum_of_inverse_errors["D"], sum_of_inverse_errors["S"]])

#     weight = {}
#     for q_acc in errors:
#         assert q_acc != target_accession and q_acc not in candidate_accessions
#         read_errors = (errors[q_acc]["I"] + errors[q_acc]["D"] + errors[q_acc]["S"])
#         if read_errors > 0:
#             q_errors_inverse =  1.0 / (errors[q_acc]["I"] + errors[q_acc]["D"] + errors[q_acc]["S"])
#         else:
#             q_errors_inverse = 0
        
#         if q_errors_inverse == 0:
#              weight[q_acc] = 1.0 / sum_of_inverse_errors
#         else:
#             weight[q_acc] = q_errors_inverse / sum_of_inverse_errors


#         # weights[q_acc] = {}
#         # for error_type in ["I", "S", "D"]:
#         #     weights[q_acc][error_type] = errors[q_acc][error_type] / float(sum_of_errors[error_type])

#     return weight



# def get_errors_per_read(target_accession, segment_length, candidate_accessions, alignment_matrix):
#     assert len(candidate_accessions) == 1
#     c_acc = list(candidate_accessions)[0]

#     errors = {}
#     target_alignment = alignment_matrix[target_accession]
#     candidate_alignment = alignment_matrix[c_acc]
    
#     for q_acc in alignment_matrix:
#         if q_acc == target_accession:
#             continue
#         if q_acc in candidate_accessions:
#             continue  

#         read_alignment = alignment_matrix[q_acc]
#         errors_to_t = sum( [1 for pos in range(len(target_alignment)) if  target_alignment[pos] != read_alignment[pos] ] )
#         errors_to_c = sum( [1 for pos in range(len(target_alignment)) if  candidate_alignment[pos] != read_alignment[pos] ] )
#         errors[q_acc] = min(errors_to_t, errors_to_c)

#     return errors


# def get_error_rates_and_lambda(target_accession, segment_length, candidate_accessions, alignment_matrix):
#     epsilon = {}
#     target_alignment = alignment_matrix[target_accession]
#     # lambda_S, lambda_D, lambda_I = 0,0,0
#     read_depth = 0
#     ed_poisson_i, ed_poisson_s, ed_poisson_d = 0, 0, 0
    
#     for q_acc in alignment_matrix:
#         if q_acc == target_accession:
#             continue
#         if q_acc in candidate_accessions:
#             continue  

#         epsilon[q_acc] = {}
#         query_alignment = alignment_matrix[q_acc]
#         ed_i, ed_s, ed_d = 0, 0, 0

#         for j in range(len(query_alignment)):
#             q_base = query_alignment[j]
#             t_base = target_alignment[j]
#             if q_base != t_base:
#                 if t_base == "-":
#                     ed_i += 1
#                 elif q_base == "-":
#                     ed_d += 1
#                 else:
#                     ed_s += 1
 

#         # get poisson counts on all positions
#         for j in range(len(query_alignment)):
#             target_alignment = alignment_matrix[target_accession]
#             # candidate_alignment = alignment_matrix[x_to_c_acc[q_acc]]
#             # if j not in forbidden:
#             q_base = query_alignment[j]
#             t_base = target_alignment[j]
#             if q_base != t_base:
#                 if t_base == "-":
#                     ed_poisson_i += 1
#                 elif q_base == "-":
#                     ed_poisson_d += 1
#                 else:
#                     ed_poisson_s += 1      



#         # here we get the probabilities for the poisson counts over each position
#         if q_acc not in candidate_accessions:
#             # lambda_S += ed_s
#             # lambda_D += ed_d
#             # lambda_I += ed_i
#             read_depth += 1


#         epsilon[q_acc]["I"] = (ed_i/float(segment_length))/4.0 
#         epsilon[q_acc]["S"] = (ed_s/float(segment_length))/3.0  
#         epsilon[q_acc]["D"] = ed_d/float(segment_length)


#     # print(ed_poisson_s, ed_poisson_d, ed_poisson_i, float(read_depth), segment_length )
#     lambda_S = max(ed_poisson_s, 2 ) / (float(segment_length) * 3.0) # divisioan by 3 because we have 3 different subs, all equally liekly under our model 
#     lambda_D = max(ed_poisson_d, 2 ) / float(segment_length)
#     lambda_I = max(ed_poisson_i, 2 ) / (float(segment_length) * 4.0)  # divisioan by 4 because we have 4 different ins, all equally liekly under our model 

#         # print(segment_length, ed_i, ed_s, ed_d, epsilon[q_acc]["I"], epsilon[q_acc]["S"], epsilon[q_acc]["D"])
#     return epsilon, lambda_S, lambda_D, lambda_I


# def arrange_query_alignments_to_target(target_accession, target_sequence, query_to_target_alignments):
#     """
#         input:
#             a dictionary query_to_target_alignments of the form
#             {query_accession1 : [(target_aligned, query_aligned, target_start, target_end)],
#              query_accession2 : [(target_aligned, query_aligned, target_start, target_end)],
#              ....
#             }
#         output: 
#             Target vector is a list of 2*len(target_seq) +1 positions
#     """

#     query_to_target_positioned_dict = {} # query_accession : [] list of length 2l_j +1 to hold insertions
#     for query_accession in query_to_target_alignments:
#         target_aligned, query_aligned, target_start, target_end = query_to_target_alignments[query_accession]

#         query_to_target_positioned, target_vector_start_position, target_vector_end_position = position_query_to_alignment(query_aligned, target_aligned, target_start)
#         query_to_target_positioned_dict[query_accession] = (query_to_target_positioned, target_vector_start_position, target_vector_end_position)
        
#     return target_accession, query_to_target_positioned_dict


# def get_non_overlapping_intervals(ranges):
#     """
#         example input:     # ranges = [(0,100,'a'),(0,75,'b'),(95,150,'c'),(120,130,'d')]
#         output: 
#     """
#     non_overlapping_parts = []
#     # ranges = [(0,100,'a'),(0,75,'b'),(95,150,'c'),(120,130,'d')]
#     endpoints = sorted(list(set([r[0] for r in ranges] + [r[1] for r in ranges])))
#     start = {}
#     end = {}
#     for e in endpoints:
#         start[e] = set()
#         end[e] = set()
#     for r in ranges:
#         start[r[0]].add(r[2])
#         end[r[1]].add(r[2])

#     current_ranges = set()
#     prev_set_size = 0
#     for e1, e2 in zip(endpoints[:-1], endpoints[1:]):
#         current_ranges.difference_update(end[e1])
#         current_ranges.update(start[e1])

#         if prev_set_size > len(current_ranges):
#             start_offset = 1
#         else:
#             start_offset = 0

#         next_current_ranges = set(current_ranges)
#         next_current_ranges.difference_update(end[e2])
#         next_current_ranges.update(start[e2])
#         next_set_size = len(next_current_ranges)
#         if next_set_size < len(current_ranges):
#             stop_offset = 0
#         else:
#             stop_offset = -1

#         if current_ranges:
#             # print('%d - %d: %s' % (e1 + start_offset, e2 + stop_offset, ','.join(current_ranges)))
#             non_overlapping_parts.append((e1 + start_offset, e2 + stop_offset, tuple(current_ranges)))
#         else:
#             pass
#             # print('%d - %d: %s' % (e1 + start_offset, e2 + stop_offset, ','.join(current_ranges)), "lol")
        
#         prev_set_size = len(current_ranges)

#     # sys.exit()
#     return(non_overlapping_parts)


# def thread_to_max_ins(max_insertion, q_ins):
#     # else, check if smaller variant can be aligned from left to right with all nucleotides matching in q_ins, e.g. say max deletion is GACG
#     # then an insertion AG may be aligned as -A-G. Take this alignment instead
#     can_be_threaded = True
#     prev_pos = -1
#     match_pos = set()
#     for q_nucl in q_ins:  # TODO: WHAT IF INSERTION AGA HERE WITH MAX_INSTERTION TAGA? match_pos will be {1,2,1} and can_be_threaded False, this is incorrect!!
#         pos = max_insertion[(prev_pos+1):].find(q_nucl) 
#         if pos < 0:
#             can_be_threaded = False
#             break
#         else:
#             match_pos.add( (pos + (prev_pos+1)) ) 
#         prev_pos = pos

#     if can_be_threaded:
#         q_insertion_modified = ""
#         for p in range(len(max_insertion)):
#             if p in match_pos:
#                 nucl = max_insertion[p]
#             else:
#                 nucl = "-"
#             q_insertion_modified = q_insertion_modified + nucl
#         # print("NEW can be threaded: q:{0}, max: {1}, new thread: {2}".format(q_ins, max_insertion, q_insertion_modified))
#         return q_insertion_modified
#     else:
#         return ""




#####################################################################################################################  
#####################################################################################################################  
##########################################   correcton_module.py  ################################################### 
#####################################################################################################################  
#####################################################################################################################  
##################################################################################################################### 




# def correct_to_nearest_neighbor(m, partition, seq_to_acc):
#     S_prime_partition = {}

#     N_t = sum([container_tuple[3] for s, container_tuple in partition.items()]) # total number of sequences in partition
    
#     if N_t == 2:
#         print("Partition has size", N_t, "no meaningful correction can be done")

#     if len(partition) > 1 and N_t > 2:
#         # all strings has not converged
#         alignment_matrix, PFM = create_position_probability_matrix(m, partition) 
        
#         # print("nearest_neighbor errors:",  math.ceil(min([ partition[s][0] for s in partition if partition[s][3] > 1 or s !=m ]) / 2.0)  )
#         # nearest_neighbor_errors = min([ partition[s][0] for s in partition if partition[s][3] > 1 or s !=m ])
#         # nearest_neighbor_errors = math.ceil(min([ partition[s][0] for s in partition if partition[s][3] > 1 or s !=m ]) / 2.0)

#         ### TEST LOG ERROR TYPES #######
#         # c = Counter()
#         # for j in range(len(PFM)):
#         #     max_v_j = max(PFM[j], key = lambda x: PFM[j][x] )
#         #     for v in PFM[j]:
#         #         if v != max_v_j:
#         #            c[v] += PFM[j][v]
#         # print("Error types:", c, "depth:", len(partition) )
#         #############################

#         ## TEST LOG ERROR TYPES #######
#         c_del = 0
#         c_ins = 0
#         c_subs = 0
#         for j in range(len(PFM)):
#             max_v_j = max(PFM[j], key = lambda x: PFM[j][x] )
#             for v in PFM[j]:
#                 if max_v_j == "-":
#                     if v != max_v_j:
#                         c_ins += PFM[j][v]
#                 else:
#                     if v != max_v_j:
#                         if v == "-":
#                             c_del += PFM[j][v]
#                         else:
#                             c_subs += PFM[j][v]

#         print("Error types:", c_del, c_ins, c_subs, "depth:", len(partition) )
#         ############################

#         for s in partition:
#             # if nearest_neighbor_errors < partition[s][0]:
#             #     nr_pos_to_correct = int( partition[s][0] - nearest_neighbor_errors ) #decide how many errors we should correct here
#             # else:
#             #     nr_pos_to_correct = int(math.ceil(partition[s][0] / 2.0)) #decide how many errors we should correct here
#             # nr_pos_to_correct = max(int( partition[s][0] - nearest_neighbor_errors ), int(math.ceil(partition[s][0] / 2.0)))
#             nr_pos_to_correct = int(math.ceil(partition[s][0] / 2.0)) #decide how many errors we should correct here

#             # print("positions to correct for sequence s:", nr_pos_to_correct, s ==m)
#             if nr_pos_to_correct  == 0:
#                 continue

#             s_alignment_in_matrix = alignment_matrix[s]
#             # find the position probabilities of the alignment of s in PFM

#             pos_freqs_for_s = []
#             for j in range(len(PFM)):
#                 pos_freqs_for_s.append( (j, PFM[j][s_alignment_in_matrix[j]]) )

#             pos_freqs_for_s.sort(key=lambda x: x[1]) # sort with respect to smallest frequencies                    
#             pos, highest_freq_of_error_to_correct = pos_freqs_for_s[ nr_pos_to_correct - 1 ]
#             end_position_in_list = nr_pos_to_correct

#             pp = nr_pos_to_correct
#             while pos_freqs_for_s[pp][1] == highest_freq_of_error_to_correct:
#                 end_position_in_list += 1
#                 pp += 1

#             J = [j for j, freq in random.sample(pos_freqs_for_s[:end_position_in_list], nr_pos_to_correct)]




#             ########### TEST WEIGHTING EACH MINORITY POSITION BY IT'S OBSERVED FREQUENCY THROUGHOUT THE ALIGNMENTS TO THE nearest_neighbor ################
#             # pos_freqs_for_s_mod = []
#             # for j in range(len(PFM)):
#             #     v_j = s_alignment_in_matrix[j]
#             #     pos_freqs_for_s_mod.append( (j, PFM[j][v_j] / float(max(c[v_j], 1) ) ))
#             # pos_freqs_for_s_mod.sort(key=lambda x: x[1]) # sort with respect to smallest frequencies                    
#             # pos, highest_freq_of_error_to_correct = pos_freqs_for_s_mod[ nr_pos_to_correct - 1 ]
#             # end_position_in_list = nr_pos_to_correct
#             # for pp in range(nr_pos_to_correct, len(pos_freqs_for_s_mod)):
#             #     # print(pos_freqs_for_s_mod[pp][1], highest_freq_of_error_to_correct)
#             #     if pos_freqs_for_s_mod[pp][1] > highest_freq_of_error_to_correct:
#             #         break
#             #     else:
#             #         end_position_in_list += 1
#             # J = [j for j, freq in random.sample(pos_freqs_for_s_mod[:end_position_in_list], nr_pos_to_correct)]
#             #############################################


#             # ####### TEST CHOOSING RANDOM SUBSET OUT OF ALL MINORITY POSITONS IN THE READ ################
#             # minority_positions_for_s = []
#             # for j in range(len(PFM)):
#             #     count_v_j = PFM[j][s_alignment_in_matrix[j]]
#             #     max_v_j = max(PFM[j], key = lambda x: PFM[j][x] )
#             #     if count_v_j < PFM[j][max_v_j]:
#             #         minority_positions_for_s.append(j)
#             # # print(len(minority_positions_for_s))

#             # if nr_pos_to_correct > len(minority_positions_for_s):
#             #     print("OMFG!!", len(minority_positions_for_s), nr_pos_to_correct)
#             #     nr_pos_to_correct = len(minority_positions_for_s)
#             # J = random.sample(minority_positions_for_s, nr_pos_to_correct)
#             ##############################################

#             ########### TEST WEIGHTING EACH MINORITY POSITION BY IT'S OBSERVED FREQUENCY THROUGHOUT THE ALIGNMENTS TO THE nearest_neighbor ################
#             pos_freqs_for_s_mod = []
#             for j in range(len(PFM)):
#                 max_v_j = max(PFM[j], key = lambda x: PFM[j][x] )
#                 v_j = s_alignment_in_matrix[j]
#                 if max_v_j == v_j:
#                     pos_freqs_for_s_mod.append( (j, PFM[j][v_j] / float(1)) )
#                 elif max_v_j == "-":
#                     pos_freqs_for_s_mod.append( (j, PFM[j][v_j] / float(max(c_ins, 1) ) ))
#                 elif v_j == "-":
#                     pos_freqs_for_s_mod.append( (j, PFM[j][v_j] / float(max(c_del, 1) ) ))
#                 else:
#                     pos_freqs_for_s_mod.append( (j, PFM[j][v_j] / float(max(c_subs, 1) ) ))

#             pos_freqs_for_s_mod.sort(key=lambda x: x[1]) # sort with respect to smallest frequencies                    
#             pos, highest_freq_of_error_to_correct = pos_freqs_for_s_mod[ nr_pos_to_correct - 1 ]
#             end_position_in_list = nr_pos_to_correct
#             for pp in range(nr_pos_to_correct, len(pos_freqs_for_s_mod)):
#                 # print(pos_freqs_for_s_mod[pp][1], highest_freq_of_error_to_correct)
#                 if pos_freqs_for_s_mod[pp][1] > highest_freq_of_error_to_correct:
#                     break
#                 else:
#                     end_position_in_list += 1
#             J = [j for j, freq in random.sample(pos_freqs_for_s_mod[:end_position_in_list], nr_pos_to_correct)]
#             #############################################

#             # J = [j for j, prob in pos_freqs_for_s[:nr_pos_to_correct]] # J is the set of the nr_pos_to_correct smallest position probabilities
#             # print(nr_pos_to_correct, end_position_in_list)
#             # print(pos_freqs_for_s[:end_position_in_list])
#             # print(J)

#             s_new = alignment_matrix[s]
#             for j in J:
#                 old_nucl = s_new[j]
#                 highest_prob_character_at_j = max(PFM[j], key=lambda k: PFM[j][k])

#                 if highest_prob_character_at_j == old_nucl: # choose the other highest on if tie (should happen only when partition consist of two sequences)
#                     pmf_j_minus_variant = copy.deepcopy(PFM[j])
#                     del pmf_j_minus_variant[old_nucl] 
#                     highest_prob_character_at_j = max(pmf_j_minus_variant, key=lambda k: pmf_j_minus_variant[k])


#                 # print("correcting", s_new[j], "to", highest_prob_character_at_j )
#                 s_new[j] = highest_prob_character_at_j
#             s_modified = "".join([nucl for nucl in s_new if nucl != "-" ])

#             # only unique strings can change in this step

#             # accession_of_s = unique_seq_to_acc[s] # this is still unique
#             # S_prime_partition[accession_of_s] = s_modified

#             accessions_of_s = seq_to_acc[s]
#             for acc in accessions_of_s:
#                 S_prime_partition[acc] = s_modified
    
#     return S_prime_partition




#####################################################################################################################  
#####################################################################################################################  
##########################################  nearest_neighbor_graph.py  #############################################
#####################################################################################################################  
#####################################################################################################################  
##################################################################################################################### 


# def histogram(data, args, name='histogram.png', x='x-axis', y='y-axis', x_cutoff=None, title=None):
#     if x_cutoff: 
#         plt.hist(data, range=[0, x_cutoff], bins = 100)
#     else:
#         plt.hist(data, bins = 100)
#     plt.xlabel(x)
#     plt.ylabel(y)
#     if title:
#         plt.title(title)

#     plt.savefig(os.path.join(args.outfolder, name))
#     plt.clf()


# def collapse(seq):
#     seen = set()
#     seen_add = seen.add
#     return [x for x in seq if not (x in seen or seen_add(x))]


# def edlib_traceback(x, y, mode="NW", task="path", k=1):
#     result = edlib.align(x, y, mode=mode, task=task, k=k)
#     ed = result["editDistance"]
#     locations =  result["locations"]
#     cigar =  result["cigar"]
#     return ed, locations, cigar

# def reverse_complement(string):
#     #rev_nuc = {'A':'T', 'C':'G', 'G':'C', 'T':'A', 'N':'N', 'X':'X'}
#     # Modified for Abyss output
#     rev_nuc = {'A':'T', 'C':'G', 'G':'C', 'T':'A', 'a':'t', 'c':'g', 'g':'c', 't':'a', 'N':'N', 'X':'X', 'n':'n', 'Y':'R', 'R':'Y', 'K':'M', 'M':'K', 'S':'S', 'W':'W', 'B':'V', 'V':'B', 'H':'D', 'D':'H', 'y':'r', 'r':'y', 'k':'m', 'm':'k', 's':'s', 'w':'w', 'b':'v', 'v':'b', 'h':'d', 'd':'h'}

#     rev_comp = ''.join([rev_nuc[nucl] for nucl in reversed(string)])
#     return(rev_comp)



# # create the parser for the "remove_barcodes" command
# remove_barcodes = subparsers.add_parser('remove_barcodes', help='Remove barcodes from sequences.')
# remove_barcodes.add_argument('-fl_reads', type=str, help='Path to the consensus fasta file')
# remove_barcodes.add_argument('outfolder', type=str, help='Output folder to results')
# remove_barcodes.add_argument('--cleanup', action='store_true', help='Remove everything except logfile.txt, candidates_converged.fa and final_candidates.fa in output folder.')
# remove_barcodes.set_defaults(which='remove_barcodes')


