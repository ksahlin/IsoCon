"""
    PFM is a list of dicts, where the inner dicts are one per column in the PFM matrix, its keys by characters A, C, G, T, -
    alignment_matrix is a representation of all alignments in a partition. this is a dictionary where sequences s_i belonging to the 
    partition as keys and the alignment of s_i with respectt to the alignment matix.
"""
import copy
import unittest
import math

from scipy.stats import poisson, binom, norm

from partitions import partition_strings_paths, partition_strings_2set_paths, partition_strings
from functions import create_position_probability_matrix, transpose, get_error_rates_and_lambda, get_difference_coordinates_for_candidates, get_supporting_reads_for_candidates
from modules import graphs
from SW_alignment_module import sw_align_sequences, sw_align_sequences_keeping_accession
from input_output import fasta_parser

from multinomial_distr import multinomial_

import networkx as nx
import matplotlib.pyplot as plt

def get_unique_seq_accessions(S):
    seq_to_acc = {}
    for acc, seq in  S.items():
        if seq in seq_to_acc:
            seq_to_acc[seq].append(acc)
        else: 
            seq_to_acc[seq] = []
            seq_to_acc[seq] = [acc]

    unique_seq_to_acc = {seq: acc_list[0] for seq, acc_list in  seq_to_acc.items() if len(acc_list) == 1 } 
    print("Unique seqs left:", len(unique_seq_to_acc))

    return unique_seq_to_acc

def get_partition_alignments(graph_partition, M, G_star):
    exact_alignments = sw_align_sequences(graph_partition, single_core = False)

    partition_alignments = {} 
    for m in M:
        indegree = 1 if m not in G_star[m] else G_star[m][m]
        partition_alignments[m] = { m : (0, m, m, indegree) }
        if m not in exact_alignments:
            continue
        else:
            for s in exact_alignments[m]:
                aln_m, aln_s, (matches, mismatches, indels) = exact_alignments[m][s]
                edit_dist = mismatches + indels
                # indegree =  1 if s not in G_star[m] else G_star[m][s]
                # if indegree > 1:
                #     print("Larger than 1!!", indegree)
                partition_alignments[m][s] = (edit_dist, aln_m, aln_s, 1)

    print("NR candidates:", len(partition_alignments))
    return partition_alignments



def find_candidate_transcripts(X):
    """
        input: a string pointing to a fasta file
        output: a string containing a path to a fasta formatted file with consensus_id_support as accession 
                    and the sequence as the read
    """ 
    S = {acc: seq for (acc, seq) in  fasta_parser.read_fasta(open(X, 'r'))}
    # S = X
    C = {}
    unique_seq_to_acc = get_unique_seq_accessions(S)
    
    G_star, graph_partition, M, converged = partition_strings_paths(S)
    partition_alignments = get_partition_alignments(graph_partition, M, G_star)       

    if converged:
        for m in M:
            N_t = sum([container_tuple[3] for s, container_tuple in partition_alignments[m].items()])
            C[m] = N_t           
        return C

    step = 1
    prev_edit_distances_2steps_ago = [2**28,2**28,2**28] # prevents 2-cycles
    prev_edit_distances = [2**28]

    while not converged:
        print("candidates:", len(M))
        edit_distances = [ partition_alignments[s1][s2][0] for s1 in partition_alignments for s2 in partition_alignments[s1]  ] 
        print("edit distances:", edit_distances) 

        ###### Different convergence criterion #########

        if prev_edit_distances_2steps_ago == edit_distances:
            # Only cyclic alignments are left, these are reads that jump between two optimal alignment of two different
            # target sequeneces. This is a product of our fast heurustics of defining a minmap score + SSW filtering to choose best alignment
            print("CYCLE!!!")
            assert len(partition_alignments) == len(M)
            break             
        if sum(edit_distances) > sum(prev_edit_distances) and  max(edit_distances) > max(prev_edit_distances) :
            #return here if there is some sequence alternating between best alignments and gets corrected and re-corrected to different candidate sequences
            assert len(partition_alignments) == len(M)
            print("exiting here!")
            break            

        has_converged = [True if ed == 0 else False for ed in edit_distances] 
        if all(has_converged):
            # we return here if tha data set contain isolated nodes.
            assert len(partition_alignments) == len(M)
            break
        #######################################################

   

        for m, partition in partition_alignments.items():
            N_t = sum([container_tuple[3] for s, container_tuple in partition.items()]) # total number of sequences in partition
            # print("cluster size:", N_t)
            # if N_t < 3:
            #     # print("skipping")
            #     print("lolling")
            #     for s in partition:
            #         print(partition_alignments[m][s][0])
            #     continue
            # print(partition)
            if len(partition) > 1:
                # all strings has not converged
                alignment_matrix, PFM = create_position_probability_matrix(m, partition) 

                for s in partition:
                    nr_pos_to_correct = int(math.ceil(partition_alignments[m][s][0] / 2.0)) #decide how many errors we should correct here
                    # print("positions to correct for sequence s:", nr_pos_to_correct, s ==m)
                    if nr_pos_to_correct  == 0:
                        continue

                    s_alignment_in_matrix = alignment_matrix[s]
                    # find the position probabilities of the alignment of s in PFM
                    pos_freqs_for_s = []
                    for j in range(len(PFM)):
                        # try:
                        pos_freqs_for_s.append( (j, PFM[j][s_alignment_in_matrix[j]]) )
                        # except KeyError:
                        #     print(j, PFM[j], s_alignment_in_matrix[j], N_t, len(partition), len(PFM), len(m) )
                        #     sys.exit()

                    pos_freqs_for_s.sort(key=lambda x: x[1]) # sort with respect to smallest frequencies                    
                    J = [j for j, prob in pos_freqs_for_s[:nr_pos_to_correct]] # J is the set of the nr_pos_to_correct smallest position probabilities
                    s_new = alignment_matrix[s]
                    for j in J:
                        old_nucl = s_new[j]
                        highest_prob_character_at_j = max(PFM[j], key=lambda k: PFM[j][k])

                        if highest_prob_character_at_j == old_nucl: # choose the other highest on if tie (should happen only when partition consist of two sequences)
                            pmf_j_minus_variant = copy.deepcopy(PFM[j])
                            del pmf_j_minus_variant[old_nucl] 
                            highest_prob_character_at_j = max(pmf_j_minus_variant, key=lambda k: pmf_j_minus_variant[k])


                        # print("correcting", s_new[j], "to", highest_prob_character_at_j )
                        s_new[j] = highest_prob_character_at_j
                    s_modified = "".join([nucl for nucl in s_new if nucl != "-" ])

                    # only unique strings can change in this step

                    accession_of_s = unique_seq_to_acc[s] # this is still unique
                    S[accession_of_s] = s_modified

        print("Tot seqs:", len(S))
        unique_seq_to_acc = get_unique_seq_accessions(S)
        # partition_alignments, partition, M, converged = partition_strings(S)

 

        G_star, graph_partition, M, converged = partition_strings_paths(S)
        partition_alignments = get_partition_alignments(graph_partition, M, G_star)  

        out_file = open("/Users/kxs624/tmp/minimizer_test_1000_step_" +  str(step) + ".fa", "w")
        for i, m in enumerate(partition_alignments):
            N_t = sum([container_tuple[3] for s, container_tuple in partition_alignments[m].items()])
            out_file.write(">{0}\n{1}\n".format("read" + str(i)+ "_support_" + str(N_t) , m))

        step += 1
        prev_edit_distances_2steps_ago = prev_edit_distances
        prev_edit_distances = edit_distances


    for m in M:
        N_t = sum([container_tuple[3] for s, container_tuple in partition_alignments[m].items()])
        C[m] = N_t    
    out_file = open("/Users/kxs624/tmp/TSPY13P_2_constant_constant_0.0001_converged.fa", "w")
    for i, m in enumerate(partition_alignments):
        out_file.write(">{0}\n{1}\n".format("read_" + str(i)+ "_support_" + str(C[m]) , m))   
    out_file.close()    
    return out_file.name

def get_partition_alignments_2set(graph_partition, C, X):
    partition_dict = {}
    for t, x_hits in graph_partition.items():
        partition_dict[t] = {}
        for x in x_hits:
            if x in X: # read assigned to alignment matrix
                partition_dict[t][x] = (C[t], X[x]) 
            else: # a candidate we incluse to test, it needs to be aligned w.r.t. the alignment matrix
                partition_dict[t][x] = (C[t], C[x]) 

    exact_alignments = sw_align_sequences_keeping_accession(partition_dict, single_core = False)
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

def filter_C_X_and_partition(X, C, G_star, partition):
    # if there are x in X that is not in best_exact_matches, these x had no (decent) alignment to any of
    # the candidates, simply skip them.

    print("total reads:", len(X), "total reads with an alignment:", len(G_star) )

    for x in X.keys():
        if x not in G_star:
            print("read missing alignment",x)
            del X[x]

    # also, filter out the candidates that did not get any alignments here.

    print("total consensus:", len(C), "total consensus with at least one alignment:", len(partition) )

    for c in C.keys():
        if c not in partition:
            print("candidate missing hit:", c)
            del C[c]
        elif len(partition[c]) == 0:
            print("candidate had edges in minimizer graph but they were all shared and it did not get chosen:", c)
            del C[c]            
            del partition[c]
        else:
            print(c, " read hits:", len(partition[c]))
    return X, C

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

def three_CO(read_file, candidate_file = ""):
    ################################### PROCESS INDATA #####################################
    if not candidate_file:
        candidate_file = find_candidate_transcripts(read_file)
    X_file = read_file
    C_file = candidate_file
    X = {acc: seq for (acc, seq) in  fasta_parser.read_fasta(open(read_file, 'r'))} 
    C = {acc: seq for (acc, seq) in  fasta_parser.read_fasta(open(candidate_file, 'r'))} 

    # C = {c: support for c, support in C.items() if support > 3}
    # TODO: eventually filter candidates with lower support than 2-3? Here?
    G_star, partition_of_X =  partition_strings_2set_paths(X, C, X_file, C_file)
    X, C = filter_C_X_and_partition(X, C, G_star, partition_of_X)
    C_seq_to_acc = {seq : acc for acc, seq in C.items()}

    #########################################################################################


    modified = True
    # changed_nodes = set(C.keys())
    step = 1
    while modified:
        modified = False
        print("NEW STEP")
        # G_star_C, alignment_graph, converged = graphs.construct_minimizer_graph(C)
        weights = { C[c_acc] : len(x_hits) for c_acc, x_hits in partition_of_X.items()} 
        G_star_C, partition_of_C, M, converged = partition_strings_paths(C, node_weights = weights, edge_creating_min_treshold = 5,  edge_creating_max_treshold = 10)
        # self edges not allowed
        print(len(C), len(partition_of_C), len(M) )


        # get all candidats that serve as null-hypothesis references and have neighbors subject to testing
        # these are all candidates that are minimizers to some other, isolated nodes is not tested
        # candidatate in G_star_C
        null_hypothesis_references_to_candidates = [c for c in partition_of_C if partition_of_C[c] ]
        print("References in testing:", len(null_hypothesis_references_to_candidates))
        nr_of_tests_this_round = len([1 for c1 in partition_of_C for c2 in  partition_of_C[c1]])

        ######### just for vizualization ###############
        # D=nx.DiGraph(G_star_C)
        D=nx.DiGraph()
        lables = {}
        for c1 in partition_of_C:
            c1_acc = C_seq_to_acc[c1]
            reads_to_c1 = [X[x_acc] for x_acc in  partition_of_X[c1_acc] ]
            w_c1 = len(reads_to_c1)
            for c2 in  partition_of_C[c1]:
                c2_acc = C_seq_to_acc[c2]
                reads_to_c2 = [X[x_acc] for x_acc in  partition_of_X[c2_acc] ]
                w_c2 = len(reads_to_c2)
                D.add_edge(c2,c1, weight = str(w_c2) + "->" + str(w_c1)  )
        labels = nx.get_edge_attributes(D, 'weight')
        # pos = nx.circular_layout(D)
        pos = dict()
        XX = [c2 for c1 in partition_of_C for c2 in partition_of_C[c1] ] # have in-edges
        YY = [c1 for c1 in partition_of_C ] # not
        pos.update( (n, (1, 4*i)) for i, n in enumerate(XX) ) # put nodes from X at x=1
        pos.update( (n, (2, 2*i)) for i, n in enumerate(YY) ) # put nodes from Y at x=2
        nx.draw_networkx_nodes(D, pos, node_size=50 )
        nx.draw_networkx_edge_labels(D, pos, arrows=True, edge_labels=labels)
        nx.draw_networkx_edges(D, pos, arrows=True, edge_labels=labels)
        plt.savefig("/Users/kxs624/tmp/Graph_bip_1000_step_" + str(step) + ".png", format="PNG")
        plt.clf()
        p_vals = []

        ####################################################

        # based on how G_star is shaped, we decide which reads to align to the candidate under the null-hypothesis
        # i.e., all reads in partition_of_X[c1] + partition_of_X[c2] should be aligned to c2. After that we estimate error rates etc.
        #or should we align al reads in S* to a cluster center c2?? that is all reads belonging to canditades that has the same minimizer?
        # we would probably get better statistical power in each test.. maybe its also simpler to implement. 



        # do one reference candidate at a time, these are all modular and this loop 
        # can easily be parallellized if we break this up to a function
        # p_values_to_t = wrapper_statistical_test()

        C_pvals = { C_seq_to_acc[c] : (len(partition_of_X[C_seq_to_acc[c]]), "not_tested", len(partition_of_X[C_seq_to_acc[c]]) ) for c in partition_of_C if not partition_of_C[c]} # initialize with transcripts not tested
        for t in null_hypothesis_references_to_candidates:
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
            
            C_pvals[t_acc] = (len(reads_to_t), -1, len(reads_to_map)) 
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
                p_vals.append(p_value)

                if p_value > 0.05/nr_of_tests_this_round or k == 0:
                    modified = True
                    del C[c_acc]
                    partition_of_X[t_acc].update(partition_of_X[c_acc])
                    del partition_of_X[c_acc]
                    print("deleting:",p_value, "k:", k, "delta size:", len(delta_t[c_acc]))
                else:
                    C_pvals[c_acc] = (k, p_value, N_t)

        print("nr candidates left:", len(C))
        print(p_vals)
        plt.hist(p_vals)
        plt.savefig("/Users/kxs624/tmp/p_vals_step_" + str(step) +".png", format="PNG")
        plt.clf()
        step += 1
        # sys.exit()
 
    out_file = open("/Users/kxs624/tmp/final_candidates_TSPY13P_pass2_2_constant_constant_0.0001.fa", "w")
    for c_acc, seq in C.items():
        support, p_value, N_t = C_pvals[c_acc] 
        out_file.write(">{0}\n{1}\n".format(c_acc + "_" + str(support) + "_" + str(p_value) + "_" + str(N_t) , seq))

    return C



class TestFunctions(unittest.TestCase):

    # def test_find_candidate_transcripts(self):
    #     self.maxDiff = None

    #     from input_output import fasta_parser
    #     try:
    #         fasta_file_name = "/Users/kxs624/Documents/data/pacbio/simulated/ISOseq_sim_n_200/simulated_pacbio_reads.fa"
    #         # fasta_file_name = "/Users/kxs624/Documents/data/pacbio/simulated/RBMY_44_-_constant_-.fa"
    #         # fasta_file_name = "/Users/kxs624/Documents/data/pacbio/simulated/DAZ2_2_exponential_constant_0.001.fa"
    #         # fasta_file_name = "/Users/kxs624/Documents/data/pacbio/simulated/TSPY13P_2_constant_constant_0.0001.fa"
    #         # fasta_file_name = "/Users/kxs624/Documents/data/pacbio/simulated/TSPY13P_4_linear_exponential_0.05.fa"
    #         # fasta_file_name = "/Users/kxs624/Documents/data/pacbio/simulated/HSFY2_2_constant_constant_0.0001.fa"

    #         S = {acc: seq for (acc, seq) in  fasta_parser.read_fasta(open(fasta_file_name, 'r'))} 
    #     except:
    #         print("test file not found:",fasta_file_name) 

    #     partition_alignments = find_candidate_transcripts(S)
    #     print(len(partition_alignments))

    def test_three_CO(self):
        self.maxDiff = None

        from input_output import fasta_parser
        try:
            # read_file_name = "/Users/kxs624/Documents/data/pacbio/simulated/RBMY_44_-_constant_-.fa"
            # read_file_name = "/Users/kxs624/Documents/data/pacbio/simulated/ISOseq_sim_n_1000/simulated_pacbio_reads.fa"
            # read_file_name = "/Users/kxs624/Documents/data/pacbio/simulated/DAZ2_2_exponential_constant_0.001.fa"
            read_file_name = "/Users/kxs624/Documents/data/pacbio/simulated/TSPY13P_2_constant_constant_0.0001.fa"
            # read_file_name = "/Users/kxs624/Documents/data/pacbio/simulated/TSPY13P_4_linear_exponential_0.05.fa"
            # read_file_name = "/Users/kxs624/Documents/data/pacbio/simulated/HSFY2_2_constant_constant_0.0001.fa"

            # X = {acc: seq for (acc, seq) in  fasta_parser.read_fasta(open(read_file_name, 'r'))} 
        except:
            print("test file not found:",read_file_name) 

        try:
            # consensus_file_name = "/Users/kxs624/tmp/minimizer_test_1000_converged.fa"
            # consensus_file_name = ""
            # consensus_file_name = "/Users/kxs624/tmp/minimizer_consensus_final_RBMY_44_-_constant_-.fa"
            # consensus_file_name = "/Users/kxs624/tmp/minimizer_consensus_DAZ2_2_exponential_constant_0.001_step10.fa"
            consensus_file_name = "/Users/kxs624/tmp/TSPY13P_2_constant_constant_0.0001_converged.fa"
            # consensus_file_name = "/Users/kxs624/tmp/final_candidates_TSPY13P_2_constant_constant_0.0001.fa"
            # consensus_file_name = "/Users/kxs624/Documents/data/pacbio/simulated/TSPY13P_4_linear_exponential_0.05.fa"
            # consensus_file_name = "/Users/kxs624/Documents/data/pacbio/simulated/HSFY2_2_constant_constant_0.0001.fa"

            # C = {acc: seq for (acc, seq) in  fasta_parser.read_fasta(open(consensus_file_name, 'r'))} 
        except:
            print("test file not found:",consensus_file_name) 

        three_CO(read_file_name, consensus_file_name)

if __name__ == '__main__':
    unittest.main()