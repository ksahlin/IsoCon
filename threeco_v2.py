"""
    PFM is a list of dicts, where the inner dicts are one per column in the PFM matrix, its keys by characters A, C, G, T, -
    alignment_matrix is a representation of all alignments in a partition. this is a dictionary where sequences s_i belonging to the 
    partition as keys and the alignment of s_i with respectt to the alignment matix.
"""
import os
import unittest
import copy
import math 

from scipy.stats import poisson

from modules.functions import transpose, create_position_probability_matrix
from modules import functions
from modules.partitions import partition_strings_paths, partition_strings_2set, partition_to_statistical_test
from modules import graphs
from modules.SW_alignment_module import sw_align_sequences, sw_align_sequences_keeping_accession
from modules.edlib_alignment_module import edlib_align_sequences, edlib_align_sequences_keeping_accession, edlib_traceback
from modules.input_output import fasta_parser, write_output
from modules import statistical_test_v2
from modules import correct_sequence_to_minimizer


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
    exact_edit_distances = edlib_align_sequences(graph_partition, single_core = False)    
    exact_alignments = sw_align_sequences(exact_edit_distances, single_core = False)

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



def find_candidate_transcripts(read_file, params):
    """
        input: a string pointing to a fasta file
        output: a string containing a path to a fasta formatted file with consensus_id_support as accession 
                    and the sequence as the read
    """ 
    S = {acc: seq for (acc, seq) in  fasta_parser.read_fasta(open(read_file, 'r'))}

    C = {}
    unique_seq_to_acc = get_unique_seq_accessions(S)

    G_star, graph_partition, M, converged = partition_strings_paths(S, params)
    partition_alignments = get_partition_alignments(graph_partition, M, G_star)       

    # if converged:
    #     for m in M:
    #         N_t = sum([container_tuple[3] for s, container_tuple in partition_alignments[m].items()])
    #         C[m] = N_t           
    #     return C

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


        # TODO: Parallelize this part over partitions: sent in the partition and return a dict s_acc : modified string
        S_prime = correct_sequence_to_minimizer.correct_strings(partition_alignments, unique_seq_to_acc, single_core = False)

        for acc, s_prime in S_prime.items():
            S[acc] = s_prime

        print("Tot seqs:", len(S))
        unique_seq_to_acc = get_unique_seq_accessions(S)
        # partition_alignments, partition, M, converged = partition_strings(S)

 

        G_star, graph_partition, M, converged = partition_strings_paths(S, params)
        partition_alignments = get_partition_alignments(graph_partition, M, G_star)  
        out_file_name = os.path.join(params.outfolder, "candidates_step_" +  str(step) + ".fa")
        out_file = open(out_file_name, "w")
        for i, m in enumerate(partition_alignments):
            N_t = sum([container_tuple[3] for s, container_tuple in partition_alignments[m].items()])
            out_file.write(">{0}\n{1}\n".format("read" + str(i)+ "_support_" + str(N_t) , m))

        step += 1
        prev_edit_distances_2steps_ago = prev_edit_distances
        prev_edit_distances = edit_distances



    for m in M:
        N_t = sum([container_tuple[3] for s, container_tuple in partition_alignments[m].items()])
        C[m] = N_t   

    original_reads = {acc: seq for (acc, seq) in  fasta_parser.read_fasta(open(read_file, 'r'))}
    reads_to_minimizers = {}
    # [for m, partition in partition_alignments.items() for s in partition]

    for read_acc, seq in original_reads.items():
        m = S[read_acc]
        # get the correspoindng minimizer for each read, if read does not belong to an m, its because it did not converge.
        if m not in C:
            print("Minimizer to read did not converge")
            continue
        elif C[m] >= params.min_candidate_support:
            reads_to_minimizers[read_acc] = { m : (original_reads[read_acc], m)}
        else:
            print("Minimizer did not pass threshold support of {0} reads.".format(C[m]))
            del C[m]

    edit_distances_of_x_to_m = edlib_align_sequences_keeping_accession(reads_to_minimizers)
    alignments_of_x_to_m = sw_align_sequences_keeping_accession(edit_distances_of_x_to_m)
    alignments_of_x_to_m_filtered, m_to_acc = filter_candidates(alignments_of_x_to_m, C, params)
    candidates_file_name = os.path.join(params.outfolder, "candidates_converged.fa")
    alignments_file_name = os.path.join(params.outfolder, "candidate_alignments.tsv")
    write_output.print_candidates_from_minimizers(candidates_file_name, alignments_file_name, alignments_of_x_to_m_filtered, C, m_to_acc, params)

    return candidates_file_name, alignments_of_x_to_m_filtered

def filter_candidates(alignments_of_x_to_c, C, params):
    alignments_of_x_to_c_transposed = transpose(alignments_of_x_to_c)   

    m_to_acc = {}

    for i, (m, support) in enumerate(list(C.items())):
        m_acc = "read_" + str(i) + "_support_" + str(support)
        m_to_acc[m] = m_acc
        #require support from at least 4 reads if not tested (consensus transcript had no close neighbors)
        # add extra constraint that the candidate has to have majority on _each_ position in c here otherwise most likely error
        if support >= params.min_candidate_support:
            # print("needs to be consensus over each base pair")
            partition_alignments_c = {m : (0, m, m, 1)}  # format: (edit_dist, aln_c, aln_x, 1)
            for x_acc in alignments_of_x_to_c_transposed[m]:
                aln_x, aln_m, (matches, mismatches, indels) = alignments_of_x_to_c_transposed[m][x_acc]
                ed = mismatches + indels
                partition_alignments_c[x_acc] = (ed, aln_m, aln_x, 1) 
                # print(ed, aln_x, aln_m)

            alignment_matrix_to_c, PFM_to_c = create_position_probability_matrix(m, partition_alignments_c)
            c_alignment = alignment_matrix_to_c[m]
            is_consensus = True
            for j in range(len(PFM_to_c)):
                c_v =  c_alignment[j]
                candidate_count = PFM_to_c[j][c_v]
                for v in PFM_to_c[j]:
                    if v != c_v and candidate_count <= PFM_to_c[j][v]: # needs to have at least one more in support than the second best as we have added c itself to the multialignment
                        # print("not consensus at:", j)
                        is_consensus = False

            if not is_consensus:
                print("Read with support {0} were not consensus".format(str(support)))
                del alignments_of_x_to_c_transposed[m]
                del C[m]
        else:
            print("deleting:")
            del alignments_of_x_to_c_transposed[m]
            del C[m]


    # we now have an accession of minimizer, change to this accession insetad of storing sequence
    alignments_of_x_to_m_filtered = transpose(alignments_of_x_to_c_transposed)
    for x_acc in list(alignments_of_x_to_m_filtered.keys()):
        for m in list(alignments_of_x_to_m_filtered[x_acc].keys()):
            m_acc = m_to_acc[m]
            aln_x, aln_m, (matches, mismatches, indels) = alignments_of_x_to_m_filtered[x_acc][m]
            del alignments_of_x_to_m_filtered[x_acc][m]
            ed =  mismatches + indels
            alignments_of_x_to_m_filtered[x_acc][m_acc] = (ed, aln_x, aln_m)

    return alignments_of_x_to_m_filtered, m_to_acc


def filter_C_X_and_partition(X, C, alignments_of_reads_to_candidates, partition):
    # if there are x in X that is not in best_exact_matches, these x had no (decent) alignment to any of
    # the candidates, simply skip them.

    print("total reads:", len(X), "total reads with an alignment:", len(alignments_of_reads_to_candidates) )
    remaining_to_align = {}
    for x in list(X.keys()):
        if x not in alignments_of_reads_to_candidates:
            print("read missing alignment",x)
            remaining_to_align[x] = X[x]
            # del X[x]

    # also, filter out the candidates that did not get any alignments here.

    print("total consensus:", len(C), "total consensus with at least one alignment:", len(partition) )

    for c_acc in list(C.keys()):
        if c_acc not in partition:
            print("candidate missing hit:", c_acc)
            del C[c_acc]
        elif len(partition[c_acc]) == 0:
            print("candidate had edges in minimizer graph but they were all shared and it did not get chosen:", c_acc)
            del C[c_acc]            
            del partition[c_acc]
        else:
            print(c_acc, " read hits:", len(partition[c_acc]))
    return remaining_to_align, C

def vizualize_test_graph(C_seq_to_acc, partition_of_X, partition_of_C):
    import networkx as nx
    import matplotlib.pyplot as plt
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
    fig_file = os.path.join(params.plotfolder, "Graph_bip_1000_step_" + str(step) + ".png")
    plt.savefig(fig_file, format="PNG")
    plt.clf()

def get_minimizer_graph(candidate_transcripts):
    best_edit_distances = {}
    isolated_nodes = set()
    for i, (c1, seq1) in enumerate(candidate_transcripts.items()):
        if i % 50 == 0:
            print("processing ", i)
        best_edit_distances[c1] = {}
        # best_cigars[c1] = {}
        best_ed = len(seq1)
        for c2, seq2 in  candidate_transcripts.items() :
            if c1 == c2:
                continue
            elif math.fabs(len(seq1) - len(seq2)) > best_ed:
                continue
            # TODO: remove task = "path" to speed up
            edit_distance, locations, cigar = edlib_traceback(seq1, seq2, mode="NW", task="path", k=min(15, best_ed))

            if 0 <= edit_distance < best_ed:
                best_ed = edit_distance
                best_edit_distances[c1] = {}
                best_edit_distances[c1][c2] = best_ed
                # best_cigars[c1] = {}
                # best_cigars[c1][c2] =  cigar
            elif edit_distance == best_ed:
                best_edit_distances[c1][c2] = best_ed
                # best_cigars[c1][c2] =  cigar

        if len(best_edit_distances[c1]) == 0: # all isolated nodes in this graph
            isolated_nodes.add(c1)

 
    minimizer_graph = transpose(best_edit_distances)
    for c_isolated in isolated_nodes:
        minimizer_graph[c_isolated] = {}
    # seen_in_test = set()
    # for m in minimizer_graph:
        # print(best_edit_distances[c1])
        # print("NEW", m, "size:", len(minimizer_graph[m]))
        # if len(best_edit_distances[c1]) > 1:
        #     print("Minimizer to multiple candidates:", m, len(minimizer_graph[m]))
        # for c in minimizer_graph[m]:
        #     ed = minimizer_graph[m][c]
            # if c in seen_in_test:
            #     print("Seen:", c)
            # seen_in_test.add(c)
            # print("ED:",ed, c )
                # print("Multiple best:", c1, len(c2), best_cigars[c1][c2])
                # print(c2)
                # print()
        # print()

    return minimizer_graph

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


def stat_filter_candidates(read_file, candidate_file, alignments_of_x_to_c, params):
    modified = True

    ############ GET READS AND CANDIDATES #################
    X = {acc: seq for (acc, seq) in  fasta_parser.read_fasta(open(read_file, 'r'))} 
    C = {acc: seq for (acc, seq) in  fasta_parser.read_fasta(open(candidate_file, 'r'))}

    ################################################################

    step = 1
    nr_of_tests = 0
    previous_round_tests = {}
    previous_candidate_p_values = {}
    realignment_to_avoid_local_max = 0
    to_realign = {}
    remaining_to_align_read_file = os.path.join(params.outfolder, "remaining_to_align.fa")


    while modified:
        modified = False
        print("STEP NR: {0}".format(step))

        ########### Write current candidates to file ##########
        temp_candidate_name = os.path.join(params.outfolder, "temp_candidates_step_{0}.fa".format(step))
        temp_candidate_file = open(temp_candidate_name, "w")

        for c_acc, c_seq in C.items():
            temp_candidate_file.write(">{0}\n{1}\n".format(c_acc, c_seq))
        temp_candidate_file.close()
        #######################################################

        # create partition
        partition_of_X = { c_acc : set() for c_acc in C.keys()}
        for x_acc in alignments_of_x_to_c:
            for c_acc in alignments_of_x_to_c[x_acc]:
                partition_of_X[c_acc].add(x_acc)

        ############ GET READ SUPORT AND ALIGNMENTS #################

        if realignment_to_avoid_local_max == 1:
            print("REALIGNING EVERYTHING FINAL STEP")
            to_realign = X      
        
        if to_realign:
            write_output.print_reads(remaining_to_align_read_file, to_realign)
            # align reads that is not yet assigned to candidate here
            G_star_rem, partition_of_remaining_X, remaining_alignments_of_x_to_c = partition_strings_2set(to_realign, C, remaining_to_align_read_file, temp_candidate_file.name)
            # add reads to best candidate given new alignments
            for c_acc in partition_of_remaining_X:
                partition_of_X[c_acc].update(partition_of_remaining_X[c_acc])

            # add the alignments to alignment structure
            for x_acc in remaining_alignments_of_x_to_c.keys():
                alignments_of_x_to_c[x_acc] = remaining_alignments_of_x_to_c[x_acc]
                # for c_acc in remaining_alignments_of_x_to_c[x_acc].keys():
                #     alignments_of_x_to_c[x_acc][c_acc] = remaining_alignments_of_x_to_c[x_acc][c_acc]

        C_seq_to_acc = {seq : acc for acc, seq in C.items()}
        ################################################################

        ############# GET THE CLOSES HIGHEST SUPPORTED REFERENCE TO TEST AGAINST FOR EACH CANDIDATE ############
        minimizer_graph = get_minimizer_graph(C)
        #####################################################################################################

        # get all candidats that serve as null-hypothesis references and have neighbors subject to testing
        # these are all candidates that are minimizers to some other, isolated nodes are not tested
        # candidatate in G_star_C
        nr_of_tests_this_round = len(minimizer_graph)
        print("NUMBER OF CANDIDATES LEFT:", len(C))
        significance_values = statistical_test_v2.do_statistical_tests(minimizer_graph, C, X, partition_of_X, single_core = params.single_core )

        ###################################### PARALLELIZING THIS SECTION ####################################
        ######################################################################################################
        ######################################################################################################

        # significance_values = {}
        # actual_tests = 0
        # # parallelize over outer for loop with function:  significance_values = do_statistical_tests(minimizer_graph, C, X, partition_of_X )
        # for t_acc in minimizer_graph:
        #     t_seq = C[t_acc]
        #     if len(minimizer_graph[t_acc]) == 0:
        #         significance_values[t_acc] = ("not_tested", len(partition_of_X[t_acc]), len(partition_of_X[t_acc]) )
        #         continue

        #     reads = set([x_acc for c_acc in minimizer_graph[t_acc] for x_acc in partition_of_X[c_acc]] )
        #     reads.update(partition_of_X[t_acc])

        #     #### Bug FIX ############
        #     total_read_in_partition = len(partition_of_X[t_acc])
        #     reads_in_partition = set(partition_of_X[t_acc])
        #     # print(partition_of_X[t_acc])
        #     for c_acc in minimizer_graph[t_acc]:
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

        #     print("N_t:", N_t, "reads in partition:", total_read_in_partition, "ref:", t_acc )
        #     print("Nr candidates:", len(minimizer_graph[t_acc]), minimizer_graph[t_acc])
        #     assert total_read_in_partition == N_t # each read should be uiquely assinged to a candidate
        #     reads_and_candidates = reads.union( [c_acc for c_acc in minimizer_graph[t_acc]]) 
        #     reads_and_candidates_and_ref = reads_and_candidates.union( [t_acc] ) 

        #     # get multialignment matrix here
        #     alignment_matrix_to_t, PFM_to_t =  arrange_alignments(t_acc, reads_and_candidates_and_ref, X, C )

        #     # get parameter estimates for statistical test
        #     candidate_accessions = set( [ c_acc for c_acc in minimizer_graph[t_acc]] )
        #     delta_t = functions.get_difference_coordinates_for_candidates(t_acc, candidate_accessions, alignment_matrix_to_t) # format: { c_acc1 : {pos:(state, char), pos2:(state, char) } , c_acc2 : {pos:(state, char), pos2:(state, char) },... }
        #     epsilon, lambda_S, lambda_D, lambda_I = functions.get_error_rates_and_lambda(t_acc, len(t_seq), candidate_accessions, alignment_matrix_to_t) 
        #     # get number of reads k supporting the given set of variants, they have to support all the variants within a candidate
        #     candidate_support = functions.get_supporting_reads_for_candidates(t_acc, candidate_accessions, alignment_matrix_to_t, delta_t, partition_of_X) # format: { c_acc1 : [x_acc1, x_acc2,.....], c_acc2 : [x_acc1, x_acc2,.....] ,... }
        #     # read_indiv_invariant_factors = adjust_probability_of_read_to_alignment_invariant(delta_t, alignment_matrix_to_t, t_acc)
        #     candidate_indiv_invariant_factors = functions.adjust_probability_of_candidate_to_alignment_invariant(delta_t, alignment_matrix_to_t, t_acc)

        #     for c_acc, c_seq in list(minimizer_graph[t_acc].items()):
        #         original_mapped_to_c = len(partition_of_X[c_acc])
        #         k = len(candidate_support[c_acc])
        #         # print("supprot:", k, "diff:", len(delta_t[c_acc]))
        #         corrected_p_value = stat_test(k, len(t_seq), epsilon, delta_t, candidate_indiv_invariant_factors, t_acc, c_acc, original_mapped_to_c)
        #         actual_tests += 1
        #         if c_acc in significance_values:
        #             if corrected_p_value > significance_values[c_acc][0]:
        #                 significance_values[c_acc] = (corrected_p_value, k, N_t)
        #         else:
        #             significance_values[c_acc] = (corrected_p_value, k, N_t)

        # print("actual tests performed this round:", actual_tests)

        ######################################################################################################
        ######################################################################################################
        ######################################################################################################

        for c_acc, (corrected_p_value, k, N_t) in list(significance_values.items()):
            if corrected_p_value == "not_tested":
                pass
            elif corrected_p_value > 0.01/nr_of_tests_this_round:
                print("removing", c_acc, "p-val:", corrected_p_value, "k", k, "N_t", N_t )
                del C[c_acc] 
                modified = True
                for x_acc in partition_of_X[c_acc]:
                    to_realign[x_acc] = X[x_acc]
                    del alignments_of_x_to_c[x_acc]

                del partition_of_X[c_acc]

        print("nr candidates left:", len(C))
        candidate_file = os.path.join(params.outfolder, "candidates_after_step_{0}.fa".format(step))
        step += 1
        print("LEN SIGN:", len(significance_values), len(C))
        write_output.print_candidates(candidate_file, alignments_of_x_to_c, C, significance_values)


        # do a last realingment to avoind local maxima of reads

        if realignment_to_avoid_local_max == 1: # we have already done a last realignment, keep going until everythin is significant never realign
            realignment_to_avoid_local_max == 2
        elif not modified and realignment_to_avoid_local_max == 0: # we have not yet done a final alignment and everythin is significant, realign to escape local maxima alignment
            realignment_to_avoid_local_max = 1
            modified = True

    final_out_file_name =  os.path.join(params.outfolder, "final_candidates.fa")
    write_output.print_candidates(final_out_file_name, alignments_of_x_to_c, C, significance_values, final = True)

    return C
