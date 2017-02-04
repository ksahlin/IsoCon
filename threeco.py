"""
    PFM is a list of dicts, where the inner dicts are one per column in the PFM matrix, its keys by characters A, C, G, T, -
    alignment_matrix is a representation of all alignments in a partition. this is a dictionary where sequences s_i belonging to the 
    partition as keys and the alignment of s_i with respectt to the alignment matix.
"""
import os
import copy
import unittest
import math


from modules.partitions import partition_strings_paths, partition_strings_2set_paths, partition_strings
from modules.functions import create_position_probability_matrix, transpose, get_error_rates_and_lambda, get_difference_coordinates_for_candidates, get_supporting_reads_for_candidates
from modules import graphs
from modules.SW_alignment_module import sw_align_sequences, sw_align_sequences_keeping_accession
from modules.input_output import fasta_parser
from modules import statistical_test

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

    out_file_name = os.path.join(params.outfolder, "candidates_converged.fa")
    out_file = open(out_file_name, "w")
    for i, m in enumerate(partition_alignments):
        if C[m] >= params.min_candidate_support:
            out_file.write(">{0}\n{1}\n".format("read_" + str(i)+ "_support_" + str(C[m]) , m))   
    out_file.close()    
    return out_file.name



def filter_C_X_and_partition(X, C, G_star, partition):
    # if there are x in X that is not in best_exact_matches, these x had no (decent) alignment to any of
    # the candidates, simply skip them.

    print("total reads:", len(X), "total reads with an alignment:", len(G_star) )

    bug_read = "m151210_031012_42146_c100926392550000001823199905121697_s1_p0/73682/0_1489_CCS"
    print("IN alignment:", bug_read in G_star)
    bug_cand = ""
    if bug_read in G_star:
        bug_cand = G_star[bug_read]
        print("to c_acc:", G_star[bug_read])


    for x in X.keys():

        if x not in G_star:
            print("read missing alignment",x)
            del X[x]

    # also, filter out the candidates that did not get any alignments here.

    print("total consensus:", len(C), "total consensus with at least one alignment:", len(partition) )

    for c_acc in C.keys():

        if c_acc == bug_cand:
            print("BUGGY TRANSCRIPT HERE:", len(partition[bug_cand]))
            if bug_read in partition[bug_cand]:
                print("in partition" )
            else:
                print("Read is not in partition" )


        if c_acc not in partition:
            print("candidate missing hit:", c_acc)
            del C[c_acc]
        elif len(partition[c_acc]) == 0:
            print("candidate had edges in minimizer graph but they were all shared and it did not get chosen:", c_acc)
            del C[c_acc]            
            del partition[c_acc]
        else:
            print(c_acc, " read hits:", len(partition[c_acc]))
    return X, C


def stat_filter_candidates(read_file, candidate_file, params):
    ################################### PROCESS INDATA #####################################
    # if not candidate_file:
    #     candidate_file = find_candidate_transcripts(read_file, params)
    X_file = read_file
    C_file = candidate_file
    X = {acc: seq for (acc, seq) in  fasta_parser.read_fasta(open(read_file, 'r'))} 
    C = {acc: seq for (acc, seq) in  fasta_parser.read_fasta(open(candidate_file, 'r'))} 

    # C = {c: support for c, support in C.items() if support > 3}
    # TODO: eventually filter candidates with lower support than 2-3? Here?
    G_star, partition_of_X, alignments_of_x_to_c =  partition_strings_2set_paths(X, C, X_file, C_file)

    X, C = filter_C_X_and_partition(X, C, G_star, partition_of_X)
    C_seq_to_acc = {seq : acc for acc, seq in C.items()}

    #########################################################################################


    modified = True
    # changed_nodes = set(C.keys())
    step = 1
    nr_of_tests = 0
    while modified:
        modified = False
        print("NEW STEP")
        weights = { C[c_acc] : len(x_hits) for c_acc, x_hits in partition_of_X.items()} 
        G_star_C, partition_of_C, M, converged = partition_strings(C, params, node_weights = weights, edge_creating_min_treshold = params.statistical_test_editdist,  edge_creating_max_treshold = 15)
        # self edges not allowed
        print(len(C), len(partition_of_C), len(M) )


        # get all candidats that serve as null-hypothesis references and have neighbors subject to testing
        # these are all candidates that are minimizers to some other, isolated nodes is not tested
        # candidatate in G_star_C
        null_hypothesis_references_to_candidates = [c for c in partition_of_C if partition_of_C[c] ]
        print("References in testing:", len(null_hypothesis_references_to_candidates))
        nr_of_tests_this_round = len([1 for c1 in partition_of_C for c2 in  partition_of_C[c1]])
        if nr_of_tests < nr_of_tests_this_round:
            nr_of_tests = nr_of_tests_this_round

        ######### just for vizualization ###############
        # D=nx.DiGraph(G_star_C)
        if params.develop_mode:
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

        ####################################################

        # based on how G_star is shaped, we decide which reads to align to the candidate under the null-hypothesis
        # i.e., all reads in partition_of_X[c1] + partition_of_X[c2] should be aligned to c2. After that we estimate error rates etc.
        #or should we align al reads in S* to a cluster center c2?? that is all reads belonging to canditades that has the same minimizer?
        # we would probably get better statistical power in each test.. maybe its also simpler to implement. 



        # do one reference candidate at a time, these are all modular and this loop 
        # can easily be parallellized if we break this up to a function
        # p_values_to_t = wrapper_statistical_test()

        C_pvals = { C_seq_to_acc[c] : (len(partition_of_X[C_seq_to_acc[c]]), "not_tested", len(partition_of_X[C_seq_to_acc[c]]) ) for c in partition_of_C if not partition_of_C[c]} # initialize with transcripts not tested
        candidate_p_values = statistical_test.do_statistical_tests(null_hypothesis_references_to_candidates, C_seq_to_acc, partition_of_X, partition_of_C, X, C, single_core = params.single_core)

        p_vals = []
        # wait for all candidate p_values to be calculated
        removed_nodes_reference_graph = {}
        print()
        print()
        print("NEW CLUSTER")
        print()
        print()
        for c_acc, (t_acc, k, p_value, N_t) in candidate_p_values.items():
            if p_value > 0.05/nr_of_tests or k == 0:
                print("Potential delete:", c_acc, k, p_value, N_t)
                if t_acc in removed_nodes_reference_graph: # t_acc is also subject to removal
                    if removed_nodes_reference_graph[t_acc][0] == c_acc:
                        print("preventing cycle!")
                        if removed_nodes_reference_graph[t_acc][1] > k:
                            del removed_nodes_reference_graph[t_acc]
                            removed_nodes_reference_graph[c_acc] = (t_acc, k)
                    else:
                        # next_t_acc, next_k  = removed_nodes_reference_graph[t_acc]
                        # removed_nodes_reference_graph[c_acc] =  (next_t_acc, next_k)
                        removed_nodes_reference_graph[c_acc] =  (t_acc, k)
                else:
                    removed_nodes_reference_graph[c_acc] = (t_acc, k)  
                # deleted.add(c_acc)
                # print(c_acc, t_acc,t_acc_transfer_reads_to )
            else:
                C_pvals[c_acc] = (k, p_value, N_t)

            p_vals.append(p_value)


        for c_acc in removed_nodes_reference_graph:
            modified = True
            t_acc, k = removed_nodes_reference_graph[c_acc]
            while t_acc in removed_nodes_reference_graph:
                t_acc, k = removed_nodes_reference_graph[t_acc]

            del C[c_acc]
            partition_of_X[t_acc].update(partition_of_X[c_acc])
            del partition_of_X[c_acc]
            print("merging:", c_acc, "into", t_acc, "k:", k)


        print("nr candidates left:", len(C))
        print(p_vals)
        if params.develop_mode:
            plt.hist(p_vals)
            fig_file = os.path.join(params.plotfolder, "p_vals_step_" + str(step) +".png")
            plt.savefig(fig_file, format="PNG")
            plt.clf()
            step += 1

    out_file_name =  os.path.join(params.outfolder, "final_candidates.fa")
    out_file = open(out_file_name, "w")
    final_candidate_count = 0
    alignments_of_x_to_c_transposed = transpose(alignments_of_x_to_c)
    for c_acc, seq in C.items():
        support, p_value, N_t = C_pvals[c_acc] 
        #require support from at least 4 reads if not tested (consensus transcript had no close neighbors)
        # add extra constraint that the candidate has to have majority on _each_ position in c here otherwise most likely error
        if support >= 4:
            if p_value == "not_tested":
                print("not tested with support", support, "needs to be consensus over each base pair")
                
                partition_alignments_c = {c_acc : (0, C[c_acc], C[c_acc], 1)}  # format: (edit_dist, aln_c, aln_x, 1)
                for x_acc in alignments_of_x_to_c_transposed[c_acc]:
                    (ed, aln_x, aln_c) = alignments_of_x_to_c_transposed[c_acc][x_acc]
                    partition_alignments_c[x_acc] = (ed, aln_c, aln_x, 1) 

                alignment_matrix_to_c, PFM_to_c = create_position_probability_matrix(C[c_acc], partition_alignments_c)
                c_alignment = alignment_matrix_to_c[c_acc]
                is_consensus = True
                for j in range(len(PFM_to_c)):
                    c_v =  c_alignment[j]
                    candidate_count = PFM_to_c[j][c_v]
                    for v in PFM_to_c[j]:
                        if v != c_v and candidate_count <= PFM_to_c[j][v]: # needs to have at least one more in support than the second best as we have added c itself to the multialignment
                            print("not consensus at:", j)
                            is_consensus = False

                if is_consensus:
                    out_file.write(">{0}\n{1}\n".format(c_acc + "_" + str(support) + "_" + str(p_value) + "_" + str(N_t) , seq))
                    final_candidate_count += 1
                else:
                    print("were not consensus")
            else:
                out_file.write(">{0}\n{1}\n".format(c_acc + "_" + str(support) + "_" + str(p_value) + "_" + str(N_t) , seq))
                final_candidate_count += 1
        else:
            print("deleting:", "support:", support, "pval:", p_value, "tot reads in partition:", N_t  )
    print("Final candidate count: ", final_candidate_count)
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

    def test_stat_filter_candidates(self):
        self.maxDiff = None

        from input_output import fasta_parser
        try:
            read_file_name = "/Users/kxs624/Documents/data/pacbio/simulated/RBMY_44_-_constant_-.fa"
            # read_file_name = "/Users/kxs624/Documents/data/pacbio/simulated/ISOseq_sim_n_1000/simulated_pacbio_reads.fa"

            # read_file_name = "/Users/kxs624/Documents/data/pacbio/simulated/DAZ2_2_exponential_constant_0.001.fa"
            # read_file_name = "/Users/kxs624/Documents/data/pacbio/simulated/TSPY13P_2_constant_constant_0.0001.fa"
            # read_file_name = "/Users/kxs624/Documents/data/pacbio/simulated/TSPY13P_4_linear_exponential_0.05.fa"
            # read_file_name = "/Users/kxs624/Documents/data/pacbio/simulated/HSFY2_2_constant_constant_0.0001.fa"

            # X = {acc: seq for (acc, seq) in  fasta_parser.read_fasta(open(read_file_name, 'r'))} 
        except:
            print("test file not found:",read_file_name) 

        try:
            # consensus_file_name = "/Users/kxs624/tmp/minimizer_test_1000_converged.fa"
            # consensus_file_name = ""
            consensus_file_name = "/Users/kxs624/tmp/initial_candidates_RBMY_44_-_constant_-.fa"
            # consensus_file_name = "/Users/kxs624/tmp/final_candidates_RBMY_44_-_constant_-pass2.fa"

            # consensus_file_name = "/Users/kxs624/tmp/minimizer_consensus_final_RBMY_44_-_constant_-.fa"
            # consensus_file_name = "/Users/kxs624/tmp/minimizer_consensus_DAZ2_2_exponential_constant_0.001_step10.fa"
            # consensus_file_name = "/Users/kxs624/tmp/TSPY13P_2_constant_constant_0.0001_converged.fa"
            # consensus_file_name = "/Users/kxs624/tmp/final_candidates_TSPY13P_2_constant_constant_0.0001.fa"
            # consensus_file_name = "/Users/kxs624/Documents/data/pacbio/simulated/TSPY13P_4_linear_exponential_0.05.fa"
            # consensus_file_name = "/Users/kxs624/Documents/data/pacbio/simulated/HSFY2_2_constant_constant_0.0001.fa"

            # C = {acc: seq for (acc, seq) in  fasta_parser.read_fasta(open(consensus_file_name, 'r'))} 
        except:
            print("test file not found:",consensus_file_name) 

        stat_filter_candidates(read_file_name, consensus_file_name, params)

if __name__ == '__main__':
    unittest.main()