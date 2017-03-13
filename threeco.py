"""
    PFM is a list of dicts, where the inner dicts are one per column in the PFM matrix, its keys by characters A, C, G, T, -
    alignment_matrix is a representation of all alignments in a partition. this is a dictionary where sequences s_i belonging to the 
    partition as keys and the alignment of s_i with respectt to the alignment matix.
"""
import os
import unittest

import copy

from modules.functions import transpose,create_position_probability_matrix
from modules.partitions import partition_strings_paths, partition_strings_2set, partition_to_statistical_test
from modules import graphs
from modules.SW_alignment_module import sw_align_sequences, sw_align_sequences_keeping_accession
from modules.edlib_alignment_module import edlib_align_sequences, edlib_align_sequences_keeping_accession
from modules.input_output import fasta_parser, write_output
from modules import statistical_test
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

def stat_filter_candidates(read_file, candidate_file, alignments_of_x_to_c, params):
    modified = True

    ############ GET READ SUPORT AND ALIGNMENTS #################
    X = {acc: seq for (acc, seq) in  fasta_parser.read_fasta(open(read_file, 'r'))} 
    ################################################################

    step = 1
    nr_of_tests = 0
    previous_round_tests = {}
    previous_candidate_p_values = {}
    realignment_to_avoid_local_max = 0

    while modified:
        modified = False
        print("NEW STEP")

        ############ GET READ SUPORT AND ALIGNMENTS #################

        C = {acc: seq for (acc, seq) in  fasta_parser.read_fasta(open(candidate_file, 'r'))}
        partition_of_X = {}

        if realignment_to_avoid_local_max != 1:
            alignments_of_x_to_c_transposed = transpose(alignments_of_x_to_c)
            # remove the alignments of candidates that didn't pass the consensus over each base pair here
            for c_acc in list(alignments_of_x_to_c_transposed.keys()):
                if c_acc not in C:
                    del alignments_of_x_to_c_transposed[c_acc]
            alignments_of_x_to_c = transpose(alignments_of_x_to_c_transposed)

            for c_acc in alignments_of_x_to_c_transposed.keys():
                partition_of_X[c_acc] = set()
                for x_acc in alignments_of_x_to_c_transposed[c_acc].keys():
                    partition_of_X[c_acc].add(x_acc)
            remaining_to_align, C = filter_C_X_and_partition(X, C, alignments_of_x_to_c, partition_of_X)
            remaining_to_align_read_file = os.path.join(params.outfolder, "remaining_to_align.fa")
            write_output.print_reads(remaining_to_align_read_file, remaining_to_align, X)

        else:
            print("REALIGNING EVERYTHING FINAL STEP")
            remaining_to_align = X
            remaining_to_align_read_file = read_file
            for c_acc in C:
                partition_of_X[c_acc] = set()

        # align reads that is not yet assigned to candidate here
        G_star_rem, partition_of_remaining_X, remaining_alignments_of_x_to_c = partition_strings_2set(remaining_to_align, C, remaining_to_align_read_file, candidate_file)

        # add reads to best candidate given new alignments
        for c_acc in partition_of_remaining_X:
            partition_of_X[c_acc].update(partition_of_remaining_X[c_acc])
            # print("adding", partition_of_remaining_X[c_acc], "to", c_acc)
        # add the alignments to alignment structure
        for x_acc in remaining_alignments_of_x_to_c.keys():
            for c_acc in remaining_alignments_of_x_to_c[x_acc].keys():
                alignments_of_x_to_c[x_acc][c_acc] = remaining_alignments_of_x_to_c[x_acc][c_acc]

        C_seq_to_acc = {seq : acc for acc, seq in C.items()}
        ################################################################

        ############# GET THE CLOSES HIGHEST SUPPORTED REFERENCE TO TEST AGAINST FOR EACH CANDIDATE ############
        weights = { C[c_acc] : len(x_hits) for c_acc, x_hits in partition_of_X.items()} 
        G_star_C, partition_of_C, M, converged = partition_to_statistical_test(C, params, node_weights = weights, edge_creating_min_treshold = params.statistical_test_editdist,  edge_creating_max_treshold = 15)
        # self edges not allowed
        print(len(C), len(partition_of_C), len(M) )
        #####################################################################################################

        # get all candidats that serve as null-hypothesis references and have neighbors subject to testing
        # these are all candidates that are minimizers to some other, isolated nodes are not tested
        # candidatate in G_star_C
        null_hypothesis_references_to_candidates = set([c for c in partition_of_C if partition_of_C[c] ])
        print("References in testing:", len(null_hypothesis_references_to_candidates))
        nr_of_tests_this_round = len([1 for c1 in partition_of_C for c2 in  partition_of_C[c1]])
        if nr_of_tests < nr_of_tests_this_round:
            nr_of_tests = nr_of_tests_this_round

        ######### just for vizualization ###############
        if params.develop_mode:
            vizualize_test_graph(C_seq_to_acc, partition_of_X, partition_of_C)
        ####################################################


        #all reads tested against a reference are aligned to that reference 
        # do one reference candidate at a time, this is therefore easily  parallellized 

        C_pvals = { C_seq_to_acc[c] : (len(partition_of_X[C_seq_to_acc[c]]), "not_tested", len(partition_of_X[C_seq_to_acc[c]]) ) for c in partition_of_C if not partition_of_C[c]} # initialize with transcripts not tested
        candidate_p_values = statistical_test.do_statistical_tests(null_hypothesis_references_to_candidates, C_seq_to_acc, partition_of_X, partition_of_C, X, C, previous_round_tests, previous_candidate_p_values, single_core = params.single_core)

        ##########################################################
        ##########################################################

        # store the performed tests to avoid identical testing in next iteration.
        previous_round_tests = {}
        for t in null_hypothesis_references_to_candidates:
            t_acc = C_seq_to_acc[t]
            N_t = len(partition_of_X[t_acc])
            for c in partition_of_C[t]:
                c_acc = C_seq_to_acc[c]
                N_t += len(partition_of_X[c_acc])
            previous_round_tests[t] = (set([c for c in partition_of_C[t]]), N_t)
            # print("PREVIOUS TEST: candidates {0}, N_t: {1}".format(previous_round_tests[t][0], previous_round_tests[t][1]))
        previous_candidate_p_values = copy.deepcopy(candidate_p_values)

        ##########################################################
        ##########################################################

        p_vals = []
        # wait for all candidate p_values to be calculated
        removed_nodes_reference_graph = {}

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

            C_pvals[c_acc] = (k, p_value, N_t)

            p_vals.append(p_value)


        for c_acc in removed_nodes_reference_graph:
            modified = True
            t_acc, k = removed_nodes_reference_graph[c_acc]
            while t_acc in removed_nodes_reference_graph:
                t_acc, k = removed_nodes_reference_graph[t_acc]

            del C[c_acc]
            # partition_of_X[t_acc].update(partition_of_X[c_acc])
            for x_acc in partition_of_X[c_acc]:
                if x_acc not in alignments_of_x_to_c:
                    print(x_acc, "not in alignments_of_x_to_c but in partition_X[c_acc]. len partition_X[c_acc]: {0}".format(len(partition_of_X[c_acc])))
                else:
                    del alignments_of_x_to_c[x_acc]

            del partition_of_X[c_acc]

            print("merging:", c_acc, "into", t_acc, "k:", k)


        print("nr candidates left:", len(C))
        print(p_vals)
        candidate_file = os.path.join(params.outfolder, "candidates_after_step_{0}.fa".format(step))
        write_output.print_candidates(candidate_file, alignments_of_x_to_c, C, C_pvals)

        if params.develop_mode:
            plt.hist(p_vals)
            fig_file = os.path.join(params.plotfolder, "p_vals_step_" + str(step) +".png")
            plt.savefig(fig_file, format="PNG")
            plt.clf()
            step += 1

        # do a last realingment to avoind local maxima of reads

        if realignment_to_avoid_local_max == 1: # we have already done a last realignment, keep going until everythin is significant never realign
            realignment_to_avoid_local_max == 2
        elif not modified and realignment_to_avoid_local_max == 0: # we have not yet done a final alignment and everythin is significant, realign to escape local maxima alignment
            realignment_to_avoid_local_max = 1
            modified = True

    final_out_file_name =  os.path.join(params.outfolder, "final_candidates.fa")
    write_output.print_candidates(final_out_file_name, alignments_of_x_to_c, C, C_pvals, final = True)
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