"""
    PFM is a list of dicts, where the inner dicts are one per column in the PFM matrix, its keys by characters A, C, G, T, -
    alignment_matrix is a representation of all alignments in a partition. this is a dictionary where sequences s_i belonging to the 
    partition as keys and the alignment of s_i with respectt to the alignment matix.
"""
import copy
import unittest
import math
from partitions import partition_strings_paths, partition_strings_2set_paths
from functions import create_position_probability_matrix, transpose
from modules import graphs
from SW_alignment_module import sw_align_sequences, sw_align_sequences_keeping_accession
from input_output import fasta_parser


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
        input: a dictionary of reads acc : sequence
        output: a string containing a path to a fasta formatted file with consensus_id_support as accession 
                    and the sequence as the read
    """ 

    S = X
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

        out_file = open("/Users/kxs624/tmp/minimizer_test_1000_step" +  str(step) + ".fa", "w")
        for i, m in enumerate(partition_alignments):
            N_t = sum([container_tuple[3] for s, container_tuple in partition_alignments[m].items()])
            out_file.write(">{0}\n{1}\n".format("read" + str(i)+ "_support_" + str(N_t) , m))

        step += 1
        prev_edit_distances_2steps_ago = prev_edit_distances
        prev_edit_distances = edit_distances


    for m in M:
        N_t = sum([container_tuple[3] for s, container_tuple in partition_alignments[m].items()])
        C[m] = N_t    
    out_file = open("/Users/kxs624/tmp/minimizer_test_1000_converged.fa", "w")
    for i, m in enumerate(partition_alignments):
        out_file.write(">{0}\n{1}\n".format("read_" + str(i)+ "_support_" + str(C[m]) , m))   
    out_file.close()    
    return out_file

def get_partition_alignments_2set(graph_partition, C, X, G_star):
    partition_dict = {}
    for c, x_hits in graph_partition.items():
        partition_dict[c] = {}
        for x in x_hits:
            partition_dict[c][x] = (C[c], X[x]) 


    exact_alignments = sw_align_sequences_keeping_accession(partition_dict, single_core = False)
    partition_alignments = {} 

    for c in exact_alignments:
        partition_alignments[c] = {}
        for x in exact_alignments[c]:
            aln_c, aln_x, (matches, mismatches, indels) = exact_alignments[c][x]
            edit_dist = mismatches + indels
            # indegree =  1 if x not in G_star[C] else G_star[C][x]
            # if indegree > 1:
            #     print("Larger than 1!!", indegree)
            partition_alignments[c][x] = (edit_dist, aln_c, aln_x, 1)

    # for c in partition_alignments:
    #     if len(partition_alignments[c]) < 1: 
    #         del partition_alignments[c]
    print("NR candidates with at least one hit:", len(partition_alignments))

    return partition_alignments

def filter_C_and_X(X, C, G_star, partition):
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
        else:
            print(c, "hits:", len(partition[c]))


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
    G_star, partition =  partition_strings_2set_paths(X, C, X_file, C_file)
    filter_C_and_X(X, C, G_star, partition)
    #########################################################################################


    modified = set(partition.keys())
    # changed_nodes = set(C.keys())

    # while modified:
    #     # based on how G_star is formed, we decide which reads to align to align to the candidate under the null-hypothesis
    #     # i.e., all reads in partition[c1] + partition[c2] should be aligned to c1. After that we estimate error rates etc.
    #     G_star_C, alignment_graph, converged = graphs.construct_minimizer_graph(C)
    #     partition_alignments = get_partition_alignments_2set(graph_partition, C, X, G_star)       

    #     modified = set()
    #     # do not recalculate significance of an edge that has not changed,
    #     # i.e., neither c1 nor c2 has gotten new reads
    #     for c1 in G_star_C.keys():
    #         for c2 in G_star_C[c1].keys():
    #             if c1 not in changed_nodes and c2 not in changed_nodes:
    #                 continue
    #             N_c2_and_c2 = len(partition[c1]) + len(partition[c2])
    #             # Identify the \Delta positions and their cordinates in the alignment between c1 and c2 here w.r.t. the coordinates in the alignment matrix
    #             # These coordinates are differenet within c1 and c2 respectively, whe need to get both for easy access
    #             # Also identify their state here so that we use proper error rates

    #             S = 0 # the sum of reads supporting m errors
    #             # calculate the individual as well as total error rates in each read here for substitutions, insertions and deletions respecively

    #             for x_i in partition[c1]:
    #                 e_s, e_i,e_d = .....
    #                 p_i = get the probability that read i has the m = |\delta| errors here given epsilons.      
    #                 # Find the number k of reads (in partition[c1] + partition[c2]) that supports the |\Delta| variants in c1.         
    #                 Z_i = a binary value 1 if read i supports m errors
    #                 S += Z_i

    #             for x_i in partition[c2]:
    #                 e_s, e_i,e_d = .....
    #                 p_i = get the probability that read i has the m = |\delta| errors here given epsilons.      
    #                 # Find the number k of reads (in partition[c1] + partition[c2]) that supports the |\Delta| variants in c1.         
    #                 Z_i = a binary value 1 if read i supports m errors
    #                 S += Z_i

    #             # We send |reads| as the total read support of c2 under the null hypothesis as well as k to the statistical test here. 
    #             p_val = significance_test(k, N_c2_and_c2, lambd)                
    #             # rearrange the alignments of reads in partition c1 to align to the consensus in partition c2 here in a smark way..
    #             # alignment_matrix, PFM = create_position_probability_matrix(m, partition) needs to be modified somehow

    #             if p_val < 0.05:
    #                 del G_star_C[c1]
    #                 # update partition_alignments, partition, M here!
    #                 # update the individual as well as total error rates in each read here for substitutions, insertions and deletions respecively
    #                 # for all the reads that has been reassigned

    #                 print("Modified!", k, N_c2_and_c2, delta, N_c1, N_c2 )
    #                 break
    #                 modified.add(c2)
    #             # else:
    #             #     changed_nodes.remove()
    #     # what happens if a node c1 is removed that is a minimizer to another sequence that has not been processed in this given step? 
    #     # we should do nothing in this step and wait for the new graph C to be generated
    # return C

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
            read_file_name = "/Users/kxs624/Documents/data/pacbio/simulated/ISOseq_sim_n_200/simulated_pacbio_reads.fa"
            # read_file_name = "/Users/kxs624/Documents/data/pacbio/simulated/DAZ2_2_exponential_constant_0.001.fa"
            # read_file_name = "/Users/kxs624/Documents/data/pacbio/simulated/TSPY13P_2_constant_constant_0.0001.fa"
            # read_file_name = "/Users/kxs624/Documents/data/pacbio/simulated/TSPY13P_4_linear_exponential_0.05.fa"
            # read_file_name = "/Users/kxs624/Documents/data/pacbio/simulated/HSFY2_2_constant_constant_0.0001.fa"

            # X = {acc: seq for (acc, seq) in  fasta_parser.read_fasta(open(read_file_name, 'r'))} 
        except:
            print("test file not found:",read_file_name) 

        try:
            consensus_file_name = "/Users/kxs624/tmp/minimizer_test_200_converged.fa"
            # consensus_file_name = "/Users/kxs624/tmp/minimizer_test_1000_converged.fa"
            # consensus_file_name = "/Users/kxs624/tmp/minimizer_consensus_DAZ2_2_exponential_constant_0.001_step10.fa"
            # consensus_file_name = "/Users/kxs624/Documents/data/pacbio/simulated/TSPY13P_2_constant_constant_0.0001.fa"
            # consensus_file_name = "/Users/kxs624/Documents/data/pacbio/simulated/TSPY13P_4_linear_exponential_0.05.fa"
            # consensus_file_name = "/Users/kxs624/Documents/data/pacbio/simulated/HSFY2_2_constant_constant_0.0001.fa"

            # C = {acc: seq for (acc, seq) in  fasta_parser.read_fasta(open(consensus_file_name, 'r'))} 
        except:
            print("test file not found:",consensus_file_name) 

        three_CO(read_file_name, consensus_file_name)

if __name__ == '__main__':
    unittest.main()