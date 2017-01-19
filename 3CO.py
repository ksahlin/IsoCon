"""
    PPM is a list of dicts, where the inner dicts are one per column in the PPM matrix, its keys by characters A, C, G, T, -
    alignment_matrix is a representation of all alignments in a partition. this is a dictionary where sequences s_i belonging to the 
    partition as keys and the alignment of s_i with respectt to the alignment matix.
"""

import unittest
import math
from graphs import partition_strings

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

def find_candidate_transcripts(X):
    S = X
    # print(len(S))
    unique_seq_to_acc = get_unique_seq_accessions(S)
    partition_alignments, partition, M, converged = partition_strings(S)

    if converged:
        return M

    while not converged:
        for m, partition in partition_alignments.items():
            if len(partition) > 1:
                # all strings has not converged
                # alignment_matrix, PPM = create_position_probability_matrix(partition)                
                for s in partition:
                    nr_pos_to_correct = int(math.ceil(partition_alignments[m][s][0] / 2.0)) #decide how many errors we should correct here
                    print(nr_pos_to_correct)
                    # s_alignment_in_matrix = alignment_matrix[s]
                    # find the position probabilities of the alignment of s in PPM
                    pos_probs_for_s = []
                    for j in range(len(PPM)):
                        pos_probs_for_s.append( (j, PPM[j][s_alignment_in_matrix[j]]) )

                    pos_probs_for_s.sort(key=lambda x: x[1]) # sort with respect to smalles probabilities                    
                    J = [j for j, prob in pos_probs_for_s[:nr_pos_to_correct]] # J is the set of the nr_pos_to_correct smallest position probabilities
                    s_new = alignment_matrix[s]
                    for j in J:
                        highest_prob_character_at_j = max(PPM[j], key=lambda k: PPM[j][k])
                        s_new[j] = highest_prob_character_at_j
                    s_modified = "".join([nucl for nucl in s_new if nucl != "-" ])

                    # only unique strings can change in this step

                    accession_of_s = seq_to_acc[s] # this is still unique
                    S[accession_of_s] = s_modified

        unique_seq_to_acc = get_unique_seq_accessions(S)
        partition_alignments, partition, M, converged = partition_strings(S)

    return M

class TestFunctions(unittest.TestCase):

    def test_find_candidate_transcripts(self):
        self.maxDiff = None

        from input_output import fasta_parser
        try:
            fasta_file_name = "/Users/kxs624/Documents/data/pacbio/simulated/ISOseq_sim_n_200/simulated_pacbio_reads.fa"
            # fasta_file_name = "/Users/kxs624/Documents/data/pacbio/simulated/DAZ2_2_exponential_constant_0.001.fa"
            # fasta_file_name = "/Users/kxs624/Documents/data/pacbio/simulated/TSPY13P_2_constant_constant_0.0001.fa"
            # fasta_file_name = "/Users/kxs624/Documents/data/pacbio/simulated/TSPY13P_4_linear_exponential_0.05.fa"
            S = {acc: seq for (acc, seq) in  fasta_parser.read_fasta(open(fasta_file_name, 'r'))} 
        except:
            print("test file not found:",fasta_file_name)  
        find_candidate_transcripts(S)

if __name__ == '__main__':
    unittest.main()