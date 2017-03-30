"""
    PFM is a list of dicts, where the inner dicts are one per column in the PFM matrix, its keys by characters A, C, G, T, -
    alignment_matrix is a representation of all alignments in a partition. this is a dictionary where sequences s_i belonging to the 
    partition as keys and the alignment of s_i with respectt to the alignment matix.
"""
import os
import unittest

import copy

from modules.functions import transpose,create_position_probability_matrix
from modules.partitions import partition_strings_paths
from modules.SW_alignment_module import sw_align_sequences, sw_align_sequences_keeping_accession
from modules.edlib_alignment_module import edlib_align_sequences, edlib_align_sequences_keeping_accession
from modules.input_output import fasta_parser, write_output
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

    lenghts = [len(seq) for seq in S.values()]
    print(sorted(lenghts))
    max_len = max(lenghts)
    min_len = min(lenghts)
    print("Max transcript length:{0}, Min transcript length:{1}".format(max_len, min_len))

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

        #################################################
        ###### temp check for isoform collapse###########
        import re
        pattern = r"[-]{8,}"
        big_edit_distances = [ partition_alignments[s1][s2] for s1 in partition_alignments for s2 in partition_alignments[s1] if re.search(pattern, partition_alignments[s1][s2][1]) or  re.search(pattern, partition_alignments[s1][s2][2]) ] 
        # for tup in big_edit_distances:
        #     print("{0}\n{1}\n{2}\n\n".format(tup[0], tup[1], tup[2]))
        cccntr = 0
        print("Barcodes:", params.barcodes )
        out_file = open(os.path.join(params.outfolder, "exon_difs.fa"), "w")
        if params.barcodes:
            for s1, s1_dict in list(partition_alignments.items()): 
                for s2, alignment_tuple in list(s1_dict.items()):
                    if re.search(pattern, alignment_tuple[1][20: -20]) or  re.search(pattern, alignment_tuple[2][20: -20]): # [20: -20] --> ignore this if-statement if missing or truncated barcode
                        del partition_alignments[s1][s2]
                        print("Deleted:", len(s2)," minimizer length:", len(s1), "length alignment:", len(alignment_tuple[2]), "edit distance:", alignment_tuple[0])
                        print(s2)
                        cccntr += 1
                        out_file.write(">{0}\n{1}\n".format(unique_seq_to_acc[s2],s2))
        else:
            for s1, s1_dict in list(partition_alignments.items()): 
                for s2, alignment_tuple in list(s1_dict.items()):
                    if re.search(pattern, alignment_tuple[1]) or  re.search(pattern, alignment_tuple[2]):
                        del partition_alignments[s1][s2]
                        cccntr += 1        
                        out_file.write(">{0}\n{1}\n".format(unique_seq_to_acc[s2],s2))

        print("Number containing exon difference and removed in this pass:", cccntr)
        # sys.exit()
        ########################################################


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
        
        # #require support from at least 4 reads if not tested (consensus transcript had no close neighbors)
        # # add extra constraint that the candidate has to have majority on _each_ position in c here otherwise most likely error
        # if support >= params.min_candidate_support:
        #     # print("needs to be consensus over each base pair")
        #     partition_alignments_c = {m : (0, m, m, 1)}  # format: (edit_dist, aln_c, aln_x, 1)
        #     for x_acc in alignments_of_x_to_c_transposed[m]:
        #         aln_x, aln_m, (matches, mismatches, indels) = alignments_of_x_to_c_transposed[m][x_acc]
        #         ed = mismatches + indels
        #         partition_alignments_c[x_acc] = (ed, aln_m, aln_x, 1) 
        #         # print(ed, aln_x, aln_m)

        #     alignment_matrix_to_c, PFM_to_c = create_position_probability_matrix(m, partition_alignments_c)
        #     c_alignment = alignment_matrix_to_c[m]
        #     is_consensus = True
        #     for j in range(len(PFM_to_c)):
        #         c_v =  c_alignment[j]
        #         candidate_count = PFM_to_c[j][c_v]
        #         for v in PFM_to_c[j]:
        #             if v != c_v and candidate_count <= PFM_to_c[j][v]: # needs to have at least one more in support than the second best as we have added c itself to the multialignment
        #                 # print("not consensus at:", j)
        #                 is_consensus = False

        #     if not is_consensus:
        #         print("Read with support {0} were not consensus".format(str(support)))
        #         del alignments_of_x_to_c_transposed[m]
        #         del C[m]
        # else:
        #     print("deleting:")
        #     del alignments_of_x_to_c_transposed[m]
        #     del C[m]


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





if __name__ == '__main__':
    unittest.main()