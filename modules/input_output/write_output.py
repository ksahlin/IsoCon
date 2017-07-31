import os
from time import time
import datetime

from modules.SW_alignment_module import sw_align_sequences_keeping_accession
from modules.edlib_alignment_module import edlib_align_sequences_keeping_accession
from modules.functions import create_position_probability_matrix, transpose

def logger(message, logfile, timestamp=True):
    if timestamp:
        currrent_time = datetime.datetime.now()
        logfile.write(str(currrent_time) + "\t" + message + "\n")
    else:
        logfile.write(message + "\n")


def check_if_consensus(c_acc, C, X, partition_of_X):

    partition_dict = {c_acc : {c_acc : (C[c_acc], C[c_acc])}}
    for x_acc in partition_of_X[c_acc]:
        partition_dict[c_acc][x_acc] = (C[c_acc], X[x_acc])

    exact_edit_distances = edlib_align_sequences_keeping_accession(partition_dict, single_core = True)    
    exact_alignments = sw_align_sequences_keeping_accession(exact_edit_distances, single_core = True)
    partition_alignments = {} 

    for c_acc in exact_alignments:
        partition_alignments[c_acc] = {}
        for x_acc in exact_alignments[c_acc]:
            aln_c, aln_x, (matches, mismatches, indels) = exact_alignments[c_acc][x_acc]
            edit_dist = mismatches + indels
            partition_alignments[c_acc][x_acc] = (edit_dist, aln_c, aln_x, 1)

    alignment_matrix_to_c, PFM_to_c = create_position_probability_matrix(C[c_acc], partition_alignments[c_acc])
    c_alignment = alignment_matrix_to_c[c_acc]
    is_consensus = True
    not_cons_positions = []
    for j in range(len(PFM_to_c)):
        c_v =  c_alignment[j]
        candidate_count = PFM_to_c[j][c_v]
        for v in PFM_to_c[j]:
            if v != c_v and candidate_count <= PFM_to_c[j][v]: # needs to have at least one more in support than the second best as we have added c itself to the multialignment
                # print("not consensus at:", j, PFM_to_c[j])
                not_cons_positions.append((j, c_v, PFM_to_c[j]))
                is_consensus = False

    if not is_consensus:
        print("Were not consensus at:", not_cons_positions)
    else:
        print("Were consensus")

    return is_consensus

def print_candidates(out_file_name, C, significance_test_values, partition_of_X, X, final = False, reads_to_consensus_tsv = "" ):
    out_file = open(out_file_name, "w")
    if final:
        reads_to_consensus_tsv_file = open(os.path.join(reads_to_consensus_tsv), "w")
        for c_acc in partition_of_X:
            for x_acc in partition_of_X[c_acc]:
                reads_to_consensus_tsv_file.write("{0}\t{1}\t{2}\t{3}\n".format(x_acc, c_acc, len(X[x_acc]), len(C[c_acc])))
        reads_to_consensus_tsv_file.close()

    final_candidate_count = 0
    for c_acc, seq in C.items():
        p_value, correction_factor, support, N_t, delta_size = significance_test_values[c_acc] 
        
        print(c_acc, "Support:", support, "P-value:", p_value, "correction factor:", correction_factor, "delta size:", delta_size)

        # add extra constraint that if the candidate cannot be statistically tested in a meaningful way (e.g., too divergent from other sequences)
        # the candidate has to have majority on _each_ position in c here otherwise most likely error
        

        # if p_value == "not_tested":
        #     is_consensus = check_if_consensus(c_acc, C, X, partition_of_X)
        #     if is_consensus:
        #         if final:
        #             if support > 2: # need at least 3 reads for meaningful consensus
        #                 out_file.write(">{0}\n{1}\n".format(c_acc + "_" + str(support) + "_" + str(p_value) + "_" + str(N_t) + "_" + str(delta_size) , seq))
        #                 final_candidate_count += 1
        #             else:
        #                 print("consensus had support: {0} and were not reported to output despite being consensus.".format(support) )
        #         else:
        #             out_file.write(">{0}\n{1}\n".format(c_acc, seq))
        #             final_candidate_count += 1

        # else:
        if final:
            out_file.write(">{0}\n{1}\n".format(c_acc + "_" + str(support) + "_" + str(p_value) + "_" + str(N_t) + "_" + str(delta_size) , seq))
            final_candidate_count += 1
        else:
            out_file.write(">{0}\n{1}\n".format(c_acc, seq))
            final_candidate_count += 1
        # else:
        #     print("Not printing candidate to file:", c_acc, "support:", support, "pval:", p_value, "tot reads in partition:", N_t  )
    print("Candidates written to file: ", final_candidate_count)
    out_file.close()


def print_candidates_from_minimizers(out_file_candidates_name, M, m_to_acc, params):
    """
        alignments_of_x_to_c has format
        {x_acc : {y_acc : (x_alignment, y_alignment, (matches, mismatches, indels)) } }
    """
    out_file_candidates = open(out_file_candidates_name, "w")
    final_candidate_count = 0
    # m_to_acc = {}
    # alignments_of_x_to_c_transposed = transpose(alignments_of_x_to_c)      
    for i, (m, support) in enumerate(M.items()):
        # m_acc = "read_" + str(i) + "_support_" + str(support)
        m_acc = m_to_acc[m] 
        out_file_candidates.write(">{0}\n{1}\n".format(m_acc, m))
        final_candidate_count += 1

    print("Final candidate count: ", final_candidate_count)
    out_file_candidates.close()

# # def print_alignments_from_reads_to_minimizers(out_file_alignments_name, alignments_of_x_to_m_filtered, params):
#     out_file_alignments = open(out_file_alignments_name, "w")
#     for x_acc in alignments_of_x_to_m_filtered.keys():
#         for m_acc in alignments_of_x_to_m_filtered[x_acc].keys():
#             ed, aln_x, aln_c = alignments_of_x_to_m_filtered[x_acc][m_acc]
#             # m_acc = m_to_acc[m]
#             out_file_alignments.write("{0}\t{1}\t{2}\t{3}\t{4}\n".format(x_acc, m_acc, ed, aln_x, aln_c) )

#     out_file_alignments.close()

def print_reads(remaining_to_align_read_file, remaining_to_align):
    read_file_align = open(remaining_to_align_read_file, "w")
    for x_acc, seq in remaining_to_align.items():
        read_file_align.write(">{0}\n{1}\n".format(x_acc, seq) )


