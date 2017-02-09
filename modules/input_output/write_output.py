import os
from time import time
import datetime

from modules.functions import create_position_probability_matrix, transpose

def logger(message, logfile, timestamp=True):
    if timestamp:
        currrent_time = datetime.datetime.now()
        logfile.write(str(currrent_time) + "\t" + message + "\n")
    else:
        logfile.write(message + "\n")


def print_candidates(out_file_name, alignments_of_x_to_c, C, C_pvals, final = False):
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
                            # print("not consensus at:", j)
                            is_consensus = False

                if is_consensus:
                    if final:
                        out_file.write(">{0}\n{1}\n".format(c_acc + "_" + str(support) + "_" + str(p_value) + "_" + str(N_t) , seq))
                    else:
                        out_file.write(">{0}\n{1}\n".format(c_acc, seq))

                    final_candidate_count += 1
                else:
                    print("were not consensus")
            else:
                if final:
                    out_file.write(">{0}\n{1}\n".format(c_acc + "_" + str(support) + "_" + str(p_value) + "_" + str(N_t) , seq))
                else:
                    out_file.write(">{0}\n{1}\n".format(c_acc, seq))

                final_candidate_count += 1
        else:
            print("deleting:", "support:", support, "pval:", p_value, "tot reads in partition:", N_t  )
    print("Final candidate count: ", final_candidate_count)
    out_file.close()


def print_candidates_from_minimizers(out_file_candidates_name, out_file_alignments_name, alignments_of_x_to_m_filtered, M, m_to_acc, params):
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

# def print_alignments_from_reads_to_minimizers(out_file_alignments_name, alignments_of_x_to_m_filtered, params):
    out_file_alignments = open(out_file_alignments_name, "w")
    for x_acc in alignments_of_x_to_m_filtered.keys():
        for m_acc in alignments_of_x_to_m_filtered[x_acc].keys():
            ed, aln_x, aln_c = alignments_of_x_to_m_filtered[x_acc][m_acc]
            # m_acc = m_to_acc[m]
            out_file_alignments.write("{0}\t{1}\t{2}\t{3}\t{4}\n".format(x_acc, m_acc, ed, aln_x, aln_c) )

    out_file_alignments.close()

def print_reads(remaining_to_align_read_file, remaining_to_align, X):
    read_file_align = open(remaining_to_align_read_file, "w")

    for x_acc, seq in X.items():
        if x_acc in remaining_to_align:
            read_file_align.write(">{0}\n{1}\n".format(x_acc, seq) )


