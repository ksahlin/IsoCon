from __future__ import print_function
import os
from time import time
import datetime

from modules.SW_alignment_module import sw_align_sequences_keeping_accession
from modules.edlib_alignment_module import edlib_align_sequences_keeping_accession

def logger(message, logfile, timestamp=True):
    if timestamp:
        currrent_time = datetime.datetime.now()
        logfile.write(str(currrent_time) + "\t" + message + "\n")
    else:
        logfile.write(message + "\n")



def print_candidates(out_file_name, C, significance_test_values, partition_of_X, X, params, final = False, reads_to_consensus_tsv = "" ):
    out_file = open(out_file_name, "w")
    if final:
        reads_to_consensus_tsv_file = open(os.path.join(reads_to_consensus_tsv), "w")
        for c_acc in partition_of_X:
            for x_acc in partition_of_X[c_acc]:
                reads_to_consensus_tsv_file.write("{0}\t{1}\t{2}\t{3}\n".format(x_acc, c_acc, len(X[x_acc]), len(C[c_acc])))
        reads_to_consensus_tsv_file.close()

    final_candidate_count = 0
    for c_acc, seq in sorted(C.items()):
        # print(significance_test_values[c_acc])
        c_acc, t_acc, p_value, correction_factor, support, N_t, delta_size = significance_test_values[c_acc] 
        if params.verbose:
            print(c_acc, "Support:", support, "P-value:", p_value, "correction factor:", correction_factor, "delta size:", delta_size, "partition size:", N_t, "tested to", t_acc)

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


def print_candidates_from_nearest_neighbors(out_file_candidates_name, C, params):
    """
        alignments_of_x_to_c has format
        {x_acc : {y_acc : (x_alignment, y_alignment, (matches, mismatches, indels)) } }
    """
    out_file_candidates = open(out_file_candidates_name, "w")
    final_candidate_count = 0
    for c_acc, c_seq in sorted(C.items()):
        out_file_candidates.write(">{0}\n{1}\n".format(c_acc, c_seq))
        final_candidate_count += 1

    print("Final candidate count: ", final_candidate_count)
    out_file_candidates.close()


def print_reads(remaining_to_align_read_file, remaining_to_align):
    read_file_align = open(remaining_to_align_read_file, "w")
    for x_acc, seq in remaining_to_align.items():
        read_file_align.write(">{0}\n{1}\n".format(x_acc, seq) )


