from __future__ import print_function
import os
import unittest
import copy
from time import time
import re
from collections import defaultdict
from collections import Counter

import pysam

from modules import functions
from modules import partitions
from modules.SW_alignment_module import sw_align_sequences, sw_align_sequences_keeping_accession
from modules.edlib_alignment_module import edlib_align_sequences, edlib_align_sequences_keeping_accession, edlib_traceback
from modules.input_output import fasta_parser, fastq_parser, write_output
from modules import correction_module
from modules import end_invariant_functions
from modules import ccs_info



def get_unique_seq_accessions(S):
    seq_to_acc = {}
    for acc, seq in  S.items():
        if seq in seq_to_acc:
            seq_to_acc[seq].append(acc)
        else: 
            seq_to_acc[seq] = []
            seq_to_acc[seq] = [acc]

    unique_seq_to_acc = {seq: acc_list[0] for seq, acc_list in  seq_to_acc.items() if len(acc_list) == 1 } 
    print("Non-converged (unique) sequences left:", len(unique_seq_to_acc))

    return seq_to_acc

def get_partition_alignments(graph_partition, M, G_star, exon_filtered, params):
    exact_edit_distances = edlib_align_sequences(graph_partition, nr_cores = params.nr_cores)    
        
    ed_temp = [ exact_edit_distances[s1][s2] for s1 in exact_edit_distances for s2 in exact_edit_distances[s1]  ] 
    ed_temp.sort()

    if params.verbose:
        print("ED from edlib:", ed_temp)
        print("number of ed calculated:", len(ed_temp))

    exact_alignments = sw_align_sequences(exact_edit_distances, nr_cores = params.nr_cores)

    ssw_temp = [ exact_alignments[s1][s2] for s1 in exact_alignments for s2 in exact_alignments[s1]  ] 
    # ssw_temp.sort()
    if params.verbose:
        print("Number of alignments returned from SSW:", len(ssw_temp))
        print("Number of alignments that were removed before correction phase -- too many mismatchas in ends (#ED-alignments - # SSW-alignments): {0} ".format(  len(ed_temp) - len(ssw_temp) ))


    filtered = functions.filter_exon_differences(exact_alignments, params.min_exon_diff, params.ignore_ends_len)
    exon_filtered.update(filtered)
    ssw_after_exon_temp = [ exact_alignments[s1][s2] for s1 in exact_alignments for s2 in exact_alignments[s1]  ] 
    print("Number of alignments that were removed before correction phase due to exon difference larger than {0}bp: {1} ".format(str(params.min_exon_diff) , len(ssw_temp) - len(ssw_after_exon_temp) ))

    partition_alignments = {} 
    for m in M:
        # selfdegree = 1 if m not in G_star[m] else G_star[m][m]
        selfdegree = G_star.node[m]["degree"]
        partition_alignments[m] = { m : (0, m, m, selfdegree) }
        if m not in exact_alignments:
            # print("nearest_neighbor did not have any anlignments, length:", M[m], "self-degree:", selfdegree)
            continue
        else:
            for s in exact_alignments[m]:
                # if s[10:len(s)-10] in m:
                #     print("LOOOOL", len(s), len(m))
                aln_m, aln_s, (matches, mismatches, indels) = exact_alignments[m][s]
                edit_dist = mismatches + indels
                # indegree =  1 if s not in G_star[m] else G_star[m][s]
                # if indegree > 1:
                #     print("Larger than 1!!", indegree)
                partition_alignments[m][s] = (edit_dist, aln_m, aln_s, 1)

    print("NR PARTITIONS :", len(partition_alignments))
    return partition_alignments



def find_candidate_transcripts(read_file, params):
    """
        input: a string pointing to a fasta file
        output: a string containing a path to a fasta formatted file with consensus_id_support as accession 
                    and the sequence as the read
    """ 
    if params.is_fastq:
        S = {acc: seq for (acc, seq, qual) in  fastq_parser.readfq(open(read_file, 'r'))}
        ccs_dict = {}
    else:
        S = {acc: seq for (acc, seq) in  fasta_parser.read_fasta(open(read_file, 'r'))}
    
        # still beta to use quaity values for correction, inactivated for now
        if False: #params.ccs:
            ccs_file = pysam.AlignmentFile(params.ccs, "rb", check_sq=False)
            ccs_dict_raw = ccs_info.get_ccs(ccs_file)
            X_ids = {  x_acc.split("/")[1] : x_acc for x_acc in S} 
            ccs_dict = ccs_info.modify_strings_and_acc(ccs_dict_raw, X_ids, S)
            for x_acc in S:
                assert S[x_acc] == ccs_dict[x_acc].seq
        else:
            ccs_dict = {}

    step = 1
    print()
    print("ITERATION:", step)
    print()

    lenghts = [len(seq) for seq in S.values()]
    C = Counter(lenghts)

    if params.verbose:
        for l in sorted(C.keys()):
            write_output.logger("seq length {0}: {1} occurances".format(l, C[l]), params.develop_logfile, timestamp=False)

    # print(sorted(lenghts))
    max_len = max(lenghts)
    min_len = min(lenghts)
    print("Max transcript length:{0}, Min transcript length:{1}".format(max_len, min_len))
    exon_filtered = set()
    
    seq_to_acc = get_unique_seq_accessions(S)

    nearest_neighbor_start = time() 
    G_star, graph_partition, M, converged = partitions.partition_strings(S, params)
    partition_alignments = get_partition_alignments(graph_partition, M, G_star, exon_filtered, params)       

    nearest_neighbor_elapsed = time() - nearest_neighbor_start
    write_output.logger('Time for nearest_neighbors and partition, step 1:{0}'.format(str(nearest_neighbor_elapsed)), params.logfile)


    prev_edit_distances_2steps_ago = [2**28,2**28,2**28] # prevents 2-cycles
    prev_edit_distances = [2**28]

    # homopolymer_mode = False

    while not converged:
        correction_start = time() 
        edit_distances = [ partition_alignments[s1][s2][0] for s1 in partition_alignments for s2 in partition_alignments[s1]  ] 
        edit_distances.sort()
        if params.verbose:
            print("edit distances from SSW:", edit_distances) 

        ###### Different convergence criterion #########

        if prev_edit_distances_2steps_ago == edit_distances:
            # Only cyclic alignments are left, these are reads that jump between two optimal alignment of two different
            # target sequeneces. This is a product of our fast heurustics of defining a minmap score + SSW filtering to choose best alignment
            print("CYCLE!!!")
            assert len(partition_alignments) == len(M)
            break
            # if homopolymer_mode:
            #     break
            # else:
            #     homopolymer_mode = True

        if sum(edit_distances) > sum(prev_edit_distances) and  max(edit_distances) > max(prev_edit_distances) :
            #return here if there is some sequence alternating between best alignments and gets corrected and re-corrected to different candidate sequences
            assert len(partition_alignments) == len(M)
            print("exiting here!")
            break
            # if homopolymer_mode:
            #     print("exiting here!")
            #     break            
            # else:
            #     homopolymer_mode = True


        has_converged = [True if ed == 0 else False for ed in edit_distances] 
        if all(has_converged):
            # we return here if tha data set contain isolated nodes.
            assert len(partition_alignments) == len(M)
            print("Normal convergence")
            break
            # if homopolymer_mode:
            #     print("Normal convergence")
            #     break
            # else: 
            #     homopolymer_mode = True
        #######################################################

        S_prime, S_prime_quality_vector = correction_module.correct_strings(partition_alignments, seq_to_acc, ccs_dict, step, nr_cores = params.nr_cores, verbose = params.verbose)

        for acc, s_prime in S_prime.items():
            S[acc] = s_prime
            if ccs_dict:
                ccs_dict[acc].qual = S_prime_quality_vector[acc]

        print("Tot seqs:", len(S))
        seq_to_acc = get_unique_seq_accessions(S)
        step += 1
        print()
        print("ITERATION:", step)
        print()
        print("Total number of unique reads put aside becauese of exon differences to all other strings:", len(exon_filtered))
        S_to_align = {acc: seq for acc, seq in S.items() if seq not in exon_filtered }
        G_star, graph_partition, M, converged = partitions.partition_strings(S_to_align, params)
        partition_alignments = get_partition_alignments(graph_partition, M, G_star, exon_filtered, params)  
        out_file_name = os.path.join(params.outfolder, "candidates_step_" +  str(step) + ".fa")
        out_file = open(out_file_name, "w")
        for i, m in enumerate(partition_alignments):
            N_t = sum([container_tuple[3] for s, container_tuple in partition_alignments[m].items()])
            out_file.write(">{0}\n{1}\n".format("read" + str(i)+ "_support_" + str(N_t) , m))


        prev_edit_distances_2steps_ago = prev_edit_distances
        prev_edit_distances = edit_distances

        correction_elapsed = time() - correction_start
        write_output.logger('Time for correction, nearest_neighbors and partition, step {0}:{1}'.format(step, str(correction_elapsed)), params.logfile)
    
        # sys.exit()
   
 
    ######################
    ###### NEW ###########
    c_seq_to_read_acc = {}
    for read_acc, seq in S.items():
        if seq in c_seq_to_read_acc:
            c_seq_to_read_acc[seq].append(read_acc)
        else:
            c_seq_to_read_acc[seq] = [read_acc]

    c_acc_to_seq = {}
    c_acc_to_support = {}            
    for i, m in enumerate(sorted(c_seq_to_read_acc)):
        if m in partition_alignments:
            N_t = partition_alignments[m][m][3] #sum([container_tuple[3] for s, container_tuple in partition_alignments[m].items()])
        else: 
            N_t = 1 # did not converge

        c_acc = "transcript_" + str(i) + "_support_" + str(N_t)
        c_acc_to_seq[c_acc] = m 
        c_acc_to_support[c_acc] = N_t

    if params.ignore_ends_len > 0:
        remaining_c_after_invariant = end_invariant_functions.collapse_candidates_under_ends_invariant(c_acc_to_seq, c_acc_to_support, params)
        # print(remaining_c_after_invariant)
        # sys.exit()
        for c_acc in remaining_c_after_invariant:
            c_seq = c_acc_to_seq[ c_acc ] 
            for removed_c_acc in remaining_c_after_invariant[c_acc]:
                removed_c_seq = c_acc_to_seq[ removed_c_acc ]
                reads_to_removed_c_acc = c_seq_to_read_acc[removed_c_seq]

                for read_acc in reads_to_removed_c_acc:
                    c_seq_to_read_acc[c_seq].append(read_acc)

                del c_acc_to_seq[ removed_c_acc ]
                del c_acc_to_support[ removed_c_acc ]
                del c_seq_to_read_acc[ removed_c_seq ]

    
    if params.is_fastq:
        original_reads = {acc: seq for (acc, seq, qual) in  fastq_parser.readfq(open(read_file, 'r'))}
    else:
        original_reads = {acc: seq for (acc, seq) in  fasta_parser.read_fasta(open(read_file, 'r'))}

    # original_reads = {acc: seq for (acc, seq) in  fasta_parser.read_fasta(open(read_file, 'r'))}
    print(len(S), len(original_reads))
    assert len(S) == len(original_reads)

    # to_realign = {}
    for c_acc in list(c_acc_to_seq.keys()):
        support = c_acc_to_support[c_acc]
        if support < params.min_candidate_support:
            c_seq = c_acc_to_seq[c_acc]
            if params.verbose:
                print("nearest_neighbor did not pass threshold. It had support of {0} reads.".format(support))
                print(c_seq_to_read_acc[c_seq])
            del c_acc_to_seq[c_acc]
            del c_seq_to_read_acc[c_seq]
            del c_acc_to_support[c_acc]
            # for read_acc in c_seq_to_read_acc[c_seq]:
            #     to_realign[read_acc] = original_reads[read_acc]

    all_reads_assigned_to_candidates = set([read_acc for c_seq in c_seq_to_read_acc for read_acc in c_seq_to_read_acc[c_seq]])
    unassigned_reads =  set(original_reads.keys()) - all_reads_assigned_to_candidates
    to_realign = {read_acc : original_reads[read_acc] for read_acc in unassigned_reads}
    print("Reads assigned to candididates:", len(all_reads_assigned_to_candidates))
    print("Reads to realign:", len(to_realign))
    print("Number of initial reads:", len(original_reads))

    candidates_file_name = os.path.join(params.outfolder, "candidates_converged.fa")
    write_output.print_candidates_from_nearest_neighbors(candidates_file_name, c_acc_to_seq, params)
    # sys.exit()
    not_converged_reads = open(os.path.join(params.outfolder, "not_converged.fa"), "w")
    not_converged_reads.close()
    assert len(to_realign) + len(all_reads_assigned_to_candidates) == len(original_reads)
    assert len(c_acc_to_seq) == len(c_seq_to_read_acc)
    c_to_reads = {}
    for c_acc, c_seq in c_acc_to_seq.items():
        c_to_reads[c_acc] = {}
        for read_acc in c_seq_to_read_acc[c_seq]:
            c_to_reads[c_acc][read_acc] = (c_seq, original_reads[read_acc])

    c_to_reads_edit_distances = edlib_align_sequences_keeping_accession(c_to_reads, nr_cores = params.nr_cores)
    print("Total reads in partition (assigned reads after edlib):", len([1 for c_acc in c_to_reads_edit_distances for read_acc in c_to_reads_edit_distances[c_acc] ]))
    read_partition = sw_align_sequences_keeping_accession(c_to_reads_edit_distances, nr_cores = params.nr_cores)
    filtered_reads = functions.filter_exon_differences(read_partition, params.min_exon_diff, params.ignore_ends_len)
    print("DEVELOP: Number of read to candidate assignments removed because of exon differences: ",len(filtered_reads))
    # sys.exit()
    for read_acc in filtered_reads:
        to_realign[read_acc] = original_reads[read_acc]


    print("Total reads in partition (assigned reads after SW):", len([1 for c_acc in read_partition for read_acc in read_partition[c_acc] ]))
    return candidates_file_name, read_partition, to_realign


    ##################################################################
    ##################################################################






if __name__ == '__main__':
    unittest.main()