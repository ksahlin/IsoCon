"""
    PFM is a list of dicts, where the inner dicts are one per column in the PFM matrix, its keys by characters A, C, G, T, -
    alignment_matrix is a representation of all alignments in a partition. this is a dictionary where sequences s_i belonging to the 
    partition as keys and the alignment of s_i with respectt to the alignment matix.
"""
from __future__ import print_function
import os
import unittest
import copy
from time import time
import re
from collections import defaultdict
from collections import Counter

import pysam

from modules.functions import transpose,create_position_probability_matrix
from modules.partitions import highest_reachable_with_edge_degrees
from modules.SW_alignment_module import sw_align_sequences, sw_align_sequences_keeping_accession
from modules.edlib_alignment_module import edlib_align_sequences, edlib_align_sequences_keeping_accession, edlib_traceback
from modules.input_output import fasta_parser, write_output
from modules import correct_sequence_to_nearest_neighbor
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

def get_partition_alignments(graph_partition, M, G_star, params):
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

    pattern = r"[-]{{{min_exon_diff},}}".format( min_exon_diff = str(params.min_exon_diff)  )  # r"[-]{20,}"

    if params.verbose:
        print(pattern)

    for s1 in list(exact_alignments.keys()): 
        for s2 in list(exact_alignments[s1].keys()):
            s1_alignment, s2_alignment, (matches, mismatches, indels) = exact_alignments[s1][s2]
            missing_exon_s1 = re.search(pattern, s1_alignment)
            missing_exon_s2 = re.search(pattern, s2_alignment)
            if missing_exon_s1:
                # print(missing_exon_s1.group(0))
                # print(s1)
                # print(s2)
                # print(len(exact_alignments[s1].keys()))
                del exact_alignments[s1][s2]
            elif missing_exon_s2:
                # print(missing_exon_s2.group(0))
                # print(s1)
                # print(s2)
                # print(len(exact_alignments[s1].keys()))
                del exact_alignments[s1][s2]

    ssw_after_exon_temp = [ exact_alignments[s1][s2] for s1 in exact_alignments for s2 in exact_alignments[s1]  ] 
    print("Number of alignments that were removed before correction phase due to exon difference larger than {0}bp: {1} ".format(str(params.min_exon_diff) , len(ssw_temp) - len(ssw_after_exon_temp) ))
    # sys.exit()


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

    seq_to_acc = get_unique_seq_accessions(S)

    nearest_neighbor_start = time() 
    G_star, graph_partition, M, converged = highest_reachable_with_edge_degrees(S, params)
    partition_alignments = get_partition_alignments(graph_partition, M, G_star, params)       

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

        S_prime, S_prime_quality_vector = correct_sequence_to_nearest_neighbor.correct_strings(partition_alignments, seq_to_acc, ccs_dict, step, nr_cores = params.nr_cores, verbose = params.verbose)

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

        # partition_alignments, partition, M, converged = partition_strings(S)

        G_star, graph_partition, M, converged = highest_reachable_with_edge_degrees(S, params)
        partition_alignments = get_partition_alignments(graph_partition, M, G_star, params)  
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
   
    C = {}
    for m in M:
        N_t = partition_alignments[m][m][3] #sum([container_tuple[3] for s, container_tuple in partition_alignments[m].items()])
        if N_t > 1: # has converged to a consensus
            C[m] = N_t   


    original_reads = {acc: seq for (acc, seq) in  fasta_parser.read_fasta(open(read_file, 'r'))}
    original_reads_seq_to_accs =  defaultdict(list)
    for (acc, seq) in original_reads.items():
        original_reads_seq_to_accs[seq].append(acc)
    # original_reads_seq_to_acc = { seq : acc for (acc, seq) in  fasta_parser.read_fasta(open(read_file, 'r'))}
    reads_to_nearest_neighbors = {}
    # [for m, partition in partition_alignments.items() for s in partition]

    not_converged_reads = open(os.path.join(params.outfolder, "not_converged.fa"), "w")
    not_converged = set()
    for read_acc, seq in original_reads.items():
        corrected_s = S[read_acc]
        if corrected_s in C:
            if C[corrected_s] >= params.min_candidate_support:
                reads_to_nearest_neighbors[read_acc] = { corrected_s : (original_reads[read_acc], corrected_s)}
            else:
                if params.verbose:
                    print("nearest_neighbor did not pass threshold. It had support of {0} reads.".format(C[corrected_s]))
                    print(read_acc)
                del C[corrected_s]
        else:
            if corrected_s in original_reads_seq_to_accs:
                if params.verbose:
                    print("Read neither converged nor was it corrected (local pair r1 <---> r2 nearest_neighbor or a isolated alignment with exon difference filtered out before each correction)")
                    print(read_acc)
                not_converged_reads.write(">{0}_not_corrected_not_converged\n{1}\n".format(read_acc, seq))
                not_converged.add(read_acc)

            else: # partially corrected but not converged
                not_converged_reads.write(">{0}\n{1}\n".format(read_acc, seq))
                not_converged.add(read_acc)
                if params.verbose:
                    print("Read partially corrected but not converged")
                    print(read_acc)
                not_converged_reads.write(">{0}_corrected_but_not_converged_version\n{1}\n".format(read_acc, corrected_s))

    not_converged_reads.close()
    edit_distances_of_x_to_m = edlib_align_sequences_keeping_accession(reads_to_nearest_neighbors, nr_cores = params.nr_cores)
    alignments_of_x_to_m = sw_align_sequences_keeping_accession(edit_distances_of_x_to_m, nr_cores = params.nr_cores)


    if params.ignore_ends_len > 0:
        C_temp_accession_to_support = {}
        C_temp_accession_to_seq = {}
        for i, (seq, support) in enumerate(sorted(C.items(), key=lambda x: x[0])): # to assure consisent naming, sort lexocographically on candidates
            c_acc = "transcript_" + str(i) + "_support_" + str(support)
            C_temp_accession_to_support[c_acc] = support
            C_temp_accession_to_seq[c_acc] = seq

        remaining_c_after_invariant = end_invariant_functions.collapse_candidates_under_ends_invariant(C_temp_accession_to_seq, C_temp_accession_to_support, params)
        alignments_of_x_to_m_transposed = transpose(alignments_of_x_to_m)   
        for c_acc in remaining_c_after_invariant:
            c_seq = C_temp_accession_to_seq[ c_acc ] 
            c_support = C_temp_accession_to_support[ c_acc ]

            for removed_c_acc in remaining_c_after_invariant[c_acc]:
                c_removed_support = C_temp_accession_to_support[ removed_c_acc ]
                c_removed_seq = C_temp_accession_to_seq[ removed_c_acc ] 
                reads_to_removed_cand = alignments_of_x_to_m_transposed[ c_removed_seq ]
                alignments_of_x_to_m_transposed[c_seq].update(reads_to_removed_cand)
                
                C[c_seq] += len(reads_to_removed_cand)
                del alignments_of_x_to_m_transposed[c_removed_seq]
                del C[c_removed_seq]
        
        alignments_of_x_to_m = transpose(alignments_of_x_to_m_transposed)

        # alignments_of_x_to_m, C  = collapse_contained_sequences(alignments_of_x_to_m, C, params)

    alignments_of_x_to_m_filtered, m_to_acc, C_filtered, partition_of_X = filter_candidates(alignments_of_x_to_m, C, params)
        
    candidates_file_name = os.path.join(params.outfolder, "candidates_converged.fa")
    write_output.print_candidates_from_nearest_neighbors(candidates_file_name, C_filtered, m_to_acc, params)

    to_realign = {}
    for acc, seq in original_reads.items():
        if acc not in alignments_of_x_to_m_filtered:
            if acc not in not_converged:
                to_realign[acc] = seq

    return candidates_file_name, partition_of_X, to_realign 



# def collapse_contained_sequences(alignments_of_x_to_m, C, params):
#     print("Number nearest_neighbors before collapsing identical super strings of a string:", len(C))

#     alignments_of_x_to_m_transposed = transpose(alignments_of_x_to_m)   
#     m_to_acc = {}
#     C_sorted_strings = sorted(C, key=lambda m: len(m))
#     for i, m in enumerate(sorted(C, key=len)):
#         print("length:",len(m))
#         if i == len(C_sorted_strings) - 1:
#             print("last sequence, skipping!") 
#             break
#         if m in C:
#             # support = C[m]
#             all_perfect_super_strings = find_all_perfect_superstrings(m, i, C_sorted_strings, params.ignore_ends_len)
#             if len(all_perfect_super_strings) > 0:
#                 print("number of superstrings:", len(all_perfect_super_strings))
#             for ss in all_perfect_super_strings: # remove super string of m and move all reads supporting ss to m
#                 if ss in C: # has not already been collapsed from shorter perfect substring
#                     reads_to_ss = alignments_of_x_to_m_transposed[ss]
#                     alignments_of_x_to_m_transposed[m].update(reads_to_ss)
#                     C[m] += len(reads_to_ss)
#                     del alignments_of_x_to_m_transposed[ss]
#                     del C[ss]
#                 else:
#                     print("Already collapsed")
#                     assert ss not in alignments_of_x_to_m_transposed

#         else:
#             print("seq was a superstring") 

#     alignments_of_x_to_m = transpose(alignments_of_x_to_m_transposed)
#     print("Number nearest_neighbors after collapsing identical super strings of a string:", len(C))
#     return alignments_of_x_to_m, C


# def find_all_perfect_superstrings(m, i, C_sorted_strings, ignore_ends_len):
#     all_super_strings = set()
#     for j in range(i+1, len(C_sorted_strings)):
#         under_diff = True
#         ss_candidate = C_sorted_strings[j]
        
#         # Stop criterion here if ss_length is longer than m + 2*ignore_ends_len we can stop searching since strings are sorteb by length
#         if len(ss_candidate) - len(m) > 2*ignore_ends_len:
#             break

#         ed, locations, cigar = edlib_traceback(m, ss_candidate, mode="HW", task="locations", k=0)
#         for start, stop in locations:
#             if start > ignore_ends_len or len(ss_candidate) - stop - 1  > ignore_ends_len:
#                 print("perfect but larger diff than ignore_ends_len:", start, len(ss_candidate) - stop - 1 )
#                 under_diff = False
#                 break

#         if ed == 0 and under_diff: # candidate was perfect super string and not too large discrepancy
#             print("perfect ss, start and stop on ss: {0}".format(locations))
#             all_super_strings.add(ss_candidate)

#     return all_super_strings



def filter_candidates(alignments_of_x_to_c, C, params):
    alignments_of_x_to_c_transposed = transpose(alignments_of_x_to_c)   

    m_to_acc = {}

    for i, (m, support) in enumerate(list(sorted(C.items(), key=lambda x: x[0]))):
        m_acc = "transcript_" + str(i) + "_support_" + str(support)
        m_to_acc[m] = m_acc
        
        #require support from at least 4 reads if not tested (consensus transcript had no close neighbors)
        # add extra constraint that the candidate has to have majority on _each_ position in c here otherwise most likely error
        if support >= params.min_candidate_support:
            if params.prefilter_candidates:
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
                    max_v_j = max(PFM_to_c[j], key = lambda x: PFM_to_c[j][x] )
                    max_count = PFM_to_c[j][max_v_j]
                    if candidate_count < max_count:
                        print("not consensus at:", j)
                        is_consensus = False                    

                    # for v in PFM_to_c[j]:
                    #     if v != c_v and candidate_count <= PFM_to_c[j][v]: # needs to have at least one more in support than the second best as we have added c itself to the multialignment
                    #         print("not consensus at:", j)
                    #         is_consensus = False

                if not is_consensus:
                    print("Read with support {0} were not consensus".format(str(support)))
                    del alignments_of_x_to_c_transposed[m]
                    del C[m]
        else:
            print("deleting:")
            del alignments_of_x_to_c_transposed[m]
            del C[m]

    partition_of_X = { m_to_acc[c] : set(alignments_of_x_to_c_transposed[c].keys()) for c in  alignments_of_x_to_c_transposed}

    # we now have an accession of nearest_neighbor, change to this accession insetad of storing sequence
    alignments_of_x_to_m_filtered = transpose(alignments_of_x_to_c_transposed)
    for x_acc in list(alignments_of_x_to_m_filtered.keys()):
        for m in list(alignments_of_x_to_m_filtered[x_acc].keys()):
            m_acc = m_to_acc[m]
            aln_x, aln_m, (matches, mismatches, indels) = alignments_of_x_to_m_filtered[x_acc][m]
            del alignments_of_x_to_m_filtered[x_acc][m]
            ed =  mismatches + indels
            alignments_of_x_to_m_filtered[x_acc][m_acc] = (ed, aln_x, aln_m)

    return alignments_of_x_to_m_filtered, m_to_acc, C, partition_of_X





if __name__ == '__main__':
    unittest.main()