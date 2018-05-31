"""
    PFM is a list of dicts, where the inner dicts are one per column in the PFM matrix, its keys by characters A, C, G, T, -
    alignment_matrix is a representation of all alignments in a partition. this is a dictionary where sequences s_i belonging to the 
    partition as keys and the alignment of s_i with respectt to the alignment matix.
"""
from __future__ import print_function
import os
import sys
import unittest
import copy
import math 
from collections import defaultdict
from itertools import combinations

import re
# from scipy.stats import poisson
from time import time
import pysam

# from modules.functions import transpose, create_position_probability_matrix
from modules import functions
from modules import partitions
from modules import graphs
from modules.SW_alignment_module import sw_align_sequences, sw_align_sequences_keeping_accession
from modules.edlib_alignment_module import edlib_align_sequences, edlib_align_sequences_keeping_accession, edlib_traceback
from modules.input_output import fasta_parser, fastq_parser, write_output
from modules import hypothesis_test_module
from modules import end_invariant_functions
from modules import ccs_info


# def vizualize_test_graph(C_seq_to_acc, read_partition, partition_of_C):
#     import networkx as nx
#     import matplotlib.pyplot as plt
#     D=nx.DiGraph()
#     lables = {}
#     for c1 in partition_of_C:
#         c1_acc = C_seq_to_acc[c1]
#         reads_to_c1 = [X[x_acc] for x_acc in  read_partition[c1_acc] ]
#         w_c1 = len(reads_to_c1)
#         for c2 in  partition_of_C[c1]:
#             c2_acc = C_seq_to_acc[c2]
#             reads_to_c2 = [X[x_acc] for x_acc in  read_partition[c2_acc] ]
#             w_c2 = len(reads_to_c2)
#             D.add_edge(c2,c1, weight = str(w_c2) + "->" + str(w_c1)  )
#     labels = nx.get_edge_attributes(D, 'weight')
#     # pos = nx.circular_layout(D)
#     pos = dict()
#     XX = [c2 for c1 in partition_of_C for c2 in partition_of_C[c1] ] # have in-edges
#     YY = [c1 for c1 in partition_of_C ] # not
#     pos.update( (n, (1, 4*i)) for i, n in enumerate(XX) ) # put nodes from X at x=1
#     pos.update( (n, (2, 2*i)) for i, n in enumerate(YY) ) # put nodes from Y at x=2
#     nx.draw_networkx_nodes(D, pos, node_size=50 )
#     nx.draw_networkx_edge_labels(D, pos, arrows=True, edge_labels=labels)
#     nx.draw_networkx_edges(D, pos, arrows=True, edge_labels=labels)
#     fig_file = os.path.join(params.plotfolder, "Graph_bip_1000_step_" + str(step) + ".png")
#     plt.savefig(fig_file, format="PNG")
#     plt.clf()




def get_nearest_neighbor_graph(candidate_transcripts):
    best_edit_distances = {}
    # no_ref_to_test_to = set()
    for i, (c1, seq1) in enumerate(candidate_transcripts.items()):
        if i % 50 == 0:
            print("processing ", i)
        best_edit_distances[c1] = {}
        # best_cigars[c1] = {}
        best_ed = len(seq1)
        for c2, seq2 in  candidate_transcripts.items() :
            # print("ALIGNING: {0} to {1}".format(c1, c2))
            if c1 == c2:
                continue
            elif math.fabs(len(seq1) - len(seq2)) > best_ed:
                # print("here", len(seq1), len(seq2), math.fabs(len(seq1) - len(seq2)), best_ed)
                continue
            # TODO: remove task = "path" to speed up
            edit_distance, locations, cigar = edlib_traceback(seq1, seq2, mode="NW", task="path", k=min(10, best_ed))
            # print(edit_distance, locations, cigar, "best:", best_ed)
            if 0 <= edit_distance < best_ed:
                best_ed = edit_distance
                best_edit_distances[c1] = {}
                best_edit_distances[c1][c2] = best_ed
                # best_cigars[c1] = {}
                # best_cigars[c1][c2] =  cigar
            elif edit_distance == best_ed:
                best_edit_distances[c1][c2] = best_ed
                # best_cigars[c1][c2] =  cigar

        # if len(best_edit_distances[c1]) == 0: # all isolated nodes in this graph
        #     no_ref_to_test_to.add(c1)

    # print(best_edit_distances["transcript_62_support_6"]) 
    # nearest_neighbor_graph = transpose(best_edit_distances)
    # print("isolated:", no_ref_to_test_to)
    # for c_isolated in no_ref_to_test_to:
    #     nearest_neighbor_graph[c_isolated] = {}

    # print(nearest_neighbor_graph["transcript_46_support_430"])
    # print(nearest_neighbor_graph["transcript_62_support_6"])

    # sys.exit()

    assert len(best_edit_distances) == len(candidate_transcripts)
    return best_edit_distances


# def check_exon_diffs(alignments_of_x_to_c, params):
#     #################################################################
#     ###### temp check for best alignment to wrong isoform ###########
#     import re
#     pattern = r"[-]{8,}"
#     cccntr = 0
#     print("Barcodes:", params.barcodes )
#     out_file = open(os.path.join(params.outfolder, "statistical_exon_difs.fa"), "w")
#     if params.barcodes:
#         for s1, s1_dict in list(alignments_of_x_to_c.items()): 
#             for s2, alignment_tuple in list(s1_dict.items()):
#                 if re.search(pattern, alignment_tuple[1][20: -20]) or  re.search(pattern, alignment_tuple[2][20: -20]): # [20: -20] --> ignore this if-statement if missing or truncated barcode
#                     # del alignments_of_x_to_c[s1][s2]
#                     # print("Deleted:", len(s1)," nearest_neighbor length:", len(s2), "length alignment:", len(alignment_tuple[2]), "edit distance:", alignment_tuple[0])
#                     # print(s2)
#                     cccntr += 1
#                     out_file.write(">{0}\n{1}\n".format(s1, alignment_tuple[1].replace("-", "") ))
#     else:
#         for s1, s1_dict in list(alignments_of_x_to_c.items()): 
#             for s2, alignment_tuple in list(s1_dict.items()):
#                 if re.search(pattern, alignment_tuple[1]) or  re.search(pattern, alignment_tuple[2]):
#                     # del alignments_of_x_to_c[s1][s2]
#                     cccntr += 1        
#                     out_file.write(">{0}\n{1}\n".format(s1, alignment_tuple[1].replace("-", "") ))

#     print("Number containing exon differences in this pass:", cccntr)
#     # sys.exit()
#     ###################################################################
#     ###################################################################


def product_with_check_overflow(p_value, mult_factor_inv):
    try:
        product = p_value*mult_factor_inv
    except OverflowError:
        print("OverflowError. p_val:{0}, mult_correction_factor:{1}".format(p_value, mult_factor_inv))
        product = 1.0
    return product

def stat_filter_candidates(read_file, candidate_file, read_partition, to_realign, params):
    modified = True

    ############ GET READS AND CANDIDATES #################
    if params.is_fastq:
        X_original = {acc: seq for (acc, seq, qual) in  fastq_parser.readfq(open(read_file, 'r'))}
    else:
        X_original = {acc: seq for (acc, seq) in  fasta_parser.read_fasta(open(read_file, 'r'))}
        # X_original = {acc: seq for (acc, seq) in  fasta_parser.read_fasta(open(read_file, 'r'))} 
    
    print("Total original reads", len(X_original))
    x_assigned_to_cluster = set([ x_acc for c_acc in read_partition for x_acc in read_partition[c_acc] ])
    X = {acc: seq for (acc, seq) in X_original.items() if acc in x_assigned_to_cluster or acc in to_realign }    # just set X to read_partition + to_realign here

    print("Original reads in fasta file:", len(X_original))
    print("Reads included in statistical testing:", len(X))
    if os.stat(candidate_file).st_size == 0:
        out_file_name = os.path.join(params.outfolder, "final_candidates.fa")
        tsv_info = os.path.join(params.outfolder, "cluster_info.tsv")
        write_output.print_candidates(out_file_name, {}, {}, {}, {}, params, final = True, reads_to_consensus_tsv = tsv_info)
        print("Candidate file is empty!")
        sys.exit(0)
    else:
        C = {acc: seq for (acc, seq) in  fasta_parser.read_fasta(open(candidate_file, 'r'))}

    ################################################################

    ### IF quality values are provided ####
    if params.is_fastq:
        ccs_dict_raw = {x_acc.split(" ")[0] : ccs_info.CCS(x_acc.split(" ")[0], seq, [ord(ascii_char) - 33 for ascii_char in qual], "NA") for (x_acc, seq, qual) in  fastq_parser.readfq(open(read_file, 'r'))}
        # int_quals = [ord(ascii_char) - 33 for ascii_char in qual] 
        X_ids = {  x_acc.split(" ")[0] : x_acc for x_acc in X} 
        print(len(X_ids), len(X), len(ccs_dict_raw))
        for x_acc in X:
            # print(ccs_dict_raw[x_acc.split(" ")[0]].qual)
            # print(ccs_dict_raw[x_acc.split(" ")[0]].seq)
            assert X_ids[x_acc.split(" ")[0]] == x_acc

        ccs_dict = ccs_info.modify_strings_and_acc_fastq(ccs_dict_raw, X_ids, X)
        for x_acc in X:
            assert X[x_acc] == ccs_dict[x_acc].seq

    elif params.ccs:
        ccs_file = pysam.AlignmentFile(params.ccs, "rb", check_sq=False)
        ccs_dict_raw = ccs_info.get_ccs(ccs_file)
        X_ids = { "/".join(x_acc.split("/")[:2]) : x_acc for x_acc in X} 
        ccs_dict = ccs_info.modify_strings_and_acc(ccs_dict_raw, X_ids, X)
        for x_acc in X:
            assert X[x_acc] == ccs_dict[x_acc].seq

    else:
        ccs_dict = {}

    ################################
    all_neighbors_graph = end_invariant_functions.get_NN_graph_ignored_ends_edlib(C, params)
    print("TOTAL EDGES G_ALL edlib:", len([1 for s in all_neighbors_graph for t in all_neighbors_graph[s]]))
    print("TOTAL Edit distances G_ALL edlib:", sum([all_neighbors_graph[s][t] for s in all_neighbors_graph for t in all_neighbors_graph[s]]))
    candidates_nn_graph_static = all_neighbors_graph

    # all_neighbors_graph_static = sw_align_sequences_keeping_accession(all_neighbors_graph, nr_cores = params.nr_cores)
    # no_alignments = set(C.keys()) - set(all_neighbors_graph_static.keys())
    # for c_acc in no_alignments:
    #     all_neighbors_graph_static[c_acc] = {}
    # # print("TOTAL EDGES G_STATIC parasail:", len([1 for s in all_neighbors_graph_static for t in all_neighbors_graph_static[s]]))
    # # print("TOTAL Edit distances G_STATIC parasail:", sum([ sum(all_neighbors_graph_static[s][t][2][1:]) for s in all_neighbors_graph_static for t in all_neighbors_graph_static[s]]))
    # candidates_nn_graph_static = {}
    # print(len(all_neighbors_graph_static), len(C), len(all_neighbors_graph))
    # for s in all_neighbors_graph_static:
    #     candidates_nn_graph_static[s] = {}
    #     for t in all_neighbors_graph_static[s]:
    #         s_aln, t_aln = all_neighbors_graph_static[s][t][0], all_neighbors_graph_static[s][t][1]
    #         mask_start, mask_end = functions.get_mask_start_and_end(s_aln, t_aln)
    #         ed = sum(all_neighbors_graph_static[s][t][2][1:]) -  min(mask_start, params.ignore_ends_len) - min(params.ignore_ends_len, (len(s_aln) - mask_end))
    #         # print(ed)
    #         if ed > 10:
    #             print(ed,"edlib:", all_neighbors_graph_static[s][t][2])
    #             print(s_aln)
    #             print(t_aln)
    #             continue
    #         else:
    #             candidates_nn_graph_static[s][t] = ed
    #         # print()
    # print("TOTAL EDGES G_STATIC parasail:", len([1 for s in candidates_nn_graph_static for t in candidates_nn_graph_static[s]]))
    # print("TOTAL Edit distances G_STATIC parasail after ignoring ends differences:", sum([ candidates_nn_graph_static[s][t] for s in candidates_nn_graph_static for t in candidates_nn_graph_static[s]]))
    print()
    # sys.exit()
    print()
    print("STARTING STATISTICAL TESTING")
    print()
    print("Number of reads to realign:", len(to_realign))
    step = 1
    prefilter = True
    previous_partition_of_X = copy.deepcopy(read_partition) 
    previous_components = { c_acc : set() for c_acc in C.keys()}
    previous_edges = { c_acc : set() for c_acc in C.keys()}
    significance_values = {}   
    realignment_to_avoid_local_max = 0
    remaining_to_align_read_file = os.path.join(params.outfolder, "remaining_to_align.fa")


    while modified:
        statistical_start = time() 

        modified = False
        print()
        print("STEP NR: {0}".format(step))
        print()
        ########### Write current candidates to file ##########
        temp_candidate_name = os.path.join(params.outfolder, "temp_candidates_step_{0}.fa".format(step))
        temp_candidate_file = open(temp_candidate_name, "w")

        for c_acc, c_seq in C.items():
            temp_candidate_file.write(">{0}\n{1}\n".format(c_acc, c_seq))
        temp_candidate_file.close()
        #######################################################

        if params.verbose:
            for c_acc in read_partition:
                print(c_acc, "has {0} reads assigned to it.".format(len(read_partition[c_acc])))

        ############ GET READ SUPORT AND ALIGNMENTS #################

        if realignment_to_avoid_local_max == 1:
            print("REALIGNING EVERYTHING FINAL STEP")
            to_realign = X      
            read_partition = { c_acc : {} for c_acc in C.keys()}



        if to_realign:
            print(len(to_realign), "reads to realign.")
            write_output.print_reads(remaining_to_align_read_file, to_realign)
            # align reads that is not yet assigned to candidate here
            G_star_rem, partition_of_realigned_reads = partitions.partition_strings_2set(to_realign, C, remaining_to_align_read_file, temp_candidate_file.name, params)
            reassigned_reads_to_candidates = {}
            for c_acc in partition_of_realigned_reads:
                reassigned_reads_to_candidates[c_acc] = {}
                for read_acc in partition_of_realigned_reads[c_acc]:
                    reassigned_reads_to_candidates[c_acc][read_acc] = (C[c_acc], X[read_acc]) 

            edit_distances_of_c_to_reads = edlib_align_sequences_keeping_accession(reassigned_reads_to_candidates, nr_cores = params.nr_cores)
            alignments_of_c_to_reads = sw_align_sequences_keeping_accession(edit_distances_of_c_to_reads, nr_cores = params.nr_cores)
            # structure: read_partition[c_acc][read_acc] = (c_aln, read_aln, (matches, mismatches, indels))

            ############## REMOVE EXON LEVEL DIFFERENCES IN ALIGNMENTS ####################
            ssw_temp = [ alignments_of_c_to_reads[c_acc][read_acc] for c_acc in alignments_of_c_to_reads for read_acc in alignments_of_c_to_reads[c_acc]  ] 
            _ = functions.filter_exon_differences(alignments_of_c_to_reads, params.min_exon_diff, params.ignore_ends_len)
            ssw_after_exon_temp = [ alignments_of_c_to_reads[c_acc][read_acc] for c_acc in alignments_of_c_to_reads for read_acc in alignments_of_c_to_reads[c_acc]  ] 
            print("Number of alignments that were removed before statistical test because best match to candidate had exon difference larger than {0}bp: {1} ".format(str(params.min_exon_diff) , len(ssw_temp) - len(ssw_after_exon_temp) ))
            #################################

            # add reads to best candidate given new alignments
            for c_acc in alignments_of_c_to_reads:
                for read_acc in alignments_of_c_to_reads[c_acc]:
                    read_partition[c_acc][read_acc] = alignments_of_c_to_reads[c_acc][read_acc]

            for c_acc in list(read_partition.keys()):
                if len(read_partition[c_acc]) == 0:
                    print(c_acc, "removed as it has no supporting reads")
                    del C[c_acc]
                    del read_partition[c_acc]
                else:
                    if params.verbose:
                        print(c_acc, "Now has {0} reads assigned to it, after aligning reads that are not assigned.".format(len(read_partition[c_acc])))

            # add the alignments to alignment structure
            # for x_acc in remaining_alignments_of_x_to_c.keys():
            #     alignments_of_x_to_c[x_acc] = remaining_alignments_of_x_to_c[x_acc]

        # C_seq_to_acc = {seq : acc for acc, seq in C.items()}
        ################################################################


        # check_exon_diffs(alignments_of_x_to_c, params)

        ############# GET THE CLOSES HIGHEST SUPPORTED REFERENCE TO TEST AGAINST FOR EACH CANDIDATE ############
        nearest_neighbor_graph = {}
        for c_acc in C.keys():
            nearest_neighbor_graph[c_acc] = {} 
            if len(candidates_nn_graph_static[c_acc]) > 0:
                candidate_edit_distances = [ed for c_nbr_acc, ed in candidates_nn_graph_static[c_acc].items() if c_nbr_acc in C]
                if candidate_edit_distances:
                    min_ed = min(candidate_edit_distances)
                    # print("new min:", min_ed)
                else:
                    print("no tests left")
                # here we get the relevant tests for the current iteration
                for c_nbr_acc in candidates_nn_graph_static[c_acc]:
                    if c_nbr_acc in C and candidates_nn_graph_static[c_acc][c_nbr_acc] == min_ed:
                        nearest_neighbor_graph[c_acc][c_nbr_acc] = min_ed

        print("Edges in NEW candidate NN graph:", len([ 1 for c_acc in nearest_neighbor_graph for t_acc in nearest_neighbor_graph[c_acc] ]) )
        print("Edit distances in NEW candidate NN graph:", sum([ nearest_neighbor_graph[c_acc][t_acc] for c_acc in nearest_neighbor_graph for t_acc in nearest_neighbor_graph[c_acc] ]) )

        # if params.ignore_ends_len > 0:
        #     nearest_neighbor_graph_old = end_invariant_functions.get_nearest_neighbors_graph_under_ignored_ends(C, params)
        # else:
        #     nearest_neighbor_graph_old = get_nearest_neighbor_graph(C)

        # for c_acc in nearest_neighbor_graph:
        #     for t_acc in nearest_neighbor_graph[c_acc]:
        #         if t_acc not in nearest_neighbor_graph_old[c_acc]:
        #             print("new test:", nearest_neighbor_graph[c_acc][t_acc])
        #             print(C[c_acc])
        #             print(C[t_acc])
        #             print(all_neighbors_graph_static[c_acc][t_acc])

        # print("Edges in candidate NN graph:", len([ 1 for c_acc in nearest_neighbor_graph_old for t_acc in nearest_neighbor_graph_old[c_acc] ]) )
        # print("Edit distances in candidate NN graph:", sum([ nearest_neighbor_graph_old[c_acc][t_acc] for c_acc in nearest_neighbor_graph_old for t_acc in nearest_neighbor_graph_old[c_acc] ]) )

        if realignment_to_avoid_local_max > 0:
            homopolymenr_invariant_graph = functions.get_homopolymer_invariants(C)
            print("Edges in candidate homopolymenr invariant graph:", len([ 1 for c_acc in homopolymenr_invariant_graph for t_acc in homopolymenr_invariant_graph[c_acc] ]) )
            for c_acc in homopolymenr_invariant_graph:
                if c_acc not in nearest_neighbor_graph:
                    print(c_acc, "not in NN_candidates graph but added now.")
                    nearest_neighbor_graph[c_acc] = {}
                for t_acc in homopolymenr_invariant_graph[c_acc]:
                    if t_acc not in nearest_neighbor_graph[c_acc]:
                        # print("Homopolymenr edge added")
                        nearest_neighbor_graph[c_acc][t_acc] = 1
            print("Total union of edges:", len([ 1 for c_acc in nearest_neighbor_graph for t_acc in nearest_neighbor_graph[c_acc] ]) ) 
        
        # print("EXTRA EDGES FROM HOMOPOLYMER IDENTICAL:", homopol_extra_added)


        ## save time if the nearest_neighbor and all cantidates in a component has identical reads assignmed to them as previous step
        # or heuristically: if candidate hase more than 2x more support than the reference itself (will give highly significant p-value anyway) to save computation time
        # Since indata is the same, the test is guaranteed to give same siginficance values as previous step        
        
        previous_significance_values = {}
        # print(nearest_neighbor_graph)
        for c_acc in list(nearest_neighbor_graph.keys()):
            # skip to test candidates with more reads than their respective references, because its redundant computation that will lead to significant values anyway..
            for t_acc in list(nearest_neighbor_graph[c_acc].keys()):
                if len(read_partition[c_acc]) >= params.min_test_ratio*len(read_partition[t_acc]):
                    if params.verbose:
                        print("skipping test for dominant candidate {0} to ref {1}".format(c_acc, t_acc))
                    del nearest_neighbor_graph[c_acc][t_acc]

            previous_significance_values[c_acc] = {}
            to_remove = set()
            for t_acc in list(nearest_neighbor_graph[c_acc].keys()):
                if (c_acc, t_acc) in previous_edges[c_acc] and ( previous_partition_of_X[t_acc] == read_partition[t_acc] ) and  (previous_partition_of_X[c_acc] == read_partition[c_acc]):
                    # print("here", (c_acc, t_acc) in previous_edges[c_acc] and ( previous_partition_of_X[t_acc] == read_partition[t_acc] ) and  (previous_partition_of_X[c_acc] == read_partition[c_acc]))
                    previous_significance_values[c_acc][t_acc] = significance_values[c_acc][t_acc]
                    to_remove.add((c_acc, t_acc))
                    if params.verbose:
                        print("TEST IDENTICAL TO PREVIOUS STEP, SKIPPING FOR", t_acc, c_acc)
                else: 
                    pass
                    # print("Modified")
            previous_edges[c_acc]  = set([(c_acc, t_acc) for t_acc in list(nearest_neighbor_graph[c_acc].keys())])
            for c_acc, t_acc in to_remove:
                del nearest_neighbor_graph[c_acc][t_acc]
        # print(nearest_neighbor_graph)
        print("Total edges after removing dominant candidates:", len([ 1 for c_acc in nearest_neighbor_graph for t_acc in nearest_neighbor_graph[c_acc] ]) ) 
        # sys.exit()
        #####################################################################################################


        # get all candidats that serve as null-hypothesis references and have neighbors subject to testing
        # these are all candidates that are nearest_neighbors to some other, isolated nodes are not tested
        # candidatate in G_star_C
        nr_of_tests_this_round = len([ 1 for c_acc in nearest_neighbor_graph for t_acc in nearest_neighbor_graph[c_acc] ] )
        print("NUMBER OF CANDIDATES LEFT:", len(C), ". Number statistical tests in this round:", nr_of_tests_this_round)
        if nr_of_tests_this_round > 0:
            new_significance_values = hypothesis_test_module.do_statistical_tests_per_edge(nearest_neighbor_graph, C, X, read_partition, ccs_dict, params )
            
            for c_acc in new_significance_values:
                for t_acc in new_significance_values[c_acc]:
                    previous_significance_values[c_acc][t_acc] = new_significance_values[c_acc][t_acc]

            # previous_significance_values.update(new_significance_values)
            significance_values = copy.deepcopy(previous_significance_values)
        else:
            significance_values = copy.deepcopy(previous_significance_values)

        assert len(significance_values) == len(C)
        highest_significance_values = {}
        for c_acc in significance_values:
            corrected_p_val_max = 0.0
            highest = (c_acc, "", "not_tested", 1.0, len(read_partition[c_acc]), len(read_partition[c_acc]), "") 
            for t_acc in significance_values[c_acc]:
                (p_value, mult_factor_inv, k, N_t, variants) = significance_values[c_acc][t_acc]
                corr_p_value = product_with_check_overflow(p_value, mult_factor_inv)
                if corr_p_value >= corrected_p_val_max:
                    corrected_p_val_max = corr_p_value
                    highest = (c_acc, t_acc, p_value, mult_factor_inv, k, N_t, variants)
            highest_significance_values[c_acc] =  highest

        if len(highest_significance_values) > 0:
            corrected_pvals = [product_with_check_overflow(p_value, mult_factor_inv)  for c_acc, (c_acc, t_acc, p_value, mult_factor_inv, k, N_t, variants) in highest_significance_values.items()  if p_value != "not_tested" ]
            if len(corrected_pvals) == 0:
                p_val_threshold = params.p_value_threshold #1.0
            else:
                corrected_pvals.sort()
                if len(corrected_pvals) % 2 == 0:
                    corrected_pvals_median = (corrected_pvals[int(len(corrected_pvals)/2)-1] + corrected_pvals[int(len(corrected_pvals)/2)]) / 2.0
                else:
                    corrected_pvals_median = corrected_pvals[int(len(corrected_pvals)/2)]
                print("Median corrected p-val:", corrected_pvals_median)
                print("Number of unique candidates tested:",  len(corrected_pvals))
                p_val_threshold = corrected_pvals_median if corrected_pvals_median > params.p_value_threshold else params.p_value_threshold
                print("Filtering threshold (p_val*mult_correction_factor):",  p_val_threshold)

        to_realign = {}
        p_value_tsv_file = open(os.path.join(params.outfolder, "p_values_{0}.tsv".format(step)), "w")

        for c_acc, (c_acc, t_acc, p_value, mult_factor_inv, k, N_t, variants) in highest_significance_values.items():
            if p_value == "not_tested":
                if params.verbose:
                    print("Did not test", c_acc)

            elif k == 0:
                if params.verbose:
                    print("Support is 0 for", c_acc) 
                print("removing", c_acc, "p-val:", p_value, "correction factor:", mult_factor_inv, "k", k, "N_t", N_t, "variants:", variants, "SUPPORT IS 0." )
                del C[c_acc] 
                modified = True
                for x_acc in read_partition[c_acc]:
                    to_realign[x_acc] = X[x_acc]
                del read_partition[c_acc]                          

            elif product_with_check_overflow(p_value, mult_factor_inv) >= p_val_threshold:
                print("removing", c_acc, "p-val:", p_value, "correction factor:", mult_factor_inv, "k", k, "N_t", N_t, "variants:", variants )
                del C[c_acc] 
                modified = True
                for x_acc in read_partition[c_acc]:
                    to_realign[x_acc] = X[x_acc]
                del read_partition[c_acc]
            
            if p_value != "not_tested": 
                p_value_tsv_file.write("{0}\t{1}\n".format( c_acc + "_" + str(k) + "_" + str(1.0 if k == 0 else min(1.0, product_with_check_overflow(p_value, mult_factor_inv))) + "_" + str(N_t) + "_" + str(len(variants)), str(p_value)))
        p_value_tsv_file.close()


        previous_partition_of_X = copy.deepcopy(read_partition)

        print("nr candidates left:", len(C))
        candidate_file = os.path.join(params.outfolder, "candidates_after_step_{0}.fa".format(step))
        step += 1

        # significance_values = new_significance_values.copy()
        
        if len(C) == 0: # no candidates were significant!
            break

        # print("LEN SIGN:", len(significance_values), len(C))
        write_output.print_candidates(candidate_file, C, highest_significance_values, read_partition, X, params)

        # do a last realingment to avoind local maxima of reads

        if realignment_to_avoid_local_max == 1: # we have already done a last realignment, keep going until everythin is significant never realign
            realignment_to_avoid_local_max = 2
        elif not modified and realignment_to_avoid_local_max == 0: # we have not yet done a final alignment and everythin is significant, realign to escape local maxima alignment
            realignment_to_avoid_local_max = 1
            modified = True
            prefilter = False

        statistical_elapsed = time() - statistical_start
        write_output.logger('Time for Statistical test, step {0}:{1}'.format(step, str(statistical_elapsed)), params.logfile)
   

    if params.ignore_ends_len > 0:
        c_acc_to_support = {c_acc : len(all_candidate_assigned_reads) for c_acc, all_candidate_assigned_reads in read_partition.items()}
        remaining_c_after_invariant = end_invariant_functions.collapse_candidates_under_ends_invariant(C, c_acc_to_support, params)
        # print(remaining_c_after_invariant)
        # sys.exit()
        for c_acc in remaining_c_after_invariant:
            c_seq = C[ c_acc ] 
            for removed_c_acc in remaining_c_after_invariant[c_acc]:
                removed_c_seq = C[ removed_c_acc ]
                reads_to_removed_c_acc = read_partition[removed_c_acc]

                for read_acc in reads_to_removed_c_acc:
                    read_partition[c_acc][read_acc] = reads_to_removed_c_acc[read_acc]

                del C[ removed_c_acc ]
                del c_acc_to_support[ removed_c_acc ]
                del read_partition[ removed_c_acc ]


    final_out_file_name =  os.path.join(params.outfolder, "final_candidates.fa")
    tsv_info = os.path.join(params.outfolder, "cluster_info.tsv")
    write_output.print_candidates(final_out_file_name, C, highest_significance_values, read_partition, X, params, final = True, reads_to_consensus_tsv = tsv_info)

    return C
