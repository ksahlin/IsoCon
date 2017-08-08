"""
    PFM is a list of dicts, where the inner dicts are one per column in the PFM matrix, its keys by characters A, C, G, T, -
    alignment_matrix is a representation of all alignments in a partition. this is a dictionary where sequences s_i belonging to the 
    partition as keys and the alignment of s_i with respectt to the alignment matix.
"""
import os
import sys
import unittest
import copy
import math 
from collections import defaultdict
from itertools import combinations


from scipy.stats import poisson
from time import time


from modules.functions import transpose, create_position_probability_matrix
from modules import functions
from modules.partitions import partition_strings_2set
from modules import graphs
from modules.SW_alignment_module import sw_align_sequences, sw_align_sequences_keeping_accession
from modules.edlib_alignment_module import edlib_align_sequences, edlib_align_sequences_keeping_accession, edlib_traceback
from modules.input_output import fasta_parser, write_output
from modules import statistical_test_v2
from modules import correct_sequence_to_minimizer
from modules import end_invariant_functions

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

def transform(read):
    transformed_seq = []
    prev_nucl = ""
    for nucl in read:
        if nucl != prev_nucl:
            transformed_seq.append(nucl)
        prev_nucl = nucl

    return "".join(transformed_seq)

def get_homopolymer_invariants(candidate_transcripts):
    seq_to_acc = { seq : acc for (acc, seq) in  candidate_transcripts.items() }
    print("Unique before compression: ", len(seq_to_acc) )

    candidate_transcripts_transformed = {}
    clusters = defaultdict(list)
    for acc in candidate_transcripts:
        seq_transformed = transform(candidate_transcripts[acc])
        candidate_transcripts_transformed[acc] = seq_transformed
        clusters[seq_transformed].append(acc)

    seq_to_acc_transformed = { seq : acc for (acc, seq) in candidate_transcripts_transformed.items()}
    print("Unique after compression: ", len(seq_to_acc_transformed) )

    edges = {}
    for seq in clusters:
        if len(clusters[seq]) > 1:
            print(clusters[seq])
            for acc in clusters[seq]:
                edges[acc] = {}
            for acc1, acc2 in combinations(clusters[seq], 2):
                edges[acc1][acc2] = 1
                edges[acc2][acc1] = 1

                # result =  edlib.align( candidate_transcripts[acc1], consensus_transcripts[acc2], task="path")
                # ed = result["editDistance"]
                # cigar = result["cigar"]
                # if ed <= 1 or acc1 == "read_3959_support_275_521_0.0_955_1" or acc2 == "read_3959_support_275_521_0.0_955_1":
                #     print("ED:", ed, cigar, acc1,acc2 )

    return edges

def get_minimizer_graph_transposed(candidate_transcripts):
    best_edit_distances = {}
    no_ref_to_test_to = set()
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

        if len(best_edit_distances[c1]) == 0: # all isolated nodes in this graph
            no_ref_to_test_to.add(c1)

    # print(best_edit_distances["transcript_62_support_6"]) 
    minimizer_graph_transposed = transpose(best_edit_distances)
    print("isolated:", no_ref_to_test_to)
    for c_isolated in no_ref_to_test_to:
        minimizer_graph_transposed[c_isolated] = {}

    # print(minimizer_graph_transposed["transcript_46_support_430"])
    # print(minimizer_graph_transposed["transcript_62_support_6"])

    # sys.exit()

    assert len(best_edit_distances) == len(candidate_transcripts)
    return minimizer_graph_transposed

def check_exon_diffs(alignments_of_x_to_c, params):
    #################################################################
    ###### temp check for best alignment to wrong isoform ###########
    import re
    pattern = r"[-]{8,}"
    cccntr = 0
    print("Barcodes:", params.barcodes )
    out_file = open(os.path.join(params.outfolder, "statistical_exon_difs.fa"), "w")
    if params.barcodes:
        for s1, s1_dict in list(alignments_of_x_to_c.items()): 
            for s2, alignment_tuple in list(s1_dict.items()):
                if re.search(pattern, alignment_tuple[1][20: -20]) or  re.search(pattern, alignment_tuple[2][20: -20]): # [20: -20] --> ignore this if-statement if missing or truncated barcode
                    # del alignments_of_x_to_c[s1][s2]
                    # print("Deleted:", len(s1)," minimizer length:", len(s2), "length alignment:", len(alignment_tuple[2]), "edit distance:", alignment_tuple[0])
                    # print(s2)
                    cccntr += 1
                    out_file.write(">{0}\n{1}\n".format(s1, alignment_tuple[1].replace("-", "") ))
    else:
        for s1, s1_dict in list(alignments_of_x_to_c.items()): 
            for s2, alignment_tuple in list(s1_dict.items()):
                if re.search(pattern, alignment_tuple[1]) or  re.search(pattern, alignment_tuple[2]):
                    # del alignments_of_x_to_c[s1][s2]
                    cccntr += 1        
                    out_file.write(">{0}\n{1}\n".format(s1, alignment_tuple[1].replace("-", "") ))

    print("Number containing exon differences in this pass:", cccntr)
    # sys.exit()
    ###################################################################
    ###################################################################

def stat_filter_candidates(read_file, candidate_file, partition_of_X, to_realign, params):
    modified = True

    ############ GET READS AND CANDIDATES #################
    X_original = {acc: seq for (acc, seq) in  fasta_parser.read_fasta(open(read_file, 'r'))} 
    print("Total original reads", len(X_original))
    x_assigned_to_cluster = set([ x_acc for c_acc in partition_of_X for x_acc in partition_of_X[c_acc] ])
    X = {acc: seq for (acc, seq) in X_original.items() if acc in x_assigned_to_cluster or acc in  to_realign }    # just set X to partition_of_X + to_realign here
    print("Reads included in statistical testing:", len(X))
    if os.stat(candidate_file).st_size == 0:
        out_file_name = os.path.join(params.outfolder, "final_candidates.fa")
        write_output.print_candidates(out_file_name, {}, {}, {}, {}, final = True)
        print("Candidate file is empty!")
        sys.exit(0)
    else:
        C = {acc: seq for (acc, seq) in  fasta_parser.read_fasta(open(candidate_file, 'r'))}

    ################################################################
    print()
    print("STARTING STATISTICAL TESTING")
    print()
    print("Number of reads to realign:", len(to_realign))
    step = 1
    
    previous_partition_of_X = copy.deepcopy(partition_of_X) #{ c_acc : set() for c_acc in C.keys()}
    previous_components = { c_acc : set() for c_acc in C.keys()}
    significance_values = {}   
    realignment_to_avoid_local_max = 0
    remaining_to_align_read_file = os.path.join(params.outfolder, "remaining_to_align.fa")


    while modified:
        statistical_start = time() 

        modified = False
        print("STEP NR: {0}".format(step))

        ########### Write current candidates to file ##########
        temp_candidate_name = os.path.join(params.outfolder, "temp_candidates_step_{0}.fa".format(step))
        temp_candidate_file = open(temp_candidate_name, "w")

        for c_acc, c_seq in C.items():
            temp_candidate_file.write(">{0}\n{1}\n".format(c_acc, c_seq))
        temp_candidate_file.close()
        #######################################################

        # # create partition
        # partition_of_X = { c_acc : set() for c_acc in C.keys()}
        # for x_acc in alignments_of_x_to_c:
        #     for c_acc in alignments_of_x_to_c[x_acc]:
        #         partition_of_X[c_acc].add(x_acc)

        for c_acc in partition_of_X:
            print(c_acc, "has {0} reads assigned to it.".format(len(partition_of_X[c_acc])))

        ############ GET READ SUPORT AND ALIGNMENTS #################

        if realignment_to_avoid_local_max == 1:
            print("REALIGNING EVERYTHING FINAL STEP")
            to_realign = X      
            partition_of_X = { c_acc : set() for c_acc in C.keys()}
            alignments_of_x_to_c = {}



        if to_realign:
            write_output.print_reads(remaining_to_align_read_file, to_realign)
            # align reads that is not yet assigned to candidate here
            G_star_rem, partition_of_remaining_X = partition_strings_2set(to_realign, C, remaining_to_align_read_file, temp_candidate_file.name, params)

            # add reads to best candidate given new alignments
            for c_acc in partition_of_remaining_X:
                partition_of_X[c_acc].update(partition_of_remaining_X[c_acc])
                print(c_acc, "Now has {0} reads assigned to it, after aligning reads that are not assigned.".format(len(partition_of_X[c_acc])))

            # add the alignments to alignment structure
            # for x_acc in remaining_alignments_of_x_to_c.keys():
            #     alignments_of_x_to_c[x_acc] = remaining_alignments_of_x_to_c[x_acc]

        C_seq_to_acc = {seq : acc for acc, seq in C.items()}
        ################################################################


        # check_exon_diffs(alignments_of_x_to_c, params)

        ############# GET THE CLOSES HIGHEST SUPPORTED REFERENCE TO TEST AGAINST FOR EACH CANDIDATE ############

        # print("NEW")
        # minimizer_graph_transposed = end_invariant_functions.get_minimizers_graph_transposed_under_ignored_ends(C, params)
        # for t in minimizer_graph_transposed:
        #     print(t, "nr candidates:", len(minimizer_graph_transposed[t]), "deltas:", [minimizer_graph_transposed[t][c] for c in minimizer_graph_transposed[t]])
        # print("OLD")
        # minimizer_graph_transposed2 = get_minimizer_graph_transposed(C)
        # for t in minimizer_graph_transposed2:
        #     print(t, "nr candidates:", len(minimizer_graph_transposed2[t]), "deltas:", [minimizer_graph_transposed2[t][c] for c in minimizer_graph_transposed2[t]])
        # sys.exit()

        if params.ignore_ends_len > 0:
            minimizer_graph_transposed = end_invariant_functions.get_minimizers_graph_transposed_under_ignored_ends(C, params)
        else:
            minimizer_graph_transposed = get_minimizer_graph_transposed(C)

        # extra_edges_from_collapsed_homopolymers = get_homopolymer_invariants(C)
        # homopol_extra_added = 0
        # for acc1 in extra_edges_from_collapsed_homopolymers:
        #     for acc2 in extra_edges_from_collapsed_homopolymers[acc1]:
        #         if acc1 in minimizer_graph_transposed:
        #             if acc2 not in minimizer_graph_transposed[acc1]:
        #                 minimizer_graph_transposed[acc1][acc2] = "homopolymer_identical"
        #                 homopol_extra_added += 1
        #         else:
        #             minimizer_graph_transposed[acc1] = {}
        #             minimizer_graph_transposed[acc1][acc2] = "homopolymer_identical"
        #             homopol_extra_added += 1

        # print("EXTRA EDGES FROM HOMOPOLYMER IDENTICAL:", homopol_extra_added)


        ## save time if the minimizer and all cantidates in a component has identical reads assignmed to them as previous step
        # Since indata is the same, the test is guaranteed to give same siginficance values as previous step

        previous_significance_values = {}
        for t_acc in list(minimizer_graph_transposed.keys()):
            accessions_in_component_to_test = set([c_acc for c_acc in minimizer_graph_transposed[t_acc].keys()] + [t_acc])
            # print(accessions_in_component_to_test)
            if (accessions_in_component_to_test == previous_components[t_acc]) and all([ len(previous_partition_of_X[acc]) == len(partition_of_X[acc]) for acc in accessions_in_component_to_test]):
                print("TEST IDENTICAL TO PREVIOUS STEP, SKIPPING FOR", t_acc)
                for acc in accessions_in_component_to_test:
                    previous_significance_values[acc] = significance_values[acc]
                del minimizer_graph_transposed[t_acc]
            else: 
                # print("Modified")
                previous_components[t_acc] = accessions_in_component_to_test

        #####################################################################################################


        # get all candidats that serve as null-hypothesis references and have neighbors subject to testing
        # these are all candidates that are minimizers to some other, isolated nodes are not tested
        # candidatate in G_star_C
        nr_of_tests_this_round = len([ 1 for t_acc in minimizer_graph_transposed for c_acc in minimizer_graph_transposed[t_acc] ] )
        print("NUMBER OF CANDIDATES LEFT:", len(C), "Number of performed statistical tests in this round:", nr_of_tests_this_round)

        if realignment_to_avoid_local_max == 1:
            new_significance_values = statistical_test_v2.do_statistical_tests_all_c_to_t(minimizer_graph_transposed, C, X, partition_of_X, params )
        else:
            new_significance_values = statistical_test_v2.do_statistical_tests_per_edge(minimizer_graph_transposed, C, X, partition_of_X, params )

        previous_partition_of_X = copy.deepcopy(partition_of_X)
        to_realign = {}
        for c_acc, (p_value, mult_factor_inv, k, N_t, delta_size) in list(new_significance_values.items()):
            if p_value == "not_tested":
                print("Did not test", c_acc)
            elif p_value * mult_factor_inv > 0.001:
                print("removing", c_acc, "p-val:", p_value, "correction factor:", mult_factor_inv, "k", k, "N_t", N_t, "delta_size:", delta_size )
                del C[c_acc] 
                modified = True
                for x_acc in partition_of_X[c_acc]:
                    to_realign[x_acc] = X[x_acc]
                    # del alignments_of_x_to_c[x_acc]

                del partition_of_X[c_acc]

        print("nr candidates left:", len(C))
        candidate_file = os.path.join(params.outfolder, "candidates_after_step_{0}.fa".format(step))
        step += 1

        significance_values = previous_significance_values.copy()
        significance_values.update(new_significance_values) # which returns None since it mutates z
        
        if len(C) == 0: # no candidates were significant!
            break

        print("LEN SIGN:", len(significance_values), len(C))
        write_output.print_candidates(candidate_file, C, significance_values, partition_of_X, X)

        # do a last realingment to avoind local maxima of reads

        if realignment_to_avoid_local_max == 1: # we have already done a last realignment, keep going until everythin is significant never realign
            realignment_to_avoid_local_max = 2
        elif not modified and realignment_to_avoid_local_max == 0: # we have not yet done a final alignment and everythin is significant, realign to escape local maxima alignment
            realignment_to_avoid_local_max = 1
            modified = True

        statistical_elapsed = time() - statistical_start
        write_output.logger('Time for Statistical test, step {0}:{1}'.format(step, str(statistical_elapsed)), params.logfile)
   

    final_out_file_name =  os.path.join(params.outfolder, "final_candidates.fa")
    tsv_info = os.path.join(params.outfolder, "cluster_info.tsv")
    write_output.print_candidates(final_out_file_name, C, significance_values, partition_of_X, X, final = True, reads_to_consensus_tsv = tsv_info)

    return C
