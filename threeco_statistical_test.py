"""
    PFM is a list of dicts, where the inner dicts are one per column in the PFM matrix, its keys by characters A, C, G, T, -
    alignment_matrix is a representation of all alignments in a partition. this is a dictionary where sequences s_i belonging to the 
    partition as keys and the alignment of s_i with respectt to the alignment matix.
"""
import os
import unittest
import copy
import math 

from scipy.stats import poisson

from modules.functions import transpose, create_position_probability_matrix
from modules import functions
from modules.partitions import partition_strings_paths, partition_strings_2set
from modules import graphs
from modules.SW_alignment_module import sw_align_sequences, sw_align_sequences_keeping_accession
from modules.edlib_alignment_module import edlib_align_sequences, edlib_align_sequences_keeping_accession, edlib_traceback
from modules.input_output import fasta_parser, write_output
from modules import statistical_test_v2
from modules import correct_sequence_to_minimizer


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

def get_minimizer_graph(candidate_transcripts):
    best_edit_distances = {}
    isolated_nodes = set()
    for i, (c1, seq1) in enumerate(candidate_transcripts.items()):
        if i % 50 == 0:
            print("processing ", i)
        best_edit_distances[c1] = {}
        # best_cigars[c1] = {}
        best_ed = len(seq1)
        for c2, seq2 in  candidate_transcripts.items() :
            if c1 == c2:
                continue
            elif math.fabs(len(seq1) - len(seq2)) > best_ed:
                continue
            # TODO: remove task = "path" to speed up
            edit_distance, locations, cigar = edlib_traceback(seq1, seq2, mode="NW", task="path", k=min(15, best_ed))

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
            isolated_nodes.add(c1)

 
    minimizer_graph = transpose(best_edit_distances)
    for c_isolated in isolated_nodes:
        minimizer_graph[c_isolated] = {}
    # seen_in_test = set()
    # for m in minimizer_graph:
        # print(best_edit_distances[c1])
        # print("NEW", m, "size:", len(minimizer_graph[m]))
        # if len(best_edit_distances[c1]) > 1:
        #     print("Minimizer to multiple candidates:", m, len(minimizer_graph[m]))
        # for c in minimizer_graph[m]:
        #     ed = minimizer_graph[m][c]
            # if c in seen_in_test:
            #     print("Seen:", c)
            # seen_in_test.add(c)
            # print("ED:",ed, c )
                # print("Multiple best:", c1, len(c2), best_cigars[c1][c2])
                # print(c2)
                # print()
        # print()

    return minimizer_graph


def stat_filter_candidates(read_file, candidate_file, alignments_of_x_to_c, params):
    modified = True

    ############ GET READS AND CANDIDATES #################
    X = {acc: seq for (acc, seq) in  fasta_parser.read_fasta(open(read_file, 'r'))} 
    C = {acc: seq for (acc, seq) in  fasta_parser.read_fasta(open(candidate_file, 'r'))}

    ################################################################

    step = 1
    
    previous_partition_of_X = { c_acc : set() for c_acc in C.keys()}
    previous_components = { c_acc : set() for c_acc in C.keys()}
    significance_values = {}   
    realignment_to_avoid_local_max = 0
    to_realign = {}
    remaining_to_align_read_file = os.path.join(params.outfolder, "remaining_to_align.fa")


    while modified:
        modified = False
        print("STEP NR: {0}".format(step))

        ########### Write current candidates to file ##########
        temp_candidate_name = os.path.join(params.outfolder, "temp_candidates_step_{0}.fa".format(step))
        temp_candidate_file = open(temp_candidate_name, "w")

        for c_acc, c_seq in C.items():
            temp_candidate_file.write(">{0}\n{1}\n".format(c_acc, c_seq))
        temp_candidate_file.close()
        #######################################################

        # create partition
        partition_of_X = { c_acc : set() for c_acc in C.keys()}
        for x_acc in alignments_of_x_to_c:
            for c_acc in alignments_of_x_to_c[x_acc]:
                partition_of_X[c_acc].add(x_acc)

        ############ GET READ SUPORT AND ALIGNMENTS #################

        if realignment_to_avoid_local_max == 1:
            print("REALIGNING EVERYTHING FINAL STEP")
            to_realign = X      
            partition_of_X = { c_acc : set() for c_acc in C.keys()}



        if to_realign:
            write_output.print_reads(remaining_to_align_read_file, to_realign)
            # align reads that is not yet assigned to candidate here
            G_star_rem, partition_of_remaining_X, remaining_alignments_of_x_to_c = partition_strings_2set(to_realign, C, remaining_to_align_read_file, temp_candidate_file.name)


            #################################################################
            ###### temp check for best alignment to wrong isoform ###########
            import re
            pattern = r"[-]{8,}"
            cccntr = 0
            print("Barcodes:", params.barcodes )
            out_file = open(os.path.join(params.outfolder, "statistical_exon_difs.fa"), "w")
            if params.barcodes:
                for s1, s1_dict in list(remaining_alignments_of_x_to_c.items()): 
                    for s2, alignment_tuple in list(s1_dict.items()):
                        if re.search(pattern, alignment_tuple[1][20: -20]) or  re.search(pattern, alignment_tuple[2][20: -20]): # [20: -20] --> ignore this if-statement if missing or truncated barcode
                            del remaining_alignments_of_x_to_c[s1][s2]
                            print("Deleted:", len(s2)," minimizer length:", len(s1), "length alignment:", len(alignment_tuple[2]), "edit distance:", alignment_tuple[0])
                            print(s2)
                            cccntr += 1
                            out_file.write(">{0}\n{1}\n".format(unique_seq_to_acc[s2],s2))
            else:
                for s1, s1_dict in list(remaining_alignments_of_x_to_c.items()): 
                    for s2, alignment_tuple in list(s1_dict.items()):
                        if re.search(pattern, alignment_tuple[1]) or  re.search(pattern, alignment_tuple[2]):
                            del remaining_alignments_of_x_to_c[s1][s2]
                            cccntr += 1        
                            out_file.write(">{0}\n{1}\n".format(unique_seq_to_acc[s2],s2))

            print("Number containing exon difference and removed in this pass:", cccntr)
            # sys.exit()
            ###################################################################
            ###################################################################


            # add reads to best candidate given new alignments
            for c_acc in partition_of_remaining_X:
                partition_of_X[c_acc].update(partition_of_remaining_X[c_acc])

            # add the alignments to alignment structure
            for x_acc in remaining_alignments_of_x_to_c.keys():
                alignments_of_x_to_c[x_acc] = remaining_alignments_of_x_to_c[x_acc]
                # for c_acc in remaining_alignments_of_x_to_c[x_acc].keys():
                #     alignments_of_x_to_c[x_acc][c_acc] = remaining_alignments_of_x_to_c[x_acc][c_acc]

        C_seq_to_acc = {seq : acc for acc, seq in C.items()}
        ################################################################


        ############# GET THE CLOSES HIGHEST SUPPORTED REFERENCE TO TEST AGAINST FOR EACH CANDIDATE ############
        minimizer_graph = get_minimizer_graph(C)

        ## save time if the minimizer and all cantidates in a component has identical reads assignmed to them as previous step
        # Since indata is the same, the test is guaranteed to give same siginficance values as previous step

        previous_significance_values = {}
        for t_acc in list(minimizer_graph.keys()):
            accessions_in_component_to_test = set([c_acc for c_acc in minimizer_graph[t_acc].keys()] + [t_acc])
            # print(accessions_in_component_to_test)
            if (accessions_in_component_to_test == previous_components[t_acc]) and all([ len(previous_partition_of_X[acc]) == len(partition_of_X[acc]) for acc in accessions_in_component_to_test]):
                print("TEST IDENTICAL TO PREVIOUS STEP, SKIPPING FOR", t_acc)
                for acc in accessions_in_component_to_test:
                    previous_significance_values[acc] = significance_values[acc]
                del minimizer_graph[t_acc]
            else: 
                # print("Modified")
                previous_components[t_acc] = accessions_in_component_to_test

        #####################################################################################################


        # get all candidats that serve as null-hypothesis references and have neighbors subject to testing
        # these are all candidates that are minimizers to some other, isolated nodes are not tested
        # candidatate in G_star_C
        nr_of_tests_this_round = len(minimizer_graph)
        print("NUMBER OF CANDIDATES LEFT:", len(C))

        if "read_2190_support_11" in minimizer_graph:
            print("DETECTED:", minimizer_graph["read_2190_support_11"])
        else:
            for t_acc in minimizer_graph:
                for c_acc in minimizer_graph[t_acc]:
                    if c_acc == "read_2190_support_11":
                        print("read_2190_support_11, A candidate:", minimizer_graph[t_acc])

        new_significance_values = statistical_test_v2.do_statistical_tests(minimizer_graph, C, X, partition_of_X, single_core = params.single_core )
        previous_partition_of_X = copy.deepcopy(partition_of_X)
        to_realign = {}
        for c_acc, (corrected_p_value, k, N_t) in list(new_significance_values.items()):

            if "read_2190_support_11" == c_acc:
                print("read_2190_support_11 is here with", (corrected_p_value, k, N_t) )

            if corrected_p_value == "not_tested":
                print("Did not test", c_acc)
            elif corrected_p_value > 0.01/nr_of_tests_this_round:
                print("removing", c_acc, "p-val:", corrected_p_value, "k", k, "N_t", N_t )
                del C[c_acc] 
                modified = True
                for x_acc in partition_of_X[c_acc]:
                    to_realign[x_acc] = X[x_acc]
                    del alignments_of_x_to_c[x_acc]

                del partition_of_X[c_acc]

        print("nr candidates left:", len(C))
        candidate_file = os.path.join(params.outfolder, "candidates_after_step_{0}.fa".format(step))
        step += 1
        significance_values = previous_significance_values.copy()
        significance_values.update(new_significance_values) # which returns None since it mutates z

        print("LEN SIGN:", len(significance_values), len(C))
        write_output.print_candidates(candidate_file, alignments_of_x_to_c, C, significance_values)

        print("GOOOOOO:", "read_2190_support_11" in C)
        # do a last realingment to avoind local maxima of reads

        if realignment_to_avoid_local_max == 1: # we have already done a last realignment, keep going until everythin is significant never realign
            realignment_to_avoid_local_max == 2
        elif not modified and realignment_to_avoid_local_max == 0: # we have not yet done a final alignment and everythin is significant, realign to escape local maxima alignment
            realignment_to_avoid_local_max = 1
            modified = True

    final_out_file_name =  os.path.join(params.outfolder, "final_candidates.fa")
    write_output.print_candidates(final_out_file_name, alignments_of_x_to_c, C, significance_values, final = True)

    return C
