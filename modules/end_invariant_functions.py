from __future__ import print_function
import os
import sys
import argparse
import re
# import numpy as np
import signal
from multiprocessing import Pool
import multiprocessing as mp
import math
from time import time
import copy

import networkx as nx
import edlib

from modules import functions
from modules.input_output import write_output
from modules.SW_alignment_module import parasail_alignment

def parasail_traceback_allow_ends(x, y, end_threshold = 0):
    (s1, s2, (s1_alignment, s2_alignment, (matches, mismatches, indels)) ) = parasail_alignment(x, y, 0, 0, opening_penalty = 3, mismatch_penalty = -3)
    ed = mismatches + indels
    p = "[-]+"
    m_start1 = re.match(p,s1_alignment)
    m_start2 = re.match(p,s2_alignment)
    if m_start1:
        cut_start_diff1 = len(m_start1.group(0))
        ed -= min(cut_start_diff1, end_threshold)
    elif m_start2:
        cut_start_diff2 = len(m_start2.group(0))
        ed -= min(cut_start_diff2, end_threshold)

    m_end1 = re.match(p,s1_alignment[::-1])
    m_end2 = re.match(p,s2_alignment[::-1])
    if m_end1:
        cut_end_diff1 = len(m_end1.group(0))
        ed -= min(cut_end_diff1, end_threshold)
    elif m_end2:
        cut_end_diff2 = len(m_end2.group(0))
        ed -= min(cut_end_diff2, end_threshold)

    return ed


def get_nearest_neighbors_parasail(batch_of_queries, global_index_in_matrix, start_index, seq_to_acc_list_sorted, neighbor_search_depth, ignore_ends_len):
    best_edit_distances = {}
    lower_target_edit_distances = {}
    # print("Processing global index:" , global_index_in_matrix)
    # error_types = {"D":0, "S": 0, "I": 0}
    for i in range(start_index, start_index + len(batch_of_queries)):
        if i % 500 == 0:
            print("processing ", i)
        seq1 = seq_to_acc_list_sorted[i][0]
        acc1 = seq_to_acc_list_sorted[i][1]
        best_edit_distances[acc1] = {}

        if acc1 in lower_target_edit_distances:
            best_ed = lower_target_edit_distances[acc1] 
            # print("already_comp", best_ed )
        else:
            best_ed = len(seq1)

        stop_up = False
        stop_down = False
        j = 1
        while True:
        # for j in range(1,len(seq_to_acc_list_sorted)):
            if i - j < 0:
                stop_down = True
            if i + j >= len(seq_to_acc_list_sorted):
                stop_up = True

            if not stop_down:
                seq2 = seq_to_acc_list_sorted[i - j][0]
                acc2 = seq_to_acc_list_sorted[i - j][1]  

                if math.fabs(len(seq1) - len(seq2)) > best_ed + 2*ignore_ends_len:
                    stop_down = True

            if not stop_up:
                seq3 = seq_to_acc_list_sorted[i + j][0]
                acc3 = seq_to_acc_list_sorted[i + j][1]  

                if math.fabs(len(seq1) - len(seq3)) > best_ed + 2*ignore_ends_len:
                    stop_up = True

            if not stop_down:
                edit_distance = parasail_traceback_allow_ends(seq1, seq2, end_threshold = ignore_ends_len)

                if 0 < edit_distance < best_ed:
                    best_ed = edit_distance
                    best_edit_distances[acc1] = {}
                    best_edit_distances[acc1][acc2] = edit_distance
                elif edit_distance == best_ed:
                    best_edit_distances[acc1][acc2] = edit_distance

                if acc2 in lower_target_edit_distances:
                    if 0 < edit_distance < lower_target_edit_distances[acc2]: 
                        lower_target_edit_distances[acc2] = edit_distance 
                else:
                    if 0 < edit_distance: 
                        lower_target_edit_distances[acc2] = edit_distance 

            if not stop_up:
                edit_distance = parasail_traceback_allow_ends(seq1, seq3, end_threshold = ignore_ends_len)

                if 0 < edit_distance < best_ed:
                    best_ed = edit_distance
                    best_edit_distances[acc1] = {}
                    best_edit_distances[acc1][acc3] = edit_distance
                elif edit_distance == best_ed:
                    best_edit_distances[acc1][acc3] = edit_distance

                if acc3 in lower_target_edit_distances:
                    if 0 < edit_distance < lower_target_edit_distances[acc3]: 
                        lower_target_edit_distances[acc3] = edit_distance 
                else:
                    if 0 < edit_distance:                 
                        lower_target_edit_distances[acc3] = edit_distance 
            
            if stop_down and stop_up:
                break

            if j >= neighbor_search_depth:
                break
            j += 1
        # print(best_edit_distances[acc1])
        # print("best ed:", best_ed)
        # if best_ed > 100:
        #     print(best_ed, "for seq with length", len(seq1), seq1)
    return best_edit_distances

def get_nearest_neighbors_under_ignored_edge_ends_parasail(seq_to_acc_list_sorted, params):
    if params.nr_cores == 1:
        best_edit_distances = get_nearest_neighbors_parasail(seq_to_acc_list_sorted, 0, 0, seq_to_acc_list_sorted, params.neighbor_search_depth, params.ignore_ends_len)

        # implement check here to se that all seqs got a nearest_neighbor, if not, print which noes that did not get a nearest_neighbor computed.!

    else:
        ####### parallelize alignment #########
        # pool = Pool(processes=mp.cpu_count())
        original_sigint_handler = signal.signal(signal.SIGINT, signal.SIG_IGN)
        signal.signal(signal.SIGINT, original_sigint_handler)
        pool = Pool(processes=mp.cpu_count())

        # here we split the input into chunks
        chunk_size = max(int(len(seq_to_acc_list_sorted) / (10*mp.cpu_count())), 20 )
        ref_seq_chunks = [ ( max(0, i - params.neighbor_search_depth -1), seq_to_acc_list_sorted[max(0, i - params.neighbor_search_depth -1) : i + chunk_size + params.neighbor_search_depth +1 ]) for i in range(0, len(seq_to_acc_list_sorted), chunk_size) ]
        chunks = [(i, seq_to_acc_list_sorted[i:i + chunk_size]) for i in range(0, len(seq_to_acc_list_sorted), chunk_size)] 

        if params.verbose:
            write_output.logger(str([j for j, ch in ref_seq_chunks]), params.develop_logfile, timestamp=False)
            write_output.logger("reference chunks:" + str([len(ch) for j,ch in ref_seq_chunks]), params.develop_logfile, timestamp=False)
            # print([j for j, ch in ref_seq_chunks])
            # print("reference chunks:", [len(ch) for j,ch in ref_seq_chunks])
            write_output.logger(str([i for i,ch in chunks]), params.develop_logfile, timestamp=False)
            write_output.logger("query chunks:" + str([len(ch) for i,ch in chunks]), params.develop_logfile, timestamp=False)

            print([i for i,ch in chunks])
            print("query chunks:", [len(ch) for i,ch in chunks])
        # get nearest_neighbors takes thre sub containers: 
        #  chunk - a container with (sequences, accesions)-tuples to be aligned (queries)
        #  ref_seq_chunks - a container with (sequences, accesions)-tuples to be aligned to (references)
        #  already_converged_chunks - a set of query sequences that has already converged 

        try:
            res = pool.map_async(get_nearest_neighbors_helper_parasail, [ ((chunks[i][1],  chunks[i][0], chunks[i][0] - ref_seq_chunks[i][0], ref_seq_chunks[i][1], params.neighbor_search_depth, params.ignore_ends_len), {}) for i in range(len(chunks))] )
            best_edit_distances_results =res.get(999999999) # Without the timeout this blocking call ignores all signals.
        except KeyboardInterrupt:
            print("Caught KeyboardInterrupt, terminating workers")
            pool.terminate()
            sys.exit()
        else:
            # print("Normal termination")
            pool.close()
        pool.join()
        best_edit_distances = {}
        for sub_graph in best_edit_distances_results:
            for seq in sub_graph:
                assert seq not in best_edit_distances
            best_edit_distances.update(sub_graph)

    return best_edit_distances

def get_nearest_neighbors_helper_parasail(arguments):
    args, kwargs = arguments
    return get_nearest_neighbors_parasail(*args, **kwargs)


def edlib_traceback_allow_ends(x, y, mode="NW", task="path", k=1, end_threshold = 0):
    result = edlib.align(x, y, mode=mode, task=task, k=k)
    ed = result["editDistance"]
    locations =  result["locations"]
    cigar =  result["cigar"]

    if cigar:
        tuples = []
        result = re.split(r'[=DXSMI]+', cigar)
        i = 0
        for length in result[:-1]:
            i += len(length)
            type_ = cigar[i]
            i += 1
            tuples.append((length, type_ ))

        ed_ignore_ends = ed
        if tuples[0][1] == "D" or  tuples[0][1] == "I":
            begin_snippet = int(tuples[0][0])
            if begin_snippet <= end_threshold:
                ed_ignore_ends -= int(begin_snippet)
        if tuples[-1][1] == "D" or  tuples[-1][1] == "I":
            end_snippet = int(tuples[-1][0])
            if end_snippet <= end_threshold:
                ed_ignore_ends -= int(end_snippet)  
        # if ed > ed_ignore_ends:          
        #     print("ed global:", ed, "ed after:", ed_ignore_ends)
        ed = ed_ignore_ends

    # if ed ==0:
    #     print("here")

    return ed, locations, cigar


def read_fasta(fasta_file):
    fasta_seqs = {}
    k = 0
    temp = ''
    accession = ''
    for line in fasta_file:
        if line[0] == '>' and k == 0:
            accession = line[1:].strip()
            fasta_seqs[accession] = ''
            k += 1
        elif line[0] == '>':
            yield accession, temp
            temp = ''
            accession = line[1:].strip()
        else:
            temp += line.strip()
    if accession:
        yield accession, temp


def get_nearest_neighbors_helper(arguments):
    args, kwargs = arguments
    return get_nearest_neighbors(*args, **kwargs)


def get_nearest_neighbors_under_ignored_edge_ends(seq_to_acc_list_sorted, params):
    if params.nr_cores == 1:
        best_edit_distances = get_nearest_neighbors(seq_to_acc_list_sorted, 0, 0, seq_to_acc_list_sorted, params.neighbor_search_depth, params.ignore_ends_len)

        # implement check here to se that all seqs got a nearest_neighbor, if not, print which noes that did not get a nearest_neighbor computed.!

    else:
        ####### parallelize alignment #########
        # pool = Pool(processes=mp.cpu_count())
        original_sigint_handler = signal.signal(signal.SIGINT, signal.SIG_IGN)
        signal.signal(signal.SIGINT, original_sigint_handler)
        pool = Pool(processes=mp.cpu_count())

        # here we split the input into chunks
        chunk_size = max(int(len(seq_to_acc_list_sorted) / (10*mp.cpu_count())), 20 )
        ref_seq_chunks = [ ( max(0, i - params.neighbor_search_depth -1), seq_to_acc_list_sorted[max(0, i - params.neighbor_search_depth -1) : i + chunk_size + params.neighbor_search_depth +1 ]) for i in range(0, len(seq_to_acc_list_sorted), chunk_size) ]
        chunks = [(i, seq_to_acc_list_sorted[i:i + chunk_size]) for i in range(0, len(seq_to_acc_list_sorted), chunk_size)] 

        if params.verbose:
            write_output.logger(str([j for j, ch in ref_seq_chunks]), params.develop_logfile, timestamp=False)
            write_output.logger("reference chunks:" + str([len(ch) for j,ch in ref_seq_chunks]), params.develop_logfile, timestamp=False)
            # print([j for j, ch in ref_seq_chunks])
            # print("reference chunks:", [len(ch) for j,ch in ref_seq_chunks])
            write_output.logger(str([i for i,ch in chunks]), params.develop_logfile, timestamp=False)
            write_output.logger("query chunks:" + str([len(ch) for i,ch in chunks]), params.develop_logfile, timestamp=False)

            print([i for i,ch in chunks])
            print("query chunks:", [len(ch) for i,ch in chunks])
        # get nearest_neighbors takes thre sub containers: 
        #  chunk - a container with (sequences, accesions)-tuples to be aligned (queries)
        #  ref_seq_chunks - a container with (sequences, accesions)-tuples to be aligned to (references)
        #  already_converged_chunks - a set of query sequences that has already converged 

        try:
            res = pool.map_async(get_nearest_neighbors_helper, [ ((chunks[i][1],  chunks[i][0], chunks[i][0] - ref_seq_chunks[i][0], ref_seq_chunks[i][1], params.neighbor_search_depth, params.ignore_ends_len), {}) for i in range(len(chunks))] )
            best_edit_distances_results =res.get(999999999) # Without the timeout this blocking call ignores all signals.
        except KeyboardInterrupt:
            print("Caught KeyboardInterrupt, terminating workers")
            pool.terminate()
            sys.exit()
        else:
            # print("Normal termination")
            pool.close()
        pool.join()
        best_edit_distances = {}
        for sub_graph in best_edit_distances_results:
            for seq in sub_graph:
                assert seq not in best_edit_distances
            best_edit_distances.update(sub_graph)

    return best_edit_distances



def get_nearest_neighbors(batch_of_queries, global_index_in_matrix, start_index, seq_to_acc_list_sorted, neighbor_search_depth, ignore_ends_threshold):
    best_edit_distances = {}
    lower_target_edit_distances = {}
    # print("Processing global index:" , global_index_in_matrix)
    # error_types = {"D":0, "S": 0, "I": 0}
    for i in range(start_index, start_index + len(batch_of_queries)):
        if i % 500 == 0:
            print("processing ", i)
        seq1 = seq_to_acc_list_sorted[i][0]
        acc1 = seq_to_acc_list_sorted[i][1]
        best_edit_distances[acc1] = {}

        if acc1 in lower_target_edit_distances:
            best_ed = lower_target_edit_distances[acc1] 
            # print("already_comp", best_ed )
        else:
            best_ed = len(seq1)

        stop_up = False
        stop_down = False
        j = 1
        while True:
        # for j in range(1,len(seq_to_acc_list_sorted)):
            if i - j < 0:
                stop_down = True
            if i + j >= len(seq_to_acc_list_sorted):
                stop_up = True

            if not stop_down:
                seq2 = seq_to_acc_list_sorted[i - j][0]
                acc2 = seq_to_acc_list_sorted[i - j][1]  

                if math.fabs(len(seq1) - len(seq2)) > best_ed + 2*ignore_ends_threshold:
                    stop_down = True

            if not stop_up:
                seq3 = seq_to_acc_list_sorted[i + j][0]
                acc3 = seq_to_acc_list_sorted[i + j][1]  

                if math.fabs(len(seq1) - len(seq3)) > best_ed + 2*ignore_ends_threshold:
                    stop_up = True

            if not stop_down:
                edit_distance_f, locations, cigar = edlib_traceback_allow_ends(seq1, seq2, mode="NW", task="path", k=best_ed+2*ignore_ends_threshold, end_threshold = ignore_ends_threshold)
                edit_distance_r, locations, cigar = edlib_traceback_allow_ends(seq1[::-1], seq2[::-1], mode="NW", task="path", k=best_ed+2*ignore_ends_threshold, end_threshold = ignore_ends_threshold)
                if edit_distance_f >= 0 and edit_distance_r >= 0:
                    edit_distance = min(edit_distance_f, edit_distance_r)
                else:
                    edit_distance = max(edit_distance_f, edit_distance_r)

                if 0 <= edit_distance < best_ed:
                    best_ed = edit_distance
                    best_edit_distances[acc1] = {}
                    best_edit_distances[acc1][acc2] = edit_distance
                elif edit_distance == best_ed:
                    best_edit_distances[acc1][acc2] = edit_distance

                if acc2 in lower_target_edit_distances:
                    if 0 < edit_distance < lower_target_edit_distances[acc2]: 
                        lower_target_edit_distances[acc2] = edit_distance 
                else:
                    if 0 < edit_distance: 
                        lower_target_edit_distances[acc2] = edit_distance 

            if not stop_up:
                edit_distance_f, locations, cigar = edlib_traceback_allow_ends(seq1, seq3, mode="NW", task="path", k=best_ed+2*ignore_ends_threshold, end_threshold = ignore_ends_threshold)
                edit_distance_r, locations, cigar = edlib_traceback_allow_ends(seq1[::-1], seq3[::-1], mode="NW", task="path", k=best_ed+2*ignore_ends_threshold, end_threshold = ignore_ends_threshold)
                if edit_distance_f >= 0 and edit_distance_r >= 0:
                    edit_distance = min(edit_distance_f, edit_distance_r)
                else:
                    edit_distance = max(edit_distance_f, edit_distance_r)

                if 0 <= edit_distance < best_ed:
                    best_ed = edit_distance
                    best_edit_distances[acc1] = {}
                    best_edit_distances[acc1][acc3] = edit_distance
                elif edit_distance == best_ed:
                    best_edit_distances[acc1][acc3] = edit_distance

                if acc3 in lower_target_edit_distances:
                    if 0 < edit_distance < lower_target_edit_distances[acc3]: 
                        lower_target_edit_distances[acc3] = edit_distance 
                else:
                    if 0 < edit_distance:                 
                        lower_target_edit_distances[acc3] = edit_distance 
            
            if stop_down and stop_up:
                break

            if j >= neighbor_search_depth:
                break
            j += 1
        # print(best_edit_distances[acc1])
        # print("best ed:", best_ed)
        # if best_ed > 100:
        #     print(best_ed, "for seq with length", len(seq1), seq1)
    return best_edit_distances



def partition_highest_reachable_with_edge_degrees(G_star, params):
    # G_star, converged = graphs.construct_exact_nearest_neighbor_graph_improved(S, params)
    unique_start_strings = set(G_star.nodes())

    # print("len G_star:", len(G_star))
    partition_sizes = []
    nr_consensus = 0
    G_transpose = nx.reverse(G_star)
    # print("len G_star_transposed (nearest_neighbors):", len(G_transpose))
    if params.verbose:
        print("Nodes in nearest_neighbor graph:", len(G_transpose))
        print("Neighbors per nodes in nearest neighbor graph", sorted([len(list(G_transpose.neighbors(n)) ) for n in G_transpose], reverse=True))

    all_weak_components = [ c for c in nx.weakly_connected_components(G_transpose)]
    M = {}
    partition = {}
    # print("here")
    for subgraph_set in sorted(all_weak_components, key=len, reverse=True):
        # print("Subgraph of size", len(subgraph.nodes()), "nr edges:", [x for x in subgraph.nodes()] )
        while subgraph_set:
            reachable_comp_sizes = []
            reachable_comp_weights = {}
            reachable_comp_nodes = []
            direct_neighbors = {}
            processed = set()

            biggest_reachable_comp_size = 0
            biggest_reachable_comp_weight = 0
            biggest_reachable_comp_nodes = set()
            biggest_reachable_comp_nearest_neighbor = "XXXXX"


            for m in subgraph_set:                
                if m in processed:
                    continue

                reachable_comp = set([m])
                reachable_comp_weight = G_transpose.node[m]["degree"]
                processed.add(m)



                ####################################################
                # take all reachable nodes
                ####################################################

                for n1,n2 in nx.dfs_edges(G_transpose, source=m): # store reachable node as processed here to avoid computation
                    if n2 == m:
                        continue
                    processed.add(n2)
                    reachable_comp.add(n2)
                    reachable_comp_weight += G_transpose.node[n2]["degree"]
                ####################################################
                ####################################################
                

                # print("total component weight:", reachable_comp_weight)

                if biggest_reachable_comp_weight == 0:
                    biggest_reachable_comp_weight = reachable_comp_weight
                    biggest_reachable_comp_nodes = set(reachable_comp)
                    biggest_reachable_comp_size = len(reachable_comp)
                    biggest_reachable_comp_nearest_neighbor = m

                elif reachable_comp_weight >= biggest_reachable_comp_weight:
                    if reachable_comp_weight > biggest_reachable_comp_weight:
                        biggest_reachable_comp_weight = reachable_comp_weight
                        biggest_reachable_comp_nodes = set(reachable_comp)
                        biggest_reachable_comp_size = len(reachable_comp)
                        biggest_reachable_comp_nearest_neighbor = m

                    elif reachable_comp_weight == biggest_reachable_comp_weight:
                        if biggest_reachable_comp_weight > 1:
                            # print("tie both in weighted partition size and total edit distance. Choosing lexographically smaller nearest_neighbor")
                            # print(" weighted partition size:", biggest_reachable_comp_weight, " total edit distance:", edit_distances_to_m[m])
                            pass
                        
                        if m < biggest_reachable_comp_nearest_neighbor:
                            biggest_reachable_comp_nodes = set(reachable_comp)
                            biggest_reachable_comp_nearest_neighbor = m
                        else:
                            pass

                    else:
                        print("BUG!")

            if biggest_reachable_comp_weight == 0: # if there were no edges! partition is nearest_neighbor itself
                M[m] = 0 
                partition[m] = set()
            else:
                nearest_neighbor = biggest_reachable_comp_nearest_neighbor # "XXXXXX" #biggest_reachable_comp_nearest_neighbor #
                max_direct_weight = 0
                # print("total nodes searched in this pass:", len(biggest_reachable_comp_nodes))
                for n in biggest_reachable_comp_nodes:
                    direct_weight = G_transpose.node[n]["degree"]                    
                    direct_weight += len(list(G_transpose.neighbors(n)))

                    if direct_weight > max_direct_weight:
                        max_direct_weight = direct_weight
                        nearest_neighbor = n
                    elif direct_weight == max_direct_weight:
                        nearest_neighbor = min(nearest_neighbor, n)
                # print("nearest_neighbor direct weight:", max_direct_weight, "nodes in reachable:", len(biggest_reachable_comp_nodes))
                M[nearest_neighbor] = biggest_reachable_comp_weight   
                partition[nearest_neighbor] = biggest_reachable_comp_nodes.difference(set([nearest_neighbor]))
                assert nearest_neighbor in biggest_reachable_comp_nodes

            G_transpose.remove_nodes_from(biggest_reachable_comp_nodes)
            subgraph_set = subgraph_set -  biggest_reachable_comp_nodes
            nr_consensus += 1

    if params.verbose:
        print("NR CONSENSUS:", nr_consensus)
        print("NR nearest_neighbors:", len(M), len(partition))
        print("partition sizes(identical strings counted once): ", sorted([len(partition[p]) +1 for p in  partition], reverse = True))

    total_strings_in_partition = sum([ len(partition[p]) +1 for p in  partition])
    partition_sequences = set()
    for m in partition:
        partition_sequences.add(m)
        # print("partition size:", len(partition[m]))
        # print(len(m))
        for s in partition[m]:
            partition_sequences.add(s)

    assert unique_start_strings == partition_sequences
    assert total_strings_in_partition == len(unique_start_strings)

    return G_star, partition, M

def get_all_NN_old(batch_of_queries, global_index_in_matrix, start_index, seq_to_acc_list_sorted, neighbor_search_depth, ignore_ends_threshold):
    all_neighbors_graph = {}
    lower_target_edit_distances = {}
    max_variants = 10
    max_ed_allowed = max_variants + 2*ignore_ends_threshold
    # print("Processing global index:" , global_index_in_matrix)
    # error_types = {"D":0, "S": 0, "I": 0}
    for i in range(start_index, start_index + len(batch_of_queries)):
        if i % 500 == 0:
            print("processing ", i)
        seq1 = seq_to_acc_list_sorted[i][0]
        acc1 = seq_to_acc_list_sorted[i][1]
        all_neighbors_graph[acc1] = {}
        stop_up = False
        stop_down = False
        j = 1
        while True:
        # for j in range(1,len(seq_to_acc_list_sorted)):
            if i - j < 0:
                stop_down = True
            if i + j >= len(seq_to_acc_list_sorted):
                stop_up = True

            if not stop_down:
                seq2 = seq_to_acc_list_sorted[i - j][0]
                acc2 = seq_to_acc_list_sorted[i - j][1]  

                if math.fabs(len(seq1) - len(seq2)) > max_variants + 2*ignore_ends_threshold:
                    stop_down = True

            if not stop_up:
                seq3 = seq_to_acc_list_sorted[i + j][0]
                acc3 = seq_to_acc_list_sorted[i + j][1]  

                if math.fabs(len(seq1) - len(seq3)) > max_variants + 2*ignore_ends_threshold:
                    stop_up = True

            if not stop_down:
                result = edlib.align(seq1, seq2, mode="NW", k=max_ed_allowed) # , task="path")
                ed = result["editDistance"]
                if 0 <= ed : #implies its smaller or equal to max_ed_allowed
                    all_neighbors_graph[acc1][acc2] = (seq1, seq2, ed)

            if not stop_up:
                result = edlib.align(seq1, seq3, mode="NW", k=max_ed_allowed) # , task="path")
                ed = result["editDistance"]
                if 0 <= ed : #implies its smaller or equal to max_ed_allowed
                    all_neighbors_graph[acc1][acc3] = (seq1, seq3, ed)
            
            if stop_down and stop_up:
                break

            if j >= neighbor_search_depth:
                break
            j += 1

    return all_neighbors_graph

def edlib_traceback(x, y, mode="HW", task="path", k=1, end_threshold = 0):
    result = edlib.align(x, y, mode=mode, task=task, k=k)
    ed = result["editDistance"]
    cigar =  result["cigar"]

    if cigar:
        start, end =  result["locations"][0]
        start_offset = start
        end_offset = len(y) - (end + 1)
        # print(start_offset, end_offset)
        ed +=  max(0, start_offset - end_threshold)
        ed +=  max(0, end_offset - end_threshold)

        tuples = []
        result = re.split(r'[=DXSMI]+', cigar)
        i = 0
        for length in result[:-1]:
            i += len(length)
            type_ = cigar[i]
            i += 1
            tuples.append((length, type_ ))

        if tuples[-1][1] == "I":
            end_snippet_length = int(tuples[-1][0])
            ed -= min(int(end_snippet_length),  end_threshold) 
        if tuples[0][1] == "I":
            begin_snippet_length = int(tuples[0][0])
            ed -= min(int(begin_snippet_length),  end_threshold) 
    return ed

def get_all_NN(batch_of_queries, global_index_in_matrix, start_index, seq_to_acc_list_sorted, neighbor_search_depth, ignore_ends_threshold):
    all_neighbors_graph = {}
    lower_target_edit_distances = {}
    max_variants = 10
    max_ed_allowed = max_variants + ignore_ends_threshold
    # print("Processing global index:" , global_index_in_matrix)
    # error_types = {"D":0, "S": 0, "I": 0}
    for i in range(start_index, start_index + len(batch_of_queries)):
        if i % 500 == 0:
            print("processing ", i)
        seq1 = seq_to_acc_list_sorted[i][0]
        acc1 = seq_to_acc_list_sorted[i][1]
        all_neighbors_graph[acc1] = {}
        stop_up = False
        stop_down = False
        j = 1
        while True:
        # for j in range(1,len(seq_to_acc_list_sorted)):
            if i - j < 0:
                stop_down = True
            if i + j >= len(seq_to_acc_list_sorted):
                stop_up = True

            if not stop_down:
                seq2 = seq_to_acc_list_sorted[i - j][0]
                acc2 = seq_to_acc_list_sorted[i - j][1]  

                if math.fabs(len(seq1) - len(seq2)) > max_variants + 2*ignore_ends_threshold:
                    stop_down = True

            if not stop_up:
                seq3 = seq_to_acc_list_sorted[i + j][0]
                acc3 = seq_to_acc_list_sorted[i + j][1]  

                if math.fabs(len(seq1) - len(seq3)) > max_variants + 2*ignore_ends_threshold:
                    stop_up = True

            if not stop_down:
                ed = edlib_traceback(seq1, seq2, mode="HW", k=max_ed_allowed, task="path", end_threshold = ignore_ends_threshold)
                # result = edlib.align(seq1, seq2, mode="HW", k=max_ed_allowed) # , task="path")
                # ed = result["editDistance"]
                if 0 <= ed <= max_variants:
                    all_neighbors_graph[acc1][acc2] = ed #(seq1, seq2, ed)

            if not stop_up:
                ed = edlib_traceback(seq1, seq3, mode="HW", k=max_ed_allowed, task="path", end_threshold = ignore_ends_threshold)
                # result = edlib.align(seq1, seq3, mode="HW", k=max_ed_allowed) # , task="path")
                # ed = result["editDistance"]
                if 0 <= ed  <= max_variants: 
                    all_neighbors_graph[acc1][acc3] = ed #(seq1, seq3, ed)
            
            if stop_down and stop_up:
                break

            if j >= neighbor_search_depth:
                break
            j += 1

    return all_neighbors_graph

def get_all_NN_helper(arguments):
    args, kwargs = arguments
    return get_all_NN(*args, **kwargs)

def get_all_NN_under_ignored_edge_ends(seq_to_acc_list_sorted, params):
    if params.nr_cores == 1:
        all_neighbors_graph = get_all_NN(seq_to_acc_list_sorted, 0, 0, seq_to_acc_list_sorted, params.neighbor_search_depth, params.ignore_ends_len)
        # all_neighbors_graph_old = get_all_NN_old(seq_to_acc_list_sorted, 0, 0, seq_to_acc_list_sorted, params.neighbor_search_depth, params.ignore_ends_len)
        # print("TOTAL EDGES G_ALL OLD edlib TMP:", len([1 for s in all_neighbors_graph_old for t in all_neighbors_graph_old[s]]))
        # all_old = set([(s,t) for s in all_neighbors_graph_old for t in all_neighbors_graph_old[s]])
        # print("TOTAL EDGES G_ALL edlib TMP:", len([1 for s in all_neighbors_graph for t in all_neighbors_graph[s]]))
        # all_ = set([(s,t) for s in all_neighbors_graph for t in all_neighbors_graph[s]] + [(t,s) for s in all_neighbors_graph for t in all_neighbors_graph[s]])
        # for s,t in  all_old - all_:
        #     print("not in old:", s,t)
        #     # print(edlib.align())
        #     res = edlib.align(all_neighbors_graph_old[s][t][0], all_neighbors_graph_old[s][t][1], mode="HW", task="path")
        #     print(res["cigar"], res["editDistance"], res["locations"], len(all_neighbors_graph_old[s][t][1]) - res["locations"][0][1] - 1, res["locations"][0][0] )
        #     res = edlib.align(all_neighbors_graph_old[s][t][1], all_neighbors_graph_old[s][t][0], mode="HW", task="path")
        #     print(res["cigar"], res["editDistance"] )
        #     res = edlib.align(all_neighbors_graph_old[s][t][0], all_neighbors_graph_old[s][t][1], mode="NW", task="path")
        #     print(res["cigar"], res["editDistance"])
        #     print(all_neighbors_graph_old[s][t][2])
        #     (s1, s2, (s1_alignment, s2_alignment, (matches, mismatches, indels)) ) = parasail_alignment(all_neighbors_graph_old[s][t][0], all_neighbors_graph_old[s][t][1], 0, 0, opening_penalty = 3)
        #     print(s1_alignment)
        #     print(s2_alignment)
        # implement check here to se that all seqs got a nearest_neighbor, if not, print which noes that did not get a nearest_neighbor computed.!

    else:
        ####### parallelize alignment #########
        # pool = Pool(processes=mp.cpu_count())
        original_sigint_handler = signal.signal(signal.SIGINT, signal.SIG_IGN)
        signal.signal(signal.SIGINT, original_sigint_handler)
        pool = Pool(processes=mp.cpu_count())

        # here we split the input into chunks
        chunk_size = max(int(len(seq_to_acc_list_sorted) / (10*mp.cpu_count())), 20 )
        ref_seq_chunks = [ ( max(0, i - params.neighbor_search_depth -1), seq_to_acc_list_sorted[max(0, i - params.neighbor_search_depth -1) : i + chunk_size + params.neighbor_search_depth +1 ]) for i in range(0, len(seq_to_acc_list_sorted), chunk_size) ]
        chunks = [(i, seq_to_acc_list_sorted[i:i + chunk_size]) for i in range(0, len(seq_to_acc_list_sorted), chunk_size)] 

        if params.verbose:
            write_output.logger(str([j for j, ch in ref_seq_chunks]), params.develop_logfile, timestamp=False)
            write_output.logger("reference chunks:" + str([len(ch) for j,ch in ref_seq_chunks]), params.develop_logfile, timestamp=False)
            # print([j for j, ch in ref_seq_chunks])
            # print("reference chunks:", [len(ch) for j,ch in ref_seq_chunks])
            write_output.logger(str([i for i,ch in chunks]), params.develop_logfile, timestamp=False)
            write_output.logger("query chunks:" + str([len(ch) for i,ch in chunks]), params.develop_logfile, timestamp=False)

            print([i for i,ch in chunks])
            print("query chunks:", [len(ch) for i,ch in chunks])
        # get nearest_neighbors takes thre sub containers: 
        #  chunk - a container with (sequences, accesions)-tuples to be aligned (queries)
        #  ref_seq_chunks - a container with (sequences, accesions)-tuples to be aligned to (references)
        #  already_converged_chunks - a set of query sequences that has already converged 

        try:
            res = pool.map_async(get_all_NN_helper, [ ((chunks[i][1],  chunks[i][0], chunks[i][0] - ref_seq_chunks[i][0], ref_seq_chunks[i][1], params.neighbor_search_depth, params.ignore_ends_len), {}) for i in range(len(chunks))] )
            best_edit_distances_results =res.get(999999999) # Without the timeout this blocking call ignores all signals.
        except KeyboardInterrupt:
            print("Caught KeyboardInterrupt, terminating workers")
            pool.terminate()
            sys.exit()
        else:
            # print("Normal termination")
            pool.close()
        pool.join()
        all_neighbors_graph = {}
        for sub_graph in best_edit_distances_results:
            for seq in sub_graph:
                assert seq not in all_neighbors_graph
            all_neighbors_graph.update(sub_graph)

    return all_neighbors_graph


def get_NN_graph_ignored_ends_edlib(candidate_transcripts, args):
    seq_to_acc = {seq: acc for (acc, seq) in candidate_transcripts.items()}
    seq_to_acc_list = list(seq_to_acc.items())
    seq_to_acc_list_sorted = sorted(seq_to_acc_list, key= lambda x: len(x[0]))
    all_neighbors_graph = get_all_NN_under_ignored_edge_ends(seq_to_acc_list_sorted, args)
    # edges needs to be symmetric
    for c1 in all_neighbors_graph:
        for c2 in all_neighbors_graph[c1]:
            ed = all_neighbors_graph[c1][c2]
            if c1 not in all_neighbors_graph[c2]:
                all_neighbors_graph[c2][c1] = ed
            else:
                min_ed = min(all_neighbors_graph[c1][c2], all_neighbors_graph[c2][c1])
                all_neighbors_graph[c2][c1] = min_ed

                # if all_neighbors_graph[c1][c2][2] < all_neighbors_graph[c2][c1][2]:
                #     all_neighbors_graph[c2][c1] = (all_neighbors_graph[c2][c1][0],all_neighbors_graph[c2][c1][1],all_neighbors_graph[c1][c2][2])

                # if all_neighbors_graph[c2][c1] != all_neighbors_graph[c1][c2]:
                #     print(all_neighbors_graph[c2][c1][2], all_neighbors_graph[c1][c2][2])
                #     res = edlib.align(all_neighbors_graph[c1][c2][0], all_neighbors_graph[c1][c2][1], mode="HW", task="path")
                #     print(res["cigar"], res["editDistance"], res["locations"], len(all_neighbors_graph[c1][c2][1]) - res["locations"][0][1] - 1, res["locations"][0][0] )
                #     res = edlib.align(all_neighbors_graph[c1][c2][1], all_neighbors_graph[c1][c2][0], mode="HW", task="path")
                #     print(res["cigar"], res["editDistance"], res["locations"], len(all_neighbors_graph[c1][c2][0]) - res["locations"][0][1] - 1, res["locations"][0][0] )

                #     (s1, s2, (s1_alignment, s2_alignment, (matches, mismatches, indels)) ) = parasail_alignment(all_neighbors_graph[c1][c2][0], all_neighbors_graph[c1][c2][1], 0, 0, opening_penalty = 3, gap_ext = 1)
                #     print(s1_alignment)
                #     print(s2_alignment)


    assert len(candidate_transcripts) == len(all_neighbors_graph)
    return all_neighbors_graph



def get_nearest_neighbors_graph_under_ignored_ends(candidate_transcripts, args):
    seq_to_acc = {seq: acc for (acc, seq) in candidate_transcripts.items()}
    seq_to_acc_list = list(seq_to_acc.items())
    seq_to_acc_list_sorted = sorted(seq_to_acc_list, key= lambda x: len(x[0]))
    nearest_neighbor_graph = get_nearest_neighbors_under_ignored_edge_ends(seq_to_acc_list_sorted, args)
    nearest_neighbor_graph_parasail = get_nearest_neighbors_under_ignored_edge_ends_parasail(seq_to_acc_list_sorted, args)
    print("edges before:", len([1 for s in nearest_neighbor_graph for r in nearest_neighbor_graph[s] ]))
    print("edges parasail:", len([1 for s in nearest_neighbor_graph_parasail for r in nearest_neighbor_graph_parasail[s] ]))

    print("EDs before:", sum([nearest_neighbor_graph[s][r] for s in nearest_neighbor_graph for r in nearest_neighbor_graph[s] ]))
    print("EDs parasail:", sum([nearest_neighbor_graph_parasail[s][r] for s in nearest_neighbor_graph_parasail for r in nearest_neighbor_graph_parasail[s] ]))


    for acc1 in list(nearest_neighbor_graph_parasail.keys()):
        for acc2 in list(nearest_neighbor_graph_parasail[acc1].keys()):
            ed = nearest_neighbor_graph_parasail[acc1][acc2]
            if ed > 10:
                del nearest_neighbor_graph_parasail[acc1][acc2]
                if args.verbose:
                    print("had ed > 10 statistical test", acc1, acc2)


    for acc1 in  nearest_neighbor_graph_parasail:
        if args.verbose:
            for acc2 in nearest_neighbor_graph_parasail[acc1]:
                if nearest_neighbor_graph_parasail[acc1][acc2] > 0:
                    print("To be tested:", acc1, acc2,  nearest_neighbor_graph_parasail[acc1][acc2])
    

    assert len(candidate_transcripts) == len(nearest_neighbor_graph_parasail)
    return nearest_neighbor_graph_parasail

def get_invariants_under_ignored_edge_ends(seq_to_acc_list_sorted, params):
    if params.nr_cores == 1:
        best_edit_distances = get_nearest_neighbors(seq_to_acc_list_sorted, 0, 0, seq_to_acc_list_sorted, params.neighbor_search_depth, params.ignore_ends_len)

        # implement check here to se that all seqs got a nearest_neighbor, if not, print which noes that did not get a nearest_neighbor computed.!

    else:
        ####### parallelize alignment #########
        # pool = Pool(processes=mp.cpu_count())
        original_sigint_handler = signal.signal(signal.SIGINT, signal.SIG_IGN)
        signal.signal(signal.SIGINT, original_sigint_handler)
        pool = Pool(processes=mp.cpu_count())

        # here we split the input into chunks
        chunk_size = max(int(len(seq_to_acc_list_sorted) / (10*mp.cpu_count())), 20 )
        ref_seq_chunks = [ ( max(0, i - params.neighbor_search_depth -1), seq_to_acc_list_sorted[max(0, i - params.neighbor_search_depth -1) : i + chunk_size + params.neighbor_search_depth +1 ]) for i in range(0, len(seq_to_acc_list_sorted), chunk_size) ]
        chunks = [(i, seq_to_acc_list_sorted[i:i + chunk_size]) for i in range(0, len(seq_to_acc_list_sorted), chunk_size)] 

        if params.verbose:
            write_output.logger(str([j for j, ch in ref_seq_chunks]), params.develop_logfile, timestamp=False)
            write_output.logger("reference chunks:" + str([len(ch) for j,ch in ref_seq_chunks]), params.develop_logfile, timestamp=False)
            # print([j for j, ch in ref_seq_chunks])
            # print("reference chunks:", [len(ch) for j,ch in ref_seq_chunks])
            write_output.logger(str([i for i,ch in chunks]), params.develop_logfile, timestamp=False)
            write_output.logger("query chunks:" + str([len(ch) for i,ch in chunks]), params.develop_logfile, timestamp=False)

            print([i for i,ch in chunks])
            print("query chunks:", [len(ch) for i,ch in chunks])

        # get nearest_neighbors takes thre sub containers: 
        #  chunk - a container with (sequences, accesions)-tuples to be aligned (queries)
        #  ref_seq_chunks - a container with (sequences, accesions)-tuples to be aligned to (references)
        #  already_converged_chunks - a set of query sequences that has already converged 

        try:
            res = pool.map_async(get_nearest_neighbors_helper, [ ((chunks[i][1],  chunks[i][0], chunks[i][0] - ref_seq_chunks[i][0], ref_seq_chunks[i][1], params.neighbor_search_depth, params.ignore_ends_len), {}) for i in range(len(chunks))] )
            best_edit_distances_results =res.get(999999999) # Without the timeout this blocking call ignores all signals.
        except KeyboardInterrupt:
            print("Caught KeyboardInterrupt, terminating workers")
            pool.terminate()
            sys.exit()
        else:
            # print("Normal termination")
            pool.close()
        pool.join()
        best_edit_distances = {}
        for sub_graph in best_edit_distances_results:
            for seq in sub_graph:
                assert seq not in best_edit_distances
            best_edit_distances.update(sub_graph)
        
    # store only invariants here, i.e., edit distance 0 when ignoring ends!
    for acc1 in list(best_edit_distances.keys()):
        for acc2 in list(best_edit_distances[acc1].keys()):
            if best_edit_distances[acc1][acc2] != 0:
                del best_edit_distances[acc1][acc2]

    return best_edit_distances


def is_overlap(text1, text2, ignore_ends_threshold):  
  # Cache the text lengths to prevent multiple calls.  
  text1_length = len(text1)  
  text2_length = len(text2)  
  # Eliminate the null case.  
  if text1_length == 0 or text2_length == 0:  
    return 0  
  # Truncate the longer string.  
  if text1_length > text2_length:  
    text1 = text1[-text2_length:]  
  elif text1_length < text2_length:  
    text2 = text2[:text1_length]  
  # Quick check for the worst case.  
  if text1 == text2:  
    return min(text1_length, text2_length)  
   
  # Start by looking for a single character match  
  # and increase length until no match is found.  
  best = 0  
  length = 1  
  while True:  
    pattern = text1[-length:]  
    found = text2.find(pattern)  
    if found == -1: 
      prefix_offset =  text1_length - best
      suffix_offset =  text2_length - best
      if prefix_offset <= ignore_ends_threshold and suffix_offset <= ignore_ends_threshold:
        return True
      else:
        return False
    length += found  
    if text1[-length:] == text2[:length]:  
      best = length  
      length += 1 


def get_invariants_under_ignored_edge_ends_speed(candidate_transcripts, candidate_support, params):
    edge_invariant_threshold = params.ignore_ends_len
    G = nx.DiGraph()
    for acc in candidate_transcripts:
        deg = candidate_support[acc]
        G.add_node(acc, degree = deg)
    sorted_lenghts = sorted(candidate_transcripts.items(), key = lambda x: len(x[1]))
    for i, (acc1, seq1) in enumerate(sorted_lenghts):
        if i % 1000 == 0:
            print(i, "candidates processed") 
        for (acc2, seq2) in sorted_lenghts:
            if acc2 == acc1:
                continue

            if len(seq2) < len(seq1) - 2*edge_invariant_threshold:
                continue
            elif len(seq1) - 2*edge_invariant_threshold <= len(seq2) <= len(seq1): # is long enough to be merged
                   
                if seq2 in seq1: # is strict substring
                    start_offset = seq1.find(seq2)
                    end_offset = len(seq1) - (start_offset + len(seq2))
                    if start_offset <= edge_invariant_threshold and end_offset <= edge_invariant_threshold:
                        G.add_edge(acc2, acc1) # directed edge seq2 --> seq1
                        G.add_edge(acc1, acc2) # directed edge seq2 --> seq1
                
                else: # check if perfect overlap within ends threshold
                    bool_overlap1 = is_overlap(seq1, seq2, edge_invariant_threshold)
                    bool_overlap2 = is_overlap(seq2, seq1, edge_invariant_threshold)
                    if bool_overlap1 or bool_overlap2:
                        G.add_edge(acc2, acc1) 
                        G.add_edge(acc1, acc2) 
            else:
                break

    return G

# def get_invariants_under_ignored_edge_ends_OLD(candidate_transcripts, candidate_support, params):

#     seq_to_acc = {seq: acc for (acc, seq) in candidate_transcripts.items()}
#     seq_to_acc_list = list(seq_to_acc.items())
#     seq_to_acc_list_sorted = sorted(seq_to_acc_list, key= lambda x: len(x[0]))
#     invariant_graph = get_invariants_under_ignored_edge_ends(seq_to_acc_list_sorted, params)

#     G = nx.DiGraph()
#     for acc in candidate_transcripts:
#         deg = candidate_support[acc]
#         G.add_node(acc, degree = deg)
#     # add edges
#     for acc1 in  invariant_graph:
#         for acc2 in invariant_graph[acc1]:
#             G.add_edge(acc1, acc2)
#     return G



def collapse_candidates_under_ends_invariant(candidate_transcripts, candidate_support, params):
    print("Final candidates before edge invariants:", len(candidate_transcripts))

    # start = time()
    # G_old = get_invariants_under_ignored_edge_ends_OLD(candidate_transcripts, candidate_support, params)
    # elapsed = time() - start
    # print("INVARIANTS OLD:", elapsed)

    start = time()
    G = get_invariants_under_ignored_edge_ends_speed(candidate_transcripts, candidate_support, params)
    elapsed = time() - start
    print("INVARIANTS NEW:", elapsed)
    # print(len(list(G_old.nodes())), len(list(G.nodes())) )

    # print(len(list(G_old.edges())), len(list(G.edges())) )

    # for edge in list(G.edges()):
    #     if not G_old.has_edge(*edge):
    #         print(edge[0])
    #         print(candidate_transcripts[edge[0]])
    #         print(edge[1])
    #         print(candidate_transcripts[edge[1]])            
    #         print()
    # print(nx.is_isomorphic(G_old, G))

    # seq_to_acc = {seq: acc for (acc, seq) in candidate_transcripts.items()}
    # seq_to_acc_list = list(seq_to_acc.items())
    # seq_to_acc_list_sorted = sorted(seq_to_acc_list, key= lambda x: len(x[0]))
    # invariant_graph = get_invariants_under_ignored_edge_ends(seq_to_acc_list_sorted, params)

    # G = nx.DiGraph()
    # for acc in candidate_transcripts:
    #     deg = candidate_support[acc]
    #     G.add_node(acc, degree = deg)
    # # add edges
    # for acc1 in  invariant_graph:
    #     for acc2 in invariant_graph[acc1]:
    #         G.add_edge(acc1, acc2)
    G_tmp = copy.deepcopy(G)
    G_star, partition, M = partition_highest_reachable_with_edge_degrees(G_tmp, params)
    # _, partition_old, _ = partition_highest_reachable_with_edge_degrees(G_old, params)

    # print("Final candidates old:", len(partition_old))
    print("Final candidates:", len(partition))
    # sys.exit()

    # # SAM_file = minimap2_alignment_module.align(targets, queries, nr_cores)


    # sorted_lenghts = sorted(candidate_transcripts.items(), key = lambda x: len(x[1]))
    # for i, (acc1, seq1) in enumerate(sorted_lenghts):
    #     if i % 1000 == 0:
    #         print(i, "candidates processed") 
    #     for (acc2, seq2) in sorted_lenghts:
    #         if acc2 == acc1:
    #             continue

    #         if len(seq2) > len(seq1):
    #             break

    #         if len(seq2) >= len(seq1) - 2*params.ignore_ends_len: # is long enough to be merged
    #             if seq2 in seq1:
    #                 # has to be within ends varinat lenghth:
    #                 start_offset = seq1.find(seq2)
    #                 end_offset = len(seq1) - (start_offset + len(seq2))
    #                 if start_offset <= params.ignore_ends_len and end_offset <= params.ignore_ends_len:
    #                     #Make sure this doesnet crach out statistical test if two candidates differ only in ends!!!
    #                     G.add_edge(acc2, acc1) # directed edge seq2 --> seq1

    
    # # sort order: largest number of neigbors, then largest degree. If still tie, sort by smallest string
    # G_transpose = nx.reverse(G)
    # largest_nr_neighbors = sorted([ (len(list(G_transpose.neighbors(n))), G_transpose.node[n]["degree"], n) for n in  sorted(G_transpose.nodes())], key= lambda x: (-x[0], -x[1], x[2]) )
    
    # marked = set()
    # partition = {}
    # for nr_nbrs, deg, c in largest_nr_neighbors:
    #     if c not in marked:
    #         nbrs = [n for n in G_transpose.neighbors(c) if n not in marked]
    #         partition[c] = set(nbrs)
    #         marked.add(c)
    #         marked.update(nbrs)

    if params.verbose:
        for t in partition:
            print(t, partition[t])
    print("Final candidates after edge invariants:", len(partition))
    print()
    return partition



def main(args):
    candidate_transcripts = {acc: seq for (acc, seq) in  read_fasta(open(args.candidate_transcripts, 'r'))}
    candidate_support = {}
    for (acc, seq) in  read_fasta(open(args.candidate_transcripts, 'r')):
        supp = acc.split("_support_")[1]
        candidate_support[acc] = int(supp)
    
    # print("Number of consensus:", len(candidate_transcripts))
    seq_to_acc = {seq: acc for (acc, seq) in  read_fasta(open(args.candidate_transcripts, 'r'))}
    seq_to_acc_list = list(seq_to_acc.items())
    seq_to_acc_list_sorted = sorted(seq_to_acc_list, key= lambda x: len(x[0]))
    collapsed_candidate_transcripts =  { acc : seq for (seq, acc) in  seq_to_acc.items() }
    # print("Number of collapsed consensus:", len(collapsed_candidate_transcripts))
    assert len(collapsed_candidate_transcripts) == len(candidate_transcripts) # all transcripts should be unique at this point
    

    nearest_neighbor_graph = get_invariants_under_ignored_edge_ends(seq_to_acc_list_sorted, args)

    outfile = open(args.outfile, "w")
    edges = 0
    tot_ed = 0
    for acc1 in  nearest_neighbor_graph:
        seq1 = candidate_transcripts[acc1]
        for acc2 in nearest_neighbor_graph[acc1]:
            seq2 = candidate_transcripts[acc2]
            edges += 1
            tot_ed += nearest_neighbor_graph[acc1][acc2]
            outfile.write("{0}\t{1}\t{2}\t{3}\t{4}\n".format(acc1, candidate_support[acc1], acc2, candidate_support[acc2],  nearest_neighbor_graph[acc1][acc2]))

    print("Number of edges:", edges)
    print("Total edit distance:", tot_ed)
    if float(edges) > 0:
        print("Avg ed (ed/edges):", tot_ed/ float(edges))

    # convert nearest_neighbor graph to nx graph object
    G = nx.DiGraph()
    # add nodes
    for acc in candidate_transcripts:
        deg = candidate_support[acc]
        G.add_node(acc, degree = deg)
    # add edges
    for acc1 in  nearest_neighbor_graph:
        for acc2 in nearest_neighbor_graph[acc1]:
            G.add_edge(acc1, acc2)
    G_star, partition, M = partition_highest_reachable_with_edge_degrees(G, params)

    print("candidates after edge invariants:", len(partition))
    # for t in partition:
    #     print(t)
    #     print(partition[t])
    #     print()




if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Print nearest_neighbor graph allowing for mismatches in ends.")
    parser.add_argument('candidate_transcripts', type=str, help='Path to the consensus fasta file')
    parser.add_argument('outfile', type=str, help='Outfile of results')
    parser.add_argument('--ignore_ends_len', type=int, default=15, help='Number of bp to ignore in ends. If two candidates are identical expept in ends of this size, they are collapses and the longest common substing is chosen to represent them. In statistical test step, nearest_neighbors are found based on ignoring the ends of this size. Also indels in ends will not be tested. [default ignore_ends_len=15].')
    parser.add_argument('--single_core', dest='single_core', action='store_true', help='Force working on single core. ')
    parser.add_argument('--neighbor_search_depth', type=int, default=2**32, help='Maximum number of pairwise alignments in search matrix to find nearest_neighbor. [default =2**32]')

    args = parser.parse_args()

    main(args)
