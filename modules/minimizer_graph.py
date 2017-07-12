from __future__ import print_function
import os,sys
import argparse
import re
import math
import numpy as np
from Bio.SubsMat import MatrixInfo as matlist
import pandas as pd
import string
import fractions
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import seaborn as sns

from collections import defaultdict
import edlib

import signal
from multiprocessing import Pool
import multiprocessing as mp
import sys

import re

import edlib

def get_minimizers_helper(arguments):
    args, kwargs = arguments
    return get_minimizers(*args, **kwargs)

def get_exact_minimizer_graph(seq_to_acc_list_sorted, has_converged, params):
    if params.single_core:
        best_edit_distances = get_minimizers(seq_to_acc_list_sorted, 0, 0, seq_to_acc_list_sorted, has_converged, params.minimizer_search_depth)

        # implement check here to se that all seqs got a minimizer, if not, print which noes that did not get a minimizer computed.!

    else:
        ####### parallelize alignment #########
        # pool = Pool(processes=mp.cpu_count())
        original_sigint_handler = signal.signal(signal.SIGINT, signal.SIG_IGN)
        signal.signal(signal.SIGINT, original_sigint_handler)
        pool = Pool(processes=mp.cpu_count())

        # here we split the input into chunks
        chunk_size = max(int(len(seq_to_acc_list_sorted) / (10*mp.cpu_count())), 20 )
        ref_seq_chunks = [ ( max(0, i - params.minimizer_search_depth -1), seq_to_acc_list_sorted[max(0, i - params.minimizer_search_depth -1) : i + chunk_size + params.minimizer_search_depth +1 ]) for i in range(0, len(seq_to_acc_list_sorted), chunk_size) ]
        print([j for j, ch in ref_seq_chunks])
        print("reference chunks:", [len(ch) for j,ch in ref_seq_chunks])

        chunks = [(i, seq_to_acc_list_sorted[i:i + chunk_size]) for i in range(0, len(seq_to_acc_list_sorted), chunk_size)] 
        print([i for i,ch in chunks])
        print("query chunks:", [len(ch) for i,ch in chunks])
        # sys.exit()
        already_converged_chunks = []
        for i, chunk in chunks:
            already_converged_chunk = set()
            for seq, acc in chunk:
                if seq in has_converged:
                    already_converged_chunk.add(seq)
            already_converged_chunks.append(already_converged_chunk)
        print("already converged chunks:", [len(ch) for ch in already_converged_chunks])

        # get minimizers takes thre sub containers: 
        #  chunk - a container with (sequences, accesions)-tuples to be aligned (queries)
        #  ref_seq_chunks - a container with (sequences, accesions)-tuples to be aligned to (references)
        #  already_converged_chunks - a set of query sequences that has already converged 

        try:
            res = pool.map_async(get_minimizers_helper, [ ((chunks[i][1],  chunks[i][0], chunks[i][0] - ref_seq_chunks[i][0], ref_seq_chunks[i][1], already_converged_chunks[i], params.minimizer_search_depth), {}) for i in range(len(chunks))] )
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

def histogram(data, args, name='histogram.png', x='x-axis', y='y-axis', x_cutoff=None, title=None):
    if x_cutoff: 
        plt.hist(data, range=[0, x_cutoff], bins = 100)
    else:
        plt.hist(data, bins = 100)
    plt.xlabel(x)
    plt.ylabel(y)
    if title:
        plt.title(title)

    plt.savefig(os.path.join(args.outfolder, name))
    plt.clf()


def collapse(seq):
    seen = set()
    seen_add = seen.add
    return [x for x in seq if not (x in seen or seen_add(x))]




def edlib_ed(x, y, mode="NW", task="distance", k=1):
    result = edlib.align(x, y, mode=mode, task=task, k=k)
    ed = result["editDistance"]
    return ed

def edlib_traceback(x, y, mode="NW", task="path", k=1):
    result = edlib.align(x, y, mode=mode, task=task, k=k)
    ed = result["editDistance"]
    locations =  result["locations"]
    cigar =  result["cigar"]
    return ed, locations, cigar


def get_minimizers(batch_of_queries, global_index_in_matrix, start_index, seq_to_acc_list_sorted, has_converged, minimizer_search_depth):
    best_edit_distances = {}
    lower_target_edit_distances = {}
    print("Processing global index:" , global_index_in_matrix)
    # error_types = {"D":0, "S": 0, "I": 0}
    for i in range(start_index, start_index + len(batch_of_queries)):
        if i % 500 == 0:
            print("processing ", i)
        seq1 = seq_to_acc_list_sorted[i][0]
        acc1 = seq_to_acc_list_sorted[i][1]
        best_edit_distances[acc1] = {}
        if seq1 in has_converged:
            # print("ctd here")
            continue

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

                if math.fabs(len(seq1) - len(seq2)) > best_ed:
                    stop_down = True

            if not stop_up:
                seq3 = seq_to_acc_list_sorted[i + j][0]
                acc3 = seq_to_acc_list_sorted[i + j][1]  

                if math.fabs(len(seq1) - len(seq3)) > best_ed:
                    stop_up = True

            if not stop_down:
                edit_distance = edlib_ed(seq1, seq2, mode="NW", task="distance", k=best_ed)
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
                edit_distance = edlib_ed(seq1, seq3, mode="NW", task="distance", k=best_ed)
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

            if j >= minimizer_search_depth:
                break
            j += 1

        # print("best ed:", best_ed)
        # if best_ed > 100:
        #     print(best_ed, "for seq with length", len(seq1), seq1)
    return best_edit_distances

def compute_2set_minimizer_graph(X, C, params):
    seq_to_acc_queries = [(seq, acc) for (acc, seq) in X.items()] #{seq: acc for (acc, seq) in  read_fasta(open(args.consensus_transcripts, 'r'))}
    # seq_to_acc_list_queries = list(seq_to_acc_queries.items())

    seq_to_acc_targets = [(seq, acc) for (acc, seq) in  C.items()] #{seq: acc for (acc, seq) in  read_fasta(open(args.target_transcripts, 'r'))}
    # seq_to_acc_list_targets = list(seq_to_acc_targets.items())
    
    seq_to_acc_list_sorted_all = sorted(seq_to_acc_queries + seq_to_acc_targets, key= lambda x: len(x[0]))

    minimizer_graph_x_to_c = get_exact_minimizer_graph_2set(seq_to_acc_list_sorted_all, set(C.keys()), params)

    # TAKE CARE OF UNALIGNED READS HERE?

    edges = 0
    tot_ed = 0
    edit_hist =[]
    neighbors = []
    for x in  minimizer_graph_x_to_c:
        for c in minimizer_graph_x_to_c[x]:
            edges += 1
            tot_ed += minimizer_graph_x_to_c[x][c]
            edit_hist.append(minimizer_graph_x_to_c[x][c])

        neighbors.append(len(minimizer_graph_x_to_c[x]))

    print("Number of edges:", edges)
    print("Total edit distance:", tot_ed)
    print("Avg ed (ed/edges):", tot_ed/ float(edges))
    # histogram(edit_hist, args, name='edit_distances.png', x='x-axis', y='y-axis', x_cutoff=100, title="Edit distances in minimizer graph")
    # histogram(edit_hist, args, name='edit_distances_zoomed.png', x='x-axis', y='y-axis', x_cutoff=5, title="Edit distances in minimizer graph")
    # histogram(neighbors, args, name='neighbours.png', x='x-axis', y='y-axis', title="Number of neighbours in minimizer graph")
    # histogram(neighbors, args, name='neighbours_zoomed.png', x='x-axis', y='y-axis', x_cutoff=20, title="Number of neighbours in minimizer graph")
    return minimizer_graph_x_to_c


def compute_minimizer_graph(S, has_converged, params):
    """
        strings S are all unique here.
    """
    # consensus_transcripts = {acc: seq for (acc, seq) in  read_fasta(open(params.consensus_transcripts, 'r'))}
    # print("Number of consensus:", len(consensus_transcripts))
    seq_to_acc = { seq : acc for (acc, seq) in S.items() }

    seq_to_acc_list = list(seq_to_acc.items())
    seq_to_acc_list_sorted = sorted(seq_to_acc_list, key= lambda x: len(x[0]))
    collapsed_consensus_transcripts =  { acc : seq for (seq, acc) in  seq_to_acc.items() }
    print("Number of collapsed consensus:", len(collapsed_consensus_transcripts))
    minimizer_graph = get_exact_minimizer_graph(seq_to_acc_list_sorted, has_converged, params)

    s1 = set()
    for acc1 in minimizer_graph:
        s1.add(S[acc1])

    s2 = set([seq for seq in seq_to_acc] )
    isolated = s2.difference(s1)
    print("isolated:", len(isolated))
    # print("isolated:", isolated)

    edges = 0
    tot_ed = 0
    edit_hist =[]
    neighbors = []
    for r1 in  minimizer_graph:
        for r2 in minimizer_graph[r1]:
            edges += 1
            tot_ed += minimizer_graph[r1][r2]
            edit_hist.append(minimizer_graph[r1][r2])

        neighbors.append(len(minimizer_graph[r1]))

    print("Number of edges:", edges)
    print("Total edit distance:", tot_ed)
    print("Avg ed (ed/edges):", tot_ed/ float(edges))
    histogram(edit_hist, params, name='edit_distances.png', x='x-axis', y='y-axis', x_cutoff=100, title="Edit distances in minimizer graph")
    histogram(neighbors, params, name='neighbours.png', x='x-axis', y='y-axis', title="Number of neighbours in minimizer graph")
    histogram(neighbors, params, name='neighbours_zoomed.png', x='x-axis', y='y-axis', x_cutoff=20, title="Number of neighbours in minimizer graph")

    return minimizer_graph, isolated



def get_exact_minimizer_graph_2set(seq_to_acc_list_sorted_all, target_accessions, params):

    if params.single_core:
        best_edit_distances = get_minimizers_2set(seq_to_acc_list_sorted_all, 0, seq_to_acc_list_sorted_all, target_accessions, params.minimizer_search_depth)

        # implement check here to se that all seqs got a minimizer, if not, print which noes that did not get a minimizer computed.!

    else:
        ####### parallelize alignment #########
        # pool = Pool(processes=mp.cpu_count())
        original_sigint_handler = signal.signal(signal.SIGINT, signal.SIG_IGN)
        signal.signal(signal.SIGINT, original_sigint_handler)
        pool = Pool(processes=mp.cpu_count())
        chunk_size = max(int(len(seq_to_acc_list_sorted_all) / (10*mp.cpu_count())), 20 )
        chunks = [(i, seq_to_acc_list_sorted_all[i:i + chunk_size]) for i in range(0, len(seq_to_acc_list_sorted_all), chunk_size)] 
        print([i for i in range(0, len(seq_to_acc_list_sorted_all), chunk_size)])
        print([len(ch) for i,ch in chunks])
        try:
            res = pool.map_async(get_minimizers_2set_helper, [ ((chunk, i , seq_to_acc_list_sorted_all, target_accessions, params.minimizer_search_depth), {}) for i,chunk in chunks] )
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


def get_minimizers_2set_helper(arguments):
    args, kwargs = arguments
    return get_minimizers_2set(*args, **kwargs)

def get_minimizers_2set(batch, start_index, seq_to_acc_list_sorted, target_accessions, minimizer_search_depth):
    best_edit_distances = {}
    error_types = {"D":0, "S": 0, "I": 0}
    for i in range(start_index, start_index + len(batch)):
        if i % 50 == 0:
            print("processing ", i)

        seq1 = seq_to_acc_list_sorted[i][0]
        acc1 = seq_to_acc_list_sorted[i][1]
        if acc1 in target_accessions:
            continue

        # reach here and we have a query read to find best alignment for
            
        best_edit_distances[acc1] = {}
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

                if math.fabs(len(seq1) - len(seq2)) > best_ed:
                    stop_down = True

            if not stop_up:
                seq3 = seq_to_acc_list_sorted[i + j][0]
                acc3 = seq_to_acc_list_sorted[i + j][1]  

                if math.fabs(len(seq1) - len(seq3)) > best_ed:
                    stop_up = True

            if not stop_down and acc2 in target_accessions:
                # if seq1 == seq2:
                #     print("ID:", acc1, acc2)
                edit_distance = edlib_ed(seq1, seq2, mode="NW", task="distance", k=best_ed)
                if 0 <= edit_distance < best_ed:
                    best_ed = edit_distance
                    best_edit_distances[acc1] = {}
                    best_edit_distances[acc1][acc2] = best_ed
                    # lower_target_edit_distances[acc2] = best_ed 

                elif edit_distance == best_ed:
                    best_edit_distances[acc1][acc2] = best_ed

            if not stop_up and acc3 in target_accessions:
                # if seq1 == seq3:
                #     print("ID:", acc1, acc3)

                edit_distance = edlib_ed(seq1, seq3, mode="NW", task="distance", k=best_ed)
                if 0 <= edit_distance < best_ed:
                    best_ed = edit_distance
                    best_edit_distances[acc1] = {}
                    best_edit_distances[acc1][acc3] = best_ed
                    # lower_target_edit_distances[acc3] = best_ed 

                elif edit_distance == best_ed:
                    best_edit_distances[acc1][acc3] = best_ed
 
            if stop_down and stop_up:
                break

            if j >= minimizer_search_depth:
                break

            j += 1
        # print("best ed:", best_ed)
        # if best_ed > 100:
        #     print(best_ed, "for seq with length", len(seq1), seq1)

    return best_edit_distances




def main(args):
    consensus_transcripts = {acc: seq for (acc, seq) in  read_fasta(open(args.consensus_transcripts, 'r'))}
    print("Number of consensus:", len(consensus_transcripts))
    seq_to_acc = {seq: acc for (acc, seq) in  read_fasta(open(args.consensus_transcripts, 'r'))}
    seq_to_acc_list = list(seq_to_acc.items())
    seq_to_acc_list_sorted = sorted(seq_to_acc_list, key= lambda x: len(x[0]))
    collapsed_consensus_transcripts =  { acc : seq for (seq, acc) in  seq_to_acc.items() }
    print("Number of collapsed consensus:", len(collapsed_consensus_transcripts))
    minimizer_graph = get_exact_minimizer_graph(seq_to_acc_list_sorted, params)

    s1 = set( [ collapsed_consensus_transcripts[acc2] for acc1 in minimizer_graph for acc2 in minimizer_graph[acc1] ])
    s2 = set([seq for seq in seq_to_acc] )
    isolated = s2.difference(s1)
    print("isolated:", len(isolated))

    edges = 0
    tot_ed = 0
    edit_hist =[]
    neighbors = []
    for r1 in  minimizer_graph:
        for r2 in minimizer_graph[r1]:
            edges += 1
            tot_ed += minimizer_graph[r1][r2]
            edit_hist.append(minimizer_graph[r1][r2])

        neighbors.append(len(minimizer_graph[r1]))

    print("Number of edges:", edges)
    print("Total edit distance:", tot_ed)
    print("Avg ed (ed/edges):", tot_ed/ float(edges))
    histogram(edit_hist, args, name='edit_distances.png', x='x-axis', y='y-axis', x_cutoff=100, title="Edit distances in minimizer graph")
    histogram(neighbors, args, name='neighbours.png', x='x-axis', y='y-axis', title="Number of neighbours in minimizer graph")
    histogram(neighbors, args, name='neighbours_zoomed.png', x='x-axis', y='y-axis', x_cutoff=20, title="Number of neighbours in minimizer graph")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Evaluate pacbio IsoSeq transcripts.")
    parser.add_argument('consensus_transcripts', type=str, help='Path to the consensus fasta file')
    parser.add_argument('outfolder', type=str, help='Output path of results')
    parser.add_argument('--single_core', dest='single_core', action='store_true', help='Force working on single core. ')
    args = parser.parse_args()

    outfolder = args.outfolder
    if not os.path.exists(outfolder):
        os.makedirs(outfolder)
    main(args)