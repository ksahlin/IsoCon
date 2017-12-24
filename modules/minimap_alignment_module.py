
import os
import sys
# import unittest
import tempfile
import subprocess
import gzip
import heapq
from collections import defaultdict
import math

import multiprocessing as mp

def paf_to_best_matches_2set(paf_file_path):
    """
        input: a PAF file
        output: A matches data structure. It is a didcionary on the following form:
                highest_paf_scores = {read: [target1, target2, ...] }
                targets as keys and a list sequences where the sequences are the X best target hits.
                for each sequence found in PAF file, store X best matches based on paf score
    """

    highest_paf_scores = defaultdict(list)


    try:
        file_object = gzip.open(paf_file_path)
        file_object.readline()
        file_object.seek(0)
    except IOError:
        file_object = open(paf_file_path)


    with file_object as paf_file:
        for line in paf_file:
            row_info = line.strip().split()
            q_acc = str(row_info[0]) #.decode('ascii'))
            q_len = int(row_info[1])
            t_acc = str(row_info[5]) #.decode('ascii'))
            t_len = int(row_info[6])
            nr_match = int(row_info[9])            
            nr_match_w_gap = int(row_info[10])            
            if q_acc == t_acc:
                # print("SELF MAPPING DETECTED")
                # print(q_acc, q_len, row_info)
                # print(t_acc, t_len, row_info)
                continue

            n_min = int(row_info[12].split(":")[-1])
            #   \frac{\min \{ |q|, |t| \} } {\max \{ |q|,|t| \} }  n^{\text{nearest_neighbors}}
            ## paf_similarity_score = math.sqrt(n_min) * min(q_len, t_len)/float(max(q_len, t_len))
            paf_similarity_score = n_min - ( max(q_len, t_len) - min(q_len, t_len))

            # if original reads noisy, this is the only approximate score that makes sense 
            # because pacbio has more insertions than deletions, thus penalizing with
            #  - ( max(q_len, t_len) - min(q_len, t_len)) would make error prone longer reads 
            # map to transcripts with an extra exon if all else identical 
            # paf_similarity_score = n_min  
            ## paf_similarity_score = nr_match - ( nr_match_w_gap - nr_match )

            if len(highest_paf_scores[q_acc]) >= 20:
                # current alignment is better than at least one of previous scores, remove the worst one so far
                if paf_similarity_score > highest_paf_scores[q_acc][0][0]: 
                    paf_score, t_acc_out = heapq.heappushpop(highest_paf_scores[q_acc], (paf_similarity_score, t_acc) )
            else:
                heapq.heappush(highest_paf_scores[q_acc], (paf_similarity_score, t_acc))


    return highest_paf_scores


def paf_to_best_matches(paf_files, acc_to_strings):
    """
        input: a list of PAF files
        output: A matches data structure. It is a didcionary on the following form:
                matches = {string1: [string2, string3, ...] }
                sequences as keys and a list sequences where the sequences are the X best hits to the sequence.
                for each sequence found in PAF file, store X best matches based on paf score
    """

    matches = {}
    highest_paf_scores = defaultdict(list)


    for paf_file_path in paf_files:
        print(paf_file_path)
        try:
            file_object = gzip.open(paf_file_path)
            file_object.readline()
            file_object.seek(0)
        except IOError:
            file_object = open(paf_file_path)


        with file_object as paf_file:
            for line in paf_file:
                row_info = line.strip().split()
                q_acc = str(row_info[0]) #.decode('ascii')
                q_len = int(row_info[1])
                t_acc = str(row_info[5]) #.decode('ascii')
                t_len = int(row_info[6])
                nr_match = int(row_info[9])            
                nr_match_w_gap = int(row_info[10])            

                if q_acc == t_acc:
                    # print("SELF MAPPING DETECTED")
                    # print(q_acc, q_len, row_info)
                    # print(t_acc, t_len, row_info)
                    continue

                n_min = int(row_info[12].split(":")[-1])
                #   \frac{\min \{ |q|, |t| \} } {\max \{ |q|,|t| \} }  n^{\text{nearest_neighbors}}
                # paf_similarity_score = math.sqrt(n_min) * min(q_len, t_len)/float(max(q_len, t_len))
                # paf_similarity_score = nr_match - ( nr_match_w_gap - nr_match )
                paf_similarity_score = n_min - ( max(q_len, t_len) - min(q_len, t_len))

                if len(highest_paf_scores[q_acc]) >= 20:
                    # current alignment is better than at least one of previous scores, remove the worst one so far
                    if paf_similarity_score > highest_paf_scores[q_acc][0][0]: 
                        paf_score, t_acc_out = heapq.heappushpop(highest_paf_scores[q_acc], (paf_similarity_score, t_acc) )
                        # print("popping:", paf_score, "pushing:", paf_similarity_score, t_acc_out )
                else:
                    heapq.heappush(highest_paf_scores[q_acc], (paf_similarity_score, t_acc))
                    # print("pushing", paf_similarity_score)

                # if len(highest_paf_scores[t_acc]) >= 5:
                #     # current alignment is better than at least one of previous scores, remove the worst one so far
                #     if paf_similarity_score > highest_paf_scores[t_acc][0][0]: 
                #         paf_score, q_acc_out = heapq.heappushpop(highest_paf_scores[t_acc], (paf_similarity_score, q_acc) )
                # else:
                #     heapq.heappush(highest_paf_scores[t_acc], (paf_similarity_score, q_acc))

    for acc1 in highest_paf_scores:
        # print(acc_to_strings)
        s1 = acc_to_strings[str(acc1)]
        matches[s1] = [] 

        matches_list = highest_paf_scores[acc1]
        for score, acc2 in matches_list:
            s2 = acc_to_strings[str(acc2)]
            matches[s1].append(s2)

    return matches


def map_with_minimap(targets, queries):
    print('Aligning with minimap.')
    sys.stdout.flush()
    # work_dir = "/tmp/" #tempfile.mkdtemp() 
    # print(work_dir)
    minimap_output = targets + ".paf"
    stderr_file = open(targets + ".minimap.stderr", 'w')
    # print('Output path: ', minimap_output)
    # print('Stderr file: ', stderr_file)
    sys.stdout.flush()
    # print(type(targets))
    processes = mp.cpu_count()
    print("minimap with {0} processes.".format(processes))
    print("minimap -f 0.00000001 -Sw5 -L40 -m0 -t {0} {1} {2}".format(processes, targets, queries))
    # print("minimap2 -x ava-pb -t {0} {1} {2}".format(processes, targets, queries))

    with open(minimap_output, "w") as minimap_file:
        sys.stdout.flush()
        subprocess.check_call([ "minimap", "-f", "0.00000001", "-Sw5", "-L40", "-m0", 
                                "-t", str(processes),
                               targets, queries ],
                                stdout=minimap_file,
                                stderr=stderr_file)
        # subprocess.check_call([ "minimap2", "-x", "ava-pb", 
        #                         "-t", str(processes),
        #                        targets, queries ],
        #                         stdout=minimap_file,
        #                         stderr=stderr_file)
        sys.stdout.flush()

    return minimap_output

def minimap_partition(targets, queries, params):
    # partition unique strings (speed optimization for minimap)
    # this works for transcript version of 3CO

    target_strings = list(targets)
    target_strings.sort(key=lambda x: len(x))

    #assign unique accessions to each string
    # queries \subset targets so we only need to have unique inexes for targets
    query_acc_to_seq = {}
    target_acc_to_seq = {}
    query_seq_to_acc = {}
    target_seq_to_acc = {}

    query_strings = []
    # query_strings = list(queries)
    # query_strings.sort(key=lambda x: len(x))
    # print([len(s) for s in query_strings])

    for i, seq in enumerate(target_strings):
        target_acc_to_seq[str(i)] = seq
        target_seq_to_acc[seq] = str(i)
        if seq in queries:
            # print(i, len(seq))
            query_acc_to_seq[str(i)] = seq
            query_seq_to_acc[seq] = str(i)
            query_strings.append(seq)



    target_bins = []
    query_bins = []
    query_bin_size =  min(len(query_strings), params.minimap_bin_size )
    # print("LEN QUERY STRINGS:", len(query_strings))
    # print("LEN TARGET STRINGS:", len(target_strings))
    # print([query_seq_to_acc[b] for b in query_strings])

    for i in range(0, len(query_strings), query_bin_size):
        # create a query bin
        bin = query_strings[i  : i+query_bin_size ]
        # print("q lengths:", len(bin[0]), len(bin[-1]))
        # print( [ len(b) for b in bin])
        # print([query_seq_to_acc[b] for b in bin])
        labeled_query_strings = {}
        for q_seq in bin:
            labeled_query_strings[ query_seq_to_acc[q_seq] ] = q_seq
        query_bins.append(labeled_query_strings)

        min_query_acc, max_query_acc = query_seq_to_acc[bin[0]], query_seq_to_acc[bin[-1]]
        print("q lengths:", len(bin[0]), len(bin[-1]), "t lengths:", len(target_acc_to_seq[ str(max(int(min_query_acc) - 50, 0)) ]) , len(target_acc_to_seq[ str(min(int(max_query_acc) + 50 +1, len(target_strings) - 1)) ]) )
        # print("Target bin size:", len(range(max(int(min_query_acc) - 50, 0) , min(int(max_query_acc) + 50 +1, len(target_strings)))), int(max_query_acc) - int(min_query_acc) )
        # print("lol", len(range(int(min_query_acc), int(max_query_acc)+1)))

        # create the respective target bin
        labeled_target_strings = {}
        for target_acc in range(max(int(min_query_acc) - 50, 0) , min(int(max_query_acc) + 50 +1, len(target_strings)) ):
            labeled_target_strings[str(target_acc)] = target_acc_to_seq[str(target_acc)]

        target_bins.append(labeled_target_strings)
        # labeled_strings_bin = [(i+j, s) for j,s in enumerate(bin)]




    work_dir = params.tempfolder 
    fasta_query_files = []
    fasta_target_files = []
    # acc_to_strings = {}
    for i in range(len(query_bins)):
        fasta_query_name = os.path.join(work_dir, "q" + str(i) + ".fa")
        fasta_query_file = open(fasta_query_name, "w")
        for acc, seq in query_bins[i].items():
            fasta_query_file.write(">{0}\n{1}\n".format(acc, seq))

        fasta_query_file.close()
        fasta_query_files.append(fasta_query_name)

        fasta_target_name = os.path.join(work_dir, "t" + str(i)+".fa")
        fasta_target_file = open(fasta_target_name, "w")
        for acc, seq in target_bins[i].items():
            fasta_target_file.write(">{0}\n{1}\n".format(acc, seq))

        fasta_target_file.close()
        fasta_target_files.append(fasta_target_name)

    # TODO: call minimap, eventually parallelize
    paf_file_names = []
    for fasta_target_file, fasta_query_file in zip(fasta_target_files, fasta_query_files):
        paf_file_name = map_with_minimap(fasta_target_file, fasta_query_file)
        paf_file_names.append(paf_file_name)

    return paf_file_names, target_acc_to_seq

