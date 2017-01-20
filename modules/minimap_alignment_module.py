
import os
import sys
# import unittest
import tempfile
import subprocess
import gzip
import heapq
from collections import defaultdict


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
            q_acc = row_info[0].decode('ascii')
            q_len = int(row_info[1])
            t_acc = row_info[5].decode('ascii')
            t_len = int(row_info[6])

            if q_acc == t_acc:
                # print("SELF MAPPING DETECTED")
                # print(q_acc, q_len, row_info)
                # print(t_acc, t_len, row_info)
                continue

            n_min = int(row_info[12].split(":")[-1])
            #   \frac{\min \{ |q|, |t| \} } {\max \{ |q|,|t| \} }  n^{\text{minimizers}}
            paf_similarity_score = n_min * min(q_len, t_len)/float(max(q_len, t_len))

            if len(highest_paf_scores[q_acc]) >= 5:
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
        try:
            file_object = gzip.open(paf_file_path)
            file_object.readline()
            file_object.seek(0)
        except IOError:
            file_object = open(paf_file_path)


        with file_object as paf_file:
            for line in paf_file:
                row_info = line.strip().split()
                q_acc = row_info[0].decode('ascii')
                q_len = int(row_info[1])
                t_acc = row_info[5].decode('ascii')
                t_len = int(row_info[6])

                if q_acc == t_acc:
                    # print("SELF MAPPING DETECTED")
                    # print(q_acc, q_len, row_info)
                    # print(t_acc, t_len, row_info)
                    continue

                n_min = int(row_info[12].split(":")[-1])
                #   \frac{\min \{ |q|, |t| \} } {\max \{ |q|,|t| \} }  n^{\text{minimizers}}
                paf_similarity_score = n_min * min(q_len, t_len)/float(max(q_len, t_len))

                if len(highest_paf_scores[q_acc]) >= 5:
                    # current alignment is better than at least one of previous scores, remove the worst one so far
                    if paf_similarity_score > highest_paf_scores[q_acc][0][0]: 
                        paf_score, t_acc_out = heapq.heappushpop(highest_paf_scores[q_acc], (paf_similarity_score, t_acc) )
                else:
                    heapq.heappush(highest_paf_scores[q_acc], (paf_similarity_score, t_acc))

                if len(highest_paf_scores[t_acc]) >= 5:
                    # current alignment is better than at least one of previous scores, remove the worst one so far
                    if paf_similarity_score > highest_paf_scores[t_acc][0][0]: 
                        paf_score, q_acc_out = heapq.heappushpop(highest_paf_scores[t_acc], (paf_similarity_score, q_acc) )
                else:
                    heapq.heappush(highest_paf_scores[t_acc], (paf_similarity_score, q_acc))

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
    with open(minimap_output, "w") as minimap_file:
        sys.stdout.flush()
        subprocess.check_call([ "minimap", "-f", "0.0000000001", "-Sw5", "-L100", "-m0",
                               targets, queries ],
                                stdout=minimap_file,
                                stderr=stderr_file)
        sys.stdout.flush()

    return minimap_output

def minimap_partition(targets, queries):
    # partition unique strings (speed optimization for minimap)
    # this works for transcript version of 3CO

    #assign unique accessions to each string
    # queries \subset targets so we only need to have unique inexes for targets
    query_accessions = {}
    target_accessions = {}

    for seq in targets:


    unique_strings = list(targets)
    unique_strings.sort(key=lambda x: len(x))
    # labeled_unique_strings.sort(key=lambda x: len(x))
    bins = []
    bin_size =  min(len(unique_strings), 200 )
    for i in range(0, len(unique_strings), bin_size):
        # create overlapping bins
        bin = unique_strings[i: i+bin_size + 50]
        labeled_strings_bin = [(i+j, s) for j,s in enumerate(bin)]
        bins.append(labeled_strings_bin)


    work_dir = "/tmp/" 
    fasta_files = []
    acc_to_strings = {}
    for i, labeled_strings_bin in enumerate(bins):
        fasta_file_name = os.path.join(work_dir,str(i)+".fa")
        fasta_file = open(fasta_file_name, "w")
        for acc, seq in labeled_strings_bin:
            fasta_file.write(">{0}\n{1}\n".format(acc, seq))
            acc_to_strings[str(acc)] = seq

        fasta_file.close()
        fasta_files.append(fasta_file_name)

    # TODO: call minimap, eventually parallelize
    paf_file_names = []
    for fa_file_name in fasta_files:
        paf_file_name = map_with_minimap(fa_file_name, fa_file_name)
        paf_file_names.append(paf_file_name)

    return paf_file_names, acc_to_strings

