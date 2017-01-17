"""
    minimizer_graph(S): creates the minimizer graph defined in..
    partition_graph(S): creates the partition of a graph defined in..

"""
import os
import sys
import unittest
import tempfile
import subprocess

import alignment_module

def paf_to_best_matches(paf_files, acc_to_strings):
    """
        input: PAF file
        output: for each sequence found in PAF file, store X best matches based on paf score
    """

    matches = {}

    for paf_file_path in paf_files:
        try:
            file_object = gzip.open(paf_file_path)
            file_object.readline()
            file_object.seek(0)
        except IOError:
            file_object = open(paf_file_path)

        highest_paf_scores = defaultdict(list)

        with file_object as paf_file:
            for line in paf_file:
                row_info = line.strip().split()
                q_acc = row_info[0].decode('ascii')
                q_len = int(row_info[1])
                t_acc = row_info[5].decode('ascii')
                t_len = int(row_info[6])

                if q_acc == t_acc:
                    print("SELF MAPPING DETECTED")
                    print(q_acc, q_len, row_info)
                    print(t_acc, t_len, row_info)
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
        s1 = acc_to_strings[acc1]
        matches[s1] = [] 

        matches_list = highest_paf_scores[acc1]
        for score, acc2 in matches_list:
            s2 = acc_to_strings[acc2]
            matches[s1].append(s2)

    return matches


def map_with_minimap(sequences_file_name):
    print('Aligning with minimap.')
    sys.stdout.flush()
    work_dir = "/tmp/" #tempfile.mkdtemp() 
    # print(work_dir)
    minimap_output = os.path.join(work_dir, "minimap.paf")
    stderr_file = open(os.path.join(work_dir,"minimap.stderr"), 'w')
    # print('Output path: ', minimap_output)
    # print('Stderr file: ', stderr_file)
    sys.stdout.flush()
    # print(type(sequences_file_name))
    with open(minimap_output, "w") as minimap_file:
        sys.stdout.flush()
        subprocess.check_call([ "minimap", "-f", "0.0000000001", "-Sw5", "-L100", "-m0",
                               sequences_file_name, sequences_file_name ],
                                stdout=minimap_file,
                                stderr=stderr_file)
        sys.stdout.flush()

    return minimap_output

def minimap_partition(unique_strings):
    # partition unique strings (speed optimization for minimap)
    # this works for transcript version of 3CO
    unique_strings.sort(key=lambda x: len(x))
    # labeled_unique_strings.sort(key=lambda x: len(x))
    bins = []
    bin_size =  min(len(unique_strings), 200 )
    for i in range(0, len(unique_strings), bin_size):
        # create overlapping bins
        bin = unique_strings[i: i+bin_size + 20]
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
            acc_to_strings[acc] = seq

        fasta_file.close()
        fasta_files.append(fasta_file_name)

    # call minimap, eventually parallelize
    paf_file_names = []
    for fa_file_name in fasta_files:
        paf_file_name = map_with_minimap(fa_file_name)
        paf_file_names.append(paf_file_name)

    return paf_file_names, acc_to_strings

def find_best_matches(approximate_matches):
    """
        input: dictionary with a string as key and a list of strings as value
        output: dictionary with a string as key and a (reduced) list of strings as value.
                Each string in the reduced list has the same (lowest) edit distance to the key 
    """
    exact_matches = alignment_module.sw_align_sequences(approximate_matches)

    # process the exact matches here

    
    # for s1 in matches:
    #     for s2 in matches[s1]:


def construct_minimizer_graph(S):

    """
        input: a dict of strings, not necesarily unique
        output: a directed graph implemented as a dict of dicts. Each edge has a weight assosiated to them.
                self edges has a weight > 1 (identical sequences) and all other edges has weight 1.
    """
    G_star = {}
    alignment_graph = {}
    # adding self edges to strings that has converged
    for acc, s in S.items():
        if s not in G_star:
            G_star[s] = {}
            alignment_graph[s] = {}
        else:
            if s in G_star[s]:
                G_star[s][s] += 1  
            else:
                G_star[s][s] = 1
                alignment_graph[s] = s


    unique_strings = set(S.values())
    paf_files, acc_to_strings = minimap_partition(unique_strings)

    approximate_matches = paf_to_best_matches(paf_files, acc_to_strings)
    find_best_matches(approximate_matches)

    # take care of isolated nodes in minimizer graph here

    #construct G_star and alignment_graph
    return G_star, alignment_graph


def construct_2set_minimizer_bipartite_graph(S, T):
    return

class TestFunctions(unittest.TestCase):

    def test_map_with_minimap(self):
        self.maxDiff = None
        temp_file_name = "/tmp/test.fa"
        temp_file = open(temp_file_name, "w")
        temp_file.write(">consensus_1139_from_2_reads_pval_1.0\nGGTCGTTTTTAAACTATTCGACACTAATTGATGGCCATCCGAATTCTTTTGGTCGCTGTCTGGCTGTCAGTAAGTATGCTAGAGTTCCGTTTCCGTTTCATTACCAACACCACGTCTCCTTGCCCAATTAGCACATTAGCCTTCTCTCCTTTCGCAAGGTTGCTCAGTTCATTTATGCTTAATGCTGGTCCATATCTCCTGTCTTCTTTGCCCAGAATGAGGAATCCTCTCAGAACTGCGGACTCAACTCCAGCTGTGCCTTCATCTGGGTCTTCAGTTAAAGGGCCAGCATCCTTTCCGAGAACTGTGAGTCTTTTAGTGGTCTTGTTGTAGTTGAATACTGGAGAATTGCCCCTTACAAGTATTCTCATTCCTGATCCCCTCACATTTATAGTCAATGAGGAGAACTGCGTTCTACTTTGCTTTGGTGGAGCGGCTGCGAAGGGAAGAAGTTTTATTATCTGAGCGGTATCAAATGTCCCAAGCACATCCCTCATTTGTTGGAACAGAGTTCTCACAAACCCCACTGTATTGGCCTCTAACGGCCTTTGGAACTAAAGACTGAAATGGCTCAAATTCCATTTTATTGTACAGCATTGTAGGATTCTGGGACCACTGAATTTTAACAGTTTCCCAGTTTCTGATGATCCACTGATAGGTATTGACCAACACTGATTCAGGACCATTAATCTCCCACATCATTGACGATGAGTAAGTTATTGTCAGTTTCTCTGTTCCCTGTGTTTCACTGATCTCCTCGGGAGACAGTAGTACATTCCCACGTTGGTCCCTAACTCTCAAAAAACGGTCAATGCTCACCACTATCTTCTCCGCGCTGGAATACTCATCTACCCCCATTTTGCTGATTCTCACTCCTCTCATTGACATCTCGGTGCTTGGAGTCATGTCGGGCAATATCCCGATCATTCCCATCACATTGTCGATGGATTCAATTCCCCAATTTTGAAAGAGCACCTTTGCATCCTTCTGAAAATGTCTCAAAAGTTGGTGCATGGGATTCAATCGCTGATTCGCCCTATTGACGAAATTCAGGTCACCTCTAACTGCTTTTATCATACAATCCTCTTGTGAAAATACCATGGCCACAATTATTGCTTCGGCAATCGACTGTTCGTCTCTCCCACTCACTATCAGCTGAATCAATCTCCTGGTTGCTTTTCTGAGTATAGCTGTTGCTCTTCTCCCAACCATTGTGAACTCTTCATATCCCTCATGCACTCTTATCTTCAATGTCTGAAGATTGCCCGTAAGCACCTCTTCCTCTCTCTTGACTGATGATCCGCTTGTTCTCTTAAATGTGAATCCACCAAAACTGAAGGATGAGCTAATTCTCAGTCCCATTGCAGCCTTGCAAATATCCACGGCTTGCTCTTCTGTTGGGTTCTGCCTAAGGATGTTTACCATCCTTATTCCACCAATCTGCGTGCTGTGGCACATCTCCAATAAAGATGCTAGTGGATCTGCTGATACTGTGGCTCTTCTTACTATGTTTCTAGCAGCAATAATTAAGCTTTGATCAACATCATCATTCCTCGCCTCCCCTCCTGGAGTGTACATCTGTTCCCAGCATGTTCCTTGGGTCAAATGCAACACTTCAATGTACACACTGCTTGTTCCACCAGCCACTGGGAGGAATCTCGTTTTGCGGACCAGTTCTCTCTCCGACATGTATGCCACCATCAGAGGAGAAATTTTGCAACCCTGGAGTTCTTCTTTCTTCTCTTTGGTTGTCGTTAGTTGCGATTCCGATGTTAGTATCCTGGCTCCCACTTCGTTAGGGAAAACAACTTCCATGATTACATCCTGTGCCTCTTTGGCACTGAGATCTGCATGACCAGGATTTATGTCAACTCTTCGACGTATTTTGACTTGGTTTCTAAAATGGACAGGGCCAAAGGTTCCATGTTTTAACCTTTCGACTTTTTCAAAATAAGTTTTGTAGATTTTTGGATAATGAACTGTACTTGTCACTGGTCCATTCCTATTCCACCATGTCACAGCCAGAGGTGATACCATCACTCGGTCTGATCCGGCGTCATTCATTTTACTCCATAAAGTTTGTCCCTGCTCATTTCTCTCAGGAATCATTTCCGTTATCCTCTTGTCTGCTGTAATTGGATATTTCATTGCCATCATCCATTTCATCCTAAGTGCTGGGTTCTTCTCCTGTCTTCCTGATGTGTACTTCTTGATTATGGCCATATGGTCCACGGTGGTTTTTGTGAGTATCTCGCGAGTGCGAGACTGCGACATTAGATTCCTTAGTTCTTTTATTCTTTCC\n>consensus_940_from_4_reads_pval_1.0\nGGTCGTTTTTAAACTATTCGACACTAATTGATGGCCATCCGAATTCTTTTGGTCGCTGTCTGGCTGTCAGTAAGTATGCTAGAGTTCCGTTTCCGTTTCATTACCAACACCACGTCTCCTTGCCCAATTAGCACATTAGCCTTCTCTCCTTTCGCAAGGTTGCTCAGTTCATTTATGCTTAATGCTGGTCCATATCTCCTGTCTTCTTTGCCCAGAATGAGGAATCCTCTCAGAACTGCGGGACTCAACTCCAGCTGTGCCTTCATCTGGGTCTTCAGTTAAAGGGCCAGCATCCTTTCCGAGAACTGTGAGTCTTTTAGTGGTCTTGTTGTAGTTGAATACTGGAGAATTGCCCCTTACAAGTATTCTCATTCCTGATCCCCTCACATTTATAGTCAATGAGGAGAACTGCGTTCTACTTTGCTTTGGTGGAGCGGCTGCGAAGGGAAGAAGTTTTATTATCTGAGCGGTATCAAATGTCCCAAGCACATCCCTCATTTGTTGGAACAGAGTTCTCACAAACCCACTGTATTGGCCTCTAACGGCCTTTGGAACTAAAGACTGAAATGGCTCAAATTCCATTTTATTGTACAGCATTGTAGGATTCTGGGACCACTGAATTTTAACAGTTTCCCAGTTTCTGATGATCCACTGATAGGTATTGACCAACACTGATTCAGGACCATTAATCTCCCACATCATTGACGATGAGTAAGTTATTGTCAGTTTCTCTGTTCCCTGTGTTTCACTGATCTCCTCGGGAGACAGTAGTACATTCCCACGTTGGTCCCTAACTCTCAAAAAACGGTCAATGCTCACCACTATCTTCTCCGCGCTGGAATACTCATCTACCCCCATTTTGCTGATTCTCACTCCTCTCATTGACATCTCGGTGCTTGGAGTCATGTCGGGCAATATCCCGATCATTCCCATCACATTGTCGATGGATTCAATTCCCCAATTTTGAAAGAGCACCTTTGCATCCTTCTGAAAATGTCTCAAAAGTTGGTGCATGGGATTCAATCGCTGATTCGCCCTATTGACGAAATTCAGGTCACCTCTAACTGCTTTTATCATACAATCCTCTTGTGAAAATACCATGGCCACAATTATTGCTTCGGCAATCGACTGTTCGTCTCTCCCACTCACTATCAGCTGAATCAATCTCCTGGTTGCTTTTCTGAGTATAGCTGTTGCTCTTCTCCCAACCATTGTGAACTCTTCATATCCCTCATGCACTCTTATCTTCAATGTCTGAAGATTGCCCGTAAGCACCTCTTCCTCTCTCTTGACTGATGATCCGCTTGTTCTCTTAAATGTGAATCCACCAAAACTGAAGGATGAGCTAATTCTCAGTCCCATTGCAGCCTTGCAAATATCCACGGCTTGCTCTTCTGTTGGGTTCTGCCTAAGGATGTTTACCATCCTTATTCCACCAATCTGCGTGCTGTGGCACATCTCCAATAAAGATGCTAGTGGATCTGCTGATACTGTGGCTCTTCTTACTATGTTTCTAGCAGCAATAATTAAGCTTTGATCAACATCATCATTCCTCGCCTCCCCTCCTGGAGTGTACATCTGTTCCCAGCATGTTCCTTGGGTCAAATGCAACACTTCAATGTACACACTGCTTGTTCCACCAGCCACTGGGAGGAATCTCGTTTTGCGGACCAGTTCTCTCTCCGACATGTATGCCACCATCAGAGGAGAAATTTTGCAACCCTGGAGTTCTTCTTTCTTCTCTTTGGTTGTCGTTAGTTGCGATTCCGATGTTAGTATCCTGGCTCCCACTTCGTTAGGGAAAACAACTTCCATGATTACATCCTGTGCCTCTTTGGCACTGAGATCTGCATGACCAGGATTTATGTCAACTCTTCGACGTATTTTGACTTGGTTTCTAAAATGGACAGGGCCAAAGGTTCCATGTTTTAACCTTTCGACTTTTTCAAAATAAGTTTTGTAGATTTTTGGATAATGAACTGTACTTGTCACTGGTCCATTCCTATTCCACCATGTCACAGCCAGAGGTGATACCATCACTCGGTCTGATCCGGCGTCATTCATTTTACTCCATAAAGTTTGTCCCTGCTCATTTCTCTCAGGAATCATTTCCGTTATCCTCTTGTCTGCTGTAATTGGATATTTCATTGCCATCATCCATTTCATCCTAAGTGCTGGGTTCTTCTCCTGTCTTCCTGATGTGTACTTCTTGATTATGGCCATATGGTCCACGGTGGTTTTTGTGAGTATCTCGCGAGTGCGAGACTGCGACATTAGATTCCTTAGTTCTTTTATTCTTTCC\n>consensus_1222_from_12_reads_pval_1.0\nGGTCGTTTTTAAACTATTCGACACTAATTGATGGCCATCCGAATTCTTTTGGTCGCTGTCTGGCTGTCAGTAAGTATGCTAGAGTTCCGTTTCCGTTTCATTACCAACACCACGTCTCCTTGCCCAATTAGCACATTAGCCTTCTCTCCTTTCGCAAGGTTGCTCAGTTCATTTATGCTTAATGCTGGTCCATATCTCCTGTCTTCTTTGCCCAGAATGAGGAATCCTCTCAGAACTGCGGACTCAACTCCAGCTGTGCCTTCATCTGGGTCTTCAGTTAAAGGGCCGGCATCCTTTCCGAGAACTGTGAGTCTTTTAGTGGTCTTGTTGTAGTTGAATACTGGAGAATTGCCCCTTACAAGTATTCTCATTCCTGATCCCCTCACATTTATAGTCAATGAGGAGAACTGCGTTCTACTTTGCTTTGGTGGAGCGGCTGCGAAGGGAAGAAGTTTTATTATCTGAGCGGTATCAAATGTCCCAAGCACATCCCTCATTTGTTGGAACAGAGTTCTCACAAACCCACTGTATTGGCCTCTAACGGCCTTTGGAACTAAAGACTGAAATGGCTCAAATTCCATTTTATTGTACAGCATTGTAGGATTCTGGGACCACTGAATTTTAACAGTTTCCCAGTTTCTGATGATCCACTGATAGGTATTGACCAACACTGATTCAGGACCATTAATCTCCCACATCATTGACGATGAGTAAGTTATTGTCAGTTTCTCTGTTCCCTGTGTTTCACTGATCTCCTCGGGAGACAGTAGTACATTCCCACGTTGGTCCCTAACTCTCAAAAAACGGTCAATGCTCACCACTATCTTCTCCGCGCTGGAATACTCATCTACCCCCATTTTGCTGATTCTCACTCCTCTCATTGACATCTCGGTGCTTGGAGTCATGTCGGGCAATATCCCGATCATTCCCATCACATTGTCGATGGATTCAATTCCCCAATTTTGAAAGAGCACCTTTGCATCCTTCTGAAAATGTCTCAAAAGTTGGTGCATGGGATTCAATCGCTGATTCGCCCTATTGACGAAATTCAGGTCACCTCTAACTGCTTTTATCATACAATCCTCTTGTGAAAATACCATGGCCACAATTATTGCTTCGGCAATCGACTGTTCGTCTCTCCCACTCACTATCATCTGAATCAATCTCCTGGTTGCTTTTCAGAGTATAGCTGTTGCTCTTCTCCCAACCATTGTGAACTCTTCATATCCCTCATGCACTCTTATCTTCAATGTCTGAAGATTGCCCGTAAGCACCTCTTCCTCTCTCTTGACTGATGATCCGCATGTTCTCTTAAATGTGAATCCACCAAAACTGAAGGATGAGCTAATTCTCAGTCCCATTGCAGCCTTGCAAATATCCACGGCTTGCTCTTCTGTTGGGTTCTGCCTAAGGATGTTTACCATCCTTATTCCACCAATCTGCGTGCTGTGGCACATCTCCAATAAAGATGCTAGTGGATCTGCTGATACTGTGGCTCTTCTTACTATGTTTCTAGCAGCAATAATTAAGCTTTGATCAACATCATCATTCCTTGCCTCCCCTCCTGGAGTGTACATCTGTTCCCAGCATGTTCCTTGGGTCAAATGCAACACTTCAATGTACACACTGCTTGTTCCACCAGCCACTGGGAGGAATCTCGTTTTGCGGACCAGTTCTCTCTCCAACATGTATGCCACCATCAGAGGAGAAATTTTGCAACCCTGGAGTTCTTCTTTCTTCTCTTTGGTTGTCGTTAGTTGCGATTCCGATGTTAGTATCCTGGCTCCCACTTCGTTAGGGAAAACAACTTCCATGATTACATCCTGTGCCTCTTTGGCACTGAGATCTGCATGACCAGGATTTATGTCAACTCTTCGACGTATTTTGACTTGGTTTCTAAAATGAACAGGGCCAAAGGTTCCATGTTTTAACCTTTCGACTTTTCAAAATAAGTTTTGTAGATTTTTGGATAATGAACTGTACTTGTCACTGGTCCATTCCTATTCCACCAGTCACAGCCAGAGGTGATACCAACACTCGGTCTGATCCGGCGTCATTCATTTTACTCCATAAAGTTTGTCCCTGCTCATTTCTCTCAGGAATCATTTCCGTTATCCTCTTGTCTGCTGTAATTGGATATTTCATTGCCATCATCCATTTCATCCTAAGTGCTGGGTTCTTCTCCTGTCTTCCTGATGTGTACTTCTTGATTATGGCCATATGGTCCACGGTGGTTTTTGTGAGTATCTCGCGAGTGCGAGACTGCGCCATTAGATTCCTTAGTTCTTTTATTCTTTCC")
        temp_file.close()
        s1,s2,s3 = "\t".join("consensus_1139_from_2_reads_pval_1.0    2302    0       2302    +       consensus_940_from_4_reads_pval_1.0     2302    0       2302    2297    2302    255     cm:i:786".split()), "\t".join("consensus_1139_from_2_reads_pval_1.0    2302    0       2302    +       consensus_1222_from_12_reads_pval_1.0   2299    0       2299    2260    2302    255     cm:i:728".split()), "\t".join("consensus_1222_from_12_reads_pval_1.0   2299    0       2299    +       consensus_940_from_4_reads_pval_1.0     2302    0       2302    2256    2302    255     cm:i:726".split())
        s1 = s1 + "\n"
        s2 = s2 + "\n"
        s3 = s3 + "\n"
        expected_result = "".join([s1,s2,s3])
        minimap_paf = map_with_minimap(temp_file_name)
        minimap_result = ""
        for line in open(minimap_paf, "r"):
            minimap_result += line
        self.assertEqual(minimap_result, expected_result)


if __name__ == '__main__':
    unittest.main()