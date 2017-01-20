"""
    minimizer_graph(S): creates the minimizer graph defined in..
    partition_graph(S): creates the partition of a graph defined in..

"""
import os
import sys
import unittest
import tempfile
import subprocess
import gzip

from collections import defaultdict
import heapq

from modules import alignment_module
from modules import functions


def paf_to_best_matches(paf_files, acc_to_strings):
    """
        input: PAF file
        output: for each sequence found in PAF file, store X best matches based on paf score
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


def map_with_minimap(sequences_file_name):
    print('Aligning with minimap.')
    sys.stdout.flush()
    # work_dir = "/tmp/" #tempfile.mkdtemp() 
    # print(work_dir)
    minimap_output = sequences_file_name + ".paf"
    stderr_file = open(sequences_file_name + ".minimap.stderr", 'w')
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

def minimap_partition(unique_strings_set):
    # partition unique strings (speed optimization for minimap)
    # this works for transcript version of 3CO
    unique_strings = list(unique_strings_set)
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
        paf_file_name = map_with_minimap(fa_file_name)
        paf_file_names.append(paf_file_name)

    return paf_file_names, acc_to_strings

def find_best_matches(approximate_matches):
    """
        input: dictionary with a string as key and a list of strings as value
        output: dictionary with a string as key and a dictionary as value. 
                The inner dict contains a tuple (edit_distance, s1_alignment, s2_alignment) as value.
                Each string in the inner dict has the same (lowest) edit distance to the key 
    """
    exact_matches = alignment_module.sw_align_sequences(approximate_matches, single_core = False )

    # process the exact matches here
    best_exact_matches = {}
    for s1 in exact_matches:
        for s2 in exact_matches[s1]:
            s1_alignment, s2_alignment, (matches, mismatches, indels) = exact_matches[s1][s2]
            edit_distance = mismatches + indels

            if s1 in best_exact_matches:
                s1_minimizer = best_exact_matches[s1].keys()[0]
                if edit_distance < best_exact_matches[s1][s1_minimizer][0]:
                    best_exact_matches[s1] = {}
                    best_exact_matches[s1][s2] = (edit_distance, s1_alignment, s2_alignment)
                elif edit_distance == best_exact_matches[s1][s1_minimizer][0]:
                    best_exact_matches[s1][s2] = (edit_distance, s1_alignment, s2_alignment)
            else:
                best_exact_matches[s1] = {}
                best_exact_matches[s1][s2] = (edit_distance, s1_alignment, s2_alignment)

            if s2 in best_exact_matches:
                s2_minimizer = best_exact_matches[s2].keys()[0]
                if edit_distance < best_exact_matches[s2][s2_minimizer][0]:
                    best_exact_matches[s2] = {}
                    best_exact_matches[s2][s1] = (edit_distance, s2_alignment, s1_alignment)
                elif edit_distance == best_exact_matches[s2][s2_minimizer][0]:
                    best_exact_matches[s2][s1] = (edit_distance, s2_alignment, s1_alignment)
            else:
                best_exact_matches[s2] = {}
                best_exact_matches[s2][s1] = (edit_distance, s2_alignment, s1_alignment)

    return best_exact_matches


def construct_minimizer_graph(S):

    """
        input: a dict of strings, not necesarily unique
        output: a directed graph implemented as a dict of dicts. Each edge has a weight assosiated to them.
                self edges has a weight > 1 (identical sequences) and all other edges has weight 1.
                Note, a node can be isolated!
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
                G_star[s][s] = 2
                alignment_graph[s][s] = (0, s, s)

    # check if converged, that is, if all nodes has self edges here, there will be no other edges added.
    converged = True
    for nbr in G_star.values():
        if len(nbr) == 0:
            converged = False

    if converged:
        return G_star, alignment_graph, converged

    unique_strings = set(S.values())
    paf_files, acc_to_strings = minimap_partition(unique_strings)

    approximate_matches = paf_to_best_matches(paf_files, acc_to_strings)
    best_exact_matches = find_best_matches(approximate_matches)

    #add remaining edges to  G_star and alignment_graph
    for s1 in best_exact_matches:
        # already have weighted self edge, i.e., identical sequence
        if s1 in G_star[s1]:
            # print("here!", G_star[s1][s1])
            continue
        for s2 in best_exact_matches[s1]:
            assert s2 not in G_star[s1]
            G_star[s1][s2] = 1
            (edit_distance, s1_alignment, s2_alignment) = best_exact_matches[s1][s2]
            alignment_graph[s1][s2] = (edit_distance, s1_alignment, s2_alignment)

    # finally, nodes here are the one where we didn't find any alignments to another sequence, point these isolated nodes to themself
    # with indegree 1
    # we also check that theu are not "leaves" in G^* 
    # that is, a sequence that has no minimizer but is a minimizer to some other sequence, this should not happen in G^*
    G_star_transposed = functions.transpose(G_star)

    for s in G_star:
        if len(G_star[s]) == 0:
            assert s not in G_star_transposed
            G_star[s][s] = 1
            alignment_graph[s][s] = (0, s, s)
            print("ISOLATED")

    return G_star, alignment_graph, converged


def partition_strings(S):
    G_star, alignment_graph, converged = construct_minimizer_graph(S)
    partition_alignments = {}
    unique_start_strings = set(G_star.keys())
    partition = {} # dict with a center as key and a set containing all sequences chosen to this partition

    if converged:
        M = {key : 1 for key in G_star.keys()}
        for m in G_star:
            partition[m] = set()
            partition_alignments[m] = {}
            indegree = G_star[m][m]
            partition_alignments[m][m] = (0, m, m , indegree)
        return partition_alignments, partition, M, converged

    marked = set()
    M = {}
    V_not_in_M = set(G_star.keys())
    partition_counter = 1

    for s in G_star:
        if s in G_star[s]:
            if  G_star[s][s] == 1: # isolate
                # isolated += 1
                # partition[n] = set()
                M[s] = partition_counter
                marked.add(s)
                partition[s] = set()
                partition_counter += 1
                partition_alignments[s] = {}
                edit_dist, a1, a_min =  alignment_graph[s][s]
                indegree = G_star[s][s]
                partition_alignments[s][s] = (edit_dist, a_min, a1 , indegree)
                print("DETECTED!!")
    isolated = set(partition.keys())

    G_star_transposed = functions.transpose(G_star)
    # check so that there are no "leaves" in G^* (a sequence that has no minimizer but is a minimizer to some other sequence, this should not happen)
    # for n in G_star_transposed:
    #     assert n not in isolated

    # do_while as long as there is a node with positive indegree of unmarked nbrs
    while True:
        # find node with biggest indegree
        max_indegree = -1
        for s in V_not_in_M:
            indegree = sum([indegree for in_nbr, indegree in G_star_transposed[s].items() if in_nbr not in marked ])
            if indegree > max_indegree:
                m, max_indegree = s, indegree
        # print(max_indegree, len(V_not_in_M), len(marked))
        if max_indegree < 1:
            break
        M[m] = partition_counter
        partition[m] = set()
        V_not_in_M.remove(m)
        partition_alignments[m] = {}
        partition_counter += 1
        marked.add(m)
        for in_nbr_to_m in G_star_transposed[m]:
            marked.add(in_nbr_to_m)

    print("Chosen minimizers:", len(M))
    for s in G_star:
        if s not in M:
            # since M covers G_star n has to have a at least one minimizer in M, choose the biggest one, i.e., lowest index
            lowest_index = len(M) + 1
            for nbr in G_star[s]:
                if nbr in M:
                    if M[nbr] < lowest_index:
                        lowest_index = M[nbr]
                        lowest_index_minimizer = nbr

            partition[lowest_index_minimizer].add(s)
            # print("1")
            edit_dist, a1, a_min =  alignment_graph[s][lowest_index_minimizer]
            indegree = G_star[s][lowest_index_minimizer]


        else:
            lowest_index_minimizer = s
            edit_dist, a1, a_min =  0, s, s
            if lowest_index_minimizer in  G_star[s]:
                indegree = G_star[s][lowest_index_minimizer]
            else:
                indegree = 1

            # print("2")

        # print(G_star[s])
        # print(lowest_index_minimizer == s)
        # print(alignment_graph[s])
        partition_alignments[lowest_index_minimizer][s] = (edit_dist, a_min, a1, indegree)


    total_strings_in_partition = sum([ len(partition[p]) +1 for p in  partition])
    partition_sequences = set()
    for m in partition:
        partition_sequences.add(m)
        # print("partition size:", len(partition[m]))
        for s in  partition[m]:
            partition_sequences.add(s)
    # if the total number of lengths in partition is equal to the original number of strings in s
    # and the number of unique strings in Partition is the same as in S, then partition is a proper partition S
    # That is, there are no bugs.
    # print(unique_start_strings == partition_sequences)
    print(total_strings_in_partition)
    # print(len(partition_sequences))
    # print(len(unique_start_strings))
    assert unique_start_strings == partition_sequences

    return partition_alignments, partition, M, converged


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

    def test_construct_minimizer_graph(self):
        S = {"1": "AAAAAAAAAAAGGGGGGGGGGAAAAAAAAAAATTTTTTTTTTTTTCCCCCCCCCCCCCCAAAAAAAAAAACCCCCCCCCCCCCGAGGAGAGAGAGAGAGAGATTTTTTTTTTTTCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC",
             "2": "AAAAATAAAAAGGGGGGGGGGAAAAAAAAAAATTTTTTTTTTTTTCCCCCCCCCCCCCCAAAAAAAAAACCCCCCCCCCCCCGAGGAGAGAGAGAGAGAGATTTTTGTTTTTTCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC",
             "3": "AAAAAAAAAAAGGGGAGGGGGAAAAAAAAAAATTTTTTTTTTTTTCCCCCCCCCCCCCAAAAAAAAAAACCCCCCCCCCCCCGAGGAGAGAGAGAGAGAGATTTTTTTCTTTTTCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC",
             "4": "AAAAAAAAAAAGGGGGGGGGGAAATAAAAAAATTTTTTTTTTTTTCCCCCCCCCCCCCAAAAAAAAAAACCCCCCCCCCCCCGAGGAGAGACAGAGAGAGATTTTTTTTTTTTCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"} 
        G_star, alignment_graph, converged = construct_minimizer_graph(S)
        # print(G_star)
        # print(alignment_graph)
        # self.assertEqual(G_star, G_star)
        # self.assertEqual(alignment_graph, alignment_graph)
        from input_output import fasta_parser
        try:
            fasta_file_name = "/Users/kxs624/Documents/data/pacbio/simulated/ISOseq_sim_n_25/simulated_pacbio_reads.fa"
            # fasta_file_name = "/Users/kxs624/Documents/data/pacbio/simulated/DAZ2_2_exponential_constant_0.001.fa"
            S = {acc: seq for (acc, seq) in  fasta_parser.read_fasta(open(fasta_file_name, 'r'))} 
        except:
            print("test file not found:",fasta_file_name)
            return
        G_star, alignment_graph, converged = construct_minimizer_graph(S)
        edit_distances = []
        nr_unique_minimizers = []
        for s1 in alignment_graph:
            # print("nr minimizers:", len(alignment_graph[s1]))
            # nr_minimizers.append(sum([ count for nbr, count in  G_star[s1].items()]))
            nr_unique_minimizers.append(len(G_star[s1].items()))
            # if len(alignment_graph[s1]) > 20:
            #     for s2 in alignment_graph[s1]:
            #         print(alignment_graph[s1][s2][0])
            #     print("------------------------")
            for s2 in alignment_graph[s1]:
                # print("edit distance:", alignment_graph[s1][s2][0])
                edit_distances.append(alignment_graph[s1][s2][0])
                # if len(alignment_graph[s1]) > 1:
                #     print(alignment_graph[s1][s2][0])
                # if alignment_graph[s1][s2][0] == 0:
                #     print("perfect match of :", G_star[s1][s2], "seqs" )
            assert len(alignment_graph[s1]) == len(G_star[s1])

        print(sorted(nr_unique_minimizers, reverse = True))
        print(sorted(edit_distances))
        # print(G_star)
        # print(alignment_graph)

    def test_partition_strings(self):
        from input_output import fasta_parser
        try:
            fasta_file_name = "/Users/kxs624/Documents/data/pacbio/simulated/ISOseq_sim_n_200/simulated_pacbio_reads.fa"
            # fasta_file_name = "/Users/kxs624/Documents/data/pacbio/simulated/DAZ2_2_exponential_constant_0.001.fa"
            # fasta_file_name = "/Users/kxs624/Documents/data/pacbio/simulated/TSPY13P_2_constant_constant_0.0001.fa"
            # fasta_file_name = "/Users/kxs624/Documents/data/pacbio/simulated/TSPY13P_4_linear_exponential_0.05.fa"
            S = {acc: seq for (acc, seq) in  fasta_parser.read_fasta(open(fasta_file_name, 'r'))} 
        except:
            print("test file not found:",fasta_file_name)    
        
        partition_alignments, partition, M, converged = partition_strings(S)
        print(len(M), len(partition), converged)

if __name__ == '__main__':
    unittest.main()