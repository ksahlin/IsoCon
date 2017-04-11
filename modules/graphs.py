"""
    minimizer_graph(S): creates the minimizer graph defined in..
    partition_graph(S): creates the partition of a graph defined in..

"""

import unittest


from modules import get_best_alignments
from modules import minimap_alignment_module
from modules import functions
from modules import minimizer_graph

def construct_exact_minimizer_graph(S, params):

    """
        input: a dict of strings, not necesarily unique
        output: a directed graph implemented as a dict of dicts. Each edge has a weight assosiated to them.
                self edges has a weight > 1 (identical sequences) and all other edges has weight 1.
                Note, a node can be isolated! An isolated node will point at itself, effectively having itself as a minimizer.

    """
    G_star = {}
    # adding self edges to strings that has converged
    for acc, s in S.items():
        if s not in G_star:
            G_star[s] = {}
        else:
            if s in G_star[s]:
                G_star[s][s] += 1  
                print(acc)
            else:
                G_star[s][s] = 2
                print(acc)

    # check if converged, that is, if all nodes has self edges here, there will be no other edges added.
    converged = False
    not_in_clusters = set()
    already_converged = set()
    for s, nbr_dict in G_star.items():
        if len(nbr_dict) == 0:
            not_in_clusters.add(s)
        else:
            already_converged.add(s)

    if len(not_in_clusters) == 0:
        converged = True
        return G_star, converged

    unique_strings = {seq : acc for acc, seq in S.items()}
    S_prime = {acc : seq for seq, acc in unique_strings.items()}
    all_internode_edges_in_minimizer_graph, isolated_nodes = minimizer_graph.compute_minimizer_graph(S_prime, params) # send in a list of nodes that already has converged, hence avoid unnnecessary computation
    # minimizer_graph, isolated_nodes = compute_minimizer_graph(S, already_converged, params) # send in a list of nodes that already has converged, hence avoid unnnecessary computation
    for s1_acc in all_internode_edges_in_minimizer_graph:
        s1 = S[s1_acc]
        if s1 in G_star[s1]: # the minimizer had already identical minimizer (ed = 0)
            continue
        else:
            for s2_acc in all_internode_edges_in_minimizer_graph[s1_acc]:
                s2 = S[s2_acc]
                # G_star[s1][s2] = 1 
                G_star[s1][s2] = all_internode_edges_in_minimizer_graph[s1_acc][s2_acc]

    for s in isolated_nodes:
        assert s in G_star
        if s not in G_star[s]:
            G_star[s][s] = 1

    return G_star, converged    

def construct_minimizer_graph_approximate(S, params, edge_creating_min_treshold = -1, edge_creating_max_treshold = 2**30):

    """
        input: a dict of strings, not necesarily unique
        output: a directed graph implemented as a dict of dicts. Each edge has a weight assosiated to them.
                self edges has a weight > 1 (identical sequences) and all other edges has weight 1.
                Note, a node can be isolated!
    """
    G_star = {}
    alignment_graph = {}
    # adding self edges to strings that has converged
    s_to_acc = {s : acc for acc, s in S.items()} 
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
    converged = False
    not_in_clusters = set()
    for s, nbr_dict in G_star.items():
        if len(nbr_dict) == 0:
            not_in_clusters.add(s)

    if len(not_in_clusters) == 0:
        converged = True
        return G_star, alignment_graph, converged

    unique_strings = set(S.values())

    G_star_exact, _ = construct_exact_minimizer_graph(S, params)

    paf_files, acc_to_strings = minimap_alignment_module.minimap_partition(unique_strings, not_in_clusters, params)
    approximate_matches = minimap_alignment_module.paf_to_best_matches(paf_files, acc_to_strings)
    best_exact_matches = get_best_alignments.find_best_matches(approximate_matches,  edge_creating_min_treshold = edge_creating_min_treshold, edge_creating_max_treshold = edge_creating_max_treshold )

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
            # print("ISOLATED")



    for s_acc in G_star:
        if s_acc in G_star[s_acc]:
            assert a_acc in G_star_exact[s_acc]
            continue
        s_approx_nbrs = len(G_star[s_acc])
        s_approx_ed = sum( [alignment_graph[s_acc][s_nbr][0] for s_nbr in alignment_graph[s_acc] ])
        s_exact_nbrs = len(G_star_exact[s_acc])
        s_exact_ed = sum( [G_star_exact[s_acc][s_nbr] for s_nbr in G_star_exact[s_acc]])
        if s_approx_nbrs != s_exact_nbrs or s_exact_ed != s_approx_ed:
            # print("appro", G_star[s_acc])
            # print("exact", G_star_exact[s_acc])
            print("Approx:", s_approx_nbrs, s_approx_ed, s_approx_ed /float(s_approx_nbrs), "exact:", s_exact_nbrs, s_exact_ed, s_exact_ed/ float(s_exact_nbrs) )
            print(s_to_acc[s_acc])
            for saae in G_star_exact[s_acc]:
                print("Exact:", s_to_acc[saae])
            for saa in G_star[s_acc]:
                print("Approx:", s_to_acc[saa])
    sys.exit()


    return G_star, converged



def construct_2set_minimizer_bipartite_graph(X, C, X_file, C_file):
    """
        X: a string pointing to a fasta file with reads  ## a dict containing original reads and their accession
        C: a string pointing to a fasta file with candidates  ## a dict containing consensus transcript candidates
    """

    # TODO: eventually filter candidates with lower support than 2-3? Here?
    paf_file_name = minimap_alignment_module.map_with_minimap(C_file, X_file)
    highest_paf_scores = minimap_alignment_module.paf_to_best_matches_2set(paf_file_name)
    best_exact_matches = get_best_alignments.find_best_matches_2set(highest_paf_scores, X, C)

    G_star = {}
    alignment_graph = {}

    for x_acc in best_exact_matches:
        G_star[x_acc] = {}
        alignment_graph[x_acc] = {}
        # if len(best_exact_matches[x_acc]) >1:
        #     print(len(best_exact_matches[x_acc]), "best matches for read to consensus", best_exact_matches[x_acc].keys())
        for c_acc in best_exact_matches[x_acc]:
            assert c_acc not in G_star[x_acc]
            G_star[x_acc][c_acc] = 1
            (edit_distance, x_alignment, c_alignment) = best_exact_matches[x_acc][c_acc]
            alignment_graph[x_acc][c_acc] = (edit_distance, x_alignment, c_alignment)

    return G_star, alignment_graph





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
        minimap_paf = map_with_minimap(temp_file_name, temp_file_name)
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