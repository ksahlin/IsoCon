"""
    nearest_neighbor_graph(S): creates the nearest_neighbor graph defined in..
    partition_graph(S): creates the partition of a graph defined in..

"""

import unittest

from collections import defaultdict
from itertools import combinations

import networkx as nx

from modules import get_best_alignments
from modules import minimap_alignment_module
from modules import functions
from modules import nearest_neighbor_graph

def transform(read):
    transformed_seq = []
    prev_nucl = ""
    for nucl in read:
        if nucl != prev_nucl:
            transformed_seq.append(nucl)
        prev_nucl = nucl

    return "".join(transformed_seq)

def construct_exact_nearest_neighbor_graph(S, params):

    """
        input: a dict of strings, not necesarily unique
        output: a directed graph implemented as a dict of dicts. A node has a weight associated as the number of identical sequences.
            An edge from a node n1 to a node n2, n1 !+n2 means that n1 has weight 1 and it's nearest_neighbor is n2. 
    """

    predicted_seq_to_acc = defaultdict(list)
    for (acc, seq) in S.items():
        predicted_seq_to_acc[seq].append(acc)
    
    converged = True
    G = nx.DiGraph()
    has_converged = set()
    for seq, list_acc in predicted_seq_to_acc.items():
        deg = len(list_acc)
        G.add_node(seq, degree = deg)
        if deg > 1:
            has_converged.add(seq)

        if deg == 1:
            converged = False
    
    if converged:
        return G, converged
    
    unique_strings = {seq : acc for acc, seq in S.items()}
    S_prime = {acc : seq for seq, acc in unique_strings.items()}
    all_internode_edges_in_nearest_neighbor_graph, isolated_nodes = nearest_neighbor_graph.compute_nearest_neighbor_graph(S_prime, has_converged, params) # send in a list of nodes that already has converged, hence avoid unnnecessary computation
 

    for s1_acc in all_internode_edges_in_nearest_neighbor_graph:
        s1 = S[s1_acc]
        if G.node[s1]["degree"] > 1:
            continue
        else:
            for s2_acc in all_internode_edges_in_nearest_neighbor_graph[s1_acc]:
                s2 = S[s2_acc]
                ed = all_internode_edges_in_nearest_neighbor_graph[s1_acc][s2_acc]
                G.add_edge(s1, s2, edit_distance=ed)

    was_assigned_to_a_nearest_neighbor = set([s for s in G.nodes() if len(list(G.neighbors(s)) ) > 0 ])
    strings_converged = set([s for s in G.nodes() if G.node[s]["degree"] > 1 ])
    isolated_nodes = set(unique_strings.keys()) - (was_assigned_to_a_nearest_neighbor | strings_converged)

    if params.verbose:
        print("{0} strings was not converged and did not find a nearest_neighbor when aligning.".format(len(isolated_nodes)))
        print("{0} edges in nearest_neighbor graph".format(len(G.edges())))

    for s in isolated_nodes:
        assert s in G

    return G, converged  


def construct_approximate_nearest_neighbor_graph(S, params, edge_creating_min_treshold = -1, edge_creating_max_treshold = 2**30):

    """
        input: a dict of strings, not necesarily unique
        output: a directed graph implemented as a networkx graph object. Each edge has a weight assosiated to them.
                self edges has a weight > 1 (identical sequences) and all other edges has weight 1.
                Note, a node can be isolated!
    """

    predicted_seq_to_acc = defaultdict(list)
    for (acc, seq) in S.items():
        predicted_seq_to_acc[seq].append(acc)
    
    graph_has_converged = True
    G = nx.DiGraph()
    has_converged = set()
    for seq, list_acc in predicted_seq_to_acc.items():
        deg = len(list_acc)
        G.add_node(seq, degree = deg)
        if deg > 1:
            has_converged.add(seq)
        if deg == 1:
            graph_has_converged = False
    
    if graph_has_converged:
        return G, graph_has_converged
    s_to_acc = {s : acc for acc, s in S.items()}   
    unique_strings = set(S.values())
    # print(G.nodes(data=True))
    not_in_clusters = set([s for s in G.nodes() if G.node[s]["degree"] == 1 ])
    print("TOTAL UNIQUE STRINGS:", len(G.nodes()))
    print("TOTAL NON CONVERGED STRINGS:", len(not_in_clusters))

    ##################################
    ##################################


    paf_files, acc_to_strings = minimap_alignment_module.minimap_partition(unique_strings, not_in_clusters, params)
    approximate_matches = minimap_alignment_module.paf_to_best_matches(paf_files, acc_to_strings)
    best_exact_matches = get_best_alignments.find_best_matches(approximate_matches, params, edge_creating_min_treshold = edge_creating_min_treshold, edge_creating_max_treshold = edge_creating_max_treshold )


    ##################################
    ##################################

    for s1 in best_exact_matches:
        if G.node[s1]["degree"] > 1:
            continue
        else:
            for s2 in best_exact_matches[s1]:
                ed = best_exact_matches[s1][s2][0]
                G.add_edge(s1, s2, edit_distance=ed)

    was_assigned_to_a_nearest_neighbor = set([s for s in G.nodes() if len(list(G.neighbors(s))) > 0 ])
    strings_converged = set([s for s in G.nodes() if G.node[s]["degree"] > 1 ])
    isolated_nodes = unique_strings - (was_assigned_to_a_nearest_neighbor | strings_converged)
    
    if params.verbose:
        print("{0} strings was not converged and did not find a nearest_neighbor when aligning.".format(len(isolated_nodes)))
        print("{0} edges in nearest_neighbor graph".format(len(G.edges())))

    for s in isolated_nodes:
        assert s in G

    return G, graph_has_converged


def construct_exact_2set_nearest_neighbor_bipartite_graph(X, C, X_file, C_file, params):

    best_exact_matches = nearest_neighbor_graph.compute_2set_nearest_neighbor_graph(X, C, params)
    read_layer =  best_exact_matches.keys()
    candidate_layer = set([cand for read in best_exact_matches for cand in best_exact_matches[read]])
    # G_star = {}
    G = nx.DiGraph()
    G.add_nodes_from(read_layer, bipartite=0)
    G.add_nodes_from(candidate_layer, bipartite=1)
    G.add_edges_from([(x,c) for x in best_exact_matches for c in best_exact_matches[x]])
    
    # for x_acc in best_exact_matches:
    #     # G_star[x_acc] = {}
    #     for c_acc in best_exact_matches[x_acc]:
    #         # assert c_acc not in G_star[x_acc]
    #         # assert c_acc not in G[x_acc]
    #         G.add_edge(x_acc, c_acc, edit_distance=best_exact_matches[x_acc][c_acc])
    #         # G_star[x_acc][c_acc] = 1
    #         # edit_distance = best_exact_matches[x_acc][c_acc]

    return G


# def construct_approximate_2set_nearest_neighbor_bipartite_graph(X, C, X_file, C_file, params):
#     """
#         X: a string pointing to a fasta file with reads  ## a dict containing original reads and their accession
#         C: a string pointing to a fasta file with candidates  ## a dict containing consensus transcript candidates
#     """

#     # TODO: eventually filter candidates with lower support than 2-3? Here?
#     paf_file_name = minimap_alignment_module.map_with_minimap(C_file, X_file)
#     highest_paf_scores = minimap_alignment_module.paf_to_best_matches_2set(paf_file_name)
#     best_exact_matches = get_best_alignments.find_best_matches_2set(highest_paf_scores, X, C, params)

#     G_star = {}
#     alignment_graph = {}

#     for x_acc in best_exact_matches:
#         G_star[x_acc] = {}
#         alignment_graph[x_acc] = {}
#         # if len(best_exact_matches[x_acc]) >1:
#         #     print(len(best_exact_matches[x_acc]), "best matches for read to consensus", best_exact_matches[x_acc].keys())
#         for c_acc in best_exact_matches[x_acc]:
#             assert c_acc not in G_star[x_acc]
#             G_star[x_acc][c_acc] = 1
#             (edit_distance, x_alignment, c_alignment) = best_exact_matches[x_acc][c_acc]
#             alignment_graph[x_acc][c_acc] = (edit_distance, x_alignment, c_alignment)

#     return G_star, alignment_graph





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

    def test_construct_nearest_neighbor_graph(self):
        S = {"1": "AAAAAAAAAAAGGGGGGGGGGAAAAAAAAAAATTTTTTTTTTTTTCCCCCCCCCCCCCCAAAAAAAAAAACCCCCCCCCCCCCGAGGAGAGAGAGAGAGAGATTTTTTTTTTTTCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC",
             "2": "AAAAATAAAAAGGGGGGGGGGAAAAAAAAAAATTTTTTTTTTTTTCCCCCCCCCCCCCCAAAAAAAAAACCCCCCCCCCCCCGAGGAGAGAGAGAGAGAGATTTTTGTTTTTTCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC",
             "3": "AAAAAAAAAAAGGGGAGGGGGAAAAAAAAAAATTTTTTTTTTTTTCCCCCCCCCCCCCAAAAAAAAAAACCCCCCCCCCCCCGAGGAGAGAGAGAGAGAGATTTTTTTCTTTTTCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC",
             "4": "AAAAAAAAAAAGGGGGGGGGGAAATAAAAAAATTTTTTTTTTTTTCCCCCCCCCCCCCAAAAAAAAAAACCCCCCCCCCCCCGAGGAGAGACAGAGAGAGATTTTTTTTTTTTCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"} 
        G_star, alignment_graph, converged = construct_nearest_neighbor_graph(S)
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
        G_star, alignment_graph, converged = construct_nearest_neighbor_graph(S)
        edit_distances = []
        nr_unique_nearest_neighbors = []
        for s1 in alignment_graph:
            # print("nr nearest_neighbors:", len(alignment_graph[s1]))
            # nr_nearest_neighbors.append(sum([ count for nbr, count in  G_star[s1].items()]))
            nr_unique_nearest_neighbors.append(len(G_star[s1].items()))
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

        print(sorted(nr_unique_nearest_neighbors, reverse = True))
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