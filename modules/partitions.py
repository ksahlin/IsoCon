import networkx as nx
import argparse, os

from operator import itemgetter
from collections import defaultdict


from modules import graphs
from modules import functions

def highest_reachable_with_edge_degrees(S, params):
    G_star, converged = graphs.construct_exact_minimizer_graph_improved(S, params)
    unique_start_strings = set(G_star.nodes())
    partition_sizes = []
    nr_consensus = 0
    G_transpose = nx.reverse(G_star)
    M = {}
    partition = {}
    print("here")
    for subgraph in sorted(nx.weakly_connected_component_subgraphs(G_transpose), key=len, reverse=True):
        # print("Subgraph of size", len(subgraph.nodes()), len(subgraph) )
        while subgraph:
            # number_connected_to = {}
            reachable_comp_sizes = []
            reachable_comp_weights = {}
            reachable_comp_nodes = []
            edit_distances = []
            processed = set()
            for m in subgraph:
                if m in processed:
                    continue
                reachable_comp = set([m])
                # print("deg m", subgraph.node[m], subgraph[m])
                reachable_comp_weight = subgraph.node[m]["degree"]
                processed.add(m)
                # print("cl size:", len([n for n in nx.dfs_postorder_nodes(subgraph, source=m)]))
                for reachable_node in nx.dfs_postorder_nodes(subgraph, source=m): # store reachable node as processed here to avoid computation
                    if reachable_node == m:
                        continue
                    
                    processed.add(reachable_node)
                    reachable_comp.add(reachable_node)
                    reachable_comp_weight += subgraph.node[reachable_node]["degree"]
                    # print(subgraph.node[reachable_node])
                    if reachable_node in subgraph[m]:
                        edit_distances.append(subgraph[m][reachable_node]["edit_distance"])
                    # print(len(G[reachable_node]), len(G_transpose[reachable_node]), reachable_node == m )
                    assert subgraph.node[reachable_node]["degree"] == 1

                reachable_comp_sizes.append(len(reachable_comp))
                reachable_comp_weights[reachable_comp_weight] = (m, reachable_comp)
                reachable_comp_nodes.append(reachable_comp)
                # number_connected_to[m] = len(reachable_comp)

            sorted_reachable_comp_sizes = sorted(reachable_comp_sizes, reverse=True)
            sorted_reachable_comp_weights = sorted(reachable_comp_weights.keys(), reverse=True)
            max_weight = max(sorted_reachable_comp_weights)
            sorted_reachable_comp_nodes = sorted(reachable_comp_nodes, key = len, reverse=True)
            biggest_comp = sorted_reachable_comp_nodes[0]
            minimizer, biggest_weighted_comp = reachable_comp_weights[max_weight]
            M[minimizer] = max_weight   
            partition[minimizer] = biggest_weighted_comp.difference(set([minimizer]))
            subgraph.remove_nodes_from(biggest_weighted_comp)

            edit_distances.sort() 
            # print("Subgraph after removal size", len(subgraph.nodes()), len(subgraph), edit_distances )
            nr_consensus += 1

    print("NR CONSENSUS:", nr_consensus)
    print("NR minimizers:", len(M), len(partition))

    print("partition sizes(identical strings counted once): ", sorted([len(partition[p]) +1 for p in  partition], reverse = True))

    total_strings_in_partition = sum([ len(partition[p]) +1 for p in  partition])
    partition_sequences = set()
    for m in partition:
        partition_sequences.add(m)
        # print("partition size:", len(partition[m]))
        # print(len(m))
        for s in  partition[m]:
            partition_sequences.add(s)
            # print(len(s))
    # if the total number of lengths in partition is equal to the original number of strings in s
    # and the number of unique strings in Partition is the same as in S, then partition is a proper partition S
    # That is, there are no bugs.
    # print(unique_start_strings == partition_sequences)
    # print(total_strings_in_partition)
    # print(len(partition_sequences))
    # print(len(unique_start_strings))
    assert unique_start_strings == partition_sequences
    assert total_strings_in_partition == len(unique_start_strings)

    return G_star, partition, M, converged


def partition_strings_paths(S, params):

    G_star, converged = graphs.construct_exact_minimizer_graph(S, params)
    # G_star, converged = graphs.construct_minimizer_graph_approximate(S, params, edge_creating_min_treshold = edge_creating_min_treshold, edge_creating_max_treshold = edge_creating_max_treshold)

    unique_start_strings = set(G_star.keys())
    partition = {} # dict with a center as key and a set containing all sequences chosen to this partition

    if converged:
        M = {key : 1 for key in G_star.keys()}
        for m in G_star:
            partition[m] = set()
            indegree = G_star[m][m]
        return G_star, partition, M, converged

    marked = set()
    M = {}
    unmarked = set(G_star.keys())
    V_G = len(G_star.keys())

    partition_counter = 1

    isolated = 0
    for s in G_star:
        if s in G_star[s]:
            if  G_star[s][s] == 1: # isolate
                isolated +=1

    print("nr isolated nodes:", isolated)

    G_star_transposed = functions.transpose(G_star)

    # do_while as long as there there are unmarked nodes
    while len(marked) < V_G:

        # find node with biggest indegree
        max_indegree = -1
        max_node_weight = -1
        indegrees = []
        for s in unmarked:
            indegree = sum([indegree for in_nbr, indegree in G_star_transposed[s].items() if in_nbr not in marked ])
            indegrees.append(indegree)
            if indegree > max_indegree:
                m, max_indegree = s, indegree
        
        # print(max_indegree, len(unmarked), len(marked))
        indegrees.sort()
        # print("PARTITION INDEGREES:", indegrees)

        # mark all nodes leading to paths to m and remove them from unmarked set
        M[m] = partition_counter
        partition[m] = set()
        m_count = 1 if m not in G_star[m] else G_star[m][m]
        unmarked.remove(m)
        assert m not in marked
        marked.add(m)
        nbrs_to_visit = set([nbr for nbr in G_star_transposed[m] if nbr in unmarked])

        while nbrs_to_visit:
            # print("ok")
            current_layer_to_visit = list(nbrs_to_visit)
            # print(len(current_layer_to_visit))
            for v in current_layer_to_visit:
                # print("semi")
                assert v not in marked
                unmarked.remove(v)
                marked.add(v)
                partition[m].add(v)
                m_count += 1
                nbrs_to_visit.remove(v)

            for v in current_layer_to_visit:
                nbrs = G_star_transposed[v]
                for u in nbrs:
                    if u in unmarked:
                        nbrs_to_visit.add(u)

        # print("Center count = ", m_count)
    print("Chosen minimizers:", len(M))


    total_strings_in_partition = sum([ len(partition[p]) +1 for p in  partition])
    partition_sequences = set()
    for m in partition:
        partition_sequences.add(m)
        # print("partition size:", len(partition[m]))
        # print(len(m))
        for s in  partition[m]:
            partition_sequences.add(s)
            # print(len(s))
    # if the total number of lengths in partition is equal to the original number of strings in s
    # and the number of unique strings in Partition is the same as in S, then partition is a proper partition S
    # That is, there are no bugs.
    # print(unique_start_strings == partition_sequences)
    # print(total_strings_in_partition)
    # print(len(partition_sequences))
    # print(len(unique_start_strings))
    assert unique_start_strings == partition_sequences

    return G_star, partition, M, converged



def partition_strings_2set(X, C, X_file, C_file, params):
    """

    """

    G_star = graphs.construct_exact_2set_minimizer_bipartite_graph(X, C, X_file, C_file, params)
    # G_star, alignment_graph = graphs.construct_2set_minimizer_bipartite_graph(X, C, X_file, C_file)

    partition = {} # dict with a center as key and a set containing all sequences chosen to this partition

    partition_counter = 1
    G_star_transposed = functions.transpose(G_star)
    nr_candidates_with_hits = len(G_star_transposed.keys())
    not_visited_candidates = set(G_star_transposed.keys())

    consensus_processed = 1
    marked = set()
    M = {}
    unmarked = set(G_star.keys())

    # do_while as long as there there are unmarked nodes
    while consensus_processed <= nr_candidates_with_hits:

        # find node with biggest indegree
        max_indegree = -1
        for c in not_visited_candidates:
            indegree = sum([indegree for in_nbr, indegree in G_star_transposed[c].items() if in_nbr not in marked ])
            if indegree > max_indegree:
                m, max_indegree = c, indegree
        # print(max_indegree, len(unmarked), len(marked))

        # mark all nodes leading to paths to m and remove them from unmarked set
        M[m] = partition_counter
        partition[m] = set()
        not_visited_candidates.remove(m)
        m_count = 1
        reads_to_visit = set([x for x in G_star_transposed[m] if x in unmarked])

        for x in reads_to_visit:
            assert x not in marked
            unmarked.remove(x)
            marked.add(x)
            partition[m].add(x)
            m_count += 1

        consensus_processed += 1
        # print("Center count = ", m_count)
        # print(len(not_visited_candidates))
    print("Chosen minimizers:", len(M))

    assert not unmarked  
    assert len(M) == nr_candidates_with_hits
    assert marked == set(G_star.keys())

    no_alignmens_at_all = set(C.keys()).difference(set(G_star_transposed.keys()))
    assert len(no_alignmens_at_all) + nr_candidates_with_hits == len(C)
    print("no_alignmens_at_all:", no_alignmens_at_all)
    for c in no_alignmens_at_all:
        assert c not in partition
        partition[c] = set()

    return G_star, partition 



