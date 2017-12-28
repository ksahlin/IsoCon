from __future__ import print_function

import networkx as nx
import argparse, os

from operator import itemgetter
from collections import defaultdict
from networkx.algorithms import bipartite

from modules import graphs
from modules import functions


def highest_reachable_with_edge_degrees(S, params):
    # if params.nontargeted:
    #     G_star, converged = graphs.construct_approximate_nearest_neighbor_graph(S, params)
    # else:
    G_star, converged = graphs.construct_exact_nearest_neighbor_graph(S, params)

    unique_start_strings = set(G_star.nodes())
    partition_sizes = []
    nr_consensus = 0
    G_transpose = nx.reverse(G_star)

    if params.verbose:
        print("Nodes in nearest_neighbor ggraph:", len(G_transpose))
        print("Neighbors per nodes in nearest neighbor graph", sorted([len(G_transpose.neighbors(n)) for n in G_transpose], reverse=True))


    M = {}
    partition = {}
    # print("here")
    for subgraph in sorted(nx.weakly_connected_component_subgraphs(G_transpose), key=len, reverse=True):
        # print("Subgraph of size", len(subgraph.nodes()), "nr edges:", len(subgraph.edges()), [len(x) for x in subgraph.nodes()] )
        while subgraph:
            reachable_comp_sizes = []
            reachable_comp_weights = {}
            reachable_comp_nodes = []
            edit_distances_to_m = {"XXXXX" : 0}
            direct_neighbors = {}
            processed = set()

            biggest_reachable_comp_size = 0
            biggest_reachable_comp_weight = 0
            biggest_reachable_comp_nodes = set()
            biggest_reachable_comp_nearest_neighbor = "XXXXX"


            for m in subgraph:
                # print(len(m))
                edit_distances_to_m[m] = 0
                
                if m in processed:
                    continue

                reachable_comp = set([m])
                reachable_comp_weight = subgraph.node[m]["degree"]
                processed.add(m)



                ####################################################
                # take all reachable nodes
                ####################################################

                for n1,n2 in nx.dfs_edges(subgraph, source=m): # store reachable node as processed here to avoid computation
                    if n2 == m:
                        continue
                    processed.add(n2)
                    reachable_comp.add(n2)
                    reachable_comp_weight += subgraph.node[n2]["degree"]
                    edit_distances_to_m[m] +=  subgraph.node[n2]["degree"] * subgraph[n1][n2]["edit_distance"]
                    assert subgraph.node[n2]["degree"] == 1
                ####################################################
                ####################################################
                

                # print("total component weight:", reachable_comp_weight)
                # print("edit distance:",  edit_distances_to_m[m])

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
                        if edit_distances_to_m[m] < edit_distances_to_m[biggest_reachable_comp_nearest_neighbor]:
                            # print("tie but smaller edit distance", edit_distances_to_m[m], edit_distances_to_m[biggest_reachable_comp_nearest_neighbor])
                            biggest_reachable_comp_nodes = set(reachable_comp)
                            biggest_reachable_comp_size = len(reachable_comp)
                            biggest_reachable_comp_nearest_neighbor = m

                        elif edit_distances_to_m[m] > edit_distances_to_m[biggest_reachable_comp_nearest_neighbor]:
                            # print("tie but bigger edit distance", edit_distances_to_m[m], edit_distances_to_m[biggest_reachable_comp_nearest_neighbor])
                            pass
                        else:
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
                    direct_weight = subgraph.node[n]["degree"]                    
                    direct_weight += len(list(subgraph.neighbors(n)))
                    # print( [ subgraph.node[nbr]["degree"] for nbr in subgraph.neighbors(n)])
                    assert all( [ subgraph.node[nbr]["degree"] == 1 for nbr in subgraph.neighbors(n)] )

                    # print("direct weight:", direct_weight)
                    if direct_weight > max_direct_weight:
                        max_direct_weight = direct_weight
                        nearest_neighbor = n
                    elif direct_weight == max_direct_weight:
                        nearest_neighbor = min(nearest_neighbor, n)
                # print("nearest_neighbor direct weight:", max_direct_weight, "nodes in reachable:", len(biggest_reachable_comp_nodes))
                M[nearest_neighbor] = biggest_reachable_comp_weight   
                partition[nearest_neighbor] = biggest_reachable_comp_nodes.difference(set([nearest_neighbor]))
                assert nearest_neighbor in biggest_reachable_comp_nodes

            # vizualize_test_graph(subgraph)
            # if len(biggest_reachable_comp_nodes) == 65:
            #     sys.exit()

            subgraph.remove_nodes_from(biggest_reachable_comp_nodes)


            nr_consensus += 1

    # for m in sorted(partition):
    #     print("min:", len(m))
    #     for p in sorted(partition[m]):
    #         print(len(p))
    nearest_neighbor_lenghts = [len(m) for m in sorted(partition)]
    
    if params.verbose:
        print("Seq lengths of centers", sorted(nearest_neighbor_lenghts))
        print("NUMBER CONSENSUS:", nr_consensus)
        print("Number of centers:", len(M), len(partition))
        print("partition sizes(identical strings are collapsed here and therefore counted as one): ", sorted([len(partition[p]) +1 for p in  partition], reverse = True))

    total_strings_in_partition = sum([ len(partition[p]) +1 for p in  partition])
    partition_sequences = set()
    for m in partition:
        partition_sequences.add(m)
        # print("partition size:", len(partition[m]))
        # print(len(m))
        for s in partition[m]:
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




def partition_strings_2set(X, C, X_file, C_file, params):
    """

    """

    G_star = graphs.construct_exact_2set_nearest_neighbor_bipartite_graph(X, C, X_file, C_file, params)
    # G_star, alignment_graph = graphs.construct_2set_nearest_neighbor_bipartite_graph(X, C, X_file, C_file)
    G_star_transposed = nx.reverse(G_star) #functions.transpose(G_star)
    partition = {} # dict with a center as key and a set containing all sequences chosen to this partition
    
    # candidate_nodes, read_nodes = bipartite.sets(G_star_transposed)
    
    read_nodes = set(n for n,d in G_star_transposed.nodes(data=True) if d['bipartite']==0)
    candidate_nodes = set(G_star_transposed) - read_nodes
    
    read_deg, cand_deg = bipartite.degrees(G_star_transposed, candidate_nodes)
    # print(len(read_nodes), len(candidate_nodes))
    # print(read_deg)
    # print(cand_deg)

    ######################
    while len(candidate_nodes) > 0:
        read_deg, cand_deg = bipartite.degrees(G_star_transposed, candidate_nodes)
        read_deg, cand_deg = dict(read_deg), dict(cand_deg)
        # print(type(read_deg), read_deg)
        # print(type(cand_deg), cand_deg)
        # print("reads left:", len(read_deg))
        # print("cands left:", len(cand_deg))
        m = max(sorted(cand_deg), key=lambda key: cand_deg[key])
        reads_supporting_m = list(G_star_transposed.neighbors(m))
        partition[m] = set(reads_supporting_m)
        G_star_transposed.remove_node(m)
        G_star_transposed.remove_nodes_from(reads_supporting_m)
        
        read_nodes = set(n for n,d in G_star_transposed.nodes(data=True) if d['bipartite']==0)
        candidate_nodes = set(G_star_transposed) - read_nodes
        # candidate_nodes, read_nodes = bipartite.sets(G_star_transposed)

        # print("total nodes left after removal:", len(G_star_transposed.nodes()), "tot candidate nodes left:", candidate_nodes)
        # print(read_nodes, [G_star[node] for node in read_nodes])
        # print(len(reads_supporting_m) , len(G_star_transposed.nodes()), G_star_transposed.nodes() )


    # print([ (m,len(partition[m])) for m in partition] )
    #####################

    return G_star, partition 



