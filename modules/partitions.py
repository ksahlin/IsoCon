from __future__ import print_function

import networkx as nx
import argparse, os
from time import time
import copy

from operator import itemgetter
from collections import defaultdict
from networkx.algorithms import bipartite

from modules import graphs
from modules import functions


# def get_partitions(G_transpose):
#     nr_consensus = 0
#     M = {}
#     partition = {}
#     # print("here")
#     all_weak_components = list(nx.weakly_connected_component_subgraphs(G_transpose))
#     for subgraph in sorted(all_weak_components, key=len, reverse=True):
#         # print("Subgraph of size", len(subgraph.nodes()), "nr edges:", len(subgraph.edges()), [len(x) for x in subgraph.nodes()] )
#         while subgraph:
 
#             # edit_distances_to_m = {"XXXXX" : 0}
#             processed = set()
#             biggest_reachable_comp_size = 0
#             biggest_reachable_comp_weight = 0
#             biggest_reachable_comp_nodes = set()
#             biggest_reachable_comp_center = "XXXXX"
#             biggest_reachable_comp_center_nr_nbrs = 0

#             for m in subgraph:
#                 # edit_distances_to_m[m] = 0
                
#                 if m in processed:
#                     continue

#                 reachable_comp = set([m])
#                 reachable_comp_weight = subgraph.node[m]["degree"]
#                 processed.add(m)



#                 ####################################################
#                 # get all reachable nodes
#                 ####################################################

#                 for n1,n2 in nx.dfs_edges(subgraph, source=m): # store reachable node as processed here to avoid computation
#                     if n2 == m:
#                         continue
#                     processed.add(n2)
#                     reachable_comp.add(n2)
#                     reachable_comp_weight += subgraph.node[n2]["degree"]
#                     # edit_distances_to_m[m] +=  subgraph.node[n2]["degree"] * subgraph[n1][n2]["edit_distance"]
#                     assert subgraph.node[n2]["degree"] == 1
#                 ####################################################
#                 ####################################################
                

#                 # print("total component weight:", reachable_comp_weight)
#                 # print("edit distance:",  edit_distances_to_m[m])

#                 if biggest_reachable_comp_weight == 0: # initialization for first processed m
#                     biggest_reachable_comp_weight = reachable_comp_weight
#                     biggest_reachable_comp_nodes = set(reachable_comp)
#                     biggest_reachable_comp_size = len(reachable_comp)
#                     biggest_reachable_comp_center = m
#                     biggest_reachable_comp_center_nr_nbrs = len(list(subgraph.neighbors(m)))

#                 # elif reachable_comp_weight >= biggest_reachable_comp_weight:
#                 elif reachable_comp_weight == biggest_reachable_comp_weight:
#                     # print("HEEERE!!",reachable_comp_weight, biggest_reachable_comp_weight)
#                     if biggest_reachable_comp_center_nr_nbrs < len(list(subgraph.neighbors(m))):
#                         biggest_reachable_comp_weight = reachable_comp_weight
#                         biggest_reachable_comp_nodes = set(reachable_comp)
#                         biggest_reachable_comp_size = len(reachable_comp)
#                         biggest_reachable_comp_center = m
#                         biggest_reachable_comp_center_nr_nbrs = len(list(subgraph.neighbors(m)))   
                    
#                     elif biggest_reachable_comp_center_nr_nbrs == len(list(subgraph.neighbors(m))):
#                         if m < biggest_reachable_comp_center: # just pick lexicographically smallest to remove non determinism
#                             biggest_reachable_comp_weight = reachable_comp_weight
#                             biggest_reachable_comp_nodes = set(reachable_comp)
#                             biggest_reachable_comp_size = len(reachable_comp)
#                             biggest_reachable_comp_center = m
#                             biggest_reachable_comp_center_nr_nbrs = len(list(subgraph.neighbors(m)))                             


#                     # if edit_distances_to_m[m] < edit_distances_to_m[biggest_reachable_comp_center]:
#                     #     # print("tie but smaller edit distance", edit_distances_to_m[m], edit_distances_to_m[biggest_reachable_comp_center])
#                     #     biggest_reachable_comp_nodes = set(reachable_comp)
#                     #     biggest_reachable_comp_size = len(reachable_comp)
#                     #     biggest_reachable_comp_center = m

#                     # elif edit_distances_to_m[m] > edit_distances_to_m[biggest_reachable_comp_center]:
#                     #     # print("tie but bigger edit distance", edit_distances_to_m[m], edit_distances_to_m[biggest_reachable_comp_center])
#                     #     pass
#                     # else:
#                     #     if biggest_reachable_comp_weight > 1:
#                     #         # print("tie both in weighted partition size and total edit distance. Choosing lexographically smaller center")
#                     #         # print(" weighted partition size:", biggest_reachable_comp_weight, " total edit distance:", edit_distances_to_m[m])
#                     #         pass
                        
#                     #     if m < biggest_reachable_comp_center:
#                     #         biggest_reachable_comp_nodes = set(reachable_comp)
#                     #         biggest_reachable_comp_center = m
#                     #     else:
#                     #         pass

#                 elif biggest_reachable_comp_weight < reachable_comp_weight:
#                     biggest_reachable_comp_weight = reachable_comp_weight
#                     biggest_reachable_comp_nodes = set(reachable_comp)
#                     biggest_reachable_comp_size = len(reachable_comp)
#                     biggest_reachable_comp_center = m
#                     biggest_reachable_comp_center_nr_nbrs = len(list(subgraph.neighbors(m)))                  


#             if biggest_reachable_comp_weight == 0: # if there were no edges! partition is center itself
#                 M[m] = 0 
#                 partition[m] = set()
#             else:
#                 center = biggest_reachable_comp_center # "XXXXXX" #biggest_reachable_comp_center #
#                 max_direct_weight = 0
#                 # print("total nodes searched in this pass:", len(biggest_reachable_comp_nodes))
#                 for n in biggest_reachable_comp_nodes:
#                     direct_weight = subgraph.node[n]["degree"]                    
#                     direct_weight += len(list(subgraph.neighbors(n)))

#                     # if len(list(subgraph.neighbors(n))) > 1 and n != center and len(list(subgraph.neighbors(center))) > 1:
#                     #     print(n in G_transpose[center], center in G_transpose[n], len(list(subgraph.neighbors(center))), len(list(subgraph.neighbors(n)))) #[n]["edit_distance"])
#                     #     if n in G_transpose[center]:
#                     #         print("ed", G_transpose[center][n]["edit_distance"])
#                     #         print(n)
#                     #         print(center)
#                     # print( [ subgraph.node[nbr]["degree"] for nbr in subgraph.neighbors(n)])

#                     assert all( [ subgraph.node[nbr]["degree"] == 1 for nbr in subgraph.neighbors(n)] )

#                     # print("direct weight:", direct_weight)
#                     if direct_weight > max_direct_weight:
#                         max_direct_weight = direct_weight
#                         center = n
#                     elif direct_weight == max_direct_weight:
#                         center = min(center, n)
#                 # print("center direct weight:", max_direct_weight, "nodes in reachable:", len(biggest_reachable_comp_nodes))
#                 M[center] = biggest_reachable_comp_weight   
#                 partition[center] = biggest_reachable_comp_nodes.difference(set([center]))
#                 assert center in biggest_reachable_comp_nodes

#             # vizualize_test_graph(subgraph)
#             # if len(biggest_reachable_comp_nodes) == 65:
#             #     sys.exit()

#             subgraph.remove_nodes_from(biggest_reachable_comp_nodes)


#             nr_consensus += 1
#     return M, partition




# def reachable(G, m):
#     reachable_new = set([n2 for n1,n2 in nx.dfs_edges(G, source=m) if n2 != m ])    
#     return reachable_new

# def weakly_connected_components(G, m):
    
#     reachable_new = set([n2 for n1,n2 in nx.dfs_edges(G, source=m) if n2 != m ])
#     # assert reachable == reachable_new
    
#     return reachable_new


# import heapq 
# def get_partitions_new2(G):
#     M_temp = {}
#     partition_temp = {}
#     G_tmp = copy.deepcopy(G)

#     reachable_for_nodes = [ ( - len(reachable(G_tmp, n)), n) for n in  sorted(G_tmp.nodes())]
#     heapq.heapify(reachable_for_nodes)
#     # highest_degree = sorted( reachable_for_nodes, key = lambda x: ( - (len(x[0]) + G.node[ x[1] ]["degree"]), x[1] ) ) # sort on partition weight first then string if tiebreakers

#     while len(reachable_for_nodes) > 0 and G_tmp:
#         (nr_reachable, m) = heapq.heappop(reachable_for_nodes)
#         nr_reachable = - nr_reachable # encoded as - because heapq sorts smallest to largest       
#         if m not in G_tmp:
#             continue
#         m_reachable = reachable(G_tmp, m)

#         ######################################################
#         # verify that m is still largest after removal of nodes
#         ######################################################
#         if reachable_for_nodes:
#             curr_min_reachable, m_tmp = min(reachable_for_nodes)
#             curr_min_reachable = - curr_min_reachable
#             if len(m_reachable) < curr_min_reachable:
#                 print("became smaller after removal of nodes!")
#                 heapq.heappush(reachable_for_nodes, (- len(m_reachable), m))
#                 continue

#         ######################################################
#         ######################################################

#         partition_temp[m] = set([node for node in m_reachable])
#         M_temp[m] = G_tmp.node[m]["degree"] + sum([G_tmp.node[node]["degree"] for node in partition_temp[m]])   
#         G_tmp.remove_nodes_from(set(partition_temp[m]))
#         G_tmp.remove_node(m)
#         # processed.update(partition_temp[m]) 
#         # processed.add(m)

#     M = {}
#     partition = {}

#     for m in list(M_temp.keys()):
#         if not partition_temp[m]:
#             # print("here")
#             M[m] =  G.node[m]["degree"]
#             partition[m] = set([])
#         else:
#             highest_degree = G.node[m]["degree"] +  len(list(G.neighbors(m))) # len(list([ nbr for nbr in G.neighbors(m) if nbr in partition_temp[m] ]))
#             center = m
#             for n in partition_temp[m]:
#                 n_degree = G.node[n]["degree"] + len(list(G.neighbors(n))) # len(list([ nbr for nbr in G.neighbors(n) if nbr in partition_temp[m] ]))
#                 if n_degree > highest_degree:
#                     center = n
#                     highest_degree = n_degree
#                 elif n_degree == highest_degree:
#                     center = min(center, n)


#             M[center] = highest_degree
#             partition_temp[m].add(m)
#             partition_temp[m].remove(center)
#             partition[center] = partition_temp[m]

#         del partition_temp[m]
#         del M_temp[m]

#     return M, partition



# def get_partitions_new(G):

#     M_temp = {}
#     partition_temp = {}
    
#     G_undirected = G.to_directed()
    
#     # subgraphs = [ ( weakly_connected_components(G_undirected, n), n) for n in  sorted(G.nodes())]
#     # highest_degree = sorted( subgraphs, key = lambda x: ( - (len(x[0]) + G.node[ x[1] ]["degree"]), x[1] ) ) # sort on partition weight first then string if tiebreakers
    
#     reachable_for_nodes = [ ( reachable(G, n), n) for n in  sorted(G.nodes())]
#     highest_degree = sorted( reachable_for_nodes, key = lambda x: ( - (len(x[0]) + G.node[ x[1] ]["degree"]), x[1] ) ) # sort on partition weight first then string if tiebreakers

#     processed = set()
#     for m_reachable, m in highest_degree: # TODO: do a while loop here          
#         if m in processed:
#             continue

#         partition_temp[m] = set([node for node in m_reachable if node not in processed])
#         M_temp[m] = G.node[m]["degree"] + sum([G.node[node]["degree"] for node in partition_temp[m]])   
#         processed.update(partition_temp[m]) 
#         processed.add(m)

#     M = {}
#     partition = {}

#     for m in list(M_temp.keys()):
#         if not partition_temp[m]:
#             # print("here")
#             M[m] =  G.node[m]["degree"]
#             partition[m] = set([])
#         else:
#             highest_degree = G.node[m]["degree"] +  len(list(G.neighbors(m))) # len(list([ nbr for nbr in G.neighbors(m) if nbr in partition_temp[m] ]))
#             center = m
#             for n in partition_temp[m]:
#                 # print(len(nx.shortest_path(G,source=m,target=n)))
#                 n_degree = G.node[n]["degree"] + len(list(G.neighbors(n))) # len(list([ nbr for nbr in G.neighbors(n) if nbr in partition_temp[m] ]))
#                 if n_degree > highest_degree:
#                     center = n
#                     highest_degree = n_degree
#                 elif n_degree == highest_degree:
#                     center = min(center, n)


#             M[center] = highest_degree
#             partition_temp[m].add(m)
#             partition_temp[m].remove(center)
#             partition[center] = partition_temp[m]

#         del partition_temp[m]
#         del M_temp[m]

#     return M, partition

def get_partitions_no_copy(G_transpose):
    nr_consensus = 0
    M = {}
    partition = {}
    # print("here")
    all_weak_components = [ c for c in nx.weakly_connected_components(G_transpose)]
    for subgraph_set in sorted(all_weak_components, key=len, reverse=True):
        while subgraph_set:
 
            processed = set()
            biggest_reachable_comp_size = 0
            biggest_reachable_comp_weight = 0
            biggest_reachable_comp_nodes = set()
            biggest_reachable_comp_center = "XXXXX"
            biggest_reachable_comp_center_nr_nbrs = 0

            # subgraph_set_nbrs_sorted = sorted( [ (m, len([n for n in G_transpose[m]])) for m in subgraph_set], key = lambda x: x[1], reverse = True)  

            for m in subgraph_set:
                # edit_distances_to_m[m] = 0
                
                if m in processed:
                    continue

                reachable_comp = set([m])
                reachable_comp_weight = G_transpose.node[m]["degree"]
                processed.add(m)



                ####################################################
                # get all reachable nodes
                ####################################################

                for n1,n2 in nx.dfs_edges(G_transpose, source=m): # store reachable node as processed here to avoid computation
                    if n2 == m:
                        continue
                    processed.add(n2)
                    reachable_comp.add(n2)
                    reachable_comp_weight += G_transpose.node[n2]["degree"]
                    assert G_transpose.node[n2]["degree"] == 1
                ####################################################
                ####################################################
                

                if biggest_reachable_comp_weight == 0: # initialization for first processed m
                    biggest_reachable_comp_weight = reachable_comp_weight
                    biggest_reachable_comp_nodes = set(reachable_comp)
                    biggest_reachable_comp_size = len(reachable_comp)
                    biggest_reachable_comp_center = m
                    biggest_reachable_comp_center_nr_nbrs = len(list(G_transpose.neighbors(m)))

                # elif reachable_comp_weight >= biggest_reachable_comp_weight:
                elif reachable_comp_weight == biggest_reachable_comp_weight:
                    # print("HEEERE!!",reachable_comp_weight, biggest_reachable_comp_weight)
                    if biggest_reachable_comp_center_nr_nbrs < len(list(G_transpose.neighbors(m))):
                        biggest_reachable_comp_weight = reachable_comp_weight
                        biggest_reachable_comp_nodes = set(reachable_comp)
                        biggest_reachable_comp_size = len(reachable_comp)
                        biggest_reachable_comp_center = m
                        biggest_reachable_comp_center_nr_nbrs = len(list(G_transpose.neighbors(m)))   
                    
                    elif biggest_reachable_comp_center_nr_nbrs == len(list(G_transpose.neighbors(m))):
                        if m < biggest_reachable_comp_center: # just pick lexicographically smallest to remove non determinism
                            biggest_reachable_comp_weight = reachable_comp_weight
                            biggest_reachable_comp_nodes = set(reachable_comp)
                            biggest_reachable_comp_size = len(reachable_comp)
                            biggest_reachable_comp_center = m
                            biggest_reachable_comp_center_nr_nbrs = len(list(G_transpose.neighbors(m)))                             



                elif biggest_reachable_comp_weight < reachable_comp_weight:
                    biggest_reachable_comp_weight = reachable_comp_weight
                    biggest_reachable_comp_nodes = set(reachable_comp)
                    biggest_reachable_comp_size = len(reachable_comp)
                    biggest_reachable_comp_center = m
                    biggest_reachable_comp_center_nr_nbrs = len(list(G_transpose.neighbors(m)))                  


            if biggest_reachable_comp_weight == 0: # if there were no edges! partition is center itself
                M[m] = 0 
                partition[m] = set()
            else:
                center = biggest_reachable_comp_center # "XXXXXX" #biggest_reachable_comp_center #
                max_direct_weight = 0
                # print("total nodes searched in this pass:", len(biggest_reachable_comp_nodes))
                for n in biggest_reachable_comp_nodes:
                    direct_weight = G_transpose.node[n]["degree"]                    
                    direct_weight += len(list(G_transpose.neighbors(n)))

                    assert all( [ G_transpose.node[nbr]["degree"] == 1 for nbr in G_transpose.neighbors(n)] )

                    # print("direct weight:", direct_weight)
                    if direct_weight > max_direct_weight:
                        max_direct_weight = direct_weight
                        center = n
                    elif direct_weight == max_direct_weight:
                        center = min(center, n)
                # print("center direct weight:", max_direct_weight, "nodes in reachable:", len(biggest_reachable_comp_nodes))
                M[center] = biggest_reachable_comp_weight   
                partition[center] = biggest_reachable_comp_nodes.difference(set([center]))
                assert center in biggest_reachable_comp_nodes

            # vizualize_test_graph(G_transpose)
            # if len(biggest_reachable_comp_nodes) == 65:
            #     sys.exit()

            G_transpose.remove_nodes_from(biggest_reachable_comp_nodes)
            subgraph_set = subgraph_set -  biggest_reachable_comp_nodes

            nr_consensus += 1
    return M, partition


def partition_strings(S, params):
    # if params.nontargeted:
    #     G_star, converged = graphs.construct_approximate_nearest_neighbor_graph(S, params)
    # else:
    G_star, converged = graphs.construct_exact_nearest_neighbor_graph(S, params)

    unique_start_strings = set(G_star.nodes())
    partition_sizes = []
    nr_consensus = 0
    G_transpose = nx.reverse(G_star)

    if params.verbose:
        print("Nodes in nearest_neighbor graph:", len(G_transpose))
        print("Neighbors per nodes in nearest neighbor graph", sorted([len(list(G_transpose.neighbors(n))) for n in G_transpose], reverse=True))


    # all_weak_components = list(nx.weakly_connected_component_subgraphs(G_transpose))
    # print([ len(subG.nodes()) for subG in sorted(all_weak_components, key = len, reverse = True) ] )
    # G_undirected = G_transpose.to_directed()
    # subgraphs = [ ( weakly_connected_components(G_undirected, n), n) for n in  sorted(G_undirected.nodes())]
    # print([len(ss[0]) for ss in sorted(subgraphs, key = lambda x: len(x[0]), reverse = True) ])
    # subgraphs = nx.weakly_connected_components(G_transpose)
    # print([len(c) for c in sorted(subgraphs,key=len, reverse=True)])
    # # print(sorted( subgraphs, key = lambda x: ( - (len(x[0]) + G.node[ x[1] ]["degree"]), x[1] ) ))
    # # print(set(all_weak_components))

    # start = time()
    # M, partition = get_partitions_new(G_transpose)
    # print("tot time new partition:", time() - start )

    # print("Number of centers:", len(M), len(partition))
    # p_lengths = sorted([len(partition[p]) +1 for p in  partition], reverse = True)
    # print("partition_new sizes(identical strings are collapsed here and therefore counted as one): ", p_lengths )
    # print("SUM PARTITIONS:", sum(sorted([len(partition[p]) +1 for p in  partition], reverse = True)))
    # print()

    # start = time()
    # M2, partition2 = get_partitions_new2(G_transpose)
    # print("tot time new2 partition:", time() - start )
    # p2_lengths = sorted([len(partition2[p]) +1 for p in  partition2], reverse = True)

    # print("Number of new2 centers:", len(M2), len(partition2))
    # print("partition_new2 sizes(identical strings are collapsed here and therefore counted as one): ", p2_lengths)
    # print("SUM PARTITIONS:", sum(sorted([len(partition2[p]) +1 for p in  partition2], reverse = True)))


    # if params.verbose:
    #     center_lenghts = [len(m) for m in sorted(partition_new)]
    #     print("Seq lengths of centers new:", sorted(center_lenghts))
    #     print("Number of centers new4:", len(M_new), len(partition_new))
    #     print("partition_new4 sizes(identical strings are collapsed here and therefore counted as one): ", sorted([len(partition_new[p]) +1 for p in  partition_new], reverse = True))
    #     print("SUM PARTITIONS:", sum(sorted([len(partition_new[p]) +1 for p in  partition_new], reverse = True)))

    # start = time()
    # M_new2, partition_new2 = get_partitions_new2(G_transpose)
    # print("tot time new2 partition:", time() - start )

    # start = time()
    # list(nx.weakly_connected_component_subgraphs(G_transpose))
    # print("tot time copy:", time() - start )

    # start = time()
    # list(nx.weakly_connected_component_subgraphs(G_transpose, copy= False))
    # print("tot time without copy:", time() - start )
    
    start = time()
    G_transpose_tmp =  copy.deepcopy(G_transpose)
    print("tot time copy G_transpose:", time() - start )

    start = time()
    M, partition = get_partitions_no_copy(G_transpose_tmp)
    print("tot time partition no copy:", time() - start )

    # start = time()
    # M_old, partition_old = get_partitions(G_transpose)
    # print("tot time old partition:", time() - start )
    # print("Number of old centers:", len(M_old), len(partition_old))
    # p_old_lengths = sorted([len(partition_old[p]) +1 for p in  partition_old], reverse = True)
    # print("partition_old sizes(identical strings are collapsed here and therefore counted as one): ", p_old_lengths)
    # print("SUM PARTITIONS:", sum(sorted([len(partition_old[p]) +1 for p in  partition_old], reverse = True)))

    # assert sorted([len(partition_old[s]) for s in partition_old]) == sorted([len(partition[s]) for s in partition])
    # partition_old_lengths = sorted([len(partition_old[p]) +1 for p in  partition_old], reverse = True)
    # partition_old_nc_lengths = sorted([len(partition_old_nc[p]) +1 for p in  partition_old_nc], reverse = True)
    # print(partition_old_lengths)
    # print(partition_old_nc_lengths)
    # print( len(set(M_old_nc.keys()) ^ set(M_old.keys())) )
    # if M_old_nc != M_old:
    #     print(len(set(M_old_nc.keys()) ^ set(M_old.keys())))

    # if partition_old_nc != partition_old:
    #     print("partitions not equal")
    # M, partition = M_old_nc, partition_old_nc


    # if list(M.keys()) == list(M2.keys()):
    #     print("M equal to M2")
    # if list(M_old.keys()) == list(M2.keys()):
    #     print("M_old equal to M2")
    # if list(M_old.keys()) == list(M.keys()):
    #     print("M_old equal to M")

    # if p_old_lengths == p_lengths:
    #     print("p_old equal to p")
    # if p_old_lengths == p2_lengths:
    #     print("p_old equal to p2")
    # if p2_lengths == p_lengths:
    #     print("p2_old equal to p")


    # print(set(M.keys()) ^ set(M_old.keys()))
    # print(set(M2.keys()) - set(M_old.keys()))
    # print(set(M_old.keys()) - set(M2.keys()))
    # print(sorted([len(partition[m]) for m in M]))
    # print(sorted([len(partition_old[m]) for m in M_old]))

    # print("deg in m:", sorted([len(partition[m]) for m in set(M) ^ set(M_old) if m in M]) )
    # print("node deg in m:", [G_transpose.node[m]["degree"] + len(list(G_transpose.neighbors(m))) for m in set(M) ^ set(M_old) if m in M] )
    # print("node nbrs in m:", [len(list(G_transpose.neighbors(m))) for m in set(M) ^ set(M_old) if m in M] )
    
    # print( "deg in m_old:", sorted([len(partition_old[m]) for m in set(M) ^ set(M_old) if m in M_old]) )
    # print( "node deg in m_old:", [G_transpose.node[m]["degree"] + len(list(G_transpose.neighbors(m))) for m in set(M) ^ set(M_old) if m in M_old] )
    # print( "node nbrs in m_old:", [len(list(G_transpose.neighbors(m))) for m in set(M) ^ set(M_old) if m in M_old] )
    # if set(M) ^ set(M_old):
    #     m_old = [m for m in set(M) ^ set(M_old) if m in M_old][0]
    #     m = [m for m in set(M) ^ set(M_old) if m in M][0]
    #     print(G_transpose.has_edge(m,m_old))
    #     print(G_transpose.has_edge(m_old,m))
    #     print(len(partition_old[m_old]))
    #     print(len(partition[m]))
    #     m_old_nns = [n for n in G_transpose.neighbors(m_old) if G_transpose.has_edge(n,m_old)]
    #     m_nns = [n for n in G_transpose.neighbors(m) if G_transpose.has_edge(n,m)]
    #     print(len( list( set(G_transpose.neighbors(m)) & set(G_transpose.neighbors(m_old))) ))
    #     print(len( list( set(G_transpose.neighbors(m)) ^ set(G_transpose.neighbors(m_old))) ))
    #     m_in_old = [len(partition_old[n]) for n in partition_old if m in partition_old[n] ]
    #     print("m in old, size:", m_in_old )
    #     m_old_in_new = [len(partition[n]) for n in partition if m_old in partition[n] ]
    #     print("m_old in new:", m_old_in_new )

    #     print("shortest_path:",len(nx.shortest_path(G_transpose,source=m,target=m_old)))
    #     print("shortest_path:",len(nx.shortest_path(G_transpose,source=m_old,target=m)))
    #     # H = G_transpose.subgraph(list(G_transpose.neighbors(m)) + list(G_transpose.neighbors(m_new)) + [m, m_new])
    #     # H = G_transpose.subgraph(list( set(G_transpose.neighbors(m)) ^ set(G_transpose.neighbors(m_new))) + [m, m_new])
    #     H = G_transpose.subgraph(m_new_nns + m_nns + [m, m_new])
    #     pos = nx.spring_layout(H)
    #     nx.draw_networkx_nodes(H, pos, node_size = 50)
    #     nx.draw_networkx_edges(H, pos) #, arrows=True)
    #     import matplotlib.pyplot as plt
    #     plt.show()
    # assert set(M.keys()) == set(M_old.keys())
    # assert partition == partition_old  

    center_lenghts = [len(m) for m in sorted(partition)]
    if params.verbose:
        print("Seq lengths of centers", sorted(center_lenghts))
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



