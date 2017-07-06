import networkx as nx
import argparse, os

from operator import itemgetter
from collections import defaultdict


from modules import graphs
from modules import functions


def vizualize_test_graph(G):
    import networkx as nx
    import matplotlib.pyplot as plt
    nx.draw_networkx(G, with_labels=False, node_size=30)
    # nx.draw_networkx_nodes(G, node_size=50 )
    # nx.draw_networkx_edge_labels(G, pos, arrows=True, edge_labels=labels)
    # nx.draw_networkx_edges(G, arrows=True, edge_labels=labels)
    fig_file = "/Users/kxs624/tmp/LOLGRAPH"
    plt.savefig(fig_file, format="PNG")
    plt.clf()

def highest_reachable_with_edge_degrees(S, params):
    G_star, converged = graphs.construct_exact_minimizer_graph_improved(S, params)
    unique_start_strings = set(G_star.nodes())

    print("len G_star:", len(G_star))
    partition_sizes = []
    nr_consensus = 0
    G_transpose = nx.reverse(G_star)
    print("len G_star_transposed (minimizers):", len(G_transpose))
    print(sorted([len(G_transpose.neighbors(n)) for n in G_transpose], reverse=True))
    M = {}
    partition = {}
    print("here")
    for subgraph in sorted(nx.weakly_connected_component_subgraphs(G_transpose), key=len, reverse=True):
        print("Subgraph of size", len(subgraph.nodes()), "nr edges:", len(subgraph.edges()), [len(x) for x in subgraph.nodes()] )
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
            biggest_reachable_comp_minimizer = "XXXXX"


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
                    biggest_reachable_comp_minimizer = m

                elif reachable_comp_weight >= biggest_reachable_comp_weight:
                    if reachable_comp_weight > biggest_reachable_comp_weight:
                        biggest_reachable_comp_weight = reachable_comp_weight
                        biggest_reachable_comp_nodes = set(reachable_comp)
                        biggest_reachable_comp_size = len(reachable_comp)
                        biggest_reachable_comp_minimizer = m

                    elif reachable_comp_weight == biggest_reachable_comp_weight:
                        if edit_distances_to_m[m] < edit_distances_to_m[biggest_reachable_comp_minimizer]:
                            # print("tie but smaller edit distance", edit_distances_to_m[m], edit_distances_to_m[biggest_reachable_comp_minimizer])
                            biggest_reachable_comp_nodes = set(reachable_comp)
                            biggest_reachable_comp_size = len(reachable_comp)
                            biggest_reachable_comp_minimizer = m

                        elif edit_distances_to_m[m] > edit_distances_to_m[biggest_reachable_comp_minimizer]:
                            # print("tie but bigger edit distance", edit_distances_to_m[m], edit_distances_to_m[biggest_reachable_comp_minimizer])
                            pass
                        else:
                            if biggest_reachable_comp_weight > 1:
                                # print("tie both in weighted partition size and total edit distance. Choosing lexographically smaller minimizer")
                                # print(" weighted partition size:", biggest_reachable_comp_weight, " total edit distance:", edit_distances_to_m[m])
                                pass
                            
                            if m < biggest_reachable_comp_minimizer:
                                biggest_reachable_comp_nodes = set(reachable_comp)
                                biggest_reachable_comp_minimizer = m
                            else:
                                pass

                    else:
                        print("BUG!")

            if biggest_reachable_comp_weight == 0: # if there were no edges! partition is minimizer itself
                M[m] = 0 
                partition[m] = set()
            else:
                minimizer = biggest_reachable_comp_minimizer # "XXXXXX" #biggest_reachable_comp_minimizer #
                max_direct_weight = 0
                print("total nodes searched in this pass:", len(biggest_reachable_comp_nodes))
                for n in biggest_reachable_comp_nodes:
                    direct_weight = subgraph.node[n]["degree"]                    
                    direct_weight += len(subgraph.neighbors(n))
                    # print( [ subgraph.node[nbr]["degree"] for nbr in subgraph.neighbors(n)])
                    assert all( [ subgraph.node[nbr]["degree"] == 1 for nbr in subgraph.neighbors(n)] )

                    # print("direct weight:", direct_weight)
                    if direct_weight > max_direct_weight:
                        max_direct_weight = direct_weight
                        minimizer = n
                    elif direct_weight == max_direct_weight:
                        minimizer = min(minimizer, n)
                print("minimizer direct weight:", max_direct_weight, "nodes in reachable:", len(biggest_reachable_comp_nodes))
                M[minimizer] = biggest_reachable_comp_weight   
                partition[minimizer] = biggest_reachable_comp_nodes.difference(set([minimizer]))
                assert minimizer in biggest_reachable_comp_nodes

            # vizualize_test_graph(subgraph)
            # if len(biggest_reachable_comp_nodes) == 65:
            #     sys.exit()

            subgraph.remove_nodes_from(biggest_reachable_comp_nodes)


            nr_consensus += 1

    # for m in sorted(partition):
    #     print("min:", len(m))
    #     for p in sorted(partition[m]):
    #         print(len(p))
    minimizer_lenghts = [len(m) for m in sorted(partition)]
    print(sorted(minimizer_lenghts))
    print("NR CONSENSUS:", nr_consensus)
    print("NR minimizers:", len(M), len(partition))
    print("partition sizes(identical strings counted once): ", sorted([len(partition[p]) +1 for p in  partition], reverse = True))
    # sys.exit()

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

# def highest_reachable_with_edge_degrees_NEW(S, params):
#     G_star, converged = graphs.construct_exact_minimizer_graph_improved(S, params)
#     unique_start_strings = set(G_star.nodes())
#     from operator import itemgetter
#     print("len G_star:", len(G_star))
#     partition_sizes = []
#     nr_consensus = 0
#     G_transpose = nx.reverse(G_star)
#     print("len G_star_transposed (minimizers):", len(G_transpose))
#     print(sorted([len(G_transpose.neighbors(n)) for n in G_transpose], reverse=True))
#     M = {}
#     partition = {}
#     print("here")
#     for subgraph in sorted(nx.weakly_connected_component_subgraphs(G_transpose), key=len, reverse=True):
#         print("Subgraph of size", len(subgraph.nodes()), "nr edges:", len(subgraph.edges()), [len(x) for x in subgraph.nodes()] )
        
#         while subgraph:
#             print([nn for string, nn in sorted(subgraph.out_degree_iter(),key=itemgetter(1),reverse=True)])
#             # sys.exit()
#             m, number_nbrs = sorted(subgraph.out_degree_iter(),key=itemgetter(1),reverse=True)[0]
#             print("NN:", number_nbrs)
#             ####################################################
#             # cut at all nodes with more than 2 minimizers
#             ####################################################
#             # minimizer_w_more_than_2_support = set()
#             nbrs_to_visit = set(subgraph.neighbors(m))        
#             print("to visit", len(nbrs_to_visit))
#             visited = set([m])
#             in_partition = set([m])
#             while nbrs_to_visit:
#                 nbr = nbrs_to_visit.pop()
#                 if nbr in visited:
#                     continue

#                 visited.add(nbr)
#                 extended_nbrs = subgraph.neighbors(nbr)
#                 next_layer_nbrs = set(extended_nbrs) - set(visited)
                
#                 if len(next_layer_nbrs) > 1: # do not add the nbr or its nbrs to reachable
#                     continue
#                 elif len(next_layer_nbrs) == 1:
#                     n_nbr = next_layer_nbrs.pop()
#                     nbrs_to_visit.add(n_nbr)
#                     in_partition.add(nbr)
#                     assert subgraph.node[nbr]["degree"] == 1
#                 elif len(next_layer_nbrs) == 0:
#                     in_partition.add(nbr)
#                     assert subgraph.node[nbr]["degree"] == 1
#                     pass

#             ####################################################
#             ####################################################

#             print("In partition of m:", len(in_partition))
#             M[m] = len(in_partition) + subgraph.node[m]["degree"]
#             partition[m] = in_partition.difference(set([m]))


#             # vizualize_test_graph(subgraph)
#             # if len(biggest_reachable_comp_nodes) == 65:
#             #     sys.exit()

#             subgraph.remove_nodes_from(in_partition)

#             nr_consensus += 1

#     # for m in sorted(partition):
#     #     print("min:", len(m))
#     #     for p in sorted(partition[m]):
#     #         print(len(p))
#     minimizer_lenghts = [len(m) for m in sorted(partition)]
#     print(sorted(minimizer_lenghts))
#     print("NR CONSENSUS:", nr_consensus)
#     print("NR minimizers:", len(M), len(partition))
#     print("partition sizes(identical strings counted once): ", sorted([len(partition[p]) +1 for p in  partition], reverse = True))
#     # sys.exit()

#     total_strings_in_partition = sum([ len(partition[p]) +1 for p in  partition])
#     partition_sequences = set()
#     for m in partition:
#         partition_sequences.add(m)
#         # print("partition size:", len(partition[m]))
#         # print(len(m))
#         for s in  partition[m]:
#             partition_sequences.add(s)
#             # print(len(s))
#     # if the total number of lengths in partition is equal to the original number of strings in s
#     # and the number of unique strings in Partition is the same as in S, then partition is a proper partition S
#     # That is, there are no bugs.
#     # print(unique_start_strings == partition_sequences)
#     # print(total_strings_in_partition)
#     # print(len(partition_sequences))
#     # print(len(unique_start_strings))
#     assert unique_start_strings == partition_sequences
#     assert total_strings_in_partition == len(unique_start_strings)

#     return G_star, partition, M, converged



# def partition_strings_paths(S, params):

#     G_star, converged = graphs.construct_exact_minimizer_graph(S, params)
#     # G_star, converged = graphs.construct_minimizer_graph_approximate(S, params, edge_creating_min_treshold = edge_creating_min_treshold, edge_creating_max_treshold = edge_creating_max_treshold)

#     unique_start_strings = set(G_star.keys())
#     partition = {} # dict with a center as key and a set containing all sequences chosen to this partition

#     if converged:
#         M = {key : 1 for key in G_star.keys()}
#         for m in G_star:
#             partition[m] = set()
#             indegree = G_star[m][m]
#         return G_star, partition, M, converged

#     marked = set()
#     M = {}
#     unmarked = set(G_star.keys())
#     V_G = len(G_star.keys())

#     partition_counter = 1

#     isolated = 0
#     for s in G_star:
#         if s in G_star[s]:
#             if  G_star[s][s] == 1: # isolate
#                 isolated +=1

#     print("nr isolated nodes:", isolated)

#     G_star_transposed = functions.transpose(G_star)

#     # do_while as long as there there are unmarked nodes
#     while len(marked) < V_G:

#         # find node with biggest indegree
#         max_indegree = -1
#         max_node_weight = -1
#         indegrees = []
#         for s in unmarked:
#             indegree = sum([indegree for in_nbr, indegree in G_star_transposed[s].items() if in_nbr not in marked ])
#             indegrees.append(indegree)
#             if indegree > max_indegree:
#                 m, max_indegree = s, indegree
        
#         # print(max_indegree, len(unmarked), len(marked))
#         indegrees.sort()
#         # print("PARTITION INDEGREES:", indegrees)

#         # mark all nodes leading to paths to m and remove them from unmarked set
#         M[m] = partition_counter
#         partition[m] = set()
#         m_count = 1 if m not in G_star[m] else G_star[m][m]
#         unmarked.remove(m)
#         assert m not in marked
#         marked.add(m)
#         nbrs_to_visit = set([nbr for nbr in G_star_transposed[m] if nbr in unmarked])

#         while nbrs_to_visit:
#             # print("ok")
#             current_layer_to_visit = list(nbrs_to_visit)
#             # print(len(current_layer_to_visit))
#             for v in current_layer_to_visit:
#                 # print("semi")
#                 assert v not in marked
#                 unmarked.remove(v)
#                 marked.add(v)
#                 partition[m].add(v)
#                 m_count += 1
#                 nbrs_to_visit.remove(v)

#             for v in current_layer_to_visit:
#                 nbrs = G_star_transposed[v]
#                 for u in nbrs:
#                     if u in unmarked:
#                         nbrs_to_visit.add(u)

#         # print("Center count = ", m_count)
#     print("Chosen minimizers:", len(M))


#     total_strings_in_partition = sum([ len(partition[p]) +1 for p in  partition])
#     partition_sequences = set()
#     for m in partition:
#         partition_sequences.add(m)
#         # print("partition size:", len(partition[m]))
#         # print(len(m))
#         for s in  partition[m]:
#             partition_sequences.add(s)
#             # print(len(s))
#     # if the total number of lengths in partition is equal to the original number of strings in s
#     # and the number of unique strings in Partition is the same as in S, then partition is a proper partition S
#     # That is, there are no bugs.
#     # print(unique_start_strings == partition_sequences)
#     # print(total_strings_in_partition)
#     # print(len(partition_sequences))
#     # print(len(unique_start_strings))
#     assert unique_start_strings == partition_sequences

#     return G_star, partition, M, converged

from networkx.algorithms import bipartite


def partition_strings_2set(X, C, X_file, C_file, params):
    """

    """

    G_star = graphs.construct_exact_2set_minimizer_bipartite_graph(X, C, X_file, C_file, params)
    # G_star, alignment_graph = graphs.construct_2set_minimizer_bipartite_graph(X, C, X_file, C_file)
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
        print("reads left:", len(read_deg))
        print("cands left:", len(cand_deg))
        m = max(cand_deg, key=lambda key: cand_deg[key])
        reads_supporting_m = G_star_transposed.neighbors(m)
        partition[m] = set(reads_supporting_m)
        G_star_transposed.remove_node(m)
        G_star_transposed.remove_nodes_from(reads_supporting_m)
        
        read_nodes = set(n for n,d in G_star_transposed.nodes(data=True) if d['bipartite']==0)
        candidate_nodes = set(G_star_transposed) - read_nodes
        # candidate_nodes, read_nodes = bipartite.sets(G_star_transposed)

        print("total nodes left after removal:", len(G_star_transposed.nodes()), "tot candidate nodes left:", candidate_nodes)
        print(read_nodes, [G_star[node] for node in read_nodes])
        # print(len(reads_supporting_m) , len(G_star_transposed.nodes()), G_star_transposed.nodes() )


    # print([ (m,len(partition[m])) for m in partition] )
    #####################



    # partition_counter = 1
    # G_star_transposed = nx.reverse(G_star) #functions.transpose(G_star)
    # nr_candidates_with_hits = len(candidate_nodes)
    # not_visited_candidates = set(candidate_nodes)

    # consensus_processed = 1
    # marked = set()
    # M = {}
    # unmarked = set(read_nodes)

    # # do_while as long as there there are unmarked nodes
    # while consensus_processed <= nr_candidates_with_hits:

    #     # find node with biggest indegree
    #     max_indegree = -1
    #     for c in not_visited_candidates:
    #         indegree = sum([indegree for in_nbr, indegree in G_star_transposed[c].items() if in_nbr not in marked ])
    #         if indegree > max_indegree:
    #             m, max_indegree = c, indegree
    #     # print(max_indegree, len(unmarked), len(marked))

    #     # mark all nodes leading to paths to m and remove them from unmarked set
    #     M[m] = partition_counter
    #     partition[m] = set()
    #     not_visited_candidates.remove(m)
    #     m_count = 1
    #     reads_to_visit = set([x for x in G_star_transposed[m] if x in unmarked])

    #     for x in reads_to_visit:
    #         assert x not in marked
    #         unmarked.remove(x)
    #         marked.add(x)
    #         partition[m].add(x)
    #         m_count += 1

    #     consensus_processed += 1
    #     # print("Center count = ", m_count)
    #     # print(len(not_visited_candidates))
    # print("Chosen minimizers:", len(M))

    # assert not unmarked  
    # assert len(M) == nr_candidates_with_hits
    # assert marked == set(G_star.keys())

    # no_alignmens_at_all = set(C.keys()).difference(set(G_star_transposed.keys()))
    # assert len(no_alignmens_at_all) + nr_candidates_with_hits == len(C)
    # print("no_alignmens_at_all:", no_alignmens_at_all)
    # for c in no_alignmens_at_all:
    #     assert c not in partition
    #     partition[c] = set()

    return G_star, partition 



