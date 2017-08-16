



# def construct_exact_minimizer_graph(S, params):

#     """
#         input: a dict of strings, not necesarily unique
#         output: a directed graph implemented as a dict of dicts. Each edge has a weight assosiated to them.
#                 self edges has a weight > 1 (identical sequences) and all other edges has weight 1.
#                 Note, a node can be isolated! An isolated node will point at itself, effectively having itself as a minimizer.

#     """



#     G_star = {}
#     # adding self edges to strings that has converged
#     for acc, s in S.items():
#         if s not in G_star:
#             G_star[s] = {}
#         else:
#             if s in G_star[s]:
#                 G_star[s][s] += 1  
#             else:
#                 G_star[s][s] = 2

#     # check if converged, that is, if all nodes has self edges here, there will be no other edges added.
#     converged = False
#     not_in_clusters = set()
#     already_converged = set()
#     for s, nbr_dict in G_star.items():
#         if len(nbr_dict) == 0:
#             not_in_clusters.add(s)
#         else:
#             already_converged.add(s)

#     if len(not_in_clusters) == 0:
#         converged = True
#         return G_star, converged


#     unique_strings = {seq : acc for acc, seq in S.items()}
#     S_prime = {acc : seq for seq, acc in unique_strings.items()}
#     all_internode_edges_in_minimizer_graph, isolated_nodes = minimizer_graph.compute_minimizer_graph(S_prime, params) # send in a list of nodes that already has converged, hence avoid unnnecessary computation
 

#     # TODO: implement already_converged to skip redundant calculations, the more important for more comverged stings we have!! 
#     # minimizer_graph, isolated_nodes = compute_minimizer_graph(S, already_converged, params) # send in a list of nodes that already has converged, hence avoid unnnecessary computation
#     for s1_acc in all_internode_edges_in_minimizer_graph:
#         s1 = S[s1_acc]
#         if s1 in G_star[s1]: # the minimizer had already identical minimizer (ed = 0)
#             continue
#         # elif len(G_star[s1]) >= 2: # already has identical homopolymer minimizers ( at least 2 for meaningful correction)
#         #     print("Homopolymer partition")
#         #     continue
#         else:
#             for s2_acc in all_internode_edges_in_minimizer_graph[s1_acc]:
#                 s2 = S[s2_acc]
#                 G_star[s1][s2] = 1 

#     for s in isolated_nodes:
#         assert s in G_star
#         if s not in G_star[s]:
#             G_star[s][s] = 1

#     return G_star, converged    


# def get_homopolymer_invariants(candidate_transcripts):
#     seq_to_acc = { seq : acc for (acc, seq) in  candidate_transcripts.items() }
#     print("Unique before compression: ", len(seq_to_acc) )

#     candidate_transcripts_transformed = {}
#     clusters = defaultdict(list)
#     for acc in candidate_transcripts:
#         seq_transformed = transform(candidate_transcripts[acc])
#         candidate_transcripts_transformed[acc] = seq_transformed
#         clusters[seq_transformed].append(acc)

#     seq_to_acc_transformed = { seq : acc for (acc, seq) in candidate_transcripts_transformed.items()}
#     print("Unique after compression: ", len(seq_to_acc_transformed) )

#     edges = {}
#     for seq in clusters:
#         if len(clusters[seq]) > 1:
#             # print(clusters[seq])
#             for acc in clusters[seq]:
#                 edges[acc] = {}
#             for acc1, acc2 in combinations(clusters[seq], 2):
#                 edges[acc1][acc2] = 1
#                 edges[acc2][acc1] = 1

#     return edges


# #############################################################
# #############################################################
# if homopolymer_compression:
#     # All converged
#     for acc, s in S.items():
#         if s not in G_star:
#             G_star[s] = {}
#         else:
#             if s in G_star[s]:
#                 G_star[s][s] += 1  
#                 # print(acc)
#             else:
#                 G_star[s][s] = 2
#                 # print(acc)

#     converged = False
#     print("ENTERING HOMOPOLYMER COMPRESSION MODE")
#     # create homopolymer equivalence class edges
#     G_homopolymer_star = {}
#     weight = {}
#     for acc, s in S.items():
#         if s not in G_star:
#             G_star[s] = {}
#             weight[s] = 1
#         else:
#             weight[s] += 1

#     homopolymer_edges = get_homopolymer_invariants(S)
#     homopol_extra_added = 0
#     for acc1 in homopolymer_edges:
#         s1 = S[acc1]
#         for acc2 in homopolymer_edges[acc1]:
#             # Do individual minimizer component graphs of the homopolymenr equivalence classes here!
#             s2 = S[acc2]
#             G_star[s1][s2] = 1
#             G_star[s2][s1] = 1
#             homopol_extra_added += 2

#     print("EDGES FROM HOMOPOLYMER IDENTICAL:", homopol_extra_added)
#     unique_strings = {transform(seq) : acc for acc, seq in S.items()}
#     S_prime_transformed = {acc : seq for seq, acc in unique_strings.items()}
#     # Send homopolymer components to this function!
#     # Keep in mind. Isolated nodes are not neccesarily isolated!
#     all_internode_edges_in_minimizer_graph, isolated_nodes = minimizer_graph.compute_minimizer_graph(S_prime_transformed, params) # send in a list of nodes that already has converged, hence avoid unnnecessary computation

#     #############################################################
#     #############################################################



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

