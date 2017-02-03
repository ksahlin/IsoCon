
from modules import graphs
from modules import functions

def partition_strings_paths(S, node_weights = {}, edge_creating_min_treshold = -1, edge_creating_max_treshold = 2**30):
    G_star, alignment_graph, converged = graphs.construct_minimizer_graph(S, edge_creating_min_treshold = edge_creating_min_treshold, edge_creating_max_treshold = edge_creating_max_treshold)
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
                # M[s] = partition_counter
                # marked.add(s)
                # partition[s] = set()
                # partition_counter += 1
                # edit_dist, a1, a_min =  alignment_graph[s][s]
                # indegree = G_star[s][s]
                # print("DETECTED!!")
    # isolated = set(partition.keys())
    print("nr isolated nodes:", isolated)

    G_star_transposed = functions.transpose(G_star)
    # check so that there are no "leaves" in G^* (a sequence that has no minimizer but is a minimizer to some other sequence, this should not happen)
    # for n in G_star_transposed:
    #     assert n not in isolated

    # do_while as long as there there are unmarked nodes
    while len(marked) < V_G:

        # find node with biggest indegree
        max_indegree = -1
        max_node_weight = -1
        for s in unmarked:
            if node_weights:
                indegree = node_weights[s]  # choose the node with the max support of reads, this is found in the node weight if specified
                if indegree >= max_indegree:
                    m, max_indegree = s, indegree        

                # for in_nbr, edge_indegree in G_star_transposed[s].items():
                #     if in_nbr not in marked:
                #         indegree += node_weights[in_nbr] # the number of reads aligned to a given candidate      

                # if indegree >= max_indegree and node_weights[s] > max_node_weight:
                #     m, max_indegree = s, indegree        
                #     max_node_weight = node_weights[s]

            else:
                indegree = sum([indegree for in_nbr, indegree in G_star_transposed[s].items() if in_nbr not in marked ])

                if indegree > max_indegree:
                    m, max_indegree = s, indegree
        # print(max_indegree, len(unmarked), len(marked))

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

        print("Center count = ", m_count)
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
    print(total_strings_in_partition)
    # print(len(partition_sequences))
    # print(len(unique_start_strings))
    assert unique_start_strings == partition_sequences

    return G_star, partition, M, converged




def partition_strings(S, node_weights = {}):
    G_star, alignment_graph, converged = graphs.construct_minimizer_graph(S)
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
                # print("DETECTED!!")
    isolated = set(partition.keys())
    print("nr isolated nodes:", len(isolated))

    G_star_transposed = functions.transpose(G_star)
    # check so that there are no "leaves" in G^* (a sequence that has no minimizer but is a minimizer to some other sequence, this should not happen)
    # for n in G_star_transposed:
    #     assert n not in isolated

    # do_while as long as there is a node with positive indegree of unmarked nbrs
    while True:
        # find node with biggest indegree
        max_indegree = -1
        for s in V_not_in_M:
            if node_weights:
                indegree = node_weights[s]  # choose the node with the max support of reads, this is found in the node weight if specified
            else:
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




def partition_strings_2set_paths(X, C, X_file, C_file):
    """

    """


    G_star, alignment_graph = graphs.construct_2set_minimizer_bipartite_graph(X, C, X_file, C_file)
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
        print("Center count = ", m_count)
        # print(len(not_visited_candidates))
    print("Chosen minimizers:", len(M))

    assert not unmarked  
    assert len(M) == nr_candidates_with_hits
    assert marked == set(G_star.keys())


    return G_star, partition, alignment_graph



