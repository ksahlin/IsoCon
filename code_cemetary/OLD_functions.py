



# def construct_exact_nearest_neighbor_graph(S, params):

#     """
#         input: a dict of strings, not necesarily unique
#         output: a directed graph implemented as a dict of dicts. Each edge has a weight assosiated to them.
#                 self edges has a weight > 1 (identical sequences) and all other edges has weight 1.
#                 Note, a node can be isolated! An isolated node will point at itself, effectively having itself as a nearest_neighbor.

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
#     all_internode_edges_in_nearest_neighbor_graph, isolated_nodes = nearest_neighbor_graph.compute_nearest_neighbor_graph(S_prime, params) # send in a list of nodes that already has converged, hence avoid unnnecessary computation
 

#     # TODO: implement already_converged to skip redundant calculations, the more important for more comverged stings we have!! 
#     # nearest_neighbor_graph, isolated_nodes = compute_nearest_neighbor_graph(S, already_converged, params) # send in a list of nodes that already has converged, hence avoid unnnecessary computation
#     for s1_acc in all_internode_edges_in_nearest_neighbor_graph:
#         s1 = S[s1_acc]
#         if s1 in G_star[s1]: # the nearest_neighbor had already identical nearest_neighbor (ed = 0)
#             continue
#         # elif len(G_star[s1]) >= 2: # already has identical homopolymer nearest_neighbors ( at least 2 for meaningful correction)
#         #     print("Homopolymer partition")
#         #     continue
#         else:
#             for s2_acc in all_internode_edges_in_nearest_neighbor_graph[s1_acc]:
#                 s2 = S[s2_acc]
#                 G_star[s1][s2] = 1 

#     for s in isolated_nodes:
#         assert s in G_star
#         if s not in G_star[s]:
#             G_star[s][s] = 1

#     return G_star, converged    


def create_multialignment_format_OLD(query_to_target_positioned_dict, start, stop):
    """
        only create multialignment format of the query sequences that cover the region [start,stop] start stop is the vector 
        coordinates where vector is of size 2*len(target) + 1
    """
    assert len(query_to_target_positioned_dict) > 0
    target_vector_length = len( list(query_to_target_positioned_dict.values())[0][0])
    assert stop < target_vector_length # vector coordinates are 0-indexed

    # get allreads alignments covering interesting segment
    # row_accessions = []
    segments = {}
    for q_acc, (q_to_t_pos, t_vector_start, t_vector_end) in query_to_target_positioned_dict.items():
        if t_vector_start <= start and t_vector_end >= stop: # read cover region
            segment = q_to_t_pos[start - t_vector_start : stop - t_vector_start +1 ]
            # row_accessions.append(q_acc)
            segments[q_acc] = segment


    # create alignment matrix of segment
    alignment_matrix = {}
    for q_acc in segments:
        alignment_matrix[q_acc] = []

    for j in range(0, stop - start + 1):
        if (start + j) % 2 == 0:  # we are between a target base pairs (need to check the longest one)
            insertions = [(segments[q_acc][j], q_acc) for q_acc in segments]
            max_insertion, q_acc_max_ins = max(insertions, key= lambda x : len(x[0]))
            max_insertion = "-" + max_insertion + "-"  # pad the max insertion

            ###### OLD WORKING CODE #########
            # for q_acc in segments:
            #     for p in range(len(max_insertion)):
            #         # all shorter insertions are left shifted -- identical indels are guaranteed to be aligned
            #         # however, no multialignment is performed among the indels
            #         if p < len(segments[q_acc][j]):
            #             alignment_matrix[q_acc].append(segments[q_acc][j][p])
            #         else:
            #             alignment_matrix[q_acc].append("-")
            #########################

            for q_acc in segments:
                # check if identical substring in biggest insertion first:
                q_ins = segments[q_acc][j]
                q_insertion_modified = ""

                if q_ins == "-":
                    # print("LOOL")
                    q_insertion_modified = "-"*len(max_insertion)

                if not q_insertion_modified:
                    pos = max_insertion.find(q_ins) 
                    q_insertion_modified = ""
                    if pos >=0:
                        if pos >= 1:
                            pass
                            # print("here perfect new!! q: {0} max: {1}, new q_ins:{2}".format(q_ins, max_insertion,  "-"*pos + max_insertion[ pos : pos + len(q_ins) ] + "-"* len(max_insertion[ pos + len(q_ins) : ])))
                        
                        q_insertion_modified = "-"*pos + max_insertion[ pos : pos + len(q_ins) ] + "-"* len(max_insertion[ pos + len(q_ins) : ])

                
                if not q_insertion_modified:
                    # else, check is smaller deletion can be aligned from left to write, e.g. say max deletion is GACG
                    # then an insertion AG may be aligned as -A-G. Take this alignment instead
                    can_be_threaded = True
                    prev_pos = -1
                    match_pos = set()
                    for q_nucl in q_ins:
                        pos = max_insertion.find(q_nucl) 
                        match_pos.add(pos)
                        if pos <= prev_pos:
                            can_be_threaded = False
                            break
                        prev_pos = pos

                    if can_be_threaded:
                        q_insertion_modified = ""
                        for p in range(len(max_insertion)):
                            if p in match_pos:
                                nucl = max_insertion[p]
                            else:
                                nucl = "-"
                            q_insertion_modified = q_insertion_modified + nucl
                        # print("NEW can be threaded: q:{0}, max: {1}, new thread: {2}".format(q_ins, max_insertion, q_insertion_modified))

                if not q_insertion_modified:
                    # otherwise just shift left
                    # print("Not solved: q:{0}, max: {1}".format(q_ins, max_insertion))
                    # check if there is at least one matching character we could align to
                    max_p = 0
                    max_matches = 0
                    for p in range(0, len(max_insertion) - len(q_ins) + 1 ):
                        nr_matches = len([1 for c1, c2 in zip(q_ins, max_insertion[p: p + len(q_ins) ] ) if c1 == c2])
                        if nr_matches > max_matches:
                            max_p = p
                            max_matches = nr_matches

                    if max_p > 0:
                        q_insertion_modified = "-"*max_p + q_ins + "-"*len(max_insertion[max_p + len(q_ins) : ])
                        # print("specially solved: q:{0} max:{1} ".format(q_insertion_modified, max_insertion) )

                if not q_insertion_modified:
                    q_insertion_modified = []
                    for p in range(len(max_insertion)):
                        # all shorter insertions are left shifted -- identical indels are guaranteed to be aligned
                        # however, no multialignment is performed among the indels
                        if p < len(q_ins):
                            q_insertion_modified.append(q_ins[p])
                        else:
                            q_insertion_modified.append("-")

                #### finally add to alignment matrix
                if len(q_insertion_modified) != len(max_insertion):
                    print(q_insertion_modified, max_insertion, q_ins)
                assert len(q_insertion_modified) == len(max_insertion)
                for p in range(len(max_insertion)):
                    alignment_matrix[q_acc].append(q_insertion_modified[p])

        else: # we are on a target base pair -- all varinats must be exactly A,C,G,T, - i.e., length 1
            for q_acc in segments:
                alignment_matrix[q_acc].append(segments[q_acc][j])
    # print(alignment_matrix)
    return alignment_matrix

def create_multialignment_format(query_to_target_positioned_dict, start, stop):
    """
        only create multialignment format of the query sequences that cover the region [start,stop] start stop is the vector 
        coordinates where vector is of size 2*len(target) + 1
    """
    assert len(query_to_target_positioned_dict) > 0
    target_vector_length = len( list(query_to_target_positioned_dict.values())[0][0])
    assert stop < target_vector_length # vector coordinates are 0-indexed

    # get allreads alignments covering interesting segment
    # row_accessions = []
    segments = {}
    for q_acc, (q_to_t_pos, t_vector_start, t_vector_end) in query_to_target_positioned_dict.items():
        if t_vector_start <= start and t_vector_end >= stop: # read cover region
            segment = q_to_t_pos[start - t_vector_start : stop - t_vector_start +1 ]
            # row_accessions.append(q_acc)
            segments[q_acc] = segment


    # create alignment matrix of segment
    alignment_matrix = {}
    for q_acc in segments:
        alignment_matrix[q_acc] = []

    for j in range(0, stop - start + 1):
        if (start + j) % 2 == 0:  # we are between a target base pairs (need to check the longest one)
            insertions = [(segments[q_acc][j], q_acc) for q_acc in segments]
            # TODO: HERE IS STOCHASTICITY! Lets just do edlib(path=True) here all the time and take leftmost best alignment if several!, that would simplify code, right?! 
            max_insertion, q_acc_max_ins = max(insertions, key= lambda x : len(x[0]))
            max_ins_len = len(max_insertion)
            all_max_ins = set(["-" + ins + "-" for (ins, acc) in insertions if len(ins) == max_ins_len])
            # if len(all_max_ins) > 1 and max_ins_len == 2:
            #     print("pos", j, all_max_ins)
            max_insertion = "-" + max_insertion + "-"  # pad the max insertion
            padded_max_ins_len = len(max_insertion)

            for q_acc in segments:
                # check if identical substring in biggest insertion first:
                q_ins = segments[q_acc][j]
                q_insertion_modified = ""

                if q_ins == "-":
                    # print("LOOL")
                    q_insertion_modified = "-"*padded_max_ins_len


                can_be_threaded = False
                if not q_insertion_modified:
                    # else, check is smaller deletion can be aligned from left to write, e.g. say max deletion is GACG
                    # then an insertion AG may be aligned as -A-G. Take this alignment instead
                    for max_ins in sorted(all_max_ins):
                        # print("max ins", max_ins)
                        q_insertion_modified = thread_to_max_ins(max_ins, q_ins)
                        if q_insertion_modified:
                            can_be_threaded = True
                            # print("mod", q_insertion_modified)
                            break


                if not q_insertion_modified:
                    q_insertion_modified = []
                    for p in range(padded_max_ins_len):
                        # all shorter insertions are left shifted -- identical indels are guaranteed to be aligned
                        # however, no multialignment is performed among the indels
                        if p < len(q_ins):
                            q_insertion_modified.append(q_ins[p])
                        else:
                            q_insertion_modified.append("-")

                #### finally add to alignment matrix
                if len(q_insertion_modified) != padded_max_ins_len:
                    print(q_insertion_modified, max_insertion, q_ins)
                assert len(q_insertion_modified) == padded_max_ins_len

                # if len(all_max_ins) > 1 and padded_max_ins_len == 4 and q_ins != "-":
                #     print(j, q_ins, "threaded:", can_be_threaded , q_insertion_modified)

                for p in range(padded_max_ins_len):
                    alignment_matrix[q_acc].append(q_insertion_modified[p])

        else: # we are on a target base pair -- all varinats must be exactly A,C,G,T, - i.e., length 1
            for q_acc in segments:
                alignment_matrix[q_acc].append(segments[q_acc][j])
    # print(alignment_matrix)
    return alignment_matrix


def adjust_probability_of_read_to_alignment_invariant(delta_t, alignment_matrix, t_acc):

    target_alignment = alignment_matrix[t_acc]
    # Get the individual invariant factors u_iv for each read and variant x_i and v in delta, v is a tuple (state, char).
    # store this 3 dimensional in dictionary u_iv = {x_acc: {pos: { (state,char) : u_iv} }} 
    u_iv = {}
    for q_acc in alignment_matrix:
        if q_acc == t_acc or q_acc in delta_t:
            continue
        u_iv[q_acc] = {}
        read_alignment = alignment_matrix[q_acc]

        for c_acc in delta_t:
            candidate_alignment = alignment_matrix[c_acc]
            for pos in delta_t[c_acc]:
                if pos not in u_iv[q_acc]:
                    u_iv[q_acc][pos] = {}
            
                state, char = delta_t[c_acc][pos]
                if state == "S": # substitutions has u_v =1 by definition
                    u_v = 1
                elif candidate_alignment[pos] != read_alignment[pos]: # read is not matching the varinat at the candidate position, this also has u_v = 1
                    u_v = 1
                else: # read matches candidate variant, now we need to check how many possible combinations an error can cause them to match due to invariant
                    u_v = get_multiplier_for_variant(state, char, pos, target_alignment, read_alignment)
                    # print("multiplier:", u_v)
                    # print(target_alignment[pos-10: pos+10],state, char)
                    # print(read_alignment[pos-10: pos+10],state, char)
                    # print(candidate_alignment[pos-10: pos+10],state, char)

                u_iv[q_acc][pos][(state, char)] = u_v

    return u_iv


def get_invariant_adjustment(delta_t, alignment_matrix, t_acc):
    invariant_factors = {}
    target_alignment = alignment_matrix[t_acc]
    stop = len(target_alignment) - 1
    for c_acc in delta_t:
        min_u_candidate_S = 1
        min_u_candidate_D = 10000
        min_u_candidate_I = 10000
        candidate_alignment = alignment_matrix[c_acc]

        for pos in delta_t[c_acc]:
            state, char = delta_t[c_acc][pos]
            u_pos = 1
            if state == "D":
                v = target_alignment[pos]
            elif state == "I":
                v = char
            else:
                # min_u_candidate = 1 # substitution by defintion has uniqueness 1 we can go to the next candidate
                continue 
            offset = 1
            upper_stop = False
            lower_stop = False
            while True:
                if pos + offset > stop:
                    upper_stop = True
                elif target_alignment[pos + offset] == candidate_alignment[pos + offset] == v:
                    u_pos += 1
                elif target_alignment[pos + offset] == candidate_alignment[pos + offset] == "-":
                    pass
                else:
                    upper_stop = True


                if pos - offset < 0:
                    lower_stop = True                    
                elif target_alignment[pos - offset] == candidate_alignment[pos - offset] == v:
                    u_pos += 1
                elif target_alignment[pos - offset] == candidate_alignment[pos - offset] == "-":
                    pass
                else:
                    lower_stop = True

                if lower_stop == upper_stop == True:
                    break


                offset += 1

            if state == "D":
                if u_pos < min_u_candidate_D:
                    min_u_candidate_D = u_pos
            elif state == "I":
                if u_pos < min_u_candidate_I:
                    min_u_candidate_I = u_pos

            # if u_pos < min_u_candidate:
            #     min_u_candidate = u_pos

        if min_u_candidate_D == 10000:
            min_u_candidate_D = 1
        if min_u_candidate_I == 10000:
            min_u_candidate_I = 1

        invariant_factors[c_acc] = (min_u_candidate_S, min_u_candidate_D, min_u_candidate_I)

    return invariant_factors

def homopolymenr_multiple_testing_correction_factor:
     #### CORRECTION FACTOR MULTIPLE TESTING W.R.T. HOMOPOLYMENR LENGHTS ##############
    # gives keyerror!!
    homopolymenr_length_numbers = functions.calculate_homopolymenr_lengths(t_seq) # dictionary { homopolymer_length : number of occurances on reference}
    print(homopolymenr_length_numbers)
    correction_factor = 1
    all_variants = {}
    for pos, (state, char) in delta_t[c_acc].items():
        print(state, char)
        u_v = candidate_indiv_invariant_factors[c_acc][pos][(state, char)]
        # print(u_v, state )
        if (u_v, state) in all_variants:
            all_variants[(u_v, state)] +=1
        else:
            all_variants[(u_v, state)] = 1

    for (u_v, state), n_v in all_variants.items():
        if state == "I":
            u_v = max(u_v - 1, 1)

        if u_v not in homopolymenr_length_numbers:
            index, u_v = min(enumerate(homopolymenr_length_numbers.keys()), key=lambda x: abs(x[1] - u_v))
            nr_hom_lengths = homopolymenr_length_numbers[u_v]
        else:
            nr_hom_lengths = homopolymenr_length_numbers[u_v]


        if u_v > 1:
            correction_factor *= functions.choose(nr_hom_lengths, n_v)

        elif state == "I":
            correction_factor *= functions.choose(4*(nr_hom_lengths + 1), n_v)
        elif state == "S":
            correction_factor *= functions.choose(3*nr_hom_lengths, n_v)

        elif state == "D":
            correction_factor *= functions.choose(nr_hom_lengths, n_v)

        else:
            print("BUG", state)
            sys.exit()
    print(correction_factor)
    #####################################################################################

def check_if_consensus(c_acc, C, X, partition_of_X):

    partition_dict = {c_acc : {c_acc : (C[c_acc], C[c_acc])}}
    for x_acc in partition_of_X[c_acc]:
        partition_dict[c_acc][x_acc] = (C[c_acc], X[x_acc])

    exact_edit_distances = edlib_align_sequences_keeping_accession(partition_dict, single_core = True)    
    exact_alignments = sw_align_sequences_keeping_accession(exact_edit_distances, single_core = True)
    partition_alignments = {} 

    for c_acc in exact_alignments:
        partition_alignments[c_acc] = {}
        for x_acc in exact_alignments[c_acc]:
            aln_c, aln_x, (matches, mismatches, indels) = exact_alignments[c_acc][x_acc]
            edit_dist = mismatches + indels
            partition_alignments[c_acc][x_acc] = (edit_dist, aln_c, aln_x, 1)

    alignment_matrix_to_c, PFM_to_c = create_position_probability_matrix(C[c_acc], partition_alignments[c_acc])
    c_alignment = alignment_matrix_to_c[c_acc]
    is_consensus = True
    not_cons_positions = []
    for j in range(len(PFM_to_c)):
        c_v =  c_alignment[j]
        candidate_count = PFM_to_c[j][c_v]
        for v in PFM_to_c[j]:
            if v != c_v and candidate_count <= PFM_to_c[j][v]: # needs to have at least one more in support than the second best as we have added c itself to the multialignment
                # print("not consensus at:", j, PFM_to_c[j])
                not_cons_positions.append((j, c_v, PFM_to_c[j]))
                is_consensus = False

    if not is_consensus:
        print("Were not consensus at:", not_cons_positions)
    else:
        print("Were consensus")

    return is_consensus

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
#             # Do individual nearest_neighbor component graphs of the homopolymenr equivalence classes here!
#             s2 = S[acc2]
#             G_star[s1][s2] = 1
#             G_star[s2][s1] = 1
#             homopol_extra_added += 2

#     print("EDGES FROM HOMOPOLYMER IDENTICAL:", homopol_extra_added)
#     unique_strings = {transform(seq) : acc for acc, seq in S.items()}
#     S_prime_transformed = {acc : seq for seq, acc in unique_strings.items()}
#     # Send homopolymer components to this function!
#     # Keep in mind. Isolated nodes are not neccesarily isolated!
#     all_internode_edges_in_nearest_neighbor_graph, isolated_nodes = nearest_neighbor_graph.compute_nearest_neighbor_graph(S_prime_transformed, params) # send in a list of nodes that already has converged, hence avoid unnnecessary computation

#     #############################################################
#     #############################################################



# def highest_reachable_with_edge_degrees_NEW(S, params):
#     G_star, converged = graphs.construct_exact_nearest_neighbor_graph_improved(S, params)
#     unique_start_strings = set(G_star.nodes())
#     from operator import itemgetter
#     print("len G_star:", len(G_star))
#     partition_sizes = []
#     nr_consensus = 0
#     G_transpose = nx.reverse(G_star)
#     print("len G_star_transposed (nearest_neighbors):", len(G_transpose))
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
#             # cut at all nodes with more than 2 nearest_neighbors
#             ####################################################
#             # nearest_neighbor_w_more_than_2_support = set()
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
#     nearest_neighbor_lenghts = [len(m) for m in sorted(partition)]
#     print(sorted(nearest_neighbor_lenghts))
#     print("NR CONSENSUS:", nr_consensus)
#     print("NR nearest_neighbors:", len(M), len(partition))
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

#     G_star, converged = graphs.construct_exact_nearest_neighbor_graph(S, params)
#     # G_star, converged = graphs.construct_nearest_neighbor_graph_approximate(S, params, edge_creating_min_treshold = edge_creating_min_treshold, edge_creating_max_treshold = edge_creating_max_treshold)

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
#     print("Chosen nearest_neighbors:", len(M))


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

