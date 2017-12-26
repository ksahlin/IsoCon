
from modules import edlib_alignment_module
from modules import SW_alignment_module

def find_best_matches(approximate_matches, params, edge_creating_min_treshold = -1, edge_creating_max_treshold = 2**30):
    """
        input: approximate_matches is a dictionary with a string as key and a list of strings as value
        output: dictionary with a string as key and a dictionary as value. 
                The inner dict contains a tuple (edit_distance, s1_alignment, s2_alignment) as value.
                Each string in the inner dict has the same (lowest) edit distance to the key 
    """

    exact_edit_distances = edlib_alignment_module.edlib_align_sequences(approximate_matches, nr_cores = params.nr_cores)
    best_exact_edit_distances = {}
    tot_ed = 0
    cntrr = 0
    for s1 in exact_edit_distances:
        for s2 in exact_edit_distances[s1]:
            edit_distance = exact_edit_distances[s1][s2]
            tot_ed += edit_distance
            cntrr += 1
            if edit_distance < edge_creating_max_treshold:
                if s1 in best_exact_edit_distances:
                    best_exact_edit_distances[s1][s2] = edit_distance
                else:
                    best_exact_edit_distances[s1] = {}
                    best_exact_edit_distances[s1][s2] = edit_distance

                if s2 in best_exact_edit_distances:
                    best_exact_edit_distances[s2][s1] = edit_distance    
                else:
                    best_exact_edit_distances[s2] = {}
                    best_exact_edit_distances[s2][s1] = edit_distance

    print("TOTAL EDIT DISTANCE:", tot_ed)
    print("TOTAL Visits:", cntrr)



    ## filter best exact edit distances here
    for s1 in list(best_exact_edit_distances.keys()):
        s1_nearest_neighbor = min(best_exact_edit_distances[s1], key = lambda x: best_exact_edit_distances[s1][x])
        # print(s1_nearest_neighbor)
        min_edit_distance = best_exact_edit_distances[s1][s1_nearest_neighbor]
        for s2 in list(best_exact_edit_distances[s1].keys()):
            ed =  best_exact_edit_distances[s1][s2]
            # print(ed, min_edit_distance)
            if ed > min_edit_distance:
                if ed > edge_creating_min_treshold:
                    del best_exact_edit_distances[s1][s2]

    cntrr = 0
    filtered_tot_ed = 0
    for s1 in best_exact_edit_distances:
        for s2 in best_exact_edit_distances[s1]:
            filtered_tot_ed += best_exact_edit_distances[s1][s2]
            cntrr += 1
    print("TOTAL FILTERED EDIT DISTANCE:", filtered_tot_ed)
    print("TOTAL EDGES:", cntrr)
    print("EDIT DISTANCE PER EDGE:", filtered_tot_ed/float(cntrr))

    exact_alignments = SW_alignment_module.sw_align_sequences(best_exact_edit_distances, nr_cores = params.nr_cores)

    # process the exact matches here
    best_exact_matches = {}
    tot_sw_ed = 0
    for s1 in exact_alignments:
        for s2 in exact_alignments[s1]:
            s1_alignment, s2_alignment, (matches, mismatches, indels) = exact_alignments[s1][s2]
            edit_distance = mismatches + indels
            tot_sw_ed += edit_distance
            if edit_distance < edge_creating_max_treshold:
                if s1 in best_exact_matches:
                    best_exact_matches[s1][s2] = (edit_distance, s1_alignment, s2_alignment)
                else:
                    best_exact_matches[s1] = {}
                    best_exact_matches[s1][s2] = (edit_distance, s1_alignment, s2_alignment)
                    # print( "lol", edit_distance, edge_creating_max_treshold)

                if s2 in best_exact_matches:
                    best_exact_matches[s2][s1] = (edit_distance, s2_alignment, s1_alignment)    
                else:
                    best_exact_matches[s2] = {}
                    best_exact_matches[s2][s1] = (edit_distance, s2_alignment, s1_alignment)


    print("TOTAL SW EDIT DISTANCE:", tot_sw_ed)
    filtered_sw_tot_ed = tot_sw_ed

    ## filter best exact matches here
    for s1 in list(best_exact_matches.keys()):
        s1_nearest_neighbor = min(best_exact_matches[s1], key = lambda x: best_exact_matches[s1][x][0])
        min_edit_distance = best_exact_matches[s1][s1_nearest_neighbor][0]
        for s2 in list(best_exact_matches[s1].keys()):
            ed =  best_exact_matches[s1][s2][0]
            if ed > min_edit_distance:
                if ed > edge_creating_min_treshold:
                    del best_exact_matches[s1][s2]
                    filtered_sw_tot_ed -= ed

    print("TOTAL FILTERED SW EDIT DISTANCE:", filtered_sw_tot_ed)



    # tot_edit = 0
    # tot_edit_edlib = 0
    # for read_acc, t_dict in best_exact_matches.items():
    #     if len(best_exact_matches[read_acc]) > 1:
    #         print(len(best_exact_matches[read_acc]), "BEST! ed:", [ed for ed, x,y in best_exact_matches[read_acc].values()] )
    #     for t_acc in best_exact_matches[read_acc]:
    #         tot_edit_edlib += edlib_alignment(read_acc, t_acc, 1,1, ends_discrepancy_threshold = 25 , x_acc = "", y_acc = "" )[2][2][2]
    #         tot_edit += best_exact_matches[read_acc][t_acc][0]
    #         break
    # print("TOT EDIT:", tot_edit)
    # print("TOT EDIT EDLIB:", tot_edit_edlib)

    # sys.exit()

    return best_exact_matches

def find_best_matches_2set(highest_paf_scores, X, C, params):
    """
        input: approximate_matches is a dictionary with a string as key and a list of strings as value
        output: dictionary with a string as key and a dictionary as value. 
                the outer dict contains the reads .
                This inner dict dict contains a tuple (edit_distance, s1_alignment, s2_alignment) as value.
                Each string in the inner dict has the same (lowest) edit distance to the key 
    """

    approximate_matches = {}
    for read_acc, best_hits in  highest_paf_scores.items():
        approximate_matches[read_acc] = {}
        best_score = 0
        for score, t_acc in best_hits:
            approximate_matches[read_acc][t_acc] = (X[read_acc], C[t_acc])


    exact_edit_distances = edlib_alignment_module.edlib_align_sequences_keeping_accession(approximate_matches, nr_cores = params.nr_cores)
    best_exact_edit_distances = {}
    for s1_acc in exact_edit_distances:
        for s2_acc in exact_edit_distances[s1_acc]:
            (s1,s2,edit_distance) = exact_edit_distances[s1_acc][s2_acc]
            # if edit_distance < edge_creating_max_treshold:
            if s1_acc in best_exact_edit_distances:
                best_exact_edit_distances[s1_acc][s2_acc] = (s1,s2,edit_distance)
            else:
                best_exact_edit_distances[s1_acc] = {}
                best_exact_edit_distances[s1_acc][s2_acc] = (s1,s2,edit_distance)


    ## filter best exact edit distances here
    for s1_acc in list(best_exact_edit_distances.keys()):
        s1_nearest_neighbor = min(best_exact_edit_distances[s1_acc], key = lambda x: best_exact_edit_distances[s1_acc][x][2])
        min_edit_distance = best_exact_edit_distances[s1_acc][s1_nearest_neighbor][2]
        # print("ed:", min_edit_distance,  s1_acc)
        for s2_acc in list(best_exact_edit_distances[s1_acc].keys()):
            ed =  best_exact_edit_distances[s1_acc][s2_acc][2]
            if ed > min_edit_distance:
                del best_exact_edit_distances[s1_acc][s2_acc]


    exact_alignments = SW_alignment_module.sw_align_sequences_keeping_accession(best_exact_edit_distances, nr_cores = params.nr_cores)

    # process the exact matches here
    best_exact_matches = {}
    for x_acc in exact_alignments:
        for c_acc in exact_alignments[x_acc]:
            x_alignment, c_alignment, (matches, mismatches, indels) = exact_alignments[x_acc][c_acc]
            edit_distance = mismatches + indels

            if x_acc in best_exact_matches:
                x_nearest_neighbor = list(best_exact_matches[x_acc].keys())[0]

                if edit_distance < best_exact_matches[x_acc][x_nearest_neighbor][0]:
                    best_exact_matches[x_acc] = {}
                    best_exact_matches[x_acc][c_acc] = (edit_distance, x_alignment, c_alignment)
                elif edit_distance == best_exact_matches[x_acc][x_nearest_neighbor][0]:
                    best_exact_matches[x_acc][c_acc] = (edit_distance, x_alignment, c_alignment)
                else:
                    pass
            else:
                best_exact_matches[x_acc] = {}
                best_exact_matches[x_acc][c_acc] = (edit_distance, x_alignment, c_alignment)




    # tot_edit = 0
    # tot_edit_edlib = 0
    # for read_acc, t_dict in best_exact_matches.items():
    #     if len(best_exact_matches[read_acc]) >1:
    #         print(len(best_exact_matches[read_acc]), "BEST! ed:", [ed for ed, x,y in best_exact_matches[read_acc].values()] )
    #     for t_acc in best_exact_matches[read_acc]:
    #         tot_edit_edlib += edlib_alignment(read_acc, t_acc, 1,1, ends_discrepancy_threshold = 25 , x_acc = "", y_acc = "" )[2][2][2]
    #         tot_edit += best_exact_matches[read_acc][t_acc][0]
    #         break
    # print("TOT EDIT:", tot_edit)
    # print("TOT EDIT EDLIB:", tot_edit_edlib)

    # sys.exit()


    return best_exact_matches
