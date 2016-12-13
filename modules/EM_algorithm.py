import sys
import operator
from scipy.stats import binom
from collections import defaultdict
import scipy
import math
import random
import signal
from multiprocessing import Pool
import multiprocessing as mp

from hitemmodules import distributions
from hitemmodules.io import write_output
from hitemmodules import align
from hitemmodules.io import calc_memory
from hitemmodules import initialize_data
from hitemmodules import misc_functions



def sample_y_by_probability(tau):
    tau_fixed = {}
    for y in tau:
        tau_fixed[y] = []
        for p in tau[y]:
            pass
    return tau_fixed

def arrange_x_alignments_to_y_helper(args):
    return arrange_x_alignments_to_y(*args)

def update_y_j_helper(args):
    return update_y_j(*args)

# to work for multiprocessing: A module-level function is a function which is defined at module level, 
# that means it is not an instance method of a class, it's not nested within another function, 
# and it is a "real" function with a name, not a lambda function.
def intdefdict():
    return defaultdict(int)
def floatdefdict():
    return defaultdict(float)
def tupledefdict():
    return defaultdict(tuple)

def arrange_x_alignments_to_y(y_j, y_j_sequence, alignment_results_y_j):

    # print(y_j_sequence)
    x_to_y_j_alignments = {} # x_i : [] list of length 2l_j +1 to hold insertions
    for x_i in alignment_results_y_j:
        y_j_aligned, x_i_aligned, stats = alignment_results_y_j[x_i]
        # print("processing alignment for: ",x_i)
        # print(y_j_aligned, equence) )
        # print(x_i_aligned)
        x_to_y_j_alignments[ x_i ] = [0]*(2*len(y_j_sequence) +1)
        y_position = 0
        temp_ins = ""
        # iterating over alignment positions
        for p in range(len(y_j_aligned)):
            if y_j_aligned[p] == "-":
                temp_ins += x_i_aligned[p]
                # x_to_y_j_alignments[x_i][2*y_position].append(x_i_aligned[p])
                # pass
            else:
                # print("temp ins", temp_ins)
                if not temp_ins:
                    x_to_y_j_alignments[x_i][2*y_position] = "-"
                else:
                    x_to_y_j_alignments[x_i][2*y_position] = temp_ins
                    temp_ins = ""

                x_to_y_j_alignments[x_i][2*y_position+1] = x_i_aligned[p]

                y_position += 1


        if not temp_ins:
            x_to_y_j_alignments[x_i][2*y_position] = "-"
        else:
            x_to_y_j_alignments[x_i][2*y_position] = temp_ins
            temp_ins = ""
        try:
            assert y_position == len(y_j_sequence)
        except AssertionError:
            print("y_position != len(y_j_sequence) !!!")
            print(y_j_aligned)
            print(x_i_aligned)
            print(y_position)
            print(len(y_j_sequence))
            print(temp_ins)
            print(y_j, y_j_sequence)
            print(x_i)
            sys.exit()

    mem_aln = calc_memory.total_size(x_to_y_j_alignments)
    # print("SIZE x_to_y_j_alignments", mem_aln)
    return y_j, x_to_y_j_alignments, mem_aln


def calculate_alpha_ij(alignment_results_dict_transposed, epsilon, epsilon_y):
    alpha_ij = {} #defaultdict(floatdefdict)
    for x_i in alignment_results_dict_transposed:
        for y_j in alignment_results_dict_transposed[x_i]:
            y_j_aligned, x_i_aligned, stats = alignment_results_dict_transposed[x_i][y_j]
            if x_i not in alpha_ij:
                alpha_ij[x_i] = {}
            matches, mismatches, indels = stats 
            n = matches + mismatches + indels
            k = mismatches + indels
            p = epsilon[x_i] + epsilon_y[y_j] - epsilon[x_i] * epsilon_y[y_j]
            # print(n, k, p, epsilon[x_i], epsilon_y[y_j])
            # print(type(p))
            alpha_ij[x_i][y_j] = binom.pmf(k, n, p)
            # if alpha_ij[x_i][y_j] == 1:
            #     print("alpha is 1 for", x_i)
    return alpha_ij


def calculate_edit_distance_ij(alignment_results_dict_transposed):
    edit_distance_ij = {} #defaultdict(floatdefdict)
    for x_i in alignment_results_dict_transposed:
        for y_j in alignment_results_dict_transposed[x_i]:
            y_j_aligned, x_i_aligned, stats = alignment_results_dict_transposed[x_i][y_j]
            if x_i not in edit_distance_ij:
                edit_distance_ij[x_i] = {}
            matches, mismatches, indels = stats 
            alignment_length = matches + mismatches + indels
            edit_distance = mismatches + indels
            edit_distance_ij[x_i][y_j] = (edit_distance, alignment_length)
            # print((edit_distance, alignment_length))

    return edit_distance_ij

def update_epsilon_min_distance_to_y(x, y, edit_distance_ij, epsilon_y, errors_y, params): 
    epsilon = {}
    errors_x = {}
    for x_i in edit_distance_ij: 
        # tmp_total_weight_ed_i = 0
        # tmp_total_weight_alignment_length_i = 0
        min_distance = len(x[x_i])
        min_y = None
        for y_j in edit_distance_ij[x_i]:
            if edit_distance_ij[x_i][y_j][0] < min_distance:
                min_distance = edit_distance_ij[x_i][y_j][0]
                min_y = y_j

        if min_distance == len(x[x_i]):
            write_output.logger("min_distance is {0} for {1}.\nIt has {2} neighbors.\nLast aligment length: {3}".format(min_distance, x_i, len(alpha_ij[x_i]), edit_distance_ij[x_i][y_j][1], alpha_ij[x_i][y_j]), params.logfile, timestamp=False)

        # epsilon[x_i] = float( min_distance - epsilon_y[min_y]*len(y[min_y]) ) /  edit_distance_ij[x_i][y_j][1]
        # epsilon[x_i] = float( min_distance - epsilon_y[min_y]*len(y[min_y]) ) / len(x[x_i])
        # epsilon[x_i] = float( min_distance ) / len(x[x_i])
        epsilon[x_i] = abs(float( min_distance - errors_y[min_y])) / edit_distance_ij[x_i][min_y][1]
        errors_x[x_i] = abs(min_distance - errors_y[min_y]) # abs(min_distance - errors_y[min_y])
        if min_distance < errors_y[min_y]:
            print("OMG X: less errors than in y -- over correction or is it error in estimate of errors in x_i_original?",min_distance, errors_y[min_y])
            print(x_i, "nr nbrs:", len(edit_distance_ij[x_i]))
        if epsilon[x_i] == 0:
            print("EPSILON estimated to 0 for ", x_i, "nr nbrs:", len(edit_distance_ij[x_i]))

    return epsilon, errors_x

# def update_epsilon(alpha_ij, edit_distance_ij, epsilon_y, epsilon_x_old, params): 
#     epsilon = {}
#     for x_i in edit_distance_ij: 
#         tmp_total_weight_ed_i = 0
#         tmp_total_weight_alignment_length_i = 0
#         for y_j in edit_distance_ij[x_i]:
#             tmp_total_weight_ed_i += alpha_ij[x_i][y_j] * edit_distance_ij[x_i][y_j][0]
#             # print(alpha_ij[x_i][y_j], edit_distance_ij[x_i][y_j][0], epsilon_y[y_j], edit_distance_ij[x_i][y_j][1])
#             # tmp_total_weight_ed_i += abs(alpha_ij[x_i][y_j] * ( edit_distance_ij[x_i][y_j][0] - epsilon_y[y_j]*edit_distance_ij[x_i][y_j][1] ))
#             tmp_total_weight_alignment_length_i +=  alpha_ij[x_i][y_j] * edit_distance_ij[x_i][y_j][1]

#         tmp_total_weight_ed_i = abs(tmp_total_weight_ed_i)
#         # assert not math.isnan(tmp_total_weight_ed_i)
#         # assert not math.isnan(tmp_total_weight_alignment_length_i)
#         if tmp_total_weight_alignment_length_i == 0:
#             write_output.logger("tmp_total_weight_alignment_length_i is 0 for {0}.\nIt has {1} neighbors.\nLast aligment length: {2}\nLast weight: {3}".format(x_i, len(alpha_ij[x_i]), edit_distance_ij[x_i][y_j][1], alpha_ij[x_i][y_j]), params.logfile, timestamp=False)
#             write_output.logger("We therefore let epsilon_x remain estimated to:{0} for {1}.".format(epsilon_x_old[x_i], x_i), params.logfile, timestamp=False)
#             epsilon[x_i] = epsilon_x_old[x_i]
#         else:
#             epsilon[x_i] = tmp_total_weight_ed_i / tmp_total_weight_alignment_length_i
#         if epsilon[x_i] == 0:
#             print("EPSILON estimated to 0 for ", x_i, "nr nbrs:", len(alpha_ij[x_i]))
#             # for y_j in edit_distance_ij[x_i]:
#             #     print("alpha: ", alpha_ij[x_i][y_j], "edit distance:", edit_distance_ij[x_i][y_j])
#         # print(epsilon[x_i])
#     return epsilon

# def update_epsilon_y(alpha_ji, edit_distance_ji, epsilon_x, epsilon_y_old, params): 
#     epsilon_y = {}
#     for y_j in edit_distance_ji: 
#         tmp_nominator_j = 0
#         tmp_demoninator_j = 0
#         # print()
#         for x_i in edit_distance_ji[y_j]:
#             tmp_nominator_j += abs(alpha_ji[y_j][x_i] * ( edit_distance_ji[y_j][x_i][0]  - epsilon_x[x_i]*edit_distance_ji[y_j][x_i][1] ))
#             tmp_demoninator_j +=  alpha_ji[y_j][x_i] * edit_distance_ji[y_j][x_i][1]

#         if tmp_demoninator_j == 0:
#             write_output.logger("tmp_denominator_j is 0 for y_j:{0}.\nIt has {1} neighbors.".format(y_j, len(alpha_ji[y_j])), params.logfile, timestamp=False)
#             write_output.logger("We therefore let epsilon_y remain estimated to:{0} for {1}.".format(epsilon_y_old[y_j], y_j), params.logfile, timestamp=False)
#             epsilon_y[y_j] = epsilon_y_old[y_j]
#         else:
#             epsilon_y[y_j] = tmp_nominator_j / tmp_demoninator_j
#         if epsilon_y[y_j] == 0:
#             print("EPSILON_Y estimated to 0 for ", y_j, "nr nbrs:", len(alpha_ji[y_j]))
#     return epsilon_y

def update_epsilon_y_y_to_y_alignments(alignment_results_dict, alignment_results_dict_transposed, y, x_to_y, y_to_x, old_epsilon_y, old_errors_y, params):
    # 1. Create a y_to_y_alignment_dict based on y_to_x_alignments, send in an empty alignment_results_dict_y_to_y = {}
    # 2. Get alignments from y to y by sending y_to_y_alignment_dict to function sw_align_sequences(y_to_y_alignment_dict, alignment_results_dict_y_to_y, y, y, params)
    # 3. Send alignment_results_dict_y_to_y to function initialize_epsilon with parameters: initialize_epsilon(y, alignment_results_dict_y_to_y, y_to_y, params)
    # output: this will give us best error estimates for y

    # 1.
    y_to_y_to_map = {}
    alignment_results_dict_y_to_y = {}
    y_j_xs_with_no_y_left = set()

    for y_j in alignment_results_dict:
        y_to_y_to_map[y_j] = {}
        alignment_results_dict_y_to_y[y_j] = {}
        for x_i in alignment_results_dict[y_j]:
            y_i = x_to_y[x_i]

            # x_i can still be here although it's corresponding starting point y_j has been removed.
            # not the other way around though, as it is implemented now. Eventually change this back. 
            if y_i not in y:
                continue 
            if y[y_j] == y[y_i]:
                del y_to_y_to_map[y_j]
                alignment_results_dict_y_to_y[y_j] = {}
                alignment_results_dict_y_to_y[y_j][y_i] = (y[y_j], y[y_i], (len(y[y_j]), 0, 0))
                print("SAVE TIME break")
                break
            else:
                # dummy initialixe data entries
                y_to_y_to_map[y_j][y_i] = 0
                alignment_results_dict_y_to_y[y_j][y_i] = 0

        if not alignment_results_dict_y_to_y[y_j]:
            y_j_xs_with_no_y_left.add(y_j)
    
    for y_j in y_j_xs_with_no_y_left:
        print("OMG, NO Y LEFT OUT OF {0} READS ALIGNED".format(len(alignment_results_dict[y_j])))
        del alignment_results_dict_y_to_y[y_j]

    # 2.
    ##  Mapping all y's that had at least one read x_i with a y left.

    print("DIM:", len(y_to_y_to_map))
    align.sw_align_sequences(y_to_y_to_map, alignment_results_dict_y_to_y, y, y, params)

    ## Some of these y_j's might not have gotten any alignments due o diverged sequences. We
    ## can then identify the best mapping x_i to y_j, named x_i_best. We take all the y_j's that x_i_best
    ## maps to and align them instead.

    y_to_y_to_map = {}
    for y_j in y_j_xs_with_no_y_left:
        alignment_results_dict_y_to_y[y_j] = {}
        (best_x_i, (x_aln,y_aln, (matches, mismatches, indels))) = sorted(alignment_results_dict[y_j].items(), key=lambda z: z[1][2][1] + z[1][2][2])[0]
        print("best x_i:", best_x_i )
        print("for y_j", y_to_x[y_j])
        best_y_matches = sorted(alignment_results_dict_transposed[best_x_i].items(), key=lambda z: z[1][2][1] + z[1][2][2])[:100]
        print(len(best_y_matches))
        for y_k, data in best_y_matches:
            if y_j == y_k:
                continue
            if y[y_j] == y[y_k]:
                alignment_results_dict_y_to_y[y_j][y_k] = (y[y_j], y[y_k], (len(y[y_j]), 0, 0))
                print("SAVE TIME break TWO 2!!")
                del y_to_y_to_map[y_j]
                break
            else:
                # dummy initialixe data entries
                y_to_y_to_map[y_j][y_k] = 0
                alignment_results_dict_y_to_y[y_j][y_k] = 0
    print("DIM2:", len(y_to_y_to_map))

    align.sw_align_sequences(y_to_y_to_map, alignment_results_dict_y_to_y, y, y, params)

    
    ## finally, there might still be some y_j that didn't have any valid alignments after this step as well, 
    ## they will have an empty alignment_results_dict_y_to_y[y_j] here
    y_j_with_no_good_alignment = set()
    for y_j in list(alignment_results_dict_y_to_y):
        if not alignment_results_dict_y_to_y[y_j]:
            y_j_with_no_good_alignment.add(y_j)
            print("No good alignments after 2nd step:", y_j)
            print("Letting number of errors be the old estimate:", old_errors_y[y_j])
    for y_j in y_j_with_no_good_alignment:
        print("NO GOOD ALIGNMENT FOR".format(y_j))
        del alignment_results_dict_y_to_y[y_j]


    # 3.
    alignment_results_dict_y_to_y_transposed = misc_functions.transpose(alignment_results_dict_y_to_y)
    y_to_y = dict([(y_j,y_j) for y_j in y])
    epsilon_y_backup, errors_y_backup = initialize_data.initialize_epsilon(y, alignment_results_dict_y_to_y, y_to_y, params)
    # sending containers in reversed order here because all keys in alignment_results_dict_y_to_y are
    # exactly all y that are left, and they have at least one alignment. This is not true if we transpose the 
    # dictionary however. The opposite holds when we want to estimate errors for x in the initialization,
    # then we do it on the transposed dictionary because we want to estimate the error on all x_i 
    epsilon_y, errors_y = initialize_data.estimate_epsilon(y, alignment_results_dict_y_to_y_transposed, alignment_results_dict_y_to_y, y_to_y, y_to_y, epsilon_y_backup, errors_y_backup, params)

    for y_j in y_j_with_no_good_alignment:
        assert y_j not in errors_y
        assert y_j not in epsilon_y
        print("Letting errors and epsilon for ", y_j, "be old estimates of: {0} and {1}".format(old_errors_y[y_j], old_epsilon_y[y_j]))
        errors_y[y_j] = old_errors_y[y_j]
        epsilon_y[y_j] = old_epsilon_y[y_j]

    return epsilon_y, errors_y


# def update_epsilon_y_min_distance_to_x_i_original(y, x, epsilon_x, errors_x, y_to_x, params): 
#     epsilon_y = {}
#     errors_y = {}
#     y_without_alignment_to_x = set()
#     if params.single_core:
#         for j, y_j in enumerate(y): 
#             x_i_original = y_to_x[y_j]
#             y_j, x_i, stats = align.ssw_alignment_helper( (x_i_original, y_j, x[x_i_original], y[y_j], j,j) )
#             if stats:
#                 matches, mismatches, indels = stats[2] 
#                 alignment_length = matches + mismatches + indels
#                 edit_distance = mismatches + indels
#                 epsilon_y[y_j] = abs( (edit_distance - epsilon_x[x_i]*alignment_length) ) / alignment_length
#                 errors_y[y_j] = abs(edit_distance - errors_x[x_i])
#                 if edit_distance < errors_x[x_i]:
#                     print("OMG, less errors than error rate of x_i!!", edit_distance, errors_x[x_i])
#                 # epsilon_y[y_j] = abs( (edit_distance - epsilon_x[x_i]*len(x[x_i])) ) / alignment_length
#                 # epsilon_y[y_j] = (edit_distance - epsilon_x[x_i]*len(x[x_i])) / alignment_length
#             else:
#                 print("No alignment between {0} and {1}. Removing y_j because it has diverged.".format(y_j, x_i))
#                 write_output.logger("No alignment between {0} and {1}. Removing y_j because it has diverged.".format(y_j, x_i), params.logfile)
#                 y_without_alignment_to_x.add(y_j)
#     else:
#         ####### parallelize alignment #########
#         # pool = Pool(processes=mp.cpu_count())
#         original_sigint_handler = signal.signal(signal.SIGINT, signal.SIG_IGN)
#         signal.signal(signal.SIGINT, original_sigint_handler)
#         pool = Pool(processes=mp.cpu_count())
#         try:
#             res = pool.map_async(align.ssw_alignment_helper, [(y_to_x[y_j], y_j, x[y_to_x[y_j]], y[y_j], j,j) for j, y_j in enumerate(y) ] )
#             alignment_results =res.get(999999999) # Without the timeout this blocking call ignores all signals.
#         except KeyboardInterrupt:
#             print("Caught KeyboardInterrupt, terminating workers")
#             pool.terminate()
#             sys.exit()
#         else:
#             print("Normal termination")
#             pool.close()
#         pool.join()
#         for y_j, x_i, stats in alignment_results:
#             if stats:
#                 matches, mismatches, indels = stats[2] 
#                 alignment_length = matches + mismatches + indels
#                 edit_distance = mismatches + indels
#                 epsilon_y[y_j] = abs( (edit_distance - epsilon_x[x_i]*len(x[x_i])) ) / alignment_length
#                 errors_y[y_j] = abs(edit_distance - errors_x[x_i])
#                 if edit_distance < errors_x[x_i]:
#                     print("OMG, less errors than error rate of x_i!!", edit_distance, errors_x[x_i])

#             else:
#                 print("No alignment between {0} and {1}. Removing y_j because it has diverged.".format(y_j, x_i))
#                 write_output.logger("No alignment between {0} and {1}. Removing y_j because it has diverged.".format(y_j, x_i), params.logfile)
#                 y_without_alignment_to_x.add(y_j)

#     return epsilon_y, errors_y, y_without_alignment_to_x


def update_y_j(j, y_j, x_i_original, y_j_seq, x_to_y_j_alignments, alpha_ji_for_y_j, epsilon, min_p, nr_errors_y_j, step):
    # p_value_threshold = min_p/(len(y_j_seq))
    # corner case, only one read connected to y
    p_values = []
    weighted_counts = []
    if nr_errors_y_j == 0:
        y_j_new = y_j_seq 
        print("Nr errors inferred to 0 ")
        return y_j, y_j_new 

    if len(x_to_y_j_alignments) <= 1:
        # to work for multiprocessing: A module-level function is a function which is defined at module level, 
        # that means it is not an instance method of a class, it's not nested within another function, 
        # and it is a "real" function with a name, not a lambda function.
        # from http://stackoverflow.com/questions/16439301/cant-pickle-defaultdict
        y_j_new = y_j_seq 
        print("HERE: update y had only one nbr ", x_to_y_j_alignments.keys())
        return y_j, y_j_new 

    if j % 20 == 0:
        print("fixing transcript", j  )

    print("FIXING", x_i_original)

    # max_significance_y_j = 1.0
    # min_epsilon = 1.0 
    # p_value_threshold = min_p/(2*len(y_j_seq)) #max(min_p/(2*len(y_j_seq)), max_significance_y_j/min_epsilon)
    p_value_threshold = min_p/10**(2*(step-1))
    print("P-value step", step, "=", p_value_threshold)
    total_depth = len(x_to_y_j_alignments.values())
    # print("depth:", total_depth)

    y_j_new = [] 
    # y_length = len(x_to_y_j_position_variants_depth_container)
    y_length = 2*len(y_j_seq) + 1
    mu_p_subs = 0
    mu_p_del = 0
    variance_p_subs = 0
    variance_p_del = 0

    # TODO: create error probability for indels longer than 1!
    # p_hat_tot_del = 0
    # p_hat_tot_subs = 0
    # for x_i in x_to_y_j_alignments:
    #     if x_i == x_i_original:
    #         print("BUG!!!!")
    #         sys.exit()
    #     p_hat_tot_del += epsilon[x_i]
    #     p_hat_tot_subs += epsilon[x_i] #/3.0

        # mu_p_subs += (epsilon[x_i]/3.0)*alpha_ji_for_y_j[x_i]
        # variance_p_subs += (epsilon[x_i]/3.0)*(1-(epsilon[x_i]/3.0))*alpha_ji_for_y_j[x_i]**2
        # mu_p_del += (epsilon[x_i]/2.0)*alpha_ji_for_y_j[x_i]
        # variance_p_del += (epsilon[x_i]/2.0)*(1-(epsilon[x_i]/2.0))*alpha_ji_for_y_j[x_i]**2

    # bin_p_hat_del = p_hat_tot_del/total_depth
    # bin_p_hat_subs = p_hat_tot_subs/total_depth
    # k_del = scipy.stats.binom.isf(p_value_threshold, total_depth, bin_p_hat_del)
    # k_subs = scipy.stats.binom.isf(p_value_threshold, total_depth, bin_p_hat_subs)

    # print("bin p del",bin_p_hat_del, k_del, 'depth:', total_depth)
    # print("bin p subs",bin_p_hat_subs, k_subs, 'depth:', total_depth)
    # weights = list(alpha_ji_for_y_j.values())
    # quantile_del = sum( [sum( random.sample(weights, min(int(k_del), total_depth -1)) ) for i in range(100)] ) / 100.0
    # quantile_subs = sum( [sum(random.sample(weights, min(int(k_subs), total_depth -1)) ) for i in range(100)] ) / 100.0
    # print("QUANTILE subs", quantile_subs, "QUANTILE del", quantile_del)


    # QUANTILE_THRESHOLD_DEL = quantile_del # max(norma_approx_quantile_del, gamma_approx_quantile_del)  #norma_approx_quantile # max(norma_approx_quantile, weight_heuristic_quantile)
    # QUANTILE_THRESHOLD_SUBS = quantile_subs # max(norma_approx_quantile_subs, gamma_approx_quantile_subs)  #norma_approx_quantile # max(norma_approx_quantile, weight_heuristic_quantile)
    nr_corrections = 0

    potential_positions_to_correct = [] # list of tuples (position, p-value/quantile, v_to_substitute )
    for p in range(y_length):
        y_j_p = y_j_seq[p//2] if p % 2 == 1 else '-'        
        observed_weight_p = 0
        # observed_count_p = 0
        # observed_epsilon = 1.0
        alternative_variants_p = defaultdict(float)
        x_to_y_j_position_variants_depth_container_p = defaultdict(int)
        for x_i in x_to_y_j_alignments:
            if x_i == x_i_original:
                print("BUG!!!!")
                sys.exit()
            x_i_p = x_to_y_j_alignments[x_i][p]
            x_to_y_j_position_variants_depth_container_p[x_i_p] += 1

            if x_i_p == y_j_p:
                observed_weight_p += alpha_ji_for_y_j[x_i]
                # observed_count_p += 1
                # if observed_epsilon > p_value_threshold:
                #     observed_epsilon *= epsilon[x_i]
            else:
                # TODO: this probability is theoretically wrong, but does it matter in practice?
                alternative_variants_p[x_i_p] += alpha_ji_for_y_j[x_i] #*(1-epsilon[x_i])

        # if y_j_p == "-" and observed_weight_p > QUANTILE_THRESHOLD_DEL: # and observed_epsilon <= p_value_threshold: 
        # # if p_value_not_error < p_value_threshold or len(alternative_variants_p) == 0: # multiple testing each transcript
        #     y_j_new.append(y_j_p)
        #     # if x_to_y_j_position_variants_depth_container_p[y_j_p] < 10 and total_depth > 80:
        #         # print("Low support not corrected", x_to_y_j_position_variants_depth_container_p[y_j_p])
        #         # print("total sum weights:", sum(map(lambda x: alpha_ji_for_y_j[x] ,alpha_ji_for_y_j)) ) 
        #         # print( "norm approx:", norma_approx_quantile_del, "gamma approx:", gamma_approx_quantile_del, "observed count:", observed_weight_p, "lower bound observed outcome pval lower than:", observed_epsilon, "pval threshold for y_j:", p_value_threshold,  "depth:", total_depth )
        #         # print( "norm approx:", norma_approx_quantile, "old normal approx:", "largest weight heur:", weight_heuristic_quantile, "largest weight heur new:", weight_heuristic_quantile_new, "sampling quantile:", sampling_quantile, "observed count:", observed_weight_p)
        #         # print( "norm approx:", norma_approx_quantile, "old normal approx:", norma_approx_quantile2, "largest weight heur:", weight_heuristic_quantile, "largest weight heur new:", weight_heuristic_quantile_new, "sampling quantile:", sampling_quantile, "observed count:", observed_weight_p)
        #         # print('pval_norm:', p_value_not_error, "quantile normal:", quantile_theshold__, "quantile emp:", quantile_emp , "quantile lower bound:", weight_heuristic_quantile, "observed count:", observed_weight_p, mu_p, variance_p)
        #         # print(map(lambda x: epsilon[x], epsilon))
        #         # print(x_to_y_j_position_variants_depth_container_p)
        #         # weighted_counts.append(observed_weight_p)
        #         # print("sampling quantile:", sampling_quantile,"largest weight heur:", weight_heuristic_quantile, "largest weight heur new:", weight_heuristic_quantile_new)

        # elif y_j_p != "-" and observed_weight_p > QUANTILE_THRESHOLD_SUBS:
        #     y_j_new.append(y_j_p)
        #     # if x_to_y_j_position_variants_depth_container_p[y_j_p] < 10 and total_depth > 80:
        #         # print("Low support not corrected", x_to_y_j_position_variants_depth_container_p[y_j_p])
        #         # print("total sum weights:", sum(map(lambda x: alpha_ji_for_y_j[x] ,alpha_ji_for_y_j)) ) 
        #         # print("norm approx:", norma_approx_quantile_subs, "gamma approx:", gamma_approx_quantile_subs, "observed count:", observed_weight_p, "lower bound observed outcome pval lower than:", observed_epsilon, "pval threshold for y_j:", p_value_threshold,  "depth:", total_depth )
        #         # print("sampling quantile:", sampling_quantile, "largest weight heur:", weight_heuristic_quantile, "largest weight heur new:", weight_heuristic_quantile_new)
        #         # print( "norm approx:", norma_approx_quantile, "old normal approx:", "largest weight heur:", weight_heuristic_quantile, "largest weight heur new:", weight_heuristic_quantile_new, "sampling quantile:", sampling_quantile, "observed count:", observed_weight_p)
        #         # print( "norm approx:", norma_approx_quantile, "old normal approx:", norma_approx_quantile2, "largest weight heur:", weight_heuristic_quantile, "largest weight heur new:", weight_heuristic_quantile_new, "sampling quantile:", sampling_quantile, "observed count:", observed_weight_p)
        #         # print('pval_norm:', p_value_not_error, "quantile normal:", quantile_theshold__, "quantile emp:", quantile_emp , "quantile lower bound:", weight_heuristic_quantile, "observed count:", observed_weight_p, mu_p, variance_p)
        #         # print(map(lambda x: epsilon[x], epsilon))
        #         # print(x_to_y_j_position_variants_depth_container_p)
        #         # weighted_counts.append(observed_weight_p)
        # else:
        #     # print(x_to_y_j_position_variants_depth_container_p)

            # print("COREECTION HERE:", p_value_not_error, x_to_y_j_position_variants_depth_container_p[y_j_p], observed_weight_p, mu_p,variance_p )
            # print("CORRECTION HERE:", x_to_y_j_position_variants_depth_container_p[y_j_p], observed_weight_p, weight_heuristic_quantile)
            # print("z_score", (observed_weight_p - mu_p)/ math.sqrt(variance_p))

            # correct to basepair of the read with largest alpha that has a different base pair than current bp in y_j
            # if len(alternative_variants_p) == 0:
            #     v_max = y_j_p
            # else:
            #     for x_i, alpha in alpha_decreasing_order:
            #         if x_to_y_j_alignments[x_i][p] != y_j_p:
            #             v_max = x_to_y_j_alignments[x_i][p]
            #             break

        if len(alternative_variants_p) > 0:
            v_max, weighted_count_max = max(alternative_variants_p.items(), key=operator.itemgetter(1))
            if observed_weight_p > weighted_count_max:
                v_max = y_j_p
            nr_corrections += 1
        else:
            v_max = y_j_p

        potential_positions_to_correct.append((p, observed_weight_p, v_max))
        y_j_new.append(v_max)            

    ### new section ###
    potential_positions_to_correct.sort(key=lambda x: x[1])
    nr_potential_to_correct = len(potential_positions_to_correct)
    nr_to_correct = int(0.9*nr_errors_y_j) if nr_errors_y_j > 1 else int(nr_errors_y_j) # int(nr_errors_y_j) #
    print(nr_errors_y_j)
    y_old = [ y_j_seq[p//2] if p % 2 == 1 else '-' for p in range(y_length)]
    for p, q_val, v_max in potential_positions_to_correct[:nr_to_correct]:
        y_old[p] = v_max
    y_j_new = y_old
    ####################
    # if nr_errors_y_j < 5 and v_max is y_j_p --> how many does not change bp?
    y_j_temp_new = filter(lambda v: v != '-', y_j_new)
    y_j_new = "".join(y_j_temp_new)
    print("INFERRED ERRORS Y:", nr_errors_y_j, "NR corrected:", len(potential_positions_to_correct[:nr_to_correct]) )
    if nr_corrections > nr_errors_y_j:
        print("OVERCORRECTION!!")

    return y_j, "".join(y_j_new) #tau_j_s_plus_one






# def calc_binomial(n_values, p_e):
#     binomial_hashtable = {}  # key is (n,k) tuple and value is probability

#     for n in n_values:
#         # infer what is the largest k we need to
#         # calculate probabilities for set threshold to 0.999
#         # it won't affect the estimations of its 0.999 or 0.9999 for example.
#         k = 0
#         while True:
#             # if n*p_e >= 10:
#             #     # fast normal approximation is good
#             #     pass

#             if p_e <= 1.0/n and n >= 20:
#                 # poisson approximation is good
#                 mu = n*p_e
#                 prob = scipy.stats.poisson.cdf(k, mu)
#                 binomial_hashtable[(n, k)] = prob
#             else:
#                 # exact binomial
#                 prob = scipy.stats.binom.cdf(k, n, p_e)
#                 binomial_hashtable[(n, k)] = prob

#             if prob > 0.999:
#                 break
#             k += 1

#     return binomial_hashtable
