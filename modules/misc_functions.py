from collections import defaultdict
from hitemmodules.io import write_output


def transpose(dct):
    d = defaultdict(dict)
    for key1, inner in dct.items():
        for key2, value in inner.items():
            d[key2][key1] = value
    return d


# def remove_starting_points_without_alignment_to_original_starting_sequence(alignment_results_dict, y_without_alignment_to_x, params):
#     """
#         clean out the y's (starting points) that do not have an x_i (a read) aligning to any y_j in the graph.
#         This will happen because of the --support_cutoff parameter --- some of the lowest quality
#         reads will not support any y.
#     """

#     removed_counter = 0
#     for y_j in list(y_without_alignment_to_x):
#         del alignment_results_dict[y_j]
#         y_without_alignment_to_x.remove(y_j)
#         removed_counter +=1

#     if removed_counter > 0:
#         write_output.logger("{0} starting points have been filtered out in this step because original sequence could not be aligned to.".format(removed_counter), params.logfile)        
#         write_output.logger("{0} starting points left.".format(len(alignment_results_dict)), params.logfile)
#         return remove_starting_points_without_alignment_to_original_starting_sequence(alignment_results_dict, y_without_alignment_to_x, params)
#         # write_output.logger("A total of {0} y's (starting points) have been filtered out at this step.".format(len(removed_ys)), params.logfile)

#     else:
#         return alignment_results_dict


def remove_starting_points_without_original_starting_sequence(alignment_results_dict, x_to_y, params):
    """
        clean out the y's (starting points) that do not have an x_i (a read) aligning to any y_j in the graph.
        This will happen because of the --support_cutoff parameter --- some of the lowest quality
        reads will not support any y.
    """

    alignment_results_dict_x_y = transpose(alignment_results_dict)
    current_supporting_reads = set(alignment_results_dict_x_y.keys())
    all_y_j_with_original_x_i_in_graph = set([x_to_y[x_i] for x_i in current_supporting_reads])
    removed_counter = 0
    for y_j in list(alignment_results_dict.keys()):
        if y_j not in all_y_j_with_original_x_i_in_graph or len(alignment_results_dict[y_j]) == 0:
            del alignment_results_dict[y_j]
            removed_counter +=1

    if removed_counter > 0:
        write_output.logger("{0} starting points have been filtered out in this step because original sequence missing.".format(removed_counter), params.logfile)        
        write_output.logger("{0} starting points left.".format(len(alignment_results_dict)), params.logfile)
        return remove_starting_points_without_original_starting_sequence(alignment_results_dict, x_to_y, params)
        # write_output.logger("A total of {0} y's (starting points) have been filtered out at this step.".format(len(removed_ys)), params.logfile)

    else:
        return alignment_results_dict


def remove_data(alignment_results_dict, x, y, x_to_y, y_to_x, params):
    x_before = len(x)
    y_before = len(y)
    write_output.logger("{0} and {1} y's and x's before.".format(y_before, x_before), params.logfile)

    filtered_alignment_results_dict = remove_starting_points_without_original_starting_sequence(alignment_results_dict, x_to_y, params)

    for y_j in list(y.keys()):
        if y_j not in filtered_alignment_results_dict:
            del y[y_j]
            del y_to_x[y_j]

    filtered_alignment_results_dict_x_y = transpose(filtered_alignment_results_dict)

    for x_i in list(x.keys()):
        if x_i not in filtered_alignment_results_dict_x_y:
            del x[x_i]
            del x_to_y[x_i]

    x_after = len(x)
    y_after = len(y)

    write_output.logger("{0} and {1} y's and x's before.".format(y_after, x_after), params.logfile)
    return filtered_alignment_results_dict


def print_graph_complexity(alignment_results_dict, params):

    y_support_histogram = defaultdict(int)
    for y in alignment_results_dict:
        write_output.logger("y-support y{0}: {1}".format(y, len(alignment_results_dict[y])), params.logfile)
        y_support_histogram[len(alignment_results_dict[y])] +=1

    alignment_results_dict_x_y = transpose(alignment_results_dict)

    x_support_histogram = defaultdict(int)
    for x in alignment_results_dict_x_y:
        write_output.logger("x-support {0}: {1}".format(x, len(alignment_results_dict_x_y[x])), params.logfile)
        x_support_histogram[len(alignment_results_dict_x_y[x])] += 1

    for support, count in sorted(x_support_histogram.items()):
        write_output.logger("reads supporting {0} y's: count {1}".format(support, count), params.logfile)

    for support, count in sorted(y_support_histogram.items()):
        write_output.logger("y linked to {0} reads: count {1}".format(support, count), params.logfile)
