from __future__ import print_function
import signal
from multiprocessing import Pool
import multiprocessing as mp
import sys
import math

import re 

from decimal import * 
getcontext().prec = 100

from modules import functions
from modules.SW_alignment_module import sw_align_sequences_keeping_accession
from modules.edlib_alignment_module import edlib_align_sequences_keeping_accession
from modules import ccs_info

def do_statistical_tests_per_edge(nearest_neighbor_graph, C, X, read_partition, ccs_dict, params):
    p_values = {}
    statistical_test_edges = {}
    for c_acc in nearest_neighbor_graph:
        statistical_test_edges[c_acc] = {}
        reads_to_c = { x_acc : X[x_acc] for x_acc in read_partition[c_acc]}
        read_alignments_to_c = read_partition[c_acc]
        for t_acc in nearest_neighbor_graph[c_acc]:
            read_alignments_to_t = read_partition[t_acc]  # { read_acc : (t_aln, c_aln, (match, mism, indel)), ...}  #{ x_acc : X[x_acc] for x_acc in read_partition[t_acc]}
            edge_specific_reads = { x_acc : X[x_acc] for acc in [c_acc, t_acc] for x_acc in read_partition[acc]}
            if ccs_dict:
                reduced_ccs_dict_for_test = {x_acc : ccs_dict[x_acc] for x_acc in list(read_partition[c_acc].keys()) + list(read_partition[t_acc].keys()) if x_acc in ccs_dict} 
            else:
                reduced_ccs_dict_for_test = {}

            statistical_test_edges[c_acc][t_acc] = (c_acc, t_acc, C[c_acc], C[t_acc], reads_to_c, read_alignments_to_t, read_alignments_to_c, params.ignore_ends_len, reduced_ccs_dict_for_test, params.max_phred_q_trusted) # ADD ALGINMENTS HERE FOR SPEEDUP

    if params.nr_cores == 1:
        for c_acc in statistical_test_edges:
            p_values[c_acc] = {}
            for t_acc in statistical_test_edges[c_acc]:
                (c_acc, t_acc, p_value, mult_factor_inv, k, N_t, variants) = statistical_test(*statistical_test_edges[c_acc][t_acc])
                p_values[c_acc][t_acc] =  (p_value, mult_factor_inv, k, N_t, variants)

    else:
        ####### parallelize statistical tests #########
        # pool = Pool(processes=mp.cpu_count())
        original_sigint_handler = signal.signal(signal.SIGINT, signal.SIG_IGN)
        signal.signal(signal.SIGINT, original_sigint_handler)
        pool = Pool(processes=params.nr_cores)

        for c_acc in list(statistical_test_edges.keys()):
            if len(statistical_test_edges[c_acc]) == 0:
                del statistical_test_edges[c_acc]

        try:
            res = pool.map_async(statistical_test_helper, [ (statistical_test_edges[c_acc][t_acc], {}) for c_acc in statistical_test_edges for t_acc in statistical_test_edges[c_acc]  ] )
            statistical_test_results =res.get(999999999) # Without the timeout this blocking call ignores all signals.
        except KeyboardInterrupt:
            print("Caught KeyboardInterrupt, terminating workers")
            pool.terminate()
            sys.exit()
        else:
            print("Normal termination")
            pool.close()
        pool.join()

        for c_acc in nearest_neighbor_graph:
            p_values[c_acc] = {}

        for (c_acc, t_acc, p_value, mult_factor_inv, k, N_t, variants) in statistical_test_results:
            p_values[c_acc][t_acc] = (p_value, mult_factor_inv, k, N_t, variants)


    print("Total number of tests performed this round:", len([ 1 for c_acc in statistical_test_edges for t_acc in statistical_test_edges[c_acc]]) )

    return p_values


def statistical_test_helper(arguments):
    args, kwargs = arguments
    return statistical_test(*args, **kwargs)


def get_mask_start_and_end(aln_t, aln_c):
    mask_start, mask_end = 0, len(aln_t)
    p = r"[-]+"
    for m in re.finditer(p, aln_t):
        # print(m.start(), m.end())
        if m.start() == 0:
            mask_start = m.end()
        if m.end() == len(aln_t):
            mask_end = m.start()

    for m in re.finditer(p, aln_c):
        # print(m.start(), m.end())
        if m.start() == 0:
            assert mask_start == 0
            mask_start = m.end()
        if m.end() == len(aln_t):
            assert mask_end == len(aln_t)
            mask_end = m.start()
    return mask_start, mask_end


def get_read_errors(read_alignments_to_c, read_alignments_to_t):
    errors = {}
    for read_acc in read_alignments_to_t:
        aln_t, aln_read, (matches, mismatches, indels) = read_alignments_to_t[read_acc]
        insertions, deletions, substitutions = functions.read_errors_from_alignment(aln_t, aln_read)
        errors[read_acc] = (insertions, deletions, substitutions)

    for read_acc in read_alignments_to_c:
        aln_c, aln_read, (matches, mismatches, indels) = read_alignments_to_c[read_acc]
        insertions, deletions, substitutions = functions.read_errors_from_alignment(aln_c, aln_read)
        errors[read_acc] = (insertions, deletions, substitutions)       
    return errors


from time import time

def merge_two_dicts(x, y):
    """Given two dicts, merge them into a new dict as a shallow copy."""
    z = x.copy()
    z.update(y)
    return z

def get_variant_coordinates(t_seq, c_seq, aln_t, aln_c, variants):
    variant_coords_t = {}
    variant_coords_c = {}
    alignment_c_to_t = {} # t-coordinate as key and alignment of c in that region as value
    alignment_t_to_c = {} # c-coordinate as key and alignment of t in that region as value

    for (i,p_t,p_c) in variants:
        t_seq_piece = "".join([n for n in aln_t[:i+1] if n != "-"])
        c_seq_piece = "".join([n for n in aln_c[:i+1] if n != "-"])

        if p_c == "-": # candidate has deletion
            v = t_seq[len(t_seq_piece)-1 ]
            p = "[{variant}]+".format(variant=v)
            m_f = re.match(p, t_seq[len(t_seq_piece)-1 +1: ])
            m_r = re.match(p, t_seq[len(t_seq_piece)-1 : : -1 ])
            if m_f and m_r:
                u_v = len(m_f.group()) + len(m_r.group())
            elif m_f:
                u_v = len(m_f.group())  
            elif m_r:
                u_v = len(m_r.group())
            else:
                u_v = 1
            # u_v = max(len(m_f.group()), len(m_r.group())) if m_f or m_r else 1
            variant_coords_t[len(t_seq_piece)-1] = ("D", "-", u_v ) 
            variant_coords_c[len(c_seq_piece)-1 +1] = ("D", "-", u_v ) # Deletion: we will get the phred base call from the pos immmediately to the right
            alignment_c_to_t[len(t_seq_piece)-1] = aln_c[i-1 : i+ u_v + 1]
            alignment_t_to_c[len(c_seq_piece)-1 +1] = aln_t[i-1 : i+ u_v + 1]

        elif p_t == "-": # candidate has insertion
            v = c_seq[len(c_seq_piece)-1 ]
            p = "[{variant}]+".format(variant=v)

            m_f = re.match(p, t_seq[len(t_seq_piece)-1 +1 : ])
            m_r = re.match(p, t_seq[len(t_seq_piece)-1 : : -1 ])
            if m_f and m_r:
                u_v = len(m_f.group()) + len(m_r.group()) + 1 # +1 because there are x+1 ways to insert the same character into a a homoplymer of x characters in order to make an insertion of length x+1 
            elif m_f:
                u_v = len(m_f.group()) + 1  # +1 because there are x+1 ways to insert the same character into a a homoplymer of x characters in order to make an insertion of length x+1 
            elif m_r:
                u_v = len(m_r.group()) + 1 # +1 because there are x+1 ways to insert the same character into a a homoplymer of x characters in order to make an insertion of length x+1 
            else:
                u_v = 1
            # m =re.match(p, t_seq[len(t_seq_piece)-1 +1 : ]) # Deletion in t: The invariant lengths start one coord to the right
            # u_v = len(m.group()) + 1 if m else 1 # +1 because there are x+1 ways to insert the same character into a a homoplymer of x characters in order to make an insertion of length x+1 

            variant_coords_t[len(t_seq_piece)-1 +1] = ("I", p_c, u_v ) # Deletion: we will get the phred base call from the pos immmediately to the right
            variant_coords_c[len(c_seq_piece)-1] = ("I", p_c, u_v )
            alignment_c_to_t[len(t_seq_piece)-1 +1] = aln_c[i-1 : i+ u_v + 1]
            alignment_t_to_c[len(c_seq_piece)-1] = aln_t[i-1 : i+ u_v + 1]

        else: # get coordinate on ref and cand if substitution
            variant_coords_t[len(t_seq_piece)-1] = ("S", p_c, 1 )
            variant_coords_c[len(c_seq_piece)-1] = ("S", p_c, 1 )
            alignment_c_to_t[len(t_seq_piece)-1] = aln_c[i-1 : i + 2]
            alignment_t_to_c[len(c_seq_piece)-1] = aln_t[i-1 : i + 2]

    return variant_coords_t, variant_coords_c, alignment_c_to_t, alignment_t_to_c



def get_support(read_alignments_to_c, variant_coords_c, read_alignments_to_t, variant_coords_t, alignment_c_to_t):
    reads_support_c = []
    for read_acc in read_alignments_to_c:
        aln_c, aln_read, (matches, mismatches, indels) = read_alignments_to_c[read_acc]
        c_seq_to_coord_in_almnt = [j for j, nucl in enumerate(aln_c) if aln_c[j] != "-"]
        support = True
        for i in variant_coords_c:
            v_type, v_nucl, u_v = variant_coords_c[i] 
            alnmt_pos = c_seq_to_coord_in_almnt[i]

            # require exact match over 3mer if S, or I,D in non-homopolymenr region.
            # If homopolymenr region of size u_v, require length match of homopolymenr
            exact_match = aln_read[alnmt_pos -1: alnmt_pos + u_v + 1] == aln_c[alnmt_pos -1: alnmt_pos + u_v +1]
            if not exact_match:
                support = False
                break
        if support:
            reads_support_c.append(read_acc)


    reads_support_c_from_t = []
    for read_acc in read_alignments_to_t:
        aln_t, aln_read, (matches, mismatches, indels) = read_alignments_to_t[read_acc]
        t_seq_to_coord_in_almnt = [j for j, nucl in enumerate(aln_t) if aln_t[j] != "-"]
        support = True
        for i in variant_coords_t:
            v_type, v_nucl, u_v = variant_coords_t[i] 
            alnmt_pos = t_seq_to_coord_in_almnt[i]

            c_align_snippet = alignment_c_to_t[i]
            if v_type == "I": # will need to check shifterd snippet on t alignment because base is indexed immediately to the right
                exact_match_to_c = aln_read[alnmt_pos -2: alnmt_pos + u_v] == c_align_snippet 
                # print(aln_read[alnmt_pos -2: alnmt_pos + u_v], c_align_snippet)
            else:
                exact_match_to_c = aln_read[alnmt_pos -1: alnmt_pos + u_v + 1] == c_align_snippet
                # print( aln_read[alnmt_pos -1: alnmt_pos + u_v + 1], c_align_snippet)
            
            # require exact match over 3mer if S, or I,D in non-homopolymenr region.
            # If homopolymenr region of size u_v, require length match of homopolymenr
            if not exact_match_to_c:
                support = False
                break
        if support:
            reads_support_c_from_t.append(read_acc)   

    print("NEW2 support", len(reads_support_c), "and", len(reads_support_c_from_t) )
    return reads_support_c + reads_support_c_from_t


def get_support_from_c(read_alignments_to_c, variant_coords_c):
    reads_support = []

    for read_acc in read_alignments_to_c:
        aln_c, aln_read, (matches, mismatches, indels) = read_alignments_to_c[read_acc]
        c_seq_to_coord_in_almnt = [j for j, nucl in enumerate(aln_c) if aln_c[j] != "-"]
        support = True
        for i in variant_coords_c:
            v_type, v_nucl, u_v = variant_coords_c[i] 
            alnmt_pos = c_seq_to_coord_in_almnt[i]

            # require exact match over 3mer if S, or I,D in non-homopolymenr region.
            # If homopolymenr region of size u_v, require length match of homopolymenr
            exact_match = aln_read[alnmt_pos -1: alnmt_pos + u_v + 1] == aln_c[alnmt_pos -1: alnmt_pos + u_v +1]
            if not exact_match:
                support = False
                break
            # if v_type == "S":
            #     if not exact_match:
            #         support = False


            # elif v_type == "I":
            #     if not exact_match:
            #         support = False

            # else:
            #     assert v_type == "D"
            #     if not exact_match:
            #         support = False

        if support:
            reads_support.append(read_acc)
    return reads_support


# def get_read_ccs_probabilities_from_reads_to_c(read_alignments_to_c, variant_coords_c, ccs_dict, errors, max_phred_q_trusted):
#     probabilities_c = {}
#     tot_errors = float(sum([ i+d+s for i, d, s in errors.values()]))
#     tot_errors = tot_errors if tot_errors > 1 else 1.0
#     subs_ratio = sum([ s for i, d, s in errors.values()]) / tot_errors
#     ins_ratio = sum([ i for i, d, s in errors.values()]) / tot_errors
#     del_ratio = sum([ d for i, d, s in errors.values()]) / tot_errors

#     for read_acc in read_alignments_to_c:
#         aln_c, aln_read, (matches, mismatches, indels) = read_alignments_to_c[read_acc]
#         c_seq_to_coord_in_almnt = [j for j, nucl in enumerate(aln_c) if aln_c[j] != "-"]
#         read_seq_to_coord_in_almnt = [j for j, nucl in enumerate(aln_read) if aln_read[j] != "-"]
#         prob = 1.0

#         for i in variant_coords_c:
#             v_type, v_nucl, u_v = variant_coords_c[i] 
#             alnmt_pos = c_seq_to_coord_in_almnt[i]

#             # require exact match over 3mer if S, or I, D in non-homopolymenr region.
#             # If homopolymenr region of size u_v, require length match of homopolymenr
#             exact_match = aln_read[alnmt_pos -1: alnmt_pos + u_v + 1] == aln_c[alnmt_pos -1: alnmt_pos + u_v +1]
#             if  exact_match:
#                 read_coord = len("".join([n for n in aln_read[ : alnmt_pos +1] if n != "-"])) - 1
#                 ccs_coord = ccs_dict[ read_acc ].read_aln_to_ccs_coord(aln_read, read_coord)

#             else:
#                 read_snippet_up_to_variant = "".join([n for n in aln_read[ : alnmt_pos +1] if n != "-"])
#                 if aln_read[alnmt_pos] == "-":
#                     alnmt_pos_after_del = read_seq_to_coord_in_almnt[len(read_snippet_up_to_variant)]
#                     read_snippet_up_to_variant = "".join([n for n in aln_read[ : alnmt_pos_after_del +1] if n != "-"])
#                 # else:
#                 #     read_snippet_up_to_variant = "".join([n for n in aln_read[ : alnmt_pos +1] if n != "-"])
#                 p = "[{variant}]+".format(variant=read_snippet_up_to_variant[-1])
#                 m = re.match(p, read_snippet_up_to_variant[ : : -1 ])
#                 read_coord = len(read_snippet_up_to_variant) - len(m.group())
#                 ccs_coord = ccs_dict[ read_acc ].read_aln_to_ccs_coord(aln_read, read_coord)


#             q_qual = ccs_dict[read_acc].qual[ccs_coord]
#             print(ccs_dict[read_acc].qual[ccs_coord -1: ccs_coord +5 ], ccs_dict[read_acc].seq[ccs_coord-1: ccs_coord +5 ])

#             # To map quality values
#             # [A, B] --> [a, b]  [3, 93] --> [3, 43]
#             # (x - A)*(b-a)/(B-A) + a
#             q_qual_mapped = (q_qual - 3)*(max_phred_q_trusted - 3.0)/(90.0) + 3

#             if u_v > 1: # probability in a homopolymenr region has all it's uncertainty attributed to the lenght of the homopolymer that 
#                 p_error =  (10**(-q_qual_mapped/10.0))
#             elif v_type == "S":
#                 p_error =  ((10**(-q_qual_mapped/10.0))*subs_ratio)/3.0 # probability that its an identical substitution error from a base call uncertainty
#             elif  v_type == "I":
#                 p_error =  ((10**(-q_qual_mapped/10.0))*ins_ratio)/4.0 # probability that its an identical insertion error from a base call uncertainty
#             elif v_type == "D":
#                 p_error =  (10**(-q_qual_mapped/10.0)) * del_ratio # probability that its a delation error from a base call uncertainty
#             else:
#                 p_error = -1 # If end up here, there is a bug

#             # print(p_error, q_qual_mapped, q_qual)
#             assert 0.0 < p_error < 1.0
#             prob *= p_error
#         probabilities_c[read_acc] = prob

#     return probabilities_c


def get_read_ccs_probabilities_c(read_alignments_to_c, variant_coords_c, alignment_t_to_c, ccs_dict, errors, max_phred_q_trusted):
    tot_errors = float(sum([ i+d+s for i, d, s in errors.values()]))
    tot_errors = tot_errors if tot_errors > 1 else 1.0
    subs_ratio = sum([ s for i, d, s in errors.values()]) / tot_errors
    ins_ratio = sum([ i for i, d, s in errors.values()]) / tot_errors
    del_ratio = sum([ d for i, d, s in errors.values()]) / tot_errors

    probabilities = {}
    reads_not_supporting_any_seq = set()

    for read_acc in read_alignments_to_c:
        aln_c, aln_read, (matches, mismatches, indels) = read_alignments_to_c[read_acc]
        c_seq_to_coord_in_almnt = [j for j, nucl in enumerate(aln_c) if aln_c[j] != "-"]
        read_seq_to_coord_in_almnt = [j for j, nucl in enumerate(aln_read) if aln_read[j] != "-"]
        prob = 1.0

        for i in variant_coords_c:
            v_type, v_nucl, u_v = variant_coords_c[i] 
            alnmt_pos = c_seq_to_coord_in_almnt[i]

            # require exact match over 3mer if S, or I, D in non-homopolymenr region.
            # If homopolymenr region of size u_v, require length match of homopolymenr
            exact_match_to_c = aln_read[alnmt_pos -1: alnmt_pos + u_v + 1] == aln_c[alnmt_pos -1: alnmt_pos + u_v +1]

            t_align_snippet = alignment_t_to_c[i]
            if v_type == "D": # will need to check shifted snippet on c alignment because base is indexed immediately to the right
                exact_match_to_t = aln_read[alnmt_pos -2: alnmt_pos + u_v] == t_align_snippet 
                # print("LOL HERE C", aln_read[alnmt_pos -2: alnmt_pos + u_v], t_align_snippet)
            else:
                exact_match_to_t = aln_read[alnmt_pos -1: alnmt_pos + u_v + 1] == t_align_snippet
                # print( aln_read[alnmt_pos -1: alnmt_pos + u_v + 1], t_align_snippet, aln_c[alnmt_pos -1: alnmt_pos + u_v +1])

            assert not (exact_match_to_c and exact_match_to_t) # they cannot both be perfect matches, that means a bug

            if  exact_match_to_c:
                read_coord = len("".join([n for n in aln_read[ : alnmt_pos +1] if n != "-"])) - 1
                ccs_coord = ccs_dict[ read_acc ].read_aln_to_ccs_coord(aln_read, read_coord)
            
            elif exact_match_to_t:
                print("HERE T!!")
                if v_type == "I":
                    read_coord = len("".join([n for n in aln_read[ : alnmt_pos +1] if n != "-"]))
                else:
                    read_coord = len("".join([n for n in aln_read[ : alnmt_pos +1] if n != "-"])) - 1
                ccs_coord = ccs_dict[ read_acc ].read_aln_to_ccs_coord(aln_read, read_coord)

            else:

                reads_not_supporting_any_seq.add(read_acc)
                prob = -1 
                break

            q_qual = ccs_dict[read_acc].qual[ccs_coord]
            # print(ccs_dict[read_acc].qual[ccs_coord -1: ccs_coord +5 ], ccs_dict[read_acc].seq[ccs_coord-1: ccs_coord +5 ], t_align_snippet)
            
            # To map quality values
            # [A, B] --> [a, b]  [3, 93] --> [3, 43]
            # (x - A)*(b-a)/(B-A) + a
            q_qual_mapped = (q_qual - 3)*(max_phred_q_trusted - 3.0)/(90.0) + 3

            if u_v > 1: # probability in a homopolymenr region has all it's uncertainty attributed to the lenght of the homopolymer that 
                p_error =  (10**(-q_qual_mapped/10.0))
            elif v_type == "S":
                p_error =  ((10**(-q_qual_mapped/10.0))*subs_ratio)/3.0 # probability that its an identical substitution error from a base call uncertainty
            elif  v_type == "I":
                p_error =  ((10**(-q_qual_mapped/10.0))*ins_ratio)/4.0 # probability that its an identical insertion error from a base call uncertainty
            elif v_type == "D":
                p_error =  (10**(-q_qual_mapped/10.0)) * del_ratio # probability that its a delation error from a base call uncertainty
            else:
                print("Bad type", v_type, u_v)
                sys.exit()

            # print(p_error, q_qual_mapped, q_qual)
            prob *= p_error

        if prob >= 0:
            probabilities[read_acc] = prob
        else:
            print("read did not support c or t in variant region")
            pass

    print("Total probs:", len(probabilities), "non-informative:", len(reads_not_supporting_any_seq))
    return probabilities, reads_not_supporting_any_seq


def get_read_ccs_probabilities_t(read_alignments_to_t, variant_coords_t, alignment_c_to_t, ccs_dict, errors, max_phred_q_trusted):
    tot_errors = float(sum([ i+d+s for i, d, s in errors.values()]))
    tot_errors = tot_errors if tot_errors > 1 else 1.0
    subs_ratio = sum([ s for i, d, s in errors.values()]) / tot_errors
    ins_ratio = sum([ i for i, d, s in errors.values()]) / tot_errors
    del_ratio = sum([ d for i, d, s in errors.values()]) / tot_errors

    probabilities = {}
    reads_not_supporting_any_seq = set()

    for read_acc in read_alignments_to_t:
        aln_t, aln_read, (matches, mismatches, indels) = read_alignments_to_t[read_acc]
        t_seq_to_coord_in_almnt = [j for j, nucl in enumerate(aln_t) if aln_t[j] != "-"]
        read_seq_to_coord_in_almnt = [j for j, nucl in enumerate(aln_read) if aln_read[j] != "-"]
        prob = 1.0

        for i in variant_coords_t:
            v_type, v_nucl, u_v = variant_coords_t[i] 
            alnmt_pos = t_seq_to_coord_in_almnt[i]

            # require exact match over 3mer if S, or I, D in non-homopolymenr region.
            # If homopolymenr region of size u_v, require length match of homopolymenr
            exact_match_to_t = aln_read[alnmt_pos -1: alnmt_pos + u_v + 1] == aln_t[alnmt_pos -1: alnmt_pos + u_v +1]

            c_align_snippet = alignment_c_to_t[i]
            if v_type == "I": # will need to check shifted snippet on t alignment because base is indexed immediately to the right
                exact_match_to_c = aln_read[alnmt_pos -2: alnmt_pos + u_v] == c_align_snippet 
                # print("LOL HERE T", aln_read[alnmt_pos -2: alnmt_pos + u_v], c_align_snippet)
            else:
                exact_match_to_c = aln_read[alnmt_pos -1: alnmt_pos + u_v + 1] == c_align_snippet
                # print( aln_read[alnmt_pos -1: alnmt_pos + u_v + 1], c_align_snippet)

            assert not (exact_match_to_c and exact_match_to_t) # they cannot both be perfect matches, that means a bug

            if  exact_match_to_t:
                read_coord = len("".join([n for n in aln_read[ : alnmt_pos +1] if n != "-"])) - 1
                ccs_coord = ccs_dict[ read_acc ].read_aln_to_ccs_coord(aln_read, read_coord)
            
            elif exact_match_to_c:
                print("HERE  C!!")
                if v_type == "D":
                    read_coord = len("".join([n for n in aln_read[ : alnmt_pos +1] if n != "-"]))  # position immediately to the right of deletion w.r.t. t
                elif v_type == "I":
                    # since alnmt_pos is immediately to the right of del in t, we need the base qual on the character before the deletion pos
                    # t: G-AACT, c: GAAACT, read = GAAACT
                    # almnt_pos is pointing to first A after the '-' in t, this has coord 1, in r this A corresponds to coord 2, but we want coord 1 (3-2 = 1)
                    read_coord = len("".join([n for n in aln_read[ : alnmt_pos +1] if n != "-"])) -2   
                else:
                    read_coord = len("".join([n for n in aln_read[ : alnmt_pos +1] if n != "-"])) - 1
                ccs_coord = ccs_dict[ read_acc ].read_aln_to_ccs_coord(aln_read, read_coord)

            else:
                reads_not_supporting_any_seq.add(read_acc)
                prob = -1 
                break

            q_qual = ccs_dict[read_acc].qual[ccs_coord]
            # print(ccs_dict[read_acc].qual[ccs_coord -1: ccs_coord +5 ], ccs_dict[read_acc].seq[ccs_coord-1: ccs_coord +5 ], c_align_snippet)

            # To map quality values
            # [A, B] --> [a, b]  [3, 93] --> [3, 43]
            # (x - A)*(b-a)/(B-A) + a
            q_qual_mapped = (q_qual - 3)*(max_phred_q_trusted - 3.0)/(90.0) + 3

            if u_v > 1: # probability in a homopolymenr region has all it's uncertainty attributed to the lenght of the homopolymer that 
                p_error =  (10**(-q_qual_mapped/10.0))
            elif v_type == "S":
                p_error =  ((10**(-q_qual_mapped/10.0))*subs_ratio)/3.0 # probability that its an identical substitution error from a base call uncertainty
            elif  v_type == "I":
                p_error =  ((10**(-q_qual_mapped/10.0))*ins_ratio)/4.0 # probability that its an identical insertion error from a base call uncertainty
            elif v_type == "D":
                p_error =  (10**(-q_qual_mapped/10.0)) * del_ratio # probability that its a delation error from a base call uncertainty
            else:
                print("Bad type", v_type, u_v)
                sys.exit()

            # print(p_error, q_qual_mapped, q_qual)
            prob *= p_error

        if prob >= 0:
            probabilities[read_acc] = prob
        else:
            print("read did not support c or t in variant region")
            pass

    print("Total probs:", len(probabilities), "non-informative:", len(reads_not_supporting_any_seq))
    return probabilities, reads_not_supporting_any_seq



def get_support_from_t(read_alignments_to_t, variant_coords_t):
    reads_support = []

    for read_acc in read_alignments_to_t:
        aln_t, aln_read, (matches, mismatches, indels) = read_alignments_to_t[read_acc]
        t_seq_to_coord_in_almnt = [j for j, nucl in enumerate(aln_t) if aln_t[j] != "-"]
        support = True
        for i in variant_coords_t:
            v_type, v_nucl, u_v = variant_coords_t[i] 
            alnmt_pos = t_seq_to_coord_in_almnt[i]


            exact_match_to_t = aln_read[alnmt_pos -1: alnmt_pos + u_v + 1] == aln_t[alnmt_pos -1: alnmt_pos + u_v +1]
            if exact_match_to_t: # matching t exactly means by definition that it doesn't match c
                support = False
                break

            elif v_type == "S":
                exact_3mer_match = aln_read[alnmt_pos - 1] == aln_t[alnmt_pos - 1] and aln_read[alnmt_pos +1] == aln_t[alnmt_pos + 1] and aln_read[alnmt_pos] == v_nucl
                if not exact_3mer_match: #aln_read[alnmt_pos] != v_nucl:
                    support = False

            elif v_type == "D":
                if aln_read[alnmt_pos] != "-": # read has deletion over the same pos as candidate
                    support = False

            else:
                assert v_type == "I"
                insertion_in_first_pos = aln_read[alnmt_pos -1] == v_nucl and  aln_t[alnmt_pos-1] == "-"
                # rare occasions insertion_in_last_pos
                # c to t: c:TTATTTTGGGGCTG, t: TTATTTT-GGGCTG
                # Exampleright shifted indel:
                # t: TTATTTTGGG--CTG
                # r: TTATTTTGGGGCCTG
                # obtain how long the homopolymenr insertion stetches in the alignment on t
                p = "[{variant}]+".format(variant=v_nucl)
                m_f = re.match(p, aln_t[alnmt_pos: ])
                u_v = len(m_f.group()) if m_f else 1 # is insertion in homopolymenr region of at least 2 bases in candidate, therefore we have at least one matching base here
                if len(aln_read) > alnmt_pos +u_v:
                    insertion_in_last_pos = aln_read[alnmt_pos +u_v] == v_nucl and aln_t[alnmt_pos + u_v] == "-" # can happen in homopolymenr
                else: # might be at the end of the t to read alignment, i.e., len(aln_read) == alnmt_pos +u_v, then there is not insertion in the read compared to t by definition
                    insertion_in_last_pos = False #aln_read[alnmt_pos +u_v -1] == v_nucl and aln_t[alnmt_pos + u_v -1] == "-" # somtimes happen in homopolymenr

                if not insertion_in_first_pos and not insertion_in_last_pos:
                # if aln_read[alnmt_pos -1] != v_nucl or aln_t[alnmt_pos] != "-": # read has to match on both position and position immediately before del to support (i.e., no indels in between) 
                    support = False


        if support:
            reads_support.append(read_acc)
    return reads_support


# def get_read_ccs_probabilities_from_reads_to_t(read_alignments_to_t, variant_coords_t, ccs_dict, errors, max_phred_q_trusted):
#     probabilities_t = {}
#     tot_errors = float(sum([ i+d+s for i, d, s in errors.values()]))
#     tot_errors = tot_errors if tot_errors > 1 else 1.0
#     subs_ratio = sum([ s for i, d, s in errors.values()]) / tot_errors
#     ins_ratio = sum([ i for i, d, s in errors.values()]) / tot_errors
#     del_ratio = sum([ d for i, d, s in errors.values()]) / tot_errors

#     for read_acc in read_alignments_to_t:
#         aln_t, aln_read, (matches, mismatches, indels) = read_alignments_to_t[read_acc]
#         t_seq_to_coord_in_almnt = [j for j, nucl in enumerate(aln_t) if aln_t[j] != "-"]
#         read_seq_to_coord_in_almnt = [j for j, nucl in enumerate(aln_read) if aln_read[j] != "-"]
#         prob = 1.0

#         for i in variant_coords_t:
#             v_type, v_nucl, u_v = variant_coords_t[i] 
#             alnmt_pos = t_seq_to_coord_in_almnt[i]

#             # require exact match over 3mer if S, or I, D in non-homopolymenr region.
#             # If homopolymenr region of size u_v, require length match of homopolymenr
#             exact_match = aln_read[alnmt_pos -1: alnmt_pos + u_v + 1] == aln_t[alnmt_pos -1: alnmt_pos + u_v +1]
#             if exact_match:
#                 read_coord = len("".join([n for n in aln_read[ : alnmt_pos +1] if n != "-"])) - 1
#                 ccs_coord = ccs_dict[ read_acc ].read_aln_to_ccs_coord(aln_read, read_coord)

#             else:
#                 read_snippet_up_to_variant = "".join([n for n in aln_read[ : alnmt_pos +1] if n != "-"])
#                 if aln_read[alnmt_pos] == "-":
#                     alnmt_pos_after_del = read_seq_to_coord_in_almnt[len(read_snippet_up_to_variant)]
#                     read_snippet_up_to_variant = "".join([n for n in aln_read[ : alnmt_pos_after_del +1] if n != "-"])
#                 # else:
#                 #     read_snippet_up_to_variant = "".join([n for n in aln_read[ : alnmt_pos +1] if n != "-"])
#                 p = "[{variant}]+".format(variant=read_snippet_up_to_variant[-1])
#                 m = re.match(p, read_snippet_up_to_variant[ : : -1 ])
#                 read_coord = len(read_snippet_up_to_variant) - len(m.group())
#                 ccs_coord = ccs_dict[ read_acc ].read_aln_to_ccs_coord(aln_read, read_coord)


#             q_qual = ccs_dict[read_acc].qual[ccs_coord]
#             # To map quality values
#             # [A, B] --> [a, b]  [3, 93] --> [3, 43]
#             # (x - A)*(b-a)/(B-A) + a
#             q_qual_mapped = (q_qual - 3)*(max_phred_q_trusted - 3.0)/(90.0) + 3

#             if u_v > 1: # probability in a homopolymenr region has all it's uncertainty attributed to the lenght of the homopolymer that 
#                 p_error =  (10**(-q_qual_mapped/10.0))
#             elif v_type == "S":
#                 p_error =  ((10**(-q_qual_mapped/10.0))*subs_ratio)/3.0 # probability that its an identical substitution error from a base call uncertainty
#             elif  v_type == "I":
#                 p_error =  ((10**(-q_qual_mapped/10.0))*ins_ratio)/4.0 # probability that its an identical insertion error from a base call uncertainty
#             elif v_type == "D":
#                 p_error =  (10**(-q_qual_mapped/10.0)) * del_ratio # probability that its a delation error from a base call uncertainty
#             else:
#                 p_error = -1 # If end up here, there is a bug

#             # print(p_error, q_qual_mapped, q_qual)
#             assert 0.0 < p_error < 1.0
#             prob *= p_error
#         probabilities_t[read_acc] = prob

#     return probabilities_t



def arrange_alignments_new_no_realign(t_acc, c_acc, t_seq, c_seq, read_alignments_to_c, read_alignments_to_t, ccs_dict, ignore_ends_len, max_phred_q_trusted):

    start_time = time()
    # pr = cProfile.Profile()
    # pr.enable()

    partition_dict = {t_acc : {}}
    partition_dict[t_acc][c_acc] = (t_seq, c_seq)
    # partition_dict[t_acc][t_acc] = (t_seq, t_seq)

    exact_edit_distances = edlib_align_sequences_keeping_accession(partition_dict, nr_cores = 1)    
    exact_alignments = sw_align_sequences_keeping_accession(exact_edit_distances, nr_cores = 1, ignore_ends_len = ignore_ends_len)
    print("saved re-aligning", len(read_alignments_to_t) + len(read_alignments_to_c), "sequences")

    # 1. Find positions differing between reference and candidate (ignoring any indel differences in ends)
    
    aln_t, aln_c, (matches, mismatches, indels) = exact_alignments[t_acc][c_acc]
    start, end = get_mask_start_and_end(aln_t, aln_c) # mask indels in ends due to legnth differences
    variants = [ (i,p_t,p_c) for i, (p_t, p_c) in  enumerate(zip(aln_t, aln_c)) if p_t != p_c and start <= i < end ]

    # 2. Get the coordinates on the candidate and reference respectively
    
    variant_coords_t, variant_coords_c, alignment_c_to_t, alignment_t_to_c = get_variant_coordinates(t_seq, c_seq, aln_t, aln_c, variants)
    print(variant_coords_t)

    # 3. Check if reads support candidate at given variant positions
    reads_support = get_support(read_alignments_to_c, variant_coords_c, read_alignments_to_t, variant_coords_t, alignment_c_to_t)

    reads_support_from_c = get_support_from_c(read_alignments_to_c, variant_coords_c)
    reads_support_from_t = get_support_from_t(read_alignments_to_t, variant_coords_t)
    
    # if set(reads_support) != set(reads_support_from_c) | set(reads_support_from_t):
    #     print("DIFF:", len(reads_support), len(reads_support_from_c), len(reads_support_from_t))
    #     sys.exit()

    # 4. get individual read error rates (again ignoring, any indel differences in ends) 
    errors = get_read_errors(read_alignments_to_c, read_alignments_to_t)

    # 5. Get position specific error rate for each variant in the reads 
    if ccs_dict:
        probability_c, non_supportive_c = get_read_ccs_probabilities_c(read_alignments_to_c, variant_coords_c, alignment_t_to_c, ccs_dict, errors, max_phred_q_trusted)
        probability_t, non_supportive_t = get_read_ccs_probabilities_t(read_alignments_to_t, variant_coords_t, alignment_c_to_t, ccs_dict, errors, max_phred_q_trusted)

        # probability_t = get_read_ccs_probabilities_from_reads_to_t(read_alignments_to_t, variant_coords_t, ccs_dict, errors, max_phred_q_trusted) 
        # probability_c = get_read_ccs_probabilities_from_reads_to_c(read_alignments_to_c, variant_coords_c, ccs_dict, errors, max_phred_q_trusted) 
        probability = merge_two_dicts(probability_c, probability_t)
    # else:
    #     probability_t = functions.get_empirical_error_probabilities(len(t_seq), errors, read_alignments_to_t, variant_coords_t) 
    #     probability_c = get_read_empirical_probabilities_from_reads_to_c(len(t_seq), errors, read_alignments_to_c, variant_coords_c) 
    #     probability = merge_two_dicts(probability_c, probability_t)

    variant_types = "".join([ str(variant_coords_c[j][0]) for j in  variant_coords_c ])
    if len(variants) == 0:
        print("{0} no difference to ref {1} after ignoring ends!".format(c_acc, t_acc))
        p_value = 0.0
    else:
        p_value = raghavan_upper_pvalue_bound(probability, reads_support)

    print("New p-value:", p_value, len(probability), len(reads_support), "non contributing:", len(non_supportive_c | non_supportive_t))

    total_elapsed = time() - start_time
    print("total new arrange:",total_elapsed)

    return reads_support_from_c, reads_support_from_t, variant_coords_t, p_value



def arrange_alignments(t_acc, reads_to_c, read_alignments_to_t, C, ignore_ends_len):
    partition_dict = {t_acc : {}}
    for read_acc in reads_to_c:
        partition_dict[t_acc][read_acc] = (C[t_acc], reads_to_c[read_acc])
    for c_acc in C:
            partition_dict[t_acc][c_acc] = (C[t_acc], C[c_acc])

    exact_edit_distances = edlib_align_sequences_keeping_accession(partition_dict, nr_cores = 1)    
    exact_alignments = sw_align_sequences_keeping_accession(exact_edit_distances, nr_cores = 1, ignore_ends_len = ignore_ends_len)
    partition_alignments = {} 

    assert len(exact_alignments) == 1
    for t_acc in exact_alignments:
        partition_alignments[t_acc] = {}
        for read_acc in exact_alignments[t_acc]:
            aln_t, aln_read, (matches, mismatches, indels) = exact_alignments[t_acc][read_acc]
            edit_dist = mismatches + indels
            partition_alignments[t_acc][read_acc] = (edit_dist, aln_t, aln_read, 1)

            x_aln_seq = "".join([n for n in aln_read if n != "-"])
            if read_acc in reads_to_c:
                assert reads_to_c[read_acc] == x_aln_seq

    for read_acc in read_alignments_to_t:
        aln_t, aln_read, (matches, mismatches, indels) = read_alignments_to_t[read_acc]
        edit_dist = mismatches + indels
        partition_alignments[t_acc][read_acc] = (edit_dist, aln_t, aln_read, 1)

    print("Re-aligned",len(reads_to_c), "sequences")
    print("Saved re-aligning",len(read_alignments_to_t), "sequences")
    alignment_matrix_to_t = functions.create_multialignment_matrix(C[t_acc], partition_alignments[t_acc])
    # PFM_to_t = functions.create_position_frequency_matrix(alignment_matrix_to_t, partition_alignments[t_acc])

    return alignment_matrix_to_t #, PFM_to_t




def statistical_test( c_acc, t_acc, c_seq, t_seq, reads_to_c, read_alignments_to_t, read_alignments_to_c, ignore_ends_len, ccs_dict, max_phred_q_trusted):
    reads = set(reads_to_c.keys()) | set(read_alignments_to_t.keys())  # dict(reads_to_c, **reads_to_t) #reads_to_c | reads_to_t
    shared_reads = set(reads_to_c.keys()) & set(read_alignments_to_t.keys()) 
    assert not shared_reads

    N_t = len(reads)

    # no reads supporting neither the candidate nor the reference t
    #  this can happen if after realignment of reads to candidates, all the reads 
    # to noth c and t got assigned to other candidates -- based on nearest_neighbor graph 
    # (should be rare) 
    if N_t == 0: 
        return c_acc, t_acc, 1.0, 1.0, 0, N_t, ""
    

    if ccs_dict:
        for x_acc in reads_to_c:
            assert reads_to_c[x_acc] == ccs_dict[x_acc].seq

    print()
    print("NEW NO REALIGN")
    x_new_from_c, x_new_from_t, delta_t_new, p_value_new = arrange_alignments_new_no_realign(t_acc, c_acc, t_seq, c_seq, read_alignments_to_c, read_alignments_to_t, ccs_dict, ignore_ends_len, max_phred_q_trusted)
    # if ignore_ends_len > 0:
    #     alignment_matrix_to_t = functions.cut_ends_of_alignment_matrix(alignment_matrix_to_t, t_acc, c_acc, ignore_ends_len)

    # get parameter estimates for statistical test
    # delta_t_new = functions.get_difference_coordinates_for_candidates(t_acc, c_acc, alignment_matrix_to_t) # format: { c_acc1 : {pos:(state, char), pos2:(state, char) } , c_acc2 : {pos:(state, char), pos2:(state, char) },... }
    # get number of reads k supporting the given set of variants, they have to support all the variants within a candidate
    # x = functions.reads_supporting_candidate(t_acc, c_acc, alignment_matrix_to_t, delta_t_new, reads) # format: { c_acc1 : [x_acc1, x_acc2,.....], c_acc2 : [x_acc1, x_acc2,.....] ,... }
    # print(len(x_new), delta_t_new, t_acc, c_acc)
    print("new support", len(x_new_from_c), "and", len(x_new_from_t) )
    print()

    # print("NEW")
    # alignment_matrix_to_t = arrange_alignments_new(t_acc, c_acc, t_seq, c_seq, reads_to_c, read_alignments_to_t, ignore_ends_len)
    # if ignore_ends_len > 0:
    #     alignment_matrix_to_t = functions.cut_ends_of_alignment_matrix(alignment_matrix_to_t, t_acc, c_acc, ignore_ends_len)

    # # get parameter estimates for statistical test
    # delta_t_new = functions.get_difference_coordinates_for_candidates(t_acc, c_acc, alignment_matrix_to_t) # format: { c_acc1 : {pos:(state, char), pos2:(state, char) } , c_acc2 : {pos:(state, char), pos2:(state, char) },... }
    # # get number of reads k supporting the given set of variants, they have to support all the variants within a candidate
    # x = functions.reads_supporting_candidate(t_acc, c_acc, alignment_matrix_to_t, delta_t_new, reads) # format: { c_acc1 : [x_acc1, x_acc2,.....], c_acc2 : [x_acc1, x_acc2,.....] ,... }
    # print(len(x), delta_t_new, t_acc, c_acc)
    # print()


    # get multialignment matrix here
    print("OLD")
    alignment_matrix_to_t =  arrange_alignments(t_acc, reads_to_c, read_alignments_to_t, {c_acc: c_seq,  t_acc: t_seq}, ignore_ends_len)
    # cut multialignment matrix first and last ignore_ends_len bases in ends of reference in the amignment matrix
    # these are bases that we disregard when testing varinats
    # We get individual cut positions depending on which candidate is being tested -- we dont want to include ends spanning over the reference or candidate
    # we cut at the start position in c or t that comes last, and the end position in c or t that comes first
    if ignore_ends_len > 0:
        alignment_matrix_to_t = functions.cut_ends_of_alignment_matrix(alignment_matrix_to_t, t_acc, c_acc, ignore_ends_len)

    # get parameter estimates for statistical test
    delta_t = functions.get_difference_coordinates_for_candidates(t_acc, c_acc, alignment_matrix_to_t) # format: { c_acc1 : {pos:(state, char), pos2:(state, char) } , c_acc2 : {pos:(state, char), pos2:(state, char) },... }
    # get number of reads k supporting the given set of variants, they have to support all the variants within a candidate
    x = functions.reads_supporting_candidate(t_acc, c_acc, alignment_matrix_to_t, delta_t, reads) # format: { c_acc1 : [x_acc1, x_acc2,.....], c_acc2 : [x_acc1, x_acc2,.....] ,... }
    invariant_factors_for_candidate = functions.get_invariant_multipliers(delta_t, alignment_matrix_to_t, t_acc)

    ############ TMP ###################
    ####################################
    ####################################
    # print(len(x), delta_t, t_acc, c_acc)
    print( (set(x_new_from_c) | set(x_new_from_t) )   ^ set(x) )
    if len(x_new_from_c) + len(x_new_from_t) != len(x):
        print()
        print("DIFFERENCE:", "new support:", len(x_new_from_c) + len(x_new_from_t), "old support:", len(x) )
        print(delta_t_new)
        # for read_acc in read_alignments_to_c:
        #     print("C")
        #     print(read_alignments_to_c[read_acc][0])
        #     print(read_alignments_to_c[read_acc][1])
        # for read_acc in read_alignments_to_t:
        #     print("T")
        #     print(read_alignments_to_t[read_acc][0])
        #     print(read_alignments_to_t[read_acc][1])

        print()
    assert len(list(delta_t.values())[0]) ==  len(delta_t_new)
    for key in delta_t:
        if len(list(delta_t[key].values())) != len(list(delta_t_new.values())) or set([ (typ, nuc) for typ, nuc, u_mult in delta_t_new.values()]) != set(list(delta_t[key].values())) :
            print(delta_t_new)
            print(delta_t)
        assert len(list(delta_t[key].values())) == len(list(delta_t_new.values())) and set([ (typ, nuc) for typ, nuc, u_mult in delta_t_new.values()]) == set(list(delta_t[key].values())) 
    new_invarinats = set([ u_mult for typ, nuc, u_mult in delta_t_new.values()])
    old_invarinats = set([tup_dict.values()[0] for cand, pos_dict in invariant_factors_for_candidate.items() for tup_dict in pos_dict.values() ])
    assert new_invarinats == old_invarinats

    ####################################
    ####################################
    # print(invariant_factors_for_candidate)


    if ccs_dict:
        insertions, deletions, substitutions = functions.get_errors_for_partitions(t_acc, len(t_seq), c_acc, alignment_matrix_to_t) 
        probability = functions.get_ccs_position_prob_per_read(t_acc, len(t_seq), alignment_matrix_to_t, invariant_factors_for_candidate, c_acc, delta_t, ccs_dict, insertions, deletions, substitutions, max_phred_q_trusted) 
    else:
        errors = functions.get_errors_per_read(t_acc, len(t_seq), c_acc, alignment_matrix_to_t) 
        # invariant_factors_for_candidate = functions.get_invariant_multipliers(delta_t, alignment_matrix_to_t, t_acc)
        probability = functions.get_prob_of_error_per_read(t_acc, len(t_seq), c_acc, errors, invariant_factors_for_candidate) 

    #TODO: do exact only if partition less than, say 200? Otherwise poisson approx
    # p_value = CLT_test(probability, weight, x)
    # p_value = poisson_approx_test(probability, weight, x)
    # p_value = exact_test(probability, weight, x)
    # print("exact p:", p_value )
    # print()
    # print(sorted(errors.values()))
    # print(sorted(probability.values()))
    # print("Weighted raghavan p:", p_value )

    delta_size = len(delta_t[c_acc])
    variant_types = "".join([ str(delta_t[c_acc][j][0]) for j in  delta_t[c_acc] ])
    if delta_size == 0:
        print("{0} no difference to ref {1} after ignoring ends!".format(c_acc, t_acc))
        p_value = 0.0
    else:
        p_value = raghavan_upper_pvalue_bound(probability, x)

    print("Old p-value:", p_value )
    if p_value_new > p_value:
        print("NEW P-VAL GREATER:", p_value_new, p_value )
    else:
        print("NEW P-VAL SMALLER:", p_value_new, p_value )

    # print("Tested", c_acc, "to ref", t_acc, "p_val:{0}, mult_factor:{1}, corrected p_val:{2} k:{3}, N_t:{4}, Delta_size:{5}".format(p_value, correction_factor, p_value * correction_factor,  len(x), N_t, delta_size) )
    # significance_values[c_acc] = (p_value, correction_factor, len(x), N_t, delta_size)
    if ccs_dict:
        # print("Tested", c_acc, "to ref", t_acc, "p_val:{0}, k:{1}, N_t:{2}, variants:{3}".format(p_value, len(x), N_t, variant_types) )
        return (c_acc, t_acc, p_value, 1.0, len(x), N_t, variant_types)
    else:
        correction_factor = calc_correction_factor(t_seq, c_acc, delta_t)
        # print("Tested", c_acc, "to ref", t_acc, "p_val:{0}, k:{1}, N_t:{2}, variants:{3}".format(p_value, len(x), N_t, variant_types) )
        return (c_acc, t_acc, p_value, correction_factor, len(x), N_t, variant_types)






def raghavan_upper_pvalue_bound(probability, x_equal_to_one):
    """ 
        Method for bounding the p-value from above based on:
        https://math.stackexchange.com/questions/1546366/distribution-of-weighted-sum-of-bernoulli-rvs
        With this paper: [1] https://www.cc.gatech.edu/~mihail/Rag88.pdf

        1. assign all weights  0 < a_i <= 1 as w_i = min(p_i) / p_i
        smallest probabilities will get weight 1.0. 
        2. Define Y = sum_i x_i*w_i calculate mean E[Y] = m = sum_i p_i*w_i
        3. Use Theorem 1 in [1] that states 
            (*) P(Y > m(1+d)) < (e^d / (1+d)^(1+d) )^m
            for d > 0, m > 0.

            (i) Get relevant d (given a value y) for us by letting d := (y/m) - 1
            (ii) Compute (e^d / (1+d)^(1+d) )^m. If d (say > 10) large and m small, avoid big and small numbers by
                (A) Use the fact that m = k/d, for some k.
                (B) Rewrite RHS in (*) as e^(d*k/d) / (1+d)^(k(1+d)/d)  =  e^k / (1+d)^(k + k/d)
        4. p_value is now bounded above (no greater than) the computed value.
    """

    # for p_i in probability.values():
    #     if p_i < 0 or p_i > 1.0:
    #         print(p_i)
    # assert max(probability.values()) <= 1.0
    # print(sorted(probability.values()))
    log_probabilities = { acc: -math.log(p_i, 10) for acc, p_i in probability.items()}
    log_p_i_max = max(log_probabilities.values())
    
    # print(log_probabilities.values())
    # print("log_p_i_max:", log_p_i_max)

    assert log_p_i_max > 0
    weight = {q_acc : log_probabilities[q_acc] / log_p_i_max  for q_acc in log_probabilities.keys()}

    # p_i_min = min(probability.values())
    # weight = {q_acc : p_i_min / probability[q_acc]  for q_acc in probability.keys()}

    m = Decimal( sum([ weight[q_acc] * probability[q_acc]  for q_acc in probability.keys()]) )
    y = Decimal( sum([weight[x_i] for x_i in x_equal_to_one ]) )
    # print(m, y, "log_p_max: ",log_p_i_max, "nr supp:", len(x_equal_to_one), sorted([(weight[x_i], x_i) for x_i in x_equal_to_one ], key = lambda x: x[0]) )
    d = y / m - 1
    k = m*d

    # if d > 10:
    if y == 0:
        raghavan_bound = 1.0
    elif d == 0:
        raghavan_bound = 0.5
    else:
        try:
            raghavan_bound = k.exp() / (d+1)**(k + k/d)
        except:
            print("Decimal computation error:")
            print("Values: m:{0}, d:{1}, y:{2}, k :{3}".format(m, d, y, k) )

        # else:
        #     raghavan_bound = (d.exp() / (d+1)**(d+1))**m
    # print("Values: m:{0}, d:{1}, y:{2}, k :{3}".format(m, d, y, k) )

    #### Our bound #####
    # will give tighter bound if all probabilities are tiny (delta is large)
    # 1. set p = max(p_i)
    # 2. Calculate probability that _any_ bernoilly event happens, i.e., Y > 0, 
    #     as P = 1 - (1 - p)^(N_t), This is less significant than k bernoilli sucesses where all p_i <= p.
    # therefore P is an upper bound on the p-value
    # p_ = Decimal(max(probability.values()))
    # our_bound = Decimal(1.0) - (Decimal(1.0) - p_ )**len(probability)
    # print("Our:{0}, Raghavan: {1}".format(our_bound, raghavan_bound))
    ##### temporary check if this is a better bound #######
    # from Section "other properties" in "https://en.wikipedia.org/wiki/Moment-generating_function
    # M_t = sum([ Decimal(1) - Decimal(p_i) + Decimal(p_i)*Decimal(y).exp() for p_i in probability.values() ])
    # wiki_bound = Decimal(1).exp()*Decimal(-y**2).exp()*y*M_t
    # # print("wiki_bound:{0}, M_t:{1}".format(wiki_bound, M_t))
    # print(round(wiki_bound, 50), round(raghavan_bound,50))
    ########################################################

    # print("m:{0}, d:{1}, y:{2}, k :{3}, p_bound={4}".format(round(m,20) , round(d,20),round(y,20), round(k,20), round(p_value_upper_bound,20) ) )
    # return min(float(our_bound), float(raghavan_bound))
    return float(raghavan_bound)



def calc_correction_factor(t_seq, c_acc, delta_t):
    m = len(t_seq)
    n_S, n_D, n_I = 0, 0, 0
    for pos, (state, char) in delta_t[c_acc].items():
        if state == "S":
            n_S += 1
        elif state == "D":
            n_D += 1
        if state == "I":
            n_I += 1
    
    correction_factor = ( (4*(m+1))**n_I ) * functions.choose(m, n_D) * functions.choose( 3*(m-n_D), n_S)
    return correction_factor
