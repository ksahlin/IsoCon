"""
    position_query_to_alignment(): arranges all query alignments onto a target sequence
    get_non_overlapping_intervals(): finds canonical intervals where a set of reads maps, 
            A canonic interval assures the same muber of reads are aligned to the target in that interval.
    create_multialignment_format(): arranges all query alignments that covers the region [start,stop] into a multialignment matrix
                    useful for statistical testing of SNVs

    FUTURE:
    identify_allele_regions(): takes a start and stop coordinate and returns positions subject to phasing (two variants have high coverage) 
    this is useful for phasing haplotypes or SNVs within viral sequences or haplotypes

"""

import unittest
from collections import defaultdict
import math
import re
import edlib

from itertools import permutations


def filter_exon_differences(pairwise_alignments, min_exon_diff, ignore_ends_len):
    pattern = r"[-]{{{min_exon_diff},}}".format( min_exon_diff = str(min_exon_diff)  )  # r"[-]{20,}"
    filtered = set()
    for s1 in list(pairwise_alignments.keys()): 
        for s2 in list(pairwise_alignments[s1].keys()):
            s1_alignment, s2_alignment, (matches, mismatches, indels) = pairwise_alignments[s1][s2]
            
            start, end = get_mask_start_and_end(s1_alignment, s2_alignment) # mask indels in ends due to legnth differences
            start = min(ignore_ends_len, start)
            end = max(len(s1_alignment) - ignore_ends_len, end )
            missing_exon_s1 = re.search(pattern, s1_alignment[start:end])
            missing_exon_s2 = re.search(pattern, s2_alignment[start:end])
            if missing_exon_s1 or missing_exon_s2:
                # print(missing_exon_s1.group(0))
                # print(s1_alignment)
                # print(s2_alignment)
                # print()
                # print(len(pairwise_alignments[s1].keys()))
                del pairwise_alignments[s1][s2]
                filtered.add(s2)
            # elif missing_exon_s2:
            #     # print(missing_exon_s2.group(0))
            #     # print(s1)
            #     # print(s2)
            #     # print(len(pairwise_alignments[s1].keys()))
            #     del pairwise_alignments[s1][s2]
            #     filtered.add(s2)
    return filtered


def transform(read):
    transformed_seq = []
    prev_nucl = ""
    for nucl in read:
        if nucl != prev_nucl:
            transformed_seq.append(nucl)
        prev_nucl = nucl

    return "".join(transformed_seq)

def get_homopolymer_invariants(candidate_transcripts):
    seq_to_acc = { seq : acc for (acc, seq) in  candidate_transcripts.items() }
    print("Unique before compression: ", len(seq_to_acc) )

    candidate_transcripts_transformed = {}
    clusters = defaultdict(list)
    for acc in candidate_transcripts:
        seq_transformed = transform(candidate_transcripts[acc])
        candidate_transcripts_transformed[acc] = seq_transformed
        clusters[seq_transformed].append(acc)

    seq_to_acc_transformed = { seq : acc for (acc, seq) in candidate_transcripts_transformed.items()}
    print("Unique after compression: ", len(seq_to_acc_transformed) )

    edges = {}
    for seq in clusters:
        if len(clusters[seq]) > 1:
            # print(clusters[seq])
            for acc in clusters[seq]:
                edges[acc] = {}
            for acc1, acc2 in permutations(clusters[seq], 2):
                edges[acc1][acc2] = 1

    return edges


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
            alignment_c_to_t[len(t_seq_piece)-1] = aln_c[max(0, i-1): i+ u_v + 1]
            alignment_t_to_c[len(c_seq_piece)-1 +1] = aln_t[max(0, i-1): i+ u_v + 1]

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
            alignment_c_to_t[len(t_seq_piece)-1 +1] = aln_c[max(0, i-1) : i+ u_v + 1]
            alignment_t_to_c[len(c_seq_piece)-1] = aln_t[max(0, i-1) : i+ u_v + 1]

        else: # get coordinate on ref and cand if substitution
            variant_coords_t[len(t_seq_piece)-1] = ("S", p_c, 1 )
            variant_coords_c[len(c_seq_piece)-1] = ("S", p_c, 1 )
            alignment_c_to_t[len(t_seq_piece)-1] = aln_c[max(0, i-1) : i + 2]
            alignment_t_to_c[len(c_seq_piece)-1] = aln_t[max(0, i-1) : i + 2]

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
            exact_match = aln_read[max(0, alnmt_pos -1) : alnmt_pos + u_v + 1] == aln_c[max(0, alnmt_pos -1): alnmt_pos + u_v +1]

            if read_acc == 'm151210_031012_42146_c100926392550000001823199905121697_s1_p0/120400/29_1460_CCS_strand=+;fiveseen=1;polyAseen=0;threeseen=1;fiveend=29;polyAend=-1;threeend=1460;primer=2;chimera=0':
                print("in get_support, c")
                print("exact_match", exact_match, "variant_coords_c", v_type, v_nucl, u_v, aln_read[max(0, alnmt_pos -1): alnmt_pos + u_v + 1], aln_c[max(0, alnmt_pos -1): alnmt_pos + u_v +1])
                print(variant_coords_c)

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
                exact_match_to_c = aln_read[max(0, alnmt_pos -2): alnmt_pos + u_v] == c_align_snippet 
                # print(aln_read[alnmt_pos -2: alnmt_pos + u_v], c_align_snippet)
            else:
                exact_match_to_c = aln_read[max(0, alnmt_pos -1): alnmt_pos + u_v + 1] == c_align_snippet
                # print( aln_read[alnmt_pos -1: alnmt_pos + u_v + 1], c_align_snippet)

            # require exact match over 3mer if S, or I,D in non-homopolymenr region.
            # If homopolymenr region of size u_v, require length match of homopolymenr
            if not exact_match_to_c:
                support = False
                break
        if support:
            reads_support_c_from_t.append(read_acc)   

    # print("NEW2 support", len(reads_support_c), "and", len(reads_support_c_from_t) )
    return reads_support_c + reads_support_c_from_t


def get_read_errors(read_alignments_to_c, read_alignments_to_t):
    errors = {}
    for read_acc in read_alignments_to_t:
        aln_t, aln_read, (matches, mismatches, indels) = read_alignments_to_t[read_acc]
        insertions, deletions, substitutions = read_errors_from_alignment(aln_t, aln_read)
        errors[read_acc] = (insertions, deletions, substitutions)

    for read_acc in read_alignments_to_c:
        aln_c, aln_read, (matches, mismatches, indels) = read_alignments_to_c[read_acc]
        insertions, deletions, substitutions = read_errors_from_alignment(aln_c, aln_read)
        errors[read_acc] = (insertions, deletions, substitutions)       
    return errors


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



def get_read_ccs_probabilities_c(read_alignments_to_c, variant_coords_c, alignment_t_to_c, ccs_dict, errors, max_phred_q_trusted):
    subs =  float(max(1.0, sum([ s for i, d, s in errors.values()])))
    ins  =  float(max(1.0, sum([ i for i, d, s in errors.values()])))
    del_ =  float(max(1.0, sum([ d for i, d, s in errors.values()])))
    tot_errors = subs + ins + del_
    subs_ratio = subs / tot_errors
    ins_ratio =  ins / tot_errors
    del_ratio =  del_ / tot_errors
    assert len(variant_coords_c) > 0
    # tot_errors = float(sum([ i+d+s for i, d, s in errors.values()]))
    # tot_errors = float(tot_errors) if tot_errors > 1 else 1.0
    # subs_ratio = sum([ s for i, d, s in errors.values()]) / tot_errors
    # ins_ratio = sum([ i for i, d, s in errors.values()]) / tot_errors
    # del_ratio = sum([ d for i, d, s in errors.values()]) / tot_errors

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
            exact_match_to_c = aln_read[max(0, alnmt_pos -1) : alnmt_pos + u_v + 1] == aln_c[max(0,alnmt_pos -1): alnmt_pos + u_v +1]

            t_align_snippet = alignment_t_to_c[i]
            if v_type == "D": # will need to check shifted snippet on c alignment because base is indexed immediately to the right
                exact_match_to_t = aln_read[max(0, alnmt_pos -2): alnmt_pos + u_v] == t_align_snippet 
                # print(i, v_type, aln_read[alnmt_pos -2: alnmt_pos + u_v], t_align_snippet, aln_c[alnmt_pos -1: alnmt_pos + u_v +1], tot_errors, subs_ratio, ins_ratio, del_ratio, read_acc)
            else:
                exact_match_to_t = aln_read[max(0, alnmt_pos -1): alnmt_pos + u_v + 1] == t_align_snippet
                # print(i, v_type, aln_read[alnmt_pos -1: alnmt_pos + u_v + 1], t_align_snippet, aln_c[alnmt_pos -1: alnmt_pos + u_v +1], tot_errors, subs_ratio, ins_ratio, del_ratio, read_acc)

            assert not (exact_match_to_c and exact_match_to_t) # they cannot both be perfect matches, that means a bug

            if  exact_match_to_c:
                read_coord = len("".join([n for n in aln_read[ : alnmt_pos +1] if n != "-"])) - 1
                ccs_coord = ccs_dict[ read_acc ].read_aln_to_ccs_coord(aln_read, read_coord)
            
            elif exact_match_to_t:
                # print("HERE T!!")
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
            assert 0.0 < prob < 1.0
            probabilities[read_acc] = prob
        else:
            assert read_acc in reads_not_supporting_any_seq
            # print("read did not support c or t in variant region")
            pass

    # print("Total probs:", len(probabilities), "non-informative:", len(reads_not_supporting_any_seq))
    return probabilities, reads_not_supporting_any_seq


def get_read_ccs_probabilities_t(read_alignments_to_t, variant_coords_t, alignment_c_to_t, ccs_dict, errors, max_phred_q_trusted):
    subs =  float(max(1.0, sum([ s for i, d, s in errors.values()])))
    ins  =  float(max(1.0, sum([ i for i, d, s in errors.values()])))
    del_ =  float(max(1.0, sum([ d for i, d, s in errors.values()])))
    tot_errors = subs + ins + del_
    subs_ratio = subs / tot_errors
    ins_ratio =  ins / tot_errors
    del_ratio =  del_ / tot_errors
    assert len(variant_coords_t) > 0

    # tot_errors = float(sum([ i+d+s for i, d, s in errors.values()]))
    # tot_errors = float(tot_errors) if tot_errors > 3 else 3.0
    # subs_ratio = max(1.0, sum([ s for i, d, s in errors.values()])) / tot_errors
    # ins_ratio =  max(1.0, sum([ i for i, d, s in errors.values()])) / tot_errors
    # del_ratio =  max(1.0, sum([ d for i, d, s in errors.values()])) / tot_errors

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
            exact_match_to_t = aln_read[max(0, alnmt_pos -1): alnmt_pos + u_v + 1] == aln_t[max(0, alnmt_pos -1): alnmt_pos + u_v +1]

            c_align_snippet = alignment_c_to_t[i]
            if v_type == "I": # will need to check shifted snippet on t alignment because base is indexed immediately to the right
                exact_match_to_c = aln_read[max(0, alnmt_pos -2): alnmt_pos + u_v] == c_align_snippet 
                # print("LOL HERE T", aln_read[max(0, alnmt_pos -2): alnmt_pos + u_v], c_align_snippet)
            else:
                exact_match_to_c = aln_read[max(0, alnmt_pos -1): alnmt_pos + u_v + 1] == c_align_snippet
                # print( aln_read[max(0, alnmt_pos -1): alnmt_pos + u_v + 1], c_align_snippet)


            assert not (exact_match_to_c and exact_match_to_t) # they cannot both be perfect matches, that means a bug

            if  exact_match_to_t:
                read_coord = len("".join([n for n in aln_read[ : alnmt_pos +1] if n != "-"])) - 1
                ccs_coord = ccs_dict[ read_acc ].read_aln_to_ccs_coord(aln_read, read_coord)
            
            elif exact_match_to_c:
                # print("HERE  C!!")
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
            assert 0.0 < prob < 1.0
            probabilities[read_acc] = prob
        else:
            assert read_acc in reads_not_supporting_any_seq
            # print("read did not support c or t in variant region")
            pass

    # print("Total probs:", len(probabilities), "non-informative:", len(reads_not_supporting_any_seq))
    return probabilities, reads_not_supporting_any_seq


def get_empirical_error_probabilities(segment_length, errors, variant_coords_t):
    probability = {}
    delta_size = float(len(variant_coords_t))
    if delta_size == 0.0:
        print(variant_coords_t, len(errors), segment_length)
        print("No variant between candidates")
    assert delta_size > 0.0
    for read_acc in errors:
        prob = 1.0
        (insertions, deletions, substitutions) = errors[read_acc]
        p_S = ( max(substitutions, delta_size) / float(segment_length) ) / 3.0   # p = 0.0 not allowed, min_p is 1/(3*len(seq))
        p_I = ( max(insertions, delta_size) / float(segment_length) ) / 4.0   # p = 0.0 not allowed, min_p is 1/(4*len(seq))
        p_D = ( max(deletions, delta_size) / float(segment_length) )         # p = 0.0 not allowed, min_p is 1/(len(seq))

        for i in variant_coords_t:
            v_type, v_nucl, u_v = variant_coords_t[i]

            if v_type == "S":
                prob *= p_S*u_v 
            elif v_type == "I":
                prob *= min(0.5, p_I*u_v) 
            elif v_type == "D":
                prob *= min(0.5, p_D*u_v)
        if prob >= 1.0:
            prob = 0.99999
        
        elif prob <= 0.0:
            print("NONPOSITIVE PVAL:",prob, insertions, deletions, substitutions, p_S, p_I, p_D, variant_coords_t)

        probability[read_acc] = prob

    return probability



def reverse_complement(string):
    #rev_nuc = {'A':'T', 'C':'G', 'G':'C', 'T':'A', 'N':'N', 'X':'X'}
    # Modified for Abyss output
    rev_nuc = {'A':'T', 'C':'G', 'G':'C', 'T':'A', 'a':'t', 'c':'g', 'g':'c', 't':'a', 'N':'N', 'X':'X', 'n':'n', 'Y':'R', 'R':'Y', 'K':'M', 'M':'K', 'S':'S', 'W':'W', 'B':'V', 'V':'B', 'H':'D', 'D':'H', 'y':'r', 'r':'y', 'k':'m', 'm':'k', 's':'s', 'w':'w', 'b':'v', 'v':'b', 'h':'d', 'd':'h'}

    rev_comp = ''.join([rev_nuc[nucl] for nucl in reversed(string)])
    return(rev_comp)


def choose(n, k):
    """
    A fast way to calculate binomial coefficients by Andrew Dalke (contrib).
    """
    if 0 <= k <= n:
        ntok = 1
        ktok = 1
        for t in range(1, min(k, n - k) + 1):
            ntok *= n
            ktok *= t
            n -= 1
        return ntok // ktok
    else:
        return 0


def read_errors_from_alignment(ref_aln, read_aln):
    p = "[-]+"
    ref_f = re.match(p, ref_aln)
    read_f = re.match(p, read_aln)
    if ref_f and read_f:
        start = max(len(ref_f.group()), len(read_f.group()))
    elif ref_f:
        start = len(ref_f.group())  
    elif read_f:
        start = len(read_f.group())
    else:
        start = 0
    
    ref_r = re.match(p, ref_aln[::-1])
    read_r = re.match(p, read_aln[::-1])
    if ref_r and read_r:
        match_length = max(len(ref_r.group()), len(read_r.group()))
    elif ref_r:
        match_length = len(ref_r.group())  
    elif read_r:
        match_length = len(read_r.group())
    else:
        match_length = 0
    stop = len(ref_aln) - match_length

    variants = ["I" if n1 == "-" else "D" if n2 == "-" else "S" for n1, n2 in zip(ref_aln[start: stop], read_aln[start: stop]) if n1 != n2 ]
    insertions, deletions, substitutions = variants.count("I"), variants.count("D"), variants.count("S")
    return insertions, deletions, substitutions



def create_position_frequency_matrix(alignment_matrix, partition):
    nr_columns = len( alignment_matrix[ list(alignment_matrix)[0] ]) # just pick a key
    PFM = [{"A": 0, "C": 0, "G": 0, "T": 0, "-": 0} for j in range(nr_columns)]
    for s in alignment_matrix:
        s_aln = alignment_matrix[s]
        indegree = partition[s][3]
        for j in range(nr_columns): # for each column
            nucl = s_aln[j]
            PFM[j][nucl] += indegree

    return PFM


# from time import time
# import io
# import cProfile, pstats #, StringIO

def create_multialignment_matrix(m, partition):
    """
        a partition is a dictionary of pairwise alignments for a given center m. "partition has the following
        structure:  partition = {s : (edit_distance, m_alignment, s_alignment, degree_of_s)}
        s can only have a degree larger than 1 if s=m, otherwise s has a degree of 1.

        This function does the following

        query_to_target_positioned_dict = {}
        for each seq in partition we call function
            position_query_to_alignment(query_aligned, target_aligned, target_start=0)
            returns the following data
            query_to_target_positioned, target_vector_start_position = 0, target_vector_end_position 
            we store all these pairwise alignments in query_to_target_positioned_dict

        where this function retuns a dictionary with query seq as key in the following form:
        query_to_target_positioned_dict[query_accession] = (query_to_target_positioned, target_vector_start_position, target_vector_end_position)

        then it calls create_multialignment_format(query_to_target_positioned_dict, start, stop)
        This function returns an alignment matric in the following smaple format:
        alignment_matrix = {"q1" : ["-", "A","-", "C", "-", "G","A","C","C","G", "G", "-", "A", "T","T","T"],
                            "q2" : ["-", "A","-", "C", "-", "G","A","G","-","-", "G", "-", "A", "T","T","T"],
                            "q3" : ["-", "A","-", "C", "-", "G","A","-","-","-", "G", "-", "A", "T","T","T"],
                            "q4" : ["-", "A","-", "C", "-", "G","C","C","-","-", "G", "-", "A", "-","-","-"],
                            "q5" : ["-", "A","-", "C", "-", "G","-","-","-","-", "G", "-", "A", "T","-","-"],
                            "q6" : ["G", "A","-", "C", "-", "G","C","-","-","-", "G", "-", "A", "-","-","-"]
                            }

        finally, we transform the alignment_matrix into a PFM by multiplying each query sequnece with the correct degree.

        PFM is a list of dicts, where the inner dicts are one per column in the PFM matrix, its keys by characters A, C, G, T, -
        alignment_matrix is a representation of all alignments in a partition. this is a dictionary where sequences s_i belonging to the 
        partition as keys and the alignment of s_i with respect to the alignment matix.
    """
    query_to_target_positioned_dict = {}
    for q_acc in partition:
        (edit_distance, m_alignment, s_alignment, degree_of_s) = partition[q_acc]
        # s_positioned_new, target_vector_start_position, target_vector_end_position = position_query_to_alignment_NEW(s_alignment, m_alignment, 0)
        s_positioned, target_vector_start_position, target_vector_end_position = position_query_to_alignment(s_alignment, m_alignment, 0)
        # assert s_positioned_new == s_positioned
        assert target_vector_start_position == 0
        assert target_vector_end_position + 1 == 2*len(m) + 1 # vector positions are 0-indexed
        query_to_target_positioned_dict[q_acc] = (s_positioned, target_vector_start_position, target_vector_end_position)

    alignment_matrix = create_multialignment_format_NEW(query_to_target_positioned_dict, 0, 2*len(m))
    return alignment_matrix


def transpose(dct):
    d = defaultdict(dict)
    for key1, inner in dct.items():
        for key2, value in inner.items():
            d[key2][key1] = value
    return d

def position_query_to_alignment(query_aligned, target_aligned, target_alignment_start_position):
    """
        input:      0-indexed target positions
        returns:    target vector is a list of 2*t_len +1 positions to keep insertions between base pairs    
                list of strings of nucleotides for each position, start position in target vector list, end position in target vector list
    """

    query_positioned = [] 
    target_position = target_alignment_start_position # 0-indexed
    temp_ins = ""
    # iterating over alignment positions
    for p in range(len(target_aligned)):
        if target_aligned[p] == "-":
            temp_ins += query_aligned[p]
        else:
            if not temp_ins:
                query_positioned.append("-") #[2*target_position] = "-"
            else:
                query_positioned.append(temp_ins)
                temp_ins = ""

            query_positioned.append(query_aligned[p])

            target_position += 1

    if not temp_ins:
        query_positioned.append("-")
    else:
        query_positioned.append(temp_ins)

    target_vector_start_position = 2*target_alignment_start_position
    target_vector_end_position = 2*(target_position-1) + 2

    return query_positioned, target_vector_start_position, target_vector_end_position



def get_best_solution(max_insertion, q_ins):

    if q_ins == "-":
        return ["-" for i in range(len(max_insertion))]

    else:
        pos = max_insertion.find(q_ins) 
        if pos >= 0:    
            q_insertion_modified = "-"*pos + max_insertion[ pos : pos + len(q_ins) ] + "-"* len(max_insertion[ pos + len(q_ins) : ])
            return [character for character in q_insertion_modified]        

        # else, check if smaller deletion can be aligned from left to write, e.g. say max deletion is GACG
        # then an insertion AG may be aligned as -A-G. Take this alignment instead
        q_insertion_modified = min_ed(max_insertion, q_ins)
        if q_insertion_modified:
            return [character for character in q_insertion_modified]

        # otherwise just shift left
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
            return [character for character in q_insertion_modified]


        q_insertion_modified = []
        for p in range(len(max_insertion)):
            # all shorter insertions are left shifted -- identical indels are guaranteed to be aligned
            # however, no multialignment is performed among the indels
            if p < len(q_ins):
                q_insertion_modified.append(q_ins[p])
            else:
                q_insertion_modified.append("-")
        return [character for character in q_insertion_modified]


def create_multialignment_format_NEW(query_to_target_positioned_dict, start, stop):
    """
            1. From segments datastructure, get a list (coordinate in MAM) of lists insertions that contains all the insertions in that position
            2. Make data structure unique_indels =  set(insertions) to get all unique insertions in a given position
            3. Make a list max_insertions containing the lengths of each max_indel in (max_indels >1bp should be padded). From this list we can get the total_matrix_length
            4. In a data structure position_solutions : {max_indel : {q_ins : solution_list, q_ins2: }, }, 
                find oplimal solution for each unique indel in unique_indels to max_indel in datastructure  (a list of dicts in each position) 
            5. for each pos in range(total_matrix_length):
                1. for q_acc in segments:
                    1. [ ] if pos odd:
                        1. [ ] trivial
                    2. [ ] elif pos even and max_indel == 1:
                        1. [ ] trivial
                    3. [ ] else:
                        1. [ ] [alignmet[q_acc].append(char) for char in solution]
    """
    assert len(query_to_target_positioned_dict) > 0
    target_vector_length = len( list(query_to_target_positioned_dict.values())[0][0])
    assert stop < target_vector_length # vector coordinates are 0-indexed

    # get allreads alignments covering interesting segment
    # row_accessions = []
    segments = {}
    segment_lists = []
    for q_acc, (q_to_t_pos, t_vector_start, t_vector_end) in query_to_target_positioned_dict.items():
        if t_vector_start <= start and t_vector_end >= stop: # read cover region
            segment = q_to_t_pos[start - t_vector_start : stop - t_vector_start +1 ]
            # row_accessions.append(q_acc)
            segments[q_acc] = segment
            segment_lists.append(segment)

    # 1
    unique_insertions = [set(i) for i in zip(*segment_lists)]
    # unique_insertions = [set() for i in range(len(segment))]
    nr_pos = len(segments[q_acc])
    # for q_acc in segments:
    #     segment = segments[q_acc]
    #     [ unique_insertions[p].add(segment[p]) for p in range(nr_pos) ]


    # 2
    # unique_insertions = []
    # for p in range(nr_pos):
    #     unique_insertions.append(set(position_insertions[p]))

    # 3
    max_insertions = []
    for p in range(nr_pos):
        max_ins_len = len(max(unique_insertions[p], key=len))
        if max_ins_len > 1:
            max_ins = sorted([ins for ins in unique_insertions[p] if len(ins) == max_ins_len ])[0]
            assert p % 2 == 0 # has to be between positions in reference
            max_ins =  "-" + max_ins + "-" 
            max_insertions.append(max_ins)    
        else:
            max_insertions.append("-")    

    # total_matrix_length = sum([len(max_ins) for max_ins in max_insertions]) # we now have all info to get total_matrix_length  

    # 4
    position_solutions = {max_ins : {} for max_ins in max_insertions}

    for nucl in ["A", "G", "C", "T", "-"]:
        position_solutions[nucl] = {"A": ["A"], "G": ["G"], "C": ["C"], "T": ["T"], "-": ["-"]}

    for p in range(nr_pos):
        max_ins = max_insertions[p]
        if len(max_ins) > 1:
            for ins in unique_insertions[p]:
                if ins not in position_solutions[max_ins]:
                    position_solutions[max_ins][ins] = get_best_solution(max_ins, ins)

    # 5 create alignment matrix
    alignment_matrix = {}
    for q_acc in segments:
        alignment_matrix[q_acc] = []
        # tmp_list = []
        seqs = segments[q_acc]
        tmp_list = [ position_solutions[max_insertions[p]][seqs[p]] for p in range(nr_pos)]

        # for p in range(nr_pos):
        #     max_ins = max_insertions[p]
        #     seq = seqs[p]
        #     tmp_list.append(position_solutions[max_ins][seq])

        alignment_matrix[q_acc] = [character for sublist in tmp_list for character in sublist]

    # assert total_matrix_length == len(alignment_matrix[q_acc])
    return alignment_matrix



def min_ed(max_insertion, q_ins):
    result = edlib.align(max_insertion, q_ins, task="path", mode= "NW")
    cigar = result["cigar"]
    tuples = []
    # do not allow deletions in max_insertion: because we do not want to alter this sequence
    if "D" in cigar:
        return ""
    matches = re.split(r'[=DXSMI]+', cigar)
    i = 0
    for length in matches[:-1]:
        i += len(length)
        type_ = cigar[i]
        i += 1
        tuples.append((int(length), type_ ))

    q_insertion_modified = ""
    q_ins_pos = 0
    for length, type_ in tuples:
        # if we reach here we are guaranteed no deletions in alignment of max_insertion
        # we therefore simply thread in the matching or mismatching characters (if type '=' or 'X')
        # or we put a "-" (if type is 'I')
        if type_ == "I":
            q_insertion_modified += "-"*length
        else:
            q_insertion_modified += q_ins[q_ins_pos : q_ins_pos + length]
            q_ins_pos += length
    return q_insertion_modified

 

class TestFunctions(unittest.TestCase):

    def test_positioning(self):
        """
            the start and end of query positions onto the target is always one position 
            before and one position after the start and end of the target coordinates.
        """
        t_seq = "ACGGA"

        q1, t, t_start = "ACGGA", "ACGGA", 0
        self.assertEqual(position_query_to_alignment(q1, t, t_start), (["-", "A","-", "C", "-","G","-", "G", "-", "A", "-"], 0, 10) )

        q, t, t_start = "TACGGA", "-ACGGA", 0
        self.assertEqual(position_query_to_alignment(q, t, t_start), (["T", "A","-", "C", "-","G","-", "G", "-", "A", "-"], 0, 10) )

        q, t, t_start = "ACGGATTT", "ACGGA---", 0
        self.assertEqual(position_query_to_alignment(q, t, t_start), (["-", "A","-", "C", "-","G","-", "G", "-", "A", "TTT"], 0, 10) )

        q, t, t_start = "ACG", "ACG", 0
        self.assertEqual(position_query_to_alignment(q, t, t_start), (["-", "A","-", "C", "-","G","-"], 0, 6) )

        q, t, t_start = "GGA", "GGA", 2
        self.assertEqual(position_query_to_alignment(q, t, t_start), (["-","G","-", "G", "-", "A", "-"], 4, 10) )

        q, t, t_start = "ACGGCC-", "ACGG--A", 0
        self.assertEqual(position_query_to_alignment(q, t, t_start), (["-", "A","-", "C", "-","G","-", "G", "CC", "-", "-"], 0, 10) )

        q, t, t_start = "ACGGCC", "ACGG--", 0
        self.assertEqual(position_query_to_alignment(q, t, t_start), (["-", "A","-", "C", "-","G","-", "G", "CC"], 0, 8) )

        q, t, t_start = "AC-GA", "ACGGA", 0
        self.assertEqual(position_query_to_alignment(q, t, t_start), (["-", "A","-", "C", "-","-","-", "G", "-", "A", "-"], 0, 10) )

    def test_intervals(self):
        ranges = [(0,100,'a'),(0,75,'b'),(95,150,'c'),(120,130,'d')]
        result = [(0,75,("a", "b")),(76, 94, ("a",)), (95,100, ("a", "c")), (101,119, ("c",)), (120,130,("c","d")), (131,150,("c",)) ]
        self.assertEqual(get_non_overlapping_intervals(ranges), result)

    def test_alignment_matrix(self):
        self.maxDiff = None
        start, stop = 0, 10
        query_to_target_positioned_dict = {"q1" : ( ["-", "A","-", "C", "-", "G", "ACCG",   "G", "-", "A", "TTT"], 0, 10),
                                            "q2" : (["-", "A","-", "C", "-", "G", "AG",     "G", "-", "A", "TTT"], 0, 10),
                                            "q3" : (["-", "A","-", "C", "-", "G", "A",      "G", "-", "A", "TTT"], 0, 10),
                                            "q4" : (["-", "A","-", "C", "-", "G", "CC",     "G", "-", "A", "-"], 0, 10),
                                            "q5" : (["-", "A","-", "C", "-", "G", "-",      "G", "-", "A", "T"], 0, 10),
                                            "q6" : (["G", "A","-", "C", "-", "G", "C",      "G", "-", "A", "-"], 0, 10)}

        alignment_matrix = {"q1" : ["-", "A","-", "C", "-", "G","A","C","C","G", "G", "-", "A", "T","T","T"],
                            "q2" : ["-", "A","-", "C", "-", "G","A","G","-","-", "G", "-", "A", "T","T","T"],
                            "q3" : ["-", "A","-", "C", "-", "G","A","-","-","-", "G", "-", "A", "T","T","T"],
                            "q4" : ["-", "A","-", "C", "-", "G","C","C","-","-", "G", "-", "A", "-","-","-"],
                            "q5" : ["-", "A","-", "C", "-", "G","-","-","-","-", "G", "-", "A", "T","-","-"],
                            "q6" : ["G", "A","-", "C", "-", "G","C","-","-","-", "G", "-", "A", "-","-","-"]
                            }

        self.assertEqual(create_multialignment_format(query_to_target_positioned_dict, start, stop), alignment_matrix )

        start, stop = 3, 7
        query_to_target_positioned_dict = {"q1" : ( ["-", "A","-", "C", "-", "G", "ACCG",   "G", "-", "A", "TTT"], 0, 10),
                                            "q2" : (["-", "A","-", "C", "-", "G", "AG",     "G", "-", "A", "TTT"], 0, 10),
                                            "q3" : (["-", "A","-", "C", "-", "G", "A",      "G", "-", "A", "TTT"], 0, 10),
                                            "q4" : (["-", "A","-", "C", "-", "G", "CC",     "G", "-", "A", "-"], 0, 10),
                                            "q5" : (["-", "A","-", "C", "-", "G", "-",      "G", "-", "A", "T"], 0, 10),
                                            "q6" : (["G", "A","-", "C", "-", "G", "C",      "G", "-", "A", "-"], 0, 10)}

        alignment_matrix = {"q1" : [ "C", "-", "G","A","C","C","G", "G"],
                            "q2" : [ "C", "-", "G","A","G","-","-", "G"],
                            "q3" : [ "C", "-", "G","A","-","-","-", "G"],
                            "q4" : [ "C", "-", "G","C","C","-","-", "G"],
                            "q5" : [ "C", "-", "G","-","-","-","-", "G"],
                            "q6" : [ "C", "-", "G","C","-","-","-", "G"]
                            }

        self.assertEqual(create_multialignment_format(query_to_target_positioned_dict, start, stop), alignment_matrix )


        start, stop = 2, 8
        query_to_target_positioned_dict = {"q1" : ( [         "-", "A","-", "C", "-", "G", "ACCG",   "G", "-", "A", "TTT"], 0, 10),
                                            "q2" : ([         "-", "A","-", "C", "-", "G", "AG",     "G", "-", "A", "TTT"], 3, 13),
                                            "q3" : (["-", "A","-", "C", "-", "G", "A",      "G", "-", "A", "TTT"], -2, 8),
                                            "q4" : ([         "-", "A","-", "C", "-", "G", "CC",     "G", "A", "A", "-"], 0, 10),
                                            "q5" : ([               "-", "A","-", "C", "-", "G", "-",      "G", "-", "A", "T"], 20, 30),
                                            "q6" : ([         "G", "A","-", "C", "-", "G", "C",      "G", "-", "A", "-"], 0, 10)}

        alignment_matrix = {"q1" : ["-", "C", "-", "G", "A", "C", "C", "G", "G", "-", "-", "-"],
                            "q3" : ["-", "G", "A", "G", "-", "-", "-", "-", "A", "T", "T", "T"],
                            "q4" : ["-", "C", "-", "G", "C", "C", "-", "-", "G", "A", "-", "-"],
                            "q6" : ["-", "C", "-", "G", "C", "-", "-", "-", "G", "-", "-", "-"]
                            }

        self.assertEqual(create_multialignment_format(query_to_target_positioned_dict, start, stop), alignment_matrix )

if __name__ == '__main__':
    unittest.main()
