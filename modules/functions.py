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

def reverse_complement(string):
    #rev_nuc = {'A':'T', 'C':'G', 'G':'C', 'T':'A', 'N':'N', 'X':'X'}
    # Modified for Abyss output
    rev_nuc = {'A':'T', 'C':'G', 'G':'C', 'T':'A', 'a':'t', 'c':'g', 'g':'c', 't':'a', 'N':'N', 'X':'X', 'n':'n', 'Y':'R', 'R':'Y', 'K':'M', 'M':'K', 'S':'S', 'W':'W', 'B':'V', 'V':'B', 'H':'D', 'D':'H', 'y':'r', 'r':'y', 'k':'m', 'm':'k', 's':'s', 'w':'w', 'b':'v', 'v':'b', 'h':'d', 'd':'h'}

    rev_comp = ''.join([rev_nuc[nucl] for nucl in reversed(string)])
    return(rev_comp)


def cut_ends_of_alignment_matrix(alignment_matrix_to_t, t_acc, c_acc, ignore_ends_len):
    target_alignment = alignment_matrix_to_t[t_acc]
    candidate_alignment = alignment_matrix_to_t[c_acc]
    cut_start = 0
    c_count = 0
    t_count = 0
    for j, (c_nucl, t_nucl) in enumerate(zip(candidate_alignment, target_alignment)):
        if c_nucl != "-":
            c_count += 1
        if t_nucl != "-":
            t_count += 1

        if math.fabs(c_count - t_count) > ignore_ends_len: # we only ignore length differences smaller or equal to ignore_ends_len
            break

        if c_count > 0 and t_count > 0: 
            cut_start = j
            break

    cut_end = 0
    c_count = 0
    t_count = 0
    for j, (c_nucl, t_nucl) in enumerate(zip( reversed(candidate_alignment), reversed(target_alignment))):
        if c_nucl != "-":
            c_count += 1
        if t_nucl != "-":
            t_count += 1

        if math.fabs(c_count - t_count) > ignore_ends_len: # we only ignore length differences smaller or equal to ignore_ends_len
            break

        if c_count > 0 and t_count > 0: 
            cut_end = j
            break
    cut_end = len(target_alignment) - cut_end


    # print("cutting from", len(target_alignment), "positions to", len(target_alignment[ cut_start : cut_end ]) )
    for acc in alignment_matrix_to_t:
        alignment_matrix_to_t[acc] = alignment_matrix_to_t[acc][ cut_start : cut_end ]
    return  alignment_matrix_to_t


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

def get_multiplier_for_variant(state, char, pos, target_alignment, candidate_alignment):
    """
        Entering this function only if state = I or D, and the read has the same character as the target at the starting position "pos".
        That is the premise.
        ## TODO: If delation in candidate, wee need to count the number of identical bases to the deleted base in the target.
        ## TODO: If insertion, the inserted base has to be identical to any of the neighboring nucleotides to count as alignment invariant 
        TODO: If delation in candidate: we need to count the number of identical bases to the deleted base in the target. In that case, the u_v = len(identical stretch in target)
        TODO: If insertion and inserted base == target base: the inserted base has to be identical to any of the neighboring nucleotides to count as alignment. In that case, the u_v = len(identical stretch in target)
    """
    stop = len(target_alignment) - 1
    u_v = 1
    if state == "D":
        v = target_alignment[pos]
    elif state == "I":
        v = char

    # print(v, state)
    offset = 1
    upper_stop = False
    lower_stop = False
    while True:

        if upper_stop:
            pass
        elif pos + offset > stop:
            upper_stop = True
        elif target_alignment[pos + offset] == v: # == candidate_alignment[pos + offset] == v:
            u_v += 1
        elif target_alignment[pos + offset] == candidate_alignment[pos + offset] == "-":
            pass
        else:
            upper_stop = True


        if lower_stop:
            pass
        elif pos - offset < 0:
            lower_stop = True                    
        elif target_alignment[pos - offset] == v: # candidate_alignment[pos - offset] == v:
            u_v += 1
        elif target_alignment[pos - offset] == candidate_alignment[pos - offset] == "-":
            pass
        else:
            lower_stop = True

        if lower_stop == upper_stop == True:
            break


        offset += 1

    return u_v


def get_invariant_multipliers(delta_t, alignment_matrix, t_acc):

    target_alignment = alignment_matrix[t_acc]
    # Get the individual invariant factors u_iv for each read and variant x_i and v in delta, v is a tuple (state, char).
    # store this 3 dimensional in dictionary u_iv = {c_acc: {pos: { (state,char) : u_iv} }} 
    u_iv = {}

    for c_acc in delta_t:
        u_iv[c_acc] = {}
        candidate_alignment = alignment_matrix[c_acc]
        for pos in delta_t[c_acc]:
            if pos not in u_iv[c_acc]:
                u_iv[c_acc][pos] = {}
        
            state, char = delta_t[c_acc][pos]
            if state == "S": # substitutions has u_v =1 by definition
                u_v = 1
            else: # read matches candidate variant, now we need to check how many possible combinations an error can cause them to match due to invariant
                u_v = get_multiplier_for_variant(state, char, pos, target_alignment, candidate_alignment)
                # print("multiplier:", u_v)
                # print(target_alignment[pos-10: pos+10],state, char)
                # print(read_alignment[pos-10: pos+10],state, char)
                # print(candidate_alignment[pos-10: pos+10],state, char)

            u_iv[c_acc][pos][(state, char)] = u_v

    return u_iv


def get_difference_coordinates_for_candidates(target_accession, candidate_accessions, alignment_matrix):
    position_differences = {}
    target_alignment = alignment_matrix[target_accession]
    
    for c_acc in candidate_accessions:
        position_differences[c_acc] = {}
        candidate_alignment = alignment_matrix[c_acc]
        for j in range(len(candidate_alignment)):
            c_base = candidate_alignment[j]
            t_base = target_alignment[j]
            if c_base != t_base:
                if t_base == "-":
                    position_differences[c_acc][j] = ("I", c_base)
                elif c_base == "-":
                    position_differences[c_acc][j] = ("D", c_base)
                else:
                    position_differences[c_acc][j] = ("S", c_base)

        # print("nr v:",len(position_differences[c_acc]))
    return position_differences


def get_ccs_position_prob_per_read(target_accession, target_length, alignment_matrix, invariant_factors_for_candidate, candidate_accessions, Delta_t, ccs_dict, insertions, deletions, substitutions, max_phred_q_trusted):
    probability = {}
    assert len(candidate_accessions) == 1
    c_acc = list(candidate_accessions)[0]
    target_alignment = alignment_matrix[target_accession]
    candidate_alignment = alignment_matrix[c_acc]
    delta_size = float(len(invariant_factors_for_candidate[c_acc]))

    tot_fraction = float(insertions + deletions + substitutions) if insertions + deletions + substitutions > 1 else 1
    subs_ratio = substitutions / tot_fraction
    ins_ratio = insertions / tot_fraction
    del_ratio = deletions / tot_fraction

    # this looping needs re-writing
    # reads can have any basepair at given position.
    # just read of the p-error at the given ccs position and that's it?
    # print(c_acc)
    tmp_weight_equal = 0
    tmp_weight_diff = 0

    # print("t", len("".join([n for n in target_alignment if n != "-"])))
    # print("c", len("".join([n for n in candidate_alignment if n != "-"])))
    for q_acc in alignment_matrix:
        if q_acc == target_accession:
            continue
        if q_acc in candidate_accessions:
            continue  

        probability[q_acc] = 1.0
        ccs_alignment = alignment_matrix[q_acc]
        for pos in Delta_t[c_acc]:
            c_state, c_base = Delta_t[c_acc][pos]
            # determine what type the read has in position
            t_nucl = target_alignment[pos]
            assert c_base != t_nucl


            u_v = invariant_factors_for_candidate[c_acc][pos][(c_state, c_base)]
            ###  base pair quality predictions ###
            ccs_coord = ccs_dict[q_acc].alignment_matrix_pos_to_ccs_coord(ccs_alignment, pos)
            # print(ccs_coord, len(ccs_alignment), pos)
            
            # p_error = ccs_dict[q_acc].get_p_error_in_base(ccs_coord)
            q_qual = ccs_dict[q_acc].qual[ccs_coord]
            # To map quality values
            # [A, B] --> [a, b]  [3, 93] --> [3, 43]
            # (x - A)*(b-a)/(B-A) + a
            q_qual_mapped = (q_qual - 3)*(max_phred_q_trusted - 3.0)/(90.0) + 3

            if u_v > 1: # probability in a homopolymenr region has all it's uncertainty attributed to the lenght of the homopolymer that 
                p_error =  (10**(-q_qual_mapped/10.0))
            elif c_state == "S":
                p_error =  ((10**(-q_qual_mapped/10.0))*subs_ratio)/3.0 # probability that its an identical substitution error from a base call uncertainty
            elif  c_state == "I":
                p_error =  ((10**(-q_qual_mapped/10.0))*ins_ratio)/4.0 # probability that its an identical insertion error from a base call uncertainty
            elif c_state == "D":
                p_error =  (10**(-q_qual_mapped/10.0)) * del_ratio # probability that its a delation error from a base call uncertainty
            else:
                p_error = -1 # If end up here, there is a bug

            # print(p_error, q_qual_mapped, q_qual)
            assert 0.0 < p_error < 1.0
            probability[q_acc] *= p_error #max(p_error, min_uncertainty)

            # probability[q_acc] *= p_error
            # print("".join([n for n in target_alignment[pos-100:pos+100]]))
            # print("".join([n for n in candidate_alignment[pos-100:pos+100]]))
            # print("".join([n for n in ccs_alignment[pos-100:pos+100]]))
            # print("acc seq:", q_acc)
            # print("ccs bam seq:", ccs_dict[q_acc].seq)
            # print(candidate_alignment[pos-10:pos+10], pos, ccs_coord)
            # print('pos:', ccs_dict[q_acc].seq.find("GTCACTGCTGGATATCA"), "pred coord:", ccs_coord)
            # print(pos, c_state, c_base, ccs_alignment[pos], ccs_dict[q_acc].seq[ccs_coord],   ccs_alignment[pos-10:pos+10],  ccs_dict[q_acc].seq[ccs_coord-6:ccs_coord +11], p_error, ccs_dict[q_acc].np, q_acc)

    return probability



def get_prob_of_error_per_read(target_accession, segment_length, candidate_accessions, errors, invariant_factors_for_candidate):
    probability = {}
    assert len(invariant_factors_for_candidate) == 1
    c_acc = list(invariant_factors_for_candidate.keys())[0]
    delta_size = float(len(invariant_factors_for_candidate[c_acc]))

    for q_acc in errors:
        probability[q_acc] = 1.0
        p_S = ( max(errors[q_acc]["S"], delta_size) / float(segment_length) ) / 3.0   # p = 0.0 not allowed, min_p is 1/(3*len(seq))
        p_I = ( max(errors[q_acc]["I"], delta_size) / float(segment_length) ) / 4.0   # p = 0.0 not allowed, min_p is 1/(4*len(seq))
        p_D = ( max(errors[q_acc]["D"], delta_size) / float(segment_length) )         # p = 0.0 not allowed, min_p is 1/(len(seq))

        for pos in invariant_factors_for_candidate[c_acc]:
            for (state, char) in invariant_factors_for_candidate[c_acc][pos]:
                u_v = invariant_factors_for_candidate[c_acc][pos][(state, char)]
                if state == "S":
                    probability[q_acc] *= p_S*u_v 
                elif state == "I":
                    probability[q_acc] *= min(0.5, p_I*u_v) 
                elif state == "D":
                    probability[q_acc] *= min(0.5, p_D*u_v)
        if probability[q_acc] >= 1.0:
            probability[q_acc] = 0.99999

    return probability


def get_errors_for_partitions(target_accession, segment_length, candidate_accessions, alignment_matrix):
    errors = {}
    target_alignment = alignment_matrix[target_accession]
    assert len(candidate_accessions) == 1
    c_acc = list(candidate_accessions)[0]

    # ed_poisson_i, ed_poisson_s, ed_poisson_d = 0, 0, 0
    candidate_alignment = alignment_matrix[c_acc]

    for q_acc in alignment_matrix:
        if q_acc == target_accession:
            continue
        if q_acc in candidate_accessions:
            continue  
        query_alignment = alignment_matrix[q_acc]

        s_min = 0
        for i, p in enumerate(query_alignment):
            if p != "-":
                s_min = i
                break
        s_max = len(query_alignment)
        for i, p in enumerate(query_alignment[::-1]):
            if p != "-":
                s_max = len(query_alignment) - i
                break

        errors[q_acc] = {}
        ed_i, ed_s, ed_d = 0.0, 0.0, 0.0
        # ed_i_to_c, ed_s_to_c, ed_d_to_c = 0.0, 0.0, 0.0

        for j in range(len(query_alignment)):
            if j < s_min or j > s_max:
                continue

            q_base = query_alignment[j]
            t_base = target_alignment[j]
            if q_base != t_base:
                if t_base == "-":
                    ed_i += 1
                elif q_base == "-":
                    ed_d += 1
                else:
                    ed_s += 1
            
            # c_base = candidate_alignment[j]
            # if q_base != c_base:
            #     if c_base == "-":
            #         ed_i_to_c += 1
            #     elif q_base == "-":
            #         ed_d_to_c += 1
            #     else:
            #         ed_s_to_c += 1

        errors[q_acc]["I"] = ed_i # min(ed_i, ed_i_to_c)
        errors[q_acc]["S"] = ed_s # min(ed_s, ed_s_to_c)
        errors[q_acc]["D"] = ed_d # min(ed_d, ed_d_to_c)
        # print(ed_i, ed_d, ed_s)
    insertions, deletions, substitutions = 0, 0, 0
    for q_acc in errors:
        insertions += errors[q_acc]["I"]
        deletions += errors[q_acc]["D"]
        substitutions += errors[q_acc]["S"]

    # print("I", insertions, "D", deletions, "S", substitutions)
    # sys.exit()
    # if insertions, deletions, substitutions == 0:
    return insertions, deletions, substitutions




def get_errors_per_read(target_accession, segment_length, candidate_accessions, alignment_matrix):
    errors = {}
    target_alignment = alignment_matrix[target_accession]
    assert len(candidate_accessions) == 1
    c_acc = list(candidate_accessions)[0]

    # ed_poisson_i, ed_poisson_s, ed_poisson_d = 0, 0, 0
    candidate_alignment = alignment_matrix[c_acc]

    for q_acc in alignment_matrix:
        if q_acc == target_accession:
            continue
        if q_acc in candidate_accessions:
            continue  

        errors[q_acc] = {}
        query_alignment = alignment_matrix[q_acc]
        ed_i, ed_s, ed_d = 0.0, 0.0, 0.0
        ed_i_to_c, ed_s_to_c, ed_d_to_c = 0.0, 0.0, 0.0

        for j in range(len(query_alignment)):
            q_base = query_alignment[j]
            
            t_base = target_alignment[j]
            if q_base != t_base:
                if t_base == "-":
                    ed_i += 1
                elif q_base == "-":
                    ed_d += 1
                else:
                    ed_s += 1
            
            c_base = candidate_alignment[j]
            if q_base != c_base:
                if c_base == "-":
                    ed_i_to_c += 1
                elif q_base == "-":
                    ed_d_to_c += 1
                else:
                    ed_s_to_c += 1

        errors[q_acc]["I"] = min(ed_i, ed_i_to_c)
        errors[q_acc]["S"] = min(ed_s, ed_s_to_c)
        errors[q_acc]["D"] = min(ed_d, ed_d_to_c)

        # get poisson counts on all positions
        # for j in range(len(query_alignment)):
        #     target_alignment = alignment_matrix[target_accession]
        #     # candidate_alignment = alignment_matrix[x_to_c_acc[q_acc]]
        #     # if j not in forbidden:
        #     q_base = query_alignment[j]
        #     t_base = target_alignment[j]
        #     if q_base != t_base:
        #         if t_base == "-":
        #             ed_poisson_i += 1
        #         elif q_base == "-":
        #             ed_poisson_d += 1
        #         else:
        #             ed_poisson_s += 1      



    # print(errors)
    return errors

def reads_supporting_candidate(target_accession, candidate_accessions, alignment_matrix, Delta_t, partition_of_X):
    assert len(candidate_accessions) == 1
    c_acc = list(candidate_accessions)[0]
    x = []
    for q_acc in partition_of_X[c_acc].union(partition_of_X[target_accession]):
        if q_acc not in alignment_matrix:
            print("READ {0} ALIGNED TO {1} BUT FAILED TO ALIGN TO {2}".format(q_acc, c_acc, target_accession) )
            continue
        query_alignment = alignment_matrix[q_acc]    
        support = 1
        for delta in Delta_t[c_acc]:
            q_base = query_alignment[delta]
            c_state, c_base = Delta_t[c_acc][delta]
            if q_base != c_base:
                support = 0

        if support:
            x.append(q_acc)
    return x


def create_position_probability_matrix(m, partition):
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
        s_positioned, target_vector_start_position, target_vector_end_position = position_query_to_alignment(s_alignment, m_alignment, 0)
        assert target_vector_start_position == 0
        assert target_vector_end_position + 1 == 2*len(m) + 1 # vector positions are 0-indexed
        query_to_target_positioned_dict[q_acc] = (s_positioned, target_vector_start_position, target_vector_end_position)

    alignment_matrix = create_multialignment_format(query_to_target_positioned_dict, 0, 2*len(m))

    # N_t = sum([container_tuple[3] for q_acc, container_tuple in partition.items()]) # total number of sequences in partition
    # print("total seq multiset:", N_t, "total seqs in set:", len(partition))
    PFM = []
    for j in range(len(alignment_matrix[q_acc])): # for each column
        PFM.append({"A": 0, "C": 0, "G": 0, "T": 0, "-": 0})
        for s in alignment_matrix:
            nucl = alignment_matrix[s][j]
            indegree = partition[s][3]
            PFM[j][nucl] += indegree
    # print( "matrix length:", len(alignment_matrix[m]))
    return alignment_matrix, PFM

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
            max_insertion, q_acc_max_ins = max(insertions, key= lambda x : len(x[0]))
            
            max_ins_len = len(max_insertion)
            all_max_ins = set([ins for (ins, acc) in insertions if len(ins) == max_ins_len])
            # if max_ins_len > 4:
            #     print(j, all_max_ins, max_insertion, q_acc_max_ins)
            # if len(all_max_ins) > 1 and max_ins_len == 2:
            #     print("pos", j, all_max_ins)
            max_insertion = sorted(all_max_ins)[0] 
            max_insertion = "-" + max_insertion + "-"  # pad the max insertion

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
                    # else, check if smaller deletion can be aligned from left to write, e.g. say max deletion is GACG
                    # then an insertion AG may be aligned as -A-G. Take this alignment instead
                    q_insertion_modified = min_ed(max_insertion, q_ins)
                    # q_insertion_modified = thread_to_max_ins(max_insertion, q_ins)
                    # if q_insertion_modified:
                    #     print("new threaded:", q_insertion_modified2)
                    #     print("threaded:", q_insertion_modified)

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
                
                # if max_ins_len > 4:
                #     print("q:{0} max:{1},q_original:{2} ".format(q_insertion_modified, max_insertion, q_ins) )

                # if q_ins != "-" and q_ins != "".join([n for n in q_insertion_modified if n != "-"]):
                #     print("BUG: q:{0} max:{1},q_original:{2} pos:{3}, new:{4} ".format(q_insertion_modified, max_insertion, q_ins, j, q_insertion_modified2) )
                #     # sys.exit()

                #### finally add to alignment matrix
                # if len(q_insertion_modified) != len(max_insertion):
                #     print(q_insertion_modified, max_insertion, q_ins)
                assert len(q_insertion_modified) == len(max_insertion)
                for p in range(len(max_insertion)):
                    alignment_matrix[q_acc].append(q_insertion_modified[p])

        else: # we are on a target base pair -- all varinats must be exactly A,C,G,T, - i.e., length 1
            for q_acc in segments:
                alignment_matrix[q_acc].append(segments[q_acc][j])
    # print(alignment_matrix)
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
