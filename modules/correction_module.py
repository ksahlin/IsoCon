import signal
from multiprocessing import Pool
import multiprocessing as mp
import sys
import math
import copy
import random
from collections import Counter

from modules.functions import create_multialignment_matrix, create_position_frequency_matrix

def correct_strings(partition_alignments, seq_to_acc, ccs_dict, step, nr_cores = 1, verbose = False):
    S_prime = {}
    S_prime_quality = {}

    partition_unique_seq_to_acc = {}
    for m, partition in partition_alignments.items():
        partition_unique_seq_to_acc[m] = {}
        partition_unique_seq_to_acc[m][m] = seq_to_acc[m]
        for s in partition:
            if s in seq_to_acc:
                s_accessions = seq_to_acc[s]
                partition_unique_seq_to_acc[m][s] = s_accessions

    if ccs_dict:
        partitioned_ccs_dict = {}
        for m, partition in partition_alignments.items():
            partitioned_ccs_dict[m] = {}
            for s in partition:
                if s in seq_to_acc:
                    s_accessions = seq_to_acc[s]
                    for s_acc in s_accessions:
                        partitioned_ccs_dict[m][s_acc] = ccs_dict[s_acc]
    else:
        partitioned_ccs_dict = {}
        for m, partition in partition_alignments.items():
            partitioned_ccs_dict[m] = {}       

    if nr_cores == 1:
        for m, partition in sorted(partition_alignments.items()):
            S_prime_partition, S_prime_quality_vectors = correct_to_consensus_helper( ((m, partition, partition_unique_seq_to_acc[m], step, verbose, partitioned_ccs_dict[m]), {}) )
            for acc, s in S_prime_partition.items():
                assert acc not in S_prime
                S_prime[acc] = s

            for acc, qual_vector in S_prime_quality_vectors.items():
                S_prime_quality[acc] = qual_vector


    else:
        ####### parallelize statistical tests #########
        # pool = Pool(processes=mp.cpu_count())
        original_sigint_handler = signal.signal(signal.SIGINT, signal.SIG_IGN)
        signal.signal(signal.SIGINT, original_sigint_handler)
        pool = Pool(processes=nr_cores)
        try:
            res = pool.map_async(correct_to_consensus_helper, [ ( (m, partition, partition_unique_seq_to_acc[m], step, verbose, partitioned_ccs_dict[m]), {}) for m, partition in partition_alignments.items() if len(partition) > 1 ] )
            S_prime_partition_dicts =res.get(999999999) # Without the timeout this blocking call ignores all signals.
        except KeyboardInterrupt:
            print("Caught KeyboardInterrupt, terminating workers")
            pool.terminate()
            sys.exit()
        else:
            # print("Normal termination")
            pool.close()
        pool.join()
        for S_prime_partition, S_prime_quality_vectors in S_prime_partition_dicts:
            for acc, s in S_prime_partition.items():
                assert acc not in S_prime
                S_prime[acc] = s

            for acc, qual_vector in S_prime_quality_vectors.items():
                S_prime_quality[acc] = qual_vector

    return S_prime, S_prime_quality


def correct_to_consensus_helper(arguments):
    args, kwargs = arguments
    if args[5]:
        print("Correction with ccs probabilities")
        return correct_to_consensus_ccs_qual(*args, **kwargs)
    else:
        return correct_to_consensus(*args[:-1], **kwargs), {}


def annotate_with_quality_values(alignment_matrix, seq_to_acc, ccs_dict):
    alignment_matrix_of_qualities = {}
    alignment_matrix_of_max_qualities = {}

    for s in alignment_matrix:
        s_accessions = seq_to_acc[s]
        all_quals = [ ccs_dict[s_acc].qual for s_acc in s_accessions ]
        sum_quals_vector = [sum(t) for t in zip(*all_quals)]
        max_quals_vector = [max(t) for t in zip(*all_quals)]
        # sum all probabilities from all reads equal to s here using all accessions in s to acc s_to_acc 
        
        list_sum_quals = []
        list_max_quals = []
        current_pos_in_s = 0

        for j in range(len(alignment_matrix[s])):
            # print(current_pos_in_s, len(sum_quals_vector) )
            current_quality = sum_quals_vector[current_pos_in_s]
            current_max_quality = max_quals_vector[current_pos_in_s]

            list_sum_quals.append(current_quality)
            list_max_quals.append(current_max_quality)

            char_at_pos = alignment_matrix[s][j]
            if char_at_pos != "-":
                if current_pos_in_s < len(sum_quals_vector) - 1:
                    current_pos_in_s += 1 

        alignment_matrix_of_qualities[s] = list_sum_quals
        alignment_matrix_of_max_qualities[s] = list_max_quals

    PFM_qualities = []
    PFM_max_qualities = []
    for j in range(len(alignment_matrix[s])): # for each column
        PFM_qualities.append({"A": 0, "C": 0, "G": 0, "T": 0, "-": 0})
        PFM_max_qualities.append({"A": 0, "C": 0, "G": 0, "T": 0, "-": 0})
        for s in alignment_matrix_of_qualities:
            nucl = alignment_matrix[s][j]
            sum_quality_at_position = alignment_matrix_of_qualities[s][j]
            PFM_qualities[j][nucl] += sum_quality_at_position

            max_quality_at_position = alignment_matrix_of_max_qualities[s][j]
            PFM_max_qualities[j][nucl] += max_quality_at_position



    # get all the differences to majority here  
    list_of_majority_nucleotides = []
    for j in range(len(PFM_qualities)):
        max_v_j = max(PFM_qualities[j], key = lambda x: PFM_qualities[j][x] )
        majority_count =  PFM_qualities[j][max_v_j]
        max_v_j_set = set([v for v in PFM_qualities[j] if PFM_qualities[j][v] == majority_count ])
        all_major = "".join(max_v_j_set)
        list_of_majority_nucleotides.append(all_major)

    assert len(list_of_majority_nucleotides) == len(PFM_qualities)
    global_all_difference_qualities = []    
    for s in alignment_matrix:
        s_aligned_vector = alignment_matrix[s]
        for j in range(len(s_aligned_vector)):
            if s_aligned_vector[j] not in list_of_majority_nucleotides[j] and len(list_of_majority_nucleotides[j]) == 1:
                global_all_difference_qualities.append(alignment_matrix_of_qualities[s][j])

    global_all_difference_qualities.sort()
    if len(global_all_difference_qualities) > 0:
        global_correction_threshold = global_all_difference_qualities[ int( math.ceil( len(global_all_difference_qualities)/2.0) ) - 1 ]
    else:
        global_correction_threshold = -1
        print("nothing to correct!")
    print("GLOBAL QUAL THRESH:", global_correction_threshold)
    return alignment_matrix_of_qualities, PFM_qualities, PFM_max_qualities, global_correction_threshold



def correct_to_consensus_ccs_qual(m, partition, seq_to_acc, step, ccs_dict):
    S_prime_partition = {}
    S_prime_quality_vector = {}
    N_t = sum([container_tuple[3] for s, container_tuple in partition.items()]) # total number of sequences in partition
    
    if len(partition) > 1:
        # all strings has not converged
        alignment_matrix = create_multialignment_matrix(m, partition) 
        PFM = create_position_frequency_matrix(alignment_matrix, partition)

        for s_before in partition:
            s_after = "".join([n for n in alignment_matrix[s_before] if n != "-"])
            assert s_before == s_after
        # print(len(partition), N_t)
        alignment_matrix_of_qualities, PFM_qualities, PFM_max_qualities, global_correction_threshold = annotate_with_quality_values(alignment_matrix, seq_to_acc, ccs_dict)

        assert len(alignment_matrix_of_qualities) == len(alignment_matrix)
        if global_correction_threshold < 0:
            return S_prime_partition, S_prime_quality_vector

        majority_vector = []
        for j in range(len(PFM_qualities)):
            max_v_j = max(PFM_qualities[j], key = lambda x: PFM_qualities[j][x] )

            majority_count =  PFM_qualities[j][max_v_j]
            max_v_j_set = set([v for v in PFM_qualities[j] if PFM_qualities[j][v] == majority_count ])
            all_major = "".join(max_v_j_set)
            majority_vector.append( all_major )
        assert len(majority_vector) == len(PFM_qualities)
        ############################
        ############################


        for s in sorted(partition):
            if partition[s][3] > 1: # at least 2 identical sequences --> its a nearest_neighbor of the partition, has converged, and should not be corrected
                print("not correcting converged sequence!")
                continue

            s_alignment_in_matrix = alignment_matrix[s]
            # s_min = 0
            for i, p in enumerate(s_alignment_in_matrix):
                if p != "-":
                    s_min = i
                    break
            for i, p in enumerate(s_alignment_in_matrix[::-1]):
                if p != "-":
                    s_max = len(s_alignment_in_matrix) - i
                    break

            print("S", s_min, s_max)
            s_quals = alignment_matrix_of_qualities[s]
            # ALL POSITIONS with sum of probabilities lower than the largest probability
            minority_positions = [ (j,majority_vector[j], s_alignment_in_matrix[j]) for j in range(len(majority_vector)) if s_alignment_in_matrix[j] not in majority_vector[j] and s_min <= j <= s_max ]
            minority_positions_correctable = [ j  for j in range(len(majority_vector)) if (len(majority_vector[j]) == 1 and majority_vector[j] != s_alignment_in_matrix[j] ) ]
            minority_positions_correctable = [ (j, PFM_qualities[j][ s_alignment_in_matrix[j] ])  for j in minority_positions_correctable ]
            nr_pos_to_correct = int(math.ceil( len(minority_positions_correctable) * 0.5)) # (step/ float(step +1)) ))
            # print("positions to correct:", nr_pos_to_correct) 

            if nr_pos_to_correct  == 0:
                print("Edit distance to nearest_neighbor:", partition[s][0], "is nearest_neighbor:", s ==m, "Minority positions:", minority_positions)
                continue
            if len(minority_positions_correctable) == 0:
                print("no unambiguous majority positions")
                continue

            minority_positions_correctable.sort(key=lambda x: x[1])
            print(len(minority_positions_correctable) ,minority_positions_correctable)
            _, quality_threshold_to_correct = minority_positions_correctable[ nr_pos_to_correct - 1 ]
            minority_positions_to_correct = [ (j, qual_j) for j, qual_j in minority_positions_correctable if qual_j <= quality_threshold_to_correct ]
            print(quality_threshold_to_correct, len(minority_positions_to_correct))
            # minority_positions_to_correct = [ (j, qual_j) for j, qual_j in minority_positions_correctable if qual_j <= global_correction_threshold ]
            # print("actual:", len(minority_positions_to_correct))
            # minority_positions_to_correct = sorted(minority_positions_correctable, key=lambda x: x[1])[:nr_pos_to_correct]  # sorted list with the smallest probabilities first

            # print(minority_positions_to_correct)
            s_new = alignment_matrix[s]
            s_qual_new = alignment_matrix_of_qualities[s]
            for j, qual_j in minority_positions_to_correct:
                highest_prob_character_at_j = majority_vector[j]
                assert len(majority_vector[j]) == 1
                s_new[j] = highest_prob_character_at_j
                s_qual_new[j] = PFM_max_qualities[j][highest_prob_character_at_j] 
            
            s_modified = "".join([nucl for nucl in s_new if nucl != "-" ])
            s_qual_modified = [s_qual_new[j] for j in range(len(s_new)) if s_new[j] != "-" ]

            # only unique strings can change in this step
            accessions_of_s = seq_to_acc[s] 
            for acc in accessions_of_s:
                S_prime_partition[acc] = s_modified
                S_prime_quality_vector[acc] = s_qual_modified
    else:
        print("Partition converged: Partition size(unique strings):{0}, partition support: {1}.".format(len(partition), N_t))

    
    return S_prime_partition, S_prime_quality_vector



def correct_to_consensus(m, partition, seq_to_acc, step, verbose):
    S_prime_partition = {}
    N_t = sum([container_tuple[3] for s, container_tuple in partition.items()]) # total number of sequences in partition
    
    # if N_t == 2:
    #     print("Partition has size", N_t, "no meaningful correction can be done")
    #     for s, container_tuple in partition.items():
    #         print(seq_to_acc[s])

    if len(partition) > 1 and N_t > 2:
        # all strings has not converged
        alignment_matrix = create_multialignment_matrix(m, partition) 
        PFM = create_position_frequency_matrix(alignment_matrix, partition)
        for s_before in partition:
            s_after = "".join([n for n in alignment_matrix[s_before] if n != "-"])
            assert s_before == s_after

        # consensus_alignment = [ max(PFM[j], key=lambda k: PFM[j][k]) for j in range(len(PFM))]
        # print("nearest_neighbor errors:",  math.ceil(min([ partition[s][0] for s in partition if partition[s][3] > 1 or s !=m ]) / 2.0)  )
        # frozen_positions = get_frozen_positions(alignment_matrix[m])
        ## TEST LOG ERROR TYPES #######
        c_del = 0
        c_ins = 0
        c_subs = 0
        majority_vector = []
        temp_majority_string = ""
        for j in range(len(PFM)):
            max_v_j = max(PFM[j], key = lambda x: PFM[j][x] )

            if max_v_j != "-":
                temp_majority_string += max_v_j

            majority_count =  PFM[j][max_v_j]
            max_v_j_set = set([v for v in PFM[j] if PFM[j][v] == majority_count ])
            all_major = "".join(max_v_j_set)
            majority_vector.append( all_major )
            if len(max_v_j_set) > 1:
                continue # DO not count errors where majority position is ambigous since we don't know the true majority here

            for v in PFM[j]:
                if max_v_j == "-":
                    if v != max_v_j:
                        c_ins += PFM[j][v]
                else:
                    if v != max_v_j:
                        if v == "-":
                            c_del += PFM[j][v]
                        else:
                            c_subs += PFM[j][v]
        if verbose:
            print("Partition error types:", c_del, c_ins, c_subs, "depth:", N_t )
        assert len(majority_vector) == len(PFM)

        ############################

        for s in sorted(partition):
            if partition[s][3] > 1: # at least 2 identical sequences --> its a nearest_neighbor of the partition, has converged, and should not be corrected
                continue

            nr_pos_to_correct2 = int(math.ceil(partition[s][0] / 2.0)) #decide how many errors we should correct here
            s_alignment_in_matrix = alignment_matrix[s]
            
            minority_positions = [ (j,majority_vector[j], s_alignment_in_matrix[j]) for j in range(len(majority_vector)) if s_alignment_in_matrix[j] not in majority_vector[j] ]
            
            nr_pos_to_correct = int(math.ceil( len([ 1 for j in range(len(majority_vector)) if (len(majority_vector[j]) == 1 and majority_vector[j] != s_alignment_in_matrix[j] ) ]) * 0.5)) # (step/ float(step +1)) ))
            # print("positions to correct:", nr_pos_to_correct)

            if verbose:
                if nr_pos_to_correct == 0:
                    print("Edit distance to nearest_neighbor:", partition[s][0], "is nearest_neighbor:", s ==m, "Minority positions:", minority_positions)

                if nr_pos_to_correct2 > 0 and nr_pos_to_correct == 0:
                    print("Edit distance to nearest_neighbor: {0}, {1} minority positions, correcting no position. Length partition (unique): {2}, total seqs: {3}".format(partition[s][0], len(minority_positions), len(partition), N_t))
                    # for s in partition:
                    #     print(s)
                    # print([ (j,majority_vector[j], s_alignment_in_matrix[j], PFM[j]) for j in range(len(majority_vector)) if majority_vector[j] != s_alignment_in_matrix[j] ])
            if nr_pos_to_correct  == 0:
                continue

            # find the position probabilities of the alignment of s in PFM


            ###################### ORIGINAL CORRECTION ######################################
            # pos_freqs_for_s = []
            # for j in range(len(PFM)):
            #     pos_freqs_for_s.append( (j, PFM[j][s_alignment_in_matrix[j]]) )

            # pos_freqs_for_s.sort(key=lambda x: x[1]) # sort with respect to smallest frequencies                    
            # pos, highest_freq_of_error_to_correct = pos_freqs_for_s[ nr_pos_to_correct - 1 ]
            # end_position_in_list = nr_pos_to_correct

            # pp = nr_pos_to_correct
            # while pos_freqs_for_s[pp][1] == highest_freq_of_error_to_correct:
            #     end_position_in_list += 1
            #     pp += 1

            # J = [j for j, freq in random.sample(pos_freqs_for_s[:end_position_in_list], nr_pos_to_correct)]

            ############################################################
            ############################################################

            ########### TEST WEIGHTING EACH MINORITY POSITION BY IT'S OBSERVED FREQUENCY THROUGHOUT THE ALIGNMENTS TO THE nearest_neighbor ################
            pos_freqs_for_s_mod = []
            for j in range(len(PFM)):
                majority_variant = majority_vector[j]
                v_j = s_alignment_in_matrix[j]
                if v_j == majority_variant or len(majority_variant) > 1: # no correction of majority position or ambigous majority!
                    continue
                else:
                    if majority_variant == "-":
                        pos_freqs_for_s_mod.append( (j, PFM[j][v_j] / float(max(c_ins, 1) ) ))
                    elif v_j == "-":
                        pos_freqs_for_s_mod.append( (j, PFM[j][v_j] / float(max(c_del, 1) ) ))
                    else:
                        pos_freqs_for_s_mod.append( (j, PFM[j][v_j] / float(max(c_subs, 1) ) ))


                # max_v_j = max(PFM[j], key = lambda x: PFM[j][x] )
                # if max_v_j == v_j:
                #     pos_freqs_for_s_mod.append( (j, PFM[j][v_j] / float(1)) )
                # elif max_v_j == "-":
                #     pos_freqs_for_s_mod.append( (j, PFM[j][v_j] / float(max(c_ins, 1) ) ))
                # elif v_j == "-":
                #     pos_freqs_for_s_mod.append( (j, PFM[j][v_j] / float(max(c_del, 1) ) ))
                # else:
                #     pos_freqs_for_s_mod.append( (j, PFM[j][v_j] / float(max(c_subs, 1) ) ))
            if len(pos_freqs_for_s_mod) == 0:
                continue

            pos_freqs_for_s_mod.sort(key=lambda x: x[1]) # sort with respect to smallest frequencies    
            if len(pos_freqs_for_s_mod) < nr_pos_to_correct:
                end_position_in_list = len(pos_freqs_for_s_mod)  
            else:   
                pos, highest_freq_of_error_to_correct = pos_freqs_for_s_mod[ nr_pos_to_correct - 1 ]
                end_position_in_list = nr_pos_to_correct
                for pp in range(nr_pos_to_correct, len(pos_freqs_for_s_mod)):
                    # print(pos_freqs_for_s_mod[pp][1], highest_freq_of_error_to_correct)
                    if pos_freqs_for_s_mod[pp][1] > highest_freq_of_error_to_correct:
                        break
                    else:
                        end_position_in_list += 1


            # J = [j for j, freq in random.sample(pos_freqs_for_s_mod[:end_position_in_list], nr_pos_to_correct)]
            J_temp = [j for j, freq in pos_freqs_for_s_mod[:end_position_in_list]] 
            # print("end pos:", end_position_in_list)
            #############################################

            s_new = alignment_matrix[s]
            for j in J_temp:
                # if j in frozen_positions:
                #     # print("tried to correct in frozen positions")
                #     continue

                old_nucl = s_new[j]
                highest_prob_character_at_j = majority_vector[j]
                assert len(majority_vector[j]) == 1
                # highest_prob_character_at_j = max(PFM[j], key=lambda k: PFM[j][k])

                if highest_prob_character_at_j == old_nucl: # choose the other highest on if tie (should happen only when partition consist of two sequences)
                    print("Highest count nucl was about to be corrected.")
                    # pmf_j_minus_variant = copy.deepcopy(PFM[j])
                    # del pmf_j_minus_variant[old_nucl] 
                    # highest_prob_character_at_j = max(pmf_j_minus_variant, key=lambda k: pmf_j_minus_variant[k])


                # print("correcting", s_new[j], "to", highest_prob_character_at_j )
                s_new[j] = highest_prob_character_at_j
            s_modified = "".join([nucl for nucl in s_new if nucl != "-" ])

            # only unique strings can change in this step

            accessions_of_s = seq_to_acc[s] 
            for acc in accessions_of_s:
                S_prime_partition[acc] = s_modified
    elif len(partition) > 1:
        ed = 0
        for s in partition:
            if  partition[s][0] > ed:
                ed = partition[s][0]
        if verbose:
            print("Partition could not be corrected: Partition size(unique strings):{0}, partition support: {1}, edit distance:{2}.".format(len(partition), N_t, ed))
    else:
        if verbose:
            print("Partition converged: Partition size(unique strings):{0}, partition support: {1}.".format(len(partition), N_t))

    return S_prime_partition

# def get_frozen_positions(m_in_alignment_matrix):
#     """
#         positions in Multialingment matrix where there are insels longer than 2 bp w.r.t. nearest_neighbor, these regions are prone to errors in alingments and
#         we wait to correct these in another partition where they are hopefully split into several true partitions.
#     """
#     frozen_pos = set()
#     ins_count = 0
#     ins_pos = set()
#     for j in range(len(m_in_alignment_matrix)):
#         if m_in_alignment_matrix[j] == "-":
#             ins_count += 1
#             ins_pos.add(j)
#         else:
#             if ins_count > 4:  # always padded, so indel of length 2 be will be of length 4 in  alignment matrixx -
#                 for k in ins_pos:
#                     frozen_pos.add(k)
#             ins_count = 0
#             ins_pos = set()
#     print("frozen:", len(frozen_pos), "tot:", len(m_in_alignment_matrix) )
#     return frozen_pos
