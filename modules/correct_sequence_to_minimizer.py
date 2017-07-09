import signal
from multiprocessing import Pool
import multiprocessing as mp
import sys
import math
import copy
import random
from collections import Counter

from modules.functions import create_position_probability_matrix

def correct_strings(partition_alignments, seq_to_acc, step, single_core = False):
    S_prime = {}

    partition_unique_seq_to_acc = {}
    for m, partition in partition_alignments.items():
        partition_unique_seq_to_acc[m] = {}
        partition_unique_seq_to_acc[m][m] = seq_to_acc[m]
        for s in partition:
            if s in seq_to_acc:
                partition_unique_seq_to_acc[m][s] = seq_to_acc[s]


    if single_core:
        for m, partition in sorted(partition_alignments.items()):
            S_prime_partition = correct_to_consensus(m, partition, partition_unique_seq_to_acc[m], step)
            for acc, s in S_prime_partition.items():
                assert acc not in S_prime
                S_prime[acc] = s

    else:
        ####### parallelize statistical tests #########
        # pool = Pool(processes=mp.cpu_count())
        original_sigint_handler = signal.signal(signal.SIGINT, signal.SIG_IGN)
        signal.signal(signal.SIGINT, original_sigint_handler)
        pool = Pool(processes=mp.cpu_count())
        try:
            res = pool.map_async(correct_to_consensus_helper, [ ( (m, partition, partition_unique_seq_to_acc[m], step), {}) for m, partition in partition_alignments.items() if len(partition) > 1 ] )
            S_prime_partition_dicts =res.get(999999999) # Without the timeout this blocking call ignores all signals.
        except KeyboardInterrupt:
            print("Caught KeyboardInterrupt, terminating workers")
            pool.terminate()
            sys.exit()
        else:
            print("Normal termination")
            pool.close()
        pool.join()
        for S_prime_partition in S_prime_partition_dicts:
            for acc, s in S_prime_partition.items():
                assert acc not in S_prime
                S_prime[acc] = s


    return S_prime

def correct_to_minimzer_helper(arguments):
    args, kwargs = arguments
    return correct_to_minimizer(*args, **kwargs)

def correct_to_consensus_helper(arguments):
    args, kwargs = arguments
    return correct_to_consensus(*args, **kwargs)

def correct_to_minimizer(m, partition, seq_to_acc):
    S_prime_partition = {}

    N_t = sum([container_tuple[3] for s, container_tuple in partition.items()]) # total number of sequences in partition
    
    if N_t == 2:
        print("Partition has size", N_t, "no meaningful correction can be done")

    if len(partition) > 1 and N_t > 2:
        # all strings has not converged
        alignment_matrix, PFM = create_position_probability_matrix(m, partition) 
        
        # print("minimizer errors:",  math.ceil(min([ partition[s][0] for s in partition if partition[s][3] > 1 or s !=m ]) / 2.0)  )
        # minimizer_errors = min([ partition[s][0] for s in partition if partition[s][3] > 1 or s !=m ])
        # minimizer_errors = math.ceil(min([ partition[s][0] for s in partition if partition[s][3] > 1 or s !=m ]) / 2.0)

        ### TEST LOG ERROR TYPES #######
        # c = Counter()
        # for j in range(len(PFM)):
        #     max_v_j = max(PFM[j], key = lambda x: PFM[j][x] )
        #     for v in PFM[j]:
        #         if v != max_v_j:
        #            c[v] += PFM[j][v]
        # print("Error types:", c, "depth:", len(partition) )
        #############################

        ## TEST LOG ERROR TYPES #######
        c_del = 0
        c_ins = 0
        c_subs = 0
        for j in range(len(PFM)):
            max_v_j = max(PFM[j], key = lambda x: PFM[j][x] )
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

        print("Error types:", c_del, c_ins, c_subs, "depth:", len(partition) )
        ############################

        for s in partition:
            # if minimizer_errors < partition[s][0]:
            #     nr_pos_to_correct = int( partition[s][0] - minimizer_errors ) #decide how many errors we should correct here
            # else:
            #     nr_pos_to_correct = int(math.ceil(partition[s][0] / 2.0)) #decide how many errors we should correct here
            # nr_pos_to_correct = max(int( partition[s][0] - minimizer_errors ), int(math.ceil(partition[s][0] / 2.0)))
            nr_pos_to_correct = int(math.ceil(partition[s][0] / 2.0)) #decide how many errors we should correct here

            # print("positions to correct for sequence s:", nr_pos_to_correct, s ==m)
            if nr_pos_to_correct  == 0:
                continue

            s_alignment_in_matrix = alignment_matrix[s]
            # find the position probabilities of the alignment of s in PFM

            pos_freqs_for_s = []
            for j in range(len(PFM)):
                pos_freqs_for_s.append( (j, PFM[j][s_alignment_in_matrix[j]]) )

            pos_freqs_for_s.sort(key=lambda x: x[1]) # sort with respect to smallest frequencies                    
            pos, highest_freq_of_error_to_correct = pos_freqs_for_s[ nr_pos_to_correct - 1 ]
            end_position_in_list = nr_pos_to_correct

            pp = nr_pos_to_correct
            while pos_freqs_for_s[pp][1] == highest_freq_of_error_to_correct:
                end_position_in_list += 1
                pp += 1

            J = [j for j, freq in random.sample(pos_freqs_for_s[:end_position_in_list], nr_pos_to_correct)]




            ########### TEST WEIGHTING EACH MINORITY POSITION BY IT'S OBSERVED FREQUENCY THROUGHOUT THE ALIGNMENTS TO THE MINIMIZER ################
            # pos_freqs_for_s_mod = []
            # for j in range(len(PFM)):
            #     v_j = s_alignment_in_matrix[j]
            #     pos_freqs_for_s_mod.append( (j, PFM[j][v_j] / float(max(c[v_j], 1) ) ))
            # pos_freqs_for_s_mod.sort(key=lambda x: x[1]) # sort with respect to smallest frequencies                    
            # pos, highest_freq_of_error_to_correct = pos_freqs_for_s_mod[ nr_pos_to_correct - 1 ]
            # end_position_in_list = nr_pos_to_correct
            # for pp in range(nr_pos_to_correct, len(pos_freqs_for_s_mod)):
            #     # print(pos_freqs_for_s_mod[pp][1], highest_freq_of_error_to_correct)
            #     if pos_freqs_for_s_mod[pp][1] > highest_freq_of_error_to_correct:
            #         break
            #     else:
            #         end_position_in_list += 1
            # J = [j for j, freq in random.sample(pos_freqs_for_s_mod[:end_position_in_list], nr_pos_to_correct)]
            #############################################


            # ####### TEST CHOOSING RANDOM SUBSET OUT OF ALL MINORITY POSITONS IN THE READ ################
            # minority_positions_for_s = []
            # for j in range(len(PFM)):
            #     count_v_j = PFM[j][s_alignment_in_matrix[j]]
            #     max_v_j = max(PFM[j], key = lambda x: PFM[j][x] )
            #     if count_v_j < PFM[j][max_v_j]:
            #         minority_positions_for_s.append(j)
            # # print(len(minority_positions_for_s))

            # if nr_pos_to_correct > len(minority_positions_for_s):
            #     print("OMFG!!", len(minority_positions_for_s), nr_pos_to_correct)
            #     nr_pos_to_correct = len(minority_positions_for_s)
            # J = random.sample(minority_positions_for_s, nr_pos_to_correct)
            ##############################################

            ########### TEST WEIGHTING EACH MINORITY POSITION BY IT'S OBSERVED FREQUENCY THROUGHOUT THE ALIGNMENTS TO THE MINIMIZER ################
            pos_freqs_for_s_mod = []
            for j in range(len(PFM)):
                max_v_j = max(PFM[j], key = lambda x: PFM[j][x] )
                v_j = s_alignment_in_matrix[j]
                if max_v_j == v_j:
                    pos_freqs_for_s_mod.append( (j, PFM[j][v_j] / float(1)) )
                elif max_v_j == "-":
                    pos_freqs_for_s_mod.append( (j, PFM[j][v_j] / float(max(c_ins, 1) ) ))
                elif v_j == "-":
                    pos_freqs_for_s_mod.append( (j, PFM[j][v_j] / float(max(c_del, 1) ) ))
                else:
                    pos_freqs_for_s_mod.append( (j, PFM[j][v_j] / float(max(c_subs, 1) ) ))

            pos_freqs_for_s_mod.sort(key=lambda x: x[1]) # sort with respect to smallest frequencies                    
            pos, highest_freq_of_error_to_correct = pos_freqs_for_s_mod[ nr_pos_to_correct - 1 ]
            end_position_in_list = nr_pos_to_correct
            for pp in range(nr_pos_to_correct, len(pos_freqs_for_s_mod)):
                # print(pos_freqs_for_s_mod[pp][1], highest_freq_of_error_to_correct)
                if pos_freqs_for_s_mod[pp][1] > highest_freq_of_error_to_correct:
                    break
                else:
                    end_position_in_list += 1
            J = [j for j, freq in random.sample(pos_freqs_for_s_mod[:end_position_in_list], nr_pos_to_correct)]
            #############################################

            # J = [j for j, prob in pos_freqs_for_s[:nr_pos_to_correct]] # J is the set of the nr_pos_to_correct smallest position probabilities
            # print(nr_pos_to_correct, end_position_in_list)
            # print(pos_freqs_for_s[:end_position_in_list])
            # print(J)

            s_new = alignment_matrix[s]
            for j in J:
                old_nucl = s_new[j]
                highest_prob_character_at_j = max(PFM[j], key=lambda k: PFM[j][k])

                if highest_prob_character_at_j == old_nucl: # choose the other highest on if tie (should happen only when partition consist of two sequences)
                    pmf_j_minus_variant = copy.deepcopy(PFM[j])
                    del pmf_j_minus_variant[old_nucl] 
                    highest_prob_character_at_j = max(pmf_j_minus_variant, key=lambda k: pmf_j_minus_variant[k])


                # print("correcting", s_new[j], "to", highest_prob_character_at_j )
                s_new[j] = highest_prob_character_at_j
            s_modified = "".join([nucl for nucl in s_new if nucl != "-" ])

            # only unique strings can change in this step

            # accession_of_s = unique_seq_to_acc[s] # this is still unique
            # S_prime_partition[accession_of_s] = s_modified

            accessions_of_s = seq_to_acc[s]
            for acc in accessions_of_s:
                S_prime_partition[acc] = s_modified
    
    return S_prime_partition


def correct_to_consensus(m, partition, seq_to_acc, step):
    S_prime_partition = {}
    N_t = sum([container_tuple[3] for s, container_tuple in partition.items()]) # total number of sequences in partition
    
    if N_t == 2:
        print("Partition has size", N_t, "no meaningful correction can be done")

    if len(partition) > 1 and N_t > 2:
        # all strings has not converged
        alignment_matrix, PFM = create_position_probability_matrix(m, partition) 
        # consensus_alignment = [ max(PFM[j], key=lambda k: PFM[j][k]) for j in range(len(PFM))]
        # print("minimizer errors:",  math.ceil(min([ partition[s][0] for s in partition if partition[s][3] > 1 or s !=m ]) / 2.0)  )
        # frozen_positions = get_frozen_positions(alignment_matrix[m])
        ## TEST LOG ERROR TYPES #######
        c_del = 0
        c_ins = 0
        c_subs = 0
        majority_vector = []
        for j in range(len(PFM)):
            max_v_j = max(PFM[j], key = lambda x: PFM[j][x] )
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

        print("Error types:", c_del, c_ins, c_subs, "depth:", N_t )
        assert len(majority_vector) == len(PFM)
        ############################

        for s in sorted(partition):
            nr_pos_to_correct2 = int(math.ceil(partition[s][0] / 2.0)) #decide how many errors we should correct here
            # print("positions to correct for sequence s:", nr_pos_to_correct, s ==m)
            s_alignment_in_matrix = alignment_matrix[s]
            nr_pos_to_correct = int(math.ceil( len([ 1 for j in range(len(majority_vector)) if (len(majority_vector[j]) == 1 and majority_vector[j] != s_alignment_in_matrix[j] ) ]) * 0.5)) # (step/ float(step +1)) ))
            # print("positions to correct:", nr_pos_to_correct)
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

            ########### TEST WEIGHTING EACH MINORITY POSITION BY IT'S OBSERVED FREQUENCY THROUGHOUT THE ALIGNMENTS TO THE MINIMIZER ################
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
    else:
        print("Partition has converged,", len(partition), N_t)
    
    return S_prime_partition

# def get_frozen_positions(m_in_alignment_matrix):
#     """
#         positions in Multialingment matrix where there are insels longer than 2 bp w.r.t. minimizer, these regions are prone to errors in alingments and
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
