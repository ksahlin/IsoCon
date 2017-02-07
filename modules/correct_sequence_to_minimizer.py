import signal
from multiprocessing import Pool
import multiprocessing as mp
import sys
import math
import copy

from modules.functions import create_position_probability_matrix

def correct_strings(partition_alignments, unique_seq_to_acc, single_core = False):
    S_prime = {}

    partition_unique_seq_to_acc = {}
    for m, partition in partition_alignments.items():
        partition_unique_seq_to_acc[m] = {}
        for s in partition:
            if s in unique_seq_to_acc:
                partition_unique_seq_to_acc[m][s] = unique_seq_to_acc[s]


    if single_core:
        for m, partition in partition_alignments.items():
            S_prime_partition = correct_to_minimizer(m, partition, unique_seq_to_acc)
        for acc, s in S_prime_partition.items():
            S_prime[acc] = s
    else:
        ####### parallelize statistical tests #########
        # pool = Pool(processes=mp.cpu_count())
        original_sigint_handler = signal.signal(signal.SIGINT, signal.SIG_IGN)
        signal.signal(signal.SIGINT, original_sigint_handler)
        pool = Pool(processes=mp.cpu_count())
        try:
            res = pool.map_async(correct_to_minimzer_helper, [ ( (m, partition, partition_unique_seq_to_acc[m]), {}) for m, partition in partition_alignments.items() if len(partition) > 1 ] )
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

def correct_to_minimizer(m, partition, unique_seq_to_acc):
    S_prime_partition = {}

    N_t = sum([container_tuple[3] for s, container_tuple in partition.items()]) # total number of sequences in partition

    if len(partition) > 1:
        # all strings has not converged
        alignment_matrix, PFM = create_position_probability_matrix(m, partition) 

        for s in partition:
            nr_pos_to_correct = int(math.ceil(partition[s][0] / 2.0)) #decide how many errors we should correct here
            # print("positions to correct for sequence s:", nr_pos_to_correct, s ==m)
            if nr_pos_to_correct  == 0:
                continue

            s_alignment_in_matrix = alignment_matrix[s]
            # find the position probabilities of the alignment of s in PFM
            pos_freqs_for_s = []
            for j in range(len(PFM)):
                # try:
                pos_freqs_for_s.append( (j, PFM[j][s_alignment_in_matrix[j]]) )
                # except KeyError:
                #     print(j, PFM[j], s_alignment_in_matrix[j], N_t, len(partition), len(PFM), len(m) )
                #     sys.exit()

            pos_freqs_for_s.sort(key=lambda x: x[1]) # sort with respect to smallest frequencies                    
            J = [j for j, prob in pos_freqs_for_s[:nr_pos_to_correct]] # J is the set of the nr_pos_to_correct smallest position probabilities
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

            accession_of_s = unique_seq_to_acc[s] # this is still unique
            S_prime_partition[accession_of_s] = s_modified
    return S_prime_partition
