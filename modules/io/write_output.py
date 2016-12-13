import os
from time import time
import datetime
from scipy.stats import poisson, binom
from hitemmodules import align

def print_filtered_data(reads, params):
    for x_i, seq in reads.items():
        params.filtered_reads.write(">{0}\n{1}\n".format(x_i, seq))

def isolated_reads(reads_not_observed_in_alignment, x, params):
    # remove the reads and corresponding y's that did not have any alignments
    isolated_reads_fasta_file = open(os.path.join(params.outfolder, "isolated_reads.fa"), 'w')
    isolated_count = 0
    for x_i in x.keys():
        if x_i in reads_not_observed_in_alignment:
            isolated_count += 1
            isolated_reads_fasta_file.write(">{0}\n{1}\n".format(x_i, x[x_i]))    
            # if x_i == "m151209_004839_42146_c100926392550000001823199905121692_s1_p0/3554/1665_60_CCS":
            #     print("XXX no alignments for m151209_004839_42146_c100926392550000001823199905121692_s1_p0/3554/1665_60_CCS")

    logger("Number isolated reads (not observed in PAF file: {0}".format(isolated_count), params.logfile)

def post_process(y_sampled, y_to_x, epsilon_y, errors_y, params, final_output = False):
    # Postprocess
    # merge identical y
    transcripts_dict = {} #seq -> [read_acc,...]
    for acc, seq in y_sampled.items():
        if seq not in transcripts_dict:
            transcripts_dict[seq] = [y_to_x[acc]]
        else:
            transcripts_dict[seq].append(y_to_x[acc])

    for seq, read_accs in transcripts_dict.items():
        print("consensus of {0} reads".format(len(read_accs)), read_accs)

    unique_consensus = [ x for x in  transcripts_dict.values() if len(x) <=1 ] #filter(lambda x: len(x)<=1, transcripts_dict.values() )
    print()
    print("Total of {0} consensus reads are formed. {1} of them converged to a unique transcript (this is a bad sign, they should be removed)".format(len(transcripts_dict),len(unique_consensus)))
    print()
    consensus_fasta_files = []
    if final_output:
        for support in range(1, params.min_consensus_support +1):
            consensus_fasta_files.append( open(os.path.join(params.outfolder, "consensus_min_support_{0}_final.fa".format(support)), 'w') )
    else:
        for support in range(1, params.min_consensus_support +1):
            consensus_fasta_files.append( open(os.path.join(params.outfolder, "consensus_min_support_{0}_step_{1}.fa".format(support, str(params.step))), 'w') )

    for number, (seq, read_accs) in enumerate(transcripts_dict.items()):
        for support in range(1, params.min_consensus_support +1):
            if len(read_accs) >= support:
                # if consensus has more support than required by the file on position support-1 in the list with file objects, it is written to that file
                consensus_fasta_files[support-1].write(">consensus_{0}_from_{1}_reads\n{2}\n".format(number, len(read_accs), seq))


    if final_output:    
        consensus_fasta_file_low_epsilon = open(os.path.join(params.outfolder, "consensus_inferred_low_epsilon_final.fa"), 'w')
    else:
        consensus_fasta_file_low_epsilon = open(os.path.join(params.outfolder, "consensus_inferred_low_epsilon_step_{0}.fa".format(str(params.step))), 'w')

    # y_errors_left = set(epsilon_y.keys())
    # y_sampled_left = set(y_sampled.keys())
    # epsilon_minus_sampled = y_errors_left.difference(y_sampled_left)
    # sampled_minus_epsilon = y_sampled_left.difference(y_errors_left)

    # for y_j in epsilon_minus_sampled:
    #     print("error minus sampled", y_j, epsilon_y[y_j],y_to_x[y_j])
    # for y_j in sampled_minus_epsilon:
    #     print("sampled minus error", y_j, y_sampled[y_j],y_to_x[y_j])


    transcripts_dict = {} #seq -> [read_acc,...]
    for y_j, error_rate in sorted(epsilon_y.items(), key=lambda x: x[1]):
        # print(y_j,epsilon_y[y_j],y_to_x[y_j])
        seq = y_sampled[y_j]
        if seq not in transcripts_dict:
            if epsilon_y[y_j] < 1.0/ len(seq):
                transcripts_dict[seq] = [y_j]
        else:
            transcripts_dict[seq].append(y_j)

    # print(transcripts_dict)
    average_error = {}
    for seq in transcripts_dict:
        error_rates = []
        for y_j in transcripts_dict[seq]:
            error_rates.append(epsilon_y[y_j])
        # print("ERROR RATES:", error_rates)
        avg_error_rate = sum(error_rates)/len(error_rates)
        # print("average:", avg_error_rate)
        average_error[seq] = avg_error_rate

    error_free_divergent_consensus = 0
    error_free_count = 0
    for number, (seq, y_js) in enumerate(transcripts_dict.items()):
        if len(y_js) > 2 and average_error[seq] < 1.0/len(seq):
            consensus_fasta_file_low_epsilon.write(">consensus_{0}_from_{1}_reads_avg_error_rate:{3}\n{2}\n".format(number, len(y_js), seq, average_error[seq]))
            error_free_count += 1
        elif len(y_js) == 1:
            error_free_divergent_consensus += 1
    print()
    print("{0} inferred low epsilon consensus reads are formed. {1} divergent consensus with low inferred error rate are formed.".format(error_free_count, error_free_divergent_consensus))
    print()


    if final_output:    
        consensus_fasta_file_errorfree = open(os.path.join(params.outfolder, "consensus_inferred_errorfree_final.fa"), 'w')
    else:
        consensus_fasta_file_errorfree = open(os.path.join(params.outfolder, "consensus_inferred_errorfree_step_{0}.fa".format(str(params.step))), 'w')

    # y_errors_left = set(errors_y.keys())
    # y_sampled_left = set(y_sampled.keys())
    # error_minus_sampled = y_errors_left.difference(y_sampled_left)
    # sampled_minus_error = y_sampled_left.difference(y_errors_left)

    # for y_j in error_minus_sampled:
    #     print("error minus sampled", y_j, errors_y[y_j],y_to_x[y_j])
    # for y_j in sampled_minus_error:
    #     print("sampled minus error", y_j, y_sampled[y_j],y_to_x[y_j])


    transcripts_dict = {} #seq -> [read_acc,...]
    for y_j, error_rate in sorted(errors_y.items(), key=lambda x: x[1]):
        # print(y_j,errors_y[y_j],y_to_x[y_j])
        seq = y_sampled[y_j]
        if seq not in transcripts_dict:
            if errors_y[y_j] == 0:
                transcripts_dict[seq] = [y_j]
        else:
            transcripts_dict[seq].append(y_j)

    # print(transcripts_dict)
    average_error = {}
    for seq in transcripts_dict:
        error_rates = []
        for y_j in transcripts_dict[seq]:
            error_rates.append(errors_y[y_j])
        # print("ERROR RATES:", error_rates)
        avg_errors = float(sum(error_rates))/len(error_rates)
        # print("average:", avg_errors)
        average_error[seq] = avg_errors

    error_free_divergent_consensus = 0
    error_free_count = 0
    for number, (seq, y_js) in enumerate(transcripts_dict.items()):
        if len(y_js) > 2 and average_error[seq] < 0.5:
            consensus_fasta_file_errorfree.write(">consensus_{0}_from_{1}_reads_avg_errors:{3}\n{2}\n".format(number, len(y_js), seq, average_error[seq]))
            error_free_count += 1
        elif len(y_js) == 1:
            error_free_divergent_consensus += 1
    print()
    print("{0} inferred errorfree consensus reads are formed. {1} divergent consensus with low inferred error rate are formed.".format(error_free_count, error_free_divergent_consensus))
    print()

def output_significant_clusters(y, x, alignment_results_dict, y_to_x, epsilon_x, params):
    """
        probabilistic method to filter out insignificant consensus. When coverage is extremely high, we 
        need a statistical method to choose if a mutation(s)/SNP(s) is true based on significant abundance.
        Fur ultra high coverage, several identical errors can uccur over same positions.
        This is a function of error rates of reads epsilon_i (estimated accurately), 
        how many reads come from the same transcript E_c (unknown but estimated approximately)
        and the consensus cluster size c_n --- supporting the same mutation/error (known).

        1. Get closest similarity between each pair of consensus (by self mappings) - factored into number of subs, ins and del

        2. Calculate probblility that n_c reads would have e_s, e_i and e_d identical subs, ins and del errors, given error rate of reads e_i + common (the hidden) errors. 
            This can be approximated by a poisson as follows.
                1.1. First consider we have exactly one deletion in common:
                    The probability of observing a deletion at a given base pair for a read x_i is epsilon_i + 1/len(x_i) (last part is the common deletion). Given x_1, ...x_n from
                    the same transcript, the probability of the number of errors over base pairs is a sum of bernoulli variables with probabilities epsilon_1 + 1/len(x_1), ..., epsilon_n_c + 1/len(x_n_c) 
                    we can approximate this distribution by a Poisson distribution with lambda = epsilon_1 + 1/len(x_1) + ... + epsilon_n_c + 1/len(x_n_c)
                    this is a good approximation if all epsilon are small http://stats.stackexchange.com/questions/93852/sum-of-bernoulli-variables-with-different-success-probabilities
                    reference to http://projecteuclid.org/euclid.aoms/1177705799. One can show that this approximates the probability
                    better than 9*max(epsilon_i + + 1/len(x_i)), i \in{1,n_c} (section 3 remark (iii) in the paper). (normal approximation can be used if all epsilon are big). 
                1.2 For an substitution (or insertion), the same method as in 1.1 is used, but all epsilon are divided by 3 (or 4) as the introduced error also need to
                    be the same.
                1.3 For several insertions,substitutions or deletions, as errors are independent, simply multiply the probabilities to get the total probability. 
        3. Output the clusters with significant values p < min(0.01, 1.0/nr_of_consensus) # also one can define set this as parameter
    """ 

    transcripts_dict = {} #seq -> [read_acc,...]
    for y_j, seq in y.items():
        if seq not in transcripts_dict:
            transcripts_dict[seq] = [y_to_x[y_j]]
        else:
            transcripts_dict[seq].append(y_to_x[y_j])

    significant_clusters = {} # format seq : (n_c, p_value)
    all_clusters = {}
    # all y's left are guaranteed to have their starting read x_j left

    y_to_y_to_map = {}
    alignment_results_dict_y_to_y = {}

    ########### GET CLOSEST SIMILARITY #################################
    ####################################################################
    ####################################################################

    for seq_1, x_i_list_1 in transcripts_dict.items():
        n_c = len(x_i_list_1)
        # already_considered_y.add(seq_1)
        if n_c < 2:
            # trivial lower bound due to method,
            # We cannot allow convergence to unique y or
            continue

        # TODO: eventually reduce number of alignments by not considering y's differing with more than
        # an inferred dynamic threshold instead of just '20', this depends on lambda and cluser size
        # prepare all relevant pairwise alignments
        print("consensus of {0} reads".format(len(x_i_list_1)), x_i_list_1)
        y_to_y_to_map[seq_1] = {}
        alignment_results_dict_y_to_y[seq_1] = {}
        for seq_2, x_i_list_2 in transcripts_dict.items():
            # we should implement "if seq_2 not in already_considered_y" to skip doing double work 
            if seq_1 != seq_2 and abs(len(seq_1) - len(seq_2)) < 20:
                # dummy initialize alignments that should be made
                y_to_y_to_map[seq_1][seq_2] = 0
                alignment_results_dict_y_to_y[seq_1][seq_2] = 0

        if not y_to_y_to_map[seq_1]:
            significant_clusters[seq_1] = (n_c, "NA")
            all_clusters[seq_1] = (n_c, -1)
            del alignment_results_dict_y_to_y[seq_1]
            # del y_to_y_to_map[seq_1]

    seq_to_seq = dict([(seq,seq) for seq in transcripts_dict])
    align.sw_align_sequences(y_to_y_to_map, alignment_results_dict_y_to_y, seq_to_seq, seq_to_seq, params)

    for seq_1 in list(alignment_results_dict_y_to_y):
        if not alignment_results_dict_y_to_y[seq_1]:
            assert seq_1 not in significant_clusters
            print("No good alignments:", seq_1)
            n_c = len(transcripts_dict[seq_1])
            significant_clusters[seq_1] = (n_c, "NA")
            all_clusters[seq_1] = (n_c, -1)
            del alignment_results_dict_y_to_y[seq_1]


    # all consensus have either a best alignment for the alignments step or they have a p-value NA
    # in case there are no similar consensus: in this case we cannot get a p-value for the cluster.

    consensus_fasta_file_significant_clusters = open(os.path.join(params.outfolder, "consensus_significant.fa"), 'w')
    consensus_fasta_file_insignificant_clusters = open(os.path.join(params.outfolder, "consensus_insignificant.fa"), 'w')
    
    ####################################################################
    ####################################################################
    ####################################################################

    # the consensus that didn't have a match
    for i, (seq,(n_c, p_val)) in enumerate(significant_clusters.items()):
        consensus_fasta_file_significant_clusters.write(">consensus_{0}_from_{1}_reads_pval_{2}\n{3}\n".format(i, n_c, p_val, seq))


    ####################################################################
    #### CALCULATE SIGNIFICANCE OF CONVERGENCE GIVEN CLOSEST MATCH #####
    ####################################################################

    nr_significant_clusters = len(significant_clusters)
    nr_insignificant_clusters = 0 
    total_nr_of_tested_clusters = len(significant_clusters) + len(alignment_results_dict_y_to_y)
    for i, y1 in enumerate(alignment_results_dict_y_to_y):
        almnts = alignment_results_dict_y_to_y[y1].items()
        # print("Y:", y1, "len total alignment_results_dict_y_to_y:", len(alignment_results_dict_y_to_y))
        # print(almnts)
        best_y_to_y1 = sorted(almnts, key = lambda z: z[1][2][1] + z[1][2][2])[0][0] 
        y1_aln, y2_aln, (matches, mismatches, indels) = alignment_results_dict_y_to_y[y1][best_y_to_y1]
        deletions = y1_aln.count("-")
        insertions = y2_aln.count("-")
        nr_not_corrected_errors_hypothesis = mismatches + deletions + insertions

        x_i_list_for_y1 = transcripts_dict[y1]
        n_c = len(x_i_list_for_y1)
        # get expected cluster size
        E_c_sum = 0
        for x_i in x_i_list_for_y1:
            E_c_sum += params.expected_support[x_i]

        # For highly abundant (targeted) sequences, a lot of the more error prone reads could 
        # be filtered away after scanning paf file, therefore, total number of remaining y after filtration 
        # might be smaller than expected support by scanning paf, whichever is the smallest serves as the 
        # esitmation of cluster size. That is, it does not make sense to have an inferred number of cluster size
        # that is larrger than size of all staring sequences, because the removed sequences will not contribute to any
        # cluster and can therefore also not contirbute to any accumulation of errors. So these should not be modeled. 

        E_c = E_c_sum / float(n_c) # min(E_c_sum / n_c, len(y))

        # get expected error rate
        e_i_sum = 0
        for x_i in x_i_list_for_y1:
            e_i = (epsilon_x[x_i]*len(x[x_i]) + nr_not_corrected_errors_hypothesis) / float(len(x[x_i]))
            e_i_sum += e_i
        e_average = e_i_sum / float(n_c)

        # simplify by assuming that insertions, deletions and substitutions have the same frequency in CCS reads
        # although this is definitely not the case in original pacbio reads, after applying the CCS protocol,
        # this is "somewhat" the case in CCS reads. Supporting data for this claim Suppl Figure2 in: 
        # http://www.nature.com/article-assets/npg/ncomms/2016/160624/ncomms11706/extref/ncomms11706-s1.pdf  
        # insertion still little bit overrepresented.
        lambda_poisson = (e_average/3.0) * E_c
        # need  -0.0001 in sf because sf is P( x > k ), therefore, the value n_c doesn't contribute 
        # for example sf(0, lambda) will not return 1 because P( x > 0 ) but sf(0 - 0.001, lambda) = 1 
        # because  P( x > -0.001 ) = 1, this holds for all integer values in the poisson distr.
        print("Expected support: {0}, cluster size: {1}, average error for each type(subs,ins,del):{2}, lambda: {3}".format(E_c, n_c, e_average/3.0, lambda_poisson))

        # the probability p that we get here from the survival function is the probability that we see the cluster size over one position
        # The probability that we see this cluster size on k occurences over n base pairs can be seen as three binomial distributions
        # with Bin(n, p_1) so the final p-value is P(x > k_mism ) where n= len(consensus) and k = nr_ins
        p_1 = poisson.sf(n_c - 0.0001, lambda_poisson/3) #**mismatches 
        p_val_mis = binom.sf(mismatches - 0.0001, len(y1), p_1)
        #get insertions here
        p_2 = poisson.sf(n_c - 0.0001, lambda_poisson/4) #**insertions
        p_val_ins = binom.sf(insertions - 0.0001, len(y1), p_2)        
        #get deletions here
        p_3 = poisson.sf(n_c - 0.0001, lambda_poisson) #**deletions
        p_val_del = binom.sf(deletions - 0.0001, len(y1), p_3)        

        p_val = p_val_mis*p_val_ins*p_val_del
        print(mismatches, insertions, deletions, p_1, p_2, p_3, p_val, len(y1))
        if p_val < 0.01/float(total_nr_of_tested_clusters):
            consensus_fasta_file_significant_clusters.write(">consensus_{0}_from_{1}_reads_pval_{2}\n{3}\n".format(i, n_c, p_val, y1))
            nr_significant_clusters += 1
            all_clusters[y1] = (n_c, p_val)
        else:
            consensus_fasta_file_insignificant_clusters.write(">consensus_{0}_from_{1}_reads_pval_{2}\n{3}\n".format(i, n_c, p_val, y1))
            nr_insignificant_clusters += 1
            all_clusters[y1] = (n_c, p_val)


    print()
    print("{0} significant clusters are formed. {1} insignificant clusters are formed.".format(nr_significant_clusters, nr_insignificant_clusters))
    print()

    # print a file with all clusters sorted from most significant. 
    # -1 are in the top because they are probably the most significant according to our
    # test (they have such a high difference to any other transcript)
    consensus_fasta_file_all_clusters = open(os.path.join(params.outfolder, "consensus_all_sorted_pval.fa"), 'w')

    for i, (seq, (n_c, p_val)) in enumerate(sorted(list(all_clusters.items()), key=lambda x: x[1][1])):
        consensus_fasta_file_all_clusters.write(">consensus_{0}_from_{1}_reads_pval_{2}\n{3}\n".format(i, n_c, p_val, seq))


def logger(message, logfile, timestamp=True):
    if timestamp:
        currrent_time = datetime.datetime.now()
        logfile.write(str(currrent_time) + "\t" + message + "\n")
    else:
        logfile.write(message + "\n")

