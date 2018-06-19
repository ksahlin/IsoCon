from __future__ import print_function
import sys
import argparse
import re
import itertools
import os
import argparse
import numpy as np
# import misc_functions
import random
from collections import defaultdict
import math
from decimal import * 
getcontext().prec = 100



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

    for read_acc, p_i in probability.items():
        if p_i < 0 or p_i > 1.0:
            print(p_i, read_acc, read_acc in x_equal_to_one)
    # print(sorted(probability.values()))
    assert max(probability.values()) <= 1.0
    assert min(probability.values()) > 0.0
    log_probabilities = { acc: -math.log(p_i, 10) for acc, p_i in probability.items()}
    log_p_i_max = max(log_probabilities.values())
    
    assert log_p_i_max > 0
    weight = {q_acc : log_probabilities[q_acc] / log_p_i_max  for q_acc in log_probabilities.keys()}

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

    return float(raghavan_bound)


def get_error_probabilities(reference_length, read_errors, variants):
    probability = {}
    delta_size = float(len(variants))
    for read_acc in read_errors:
        prob = 1.0
        # print(read_errors)
        (insertions, deletions, substitutions) = read_errors[read_acc]
        p_S = ( max(substitutions, delta_size) / float(reference_length) ) / 3.0   # p = 0.0 not allowed, min_p is 1/(3*len(seq))
        p_I = ( max(insertions, delta_size) / float(reference_length) ) / 4.0   # p = 0.0 not allowed, min_p is 1/(4*len(seq))
        p_D = ( max(deletions, delta_size) / float(reference_length) )         # p = 0.0 not allowed, min_p is 1/(len(seq))

        for i in variants:
            v_type, v_nucl, u_v = variants[i]

            if v_type == "S":
                prob *= p_S*u_v 
            elif v_type == "I":
                prob *= min(0.5, p_I*u_v) 
            elif v_type == "D":
                prob *= min(0.5, p_D*u_v)
        if prob >= 1.0:
            prob = 0.99999
        
        elif prob <= 0.0:
            print("NONPOSITIVE PVAL:",prob, insertions, deletions, substitutions, p_S, p_I, p_D, variants)

        probability[read_acc] = prob

    return probability


def mutate(sequence, positions):
    variants = {p : "" for p in positions}
    new_sequence = list(sequence)
    choices = [("m", 0.3333), ("i", 0.3333), ("d", 0.33333)]
    for p in positions:
        r = random.uniform(0,1)
        if r < 0.3333:
            new_sequence[p] = random.choice( list(set(["A","G","C","T"]) - set(new_sequence[p])) )
            variants[p] = ("S", new_sequence[p], 1)
            assert new_sequence[p] != sequence[p]
        elif 0.3333 <= r < 0.6666:
            new_sequence[p] = new_sequence[p] + random.choice(["A","G","C","T"])
            variants[p] = ("I", new_sequence[p][1], 1)
            assert new_sequence[p] != sequence[p]
        else:
            new_sequence[p] = ""
            variants[p] = ("D", "-", 1)

    return "".join([n for n in new_sequence]), variants


def generate_reads(ref_seq, mutated_ref_seq, params):
    sequence_transcripts = {"ref_seq": ref_seq, "mutated_ref_seq" : mutated_ref_seq}
    reads_supports, read_errors, pass_numbers = set(), {}, []
    # read lengths ~ according to P6-C4 chemistry histogram from here 
    # http://www.slideshare.net/GenomeInABottle/jan2016-pac-bio-giab   slide 13
    # this looks like it can be well approximated by triangiular distributions with parameters
    # 0 (base start), 10000 (peak), ~45000 (base end)
    # http://docs.scipy.org/doc/numpy/reference/generated/numpy.random.triangular.html
    # pacbios own distribution is here:
    # http://www.pacb.com/blog/new-chemistry-boosts-average-read/
    
    # read_lengths = np.random.triangular(0, 10000, 45000, config["read_count"])

    # Get average quality based on subread length and the length of the transcript
    # while read count is less than the red count we want:
    #   1. Draw read length from distribution for each read length in simulated read lengths
    #   2. Randomly select a transcript from pool
    #   3. Get average quality based on the number of passes = floor(read_length/transcript_length)
    #       Avg quality is derived from this plot: https://speakerdeck.com/pacbio/specifics-of-smrt-sequencing-data
    #          slide 21, the P4 C2 chemistry, in case of one pass we chose 13 percent error rate from here : http://www.sciencedirect.com/science/article/pii/S1672022915001345.
    #          Out of the errors we follow the data here: http://bib.oxfordjournals.org/content/17/1/154.full.pdf
    #           and here http://www.homolog.us/Tutorials/index.php?p=2.8&s=2
    #           that indicates that for a pacbio genomic read, we have roughly 13-15 percent error rate (older chemistry)
    #           we choose 13. Out of the total errors, we let 11/16 = 68.75 be insertions
    #           4/16= 25% be deletions and 1/16 = 6.25% be substitutions  (given here http://www.homolog.us/Tutorials/index.php?p=2.8&s=2 and http://bib.oxfordjournals.org/content/17/1/154.full.pdf)
    #   4. generate the read

    quality_function = {1: 0.87, 2:0.95, 3: 0.957, 4:0.969, 5:0.981, 6:0.985, 7:0.99, 
                8:0.992, 9:0.994, 10: 0.995,11: 0.995,12: 0.995, 13: 0.996, 14: 0.996, 15: 0.996,
                16: 0.999, 17: 0.999, 18: 0.999} # passes : quality
    read_count = 1
    # just generate all numbers at once and draw from this 5x should be enough
    it = 0
    lengths = np.random.triangular(0, 10000, 45000, 5*params.read_count)
    ref_choice = np.random.choice(["ref_seq", "mutated_ref_seq" ], size = 5*params.read_count, p = [1.0-params.abundance_ratio, params.abundance_ratio])
    # pacbio_reads = {}
    # reads_generated_log = defaultdict(int)
    # errors = []
    while read_count <= params.read_count:
        if it >= len(lengths):
            lengths = np.random.triangular(0, 10000, 45000, 5*params.read_count)
            it = 0

        read_len = lengths[it]
        acc = ref_choice[it]
        transcript = sequence_transcripts[acc]
        passes =  int(read_len/ len(transcript))
        # print(passes, read_len, len(transcript))
        if passes > 0:
            if passes < 18:
                quality = quality_function[passes]
            else:
                quality = 0.999
            pass_numbers.append(passes)
            subs_rate = (1.0 - quality)*0.0625
            ins_rate = (1.0 - quality)*0.6875
            del_rate = (1.0 - quality)*0.25
            read_errors[it] = (params.gene_length*ins_rate, params.gene_length*del_rate, params.gene_length*subs_rate) 

            if acc == "mutated_ref_seq":
                tot_error_rate =  float(sum(read_errors[it])) / float(params.gene_length)
                # print(np.random.choice([0, 1 ], size =1, p = [tot_error_rate, 1.0-tot_error_rate]))
                is_supporting = np.random.choice([0, 1 ], size =1, p = [tot_error_rate, 1.0-tot_error_rate])[0]
                # print(is_supporting)
                if is_supporting:
                    reads_supports.add(it)
            # read, error_log, total_error_length, total_indel_length, total_del_length = misc_functions.mutate_sequence(transcript, subs_rate, ins_rate, del_rate)
            # read_acc = "{0}_read_{1}_error_rate_{2}_total_errors_{3}".format(acc, str(read_count), total_error_length/float(len(read) + total_del_length), total_error_length)
            # reads_generated_log[acc.split(":copy")[0]] += 1
            # errors.append(total_error_length)
            # pacbio_reads[read_acc] = read

            read_count += 1

        it += 1


    # for acc, abundance in misc_functions.iteritems(reads_generated_log):
    #     params.logfile.write("{0}\t{1}\n".format(acc, abundance))
    # print(read_errors)
    return read_errors, reads_supports, pass_numbers



def main(params):

    # simulate transcript
    # reference = {acc : seq for acc, seq in misc_functions.read_fasta(open(params.transcript,"r"))} 
    reference = { "sim" : "".join([random.choice(["A","G","C","T"]) for i in range(params.gene_length)])}
    ref_seq = list(reference.values())[0]
    # print(params.nr_differences)

    p_values = {}
    divergences = []
    for read_depth in params.read_depths:
        print("RUNNING EXPERIMENTS FOR READ DEPTH:", read_depth )
        p_values[read_depth] = []
        params.read_count = read_depth
        for exp_id in range(params.replicates):
            # simulate  number of mutations
            # given gene_length, similarity, calculate average edit distance (conservatively by rounding down)
            if params.nr_differences:
                nr_mut = params.nr_differences
            else:
                nr_mut = max(1, sum(np.random.choice([1, 0], size = params.gene_length, p = [1.0-params.similarity, params.similarity])))
                # params.nr_differences = max(1, int(params.gene_length * (1.0 -params.similarity) ) )
            divergences.append(nr_mut)
            # print("mutations:", nr_mut)

            while True:
                random_positions = set([random.randint(0, len(ref_seq) - 1) for i in range(nr_mut)] )
                if len(random_positions) == nr_mut:
                    break

            mutated_ref_seq, variants = mutate(ref_seq, random_positions)

            # simulate reads
            read_errors, reads_supports, pass_numbers = generate_reads(ref_seq, mutated_ref_seq, params)
            if not reads_supports:
                print("no support for read depth", read_depth)
                p_value = 1.0
            else:
                probabilities =  get_error_probabilities(params.gene_length, read_errors, variants)
                # print("supports", len(reads_supports))
                p_value = raghavan_upper_pvalue_bound(probabilities, reads_supports)

            p_values[read_depth].append(p_value)
            # print("p-val:", p_value)

    # print output
    print("Mean divergence of simulated gene copies: {0}nt.".format( sum(divergences)/len(divergences) ) )
    print("Median divergence of simulated gene copies: {0}nt. ".format(sorted(divergences)[len(divergences)/2]) )
    print("Estimated mean passes of reads: {0}.".format( sum(pass_numbers)/len(pass_numbers) ) )
    print("Estimated median passes of reads: {0}.".format(sorted(pass_numbers)[len(pass_numbers)/2]) )

    # recall rates
    for read_depth in params.read_depths:
        print("Median p-value for the transcript: {0}, at read depth {1}. ".format(read_depth, sorted(p_values[read_depth])[len(p_values[read_depth])/2]) )
        print("Mean p-value for the transcript: {0}, at read depth {1}. ".format(read_depth, sum(p_values[read_depth])/len(p_values[read_depth]) ) )
        print("90% upper quantile p-value for the transcript: {0}, at read depth {1}. ".format(read_depth, sorted(p_values[read_depth])[int(len(p_values[read_depth])*0.9)] ) )
        print()


def mkdir_p(path):
    try:
        os.makedirs(path)
        print("creating", path)
    except OSError as exc:  # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Generate a gene family and isoforms from a set of original exons.")
    # group = parser.add_mutually_exclusive_group(required=True)
    # group.add_argument('--alpha', type=float, default = 0.01, help='The p-value cutoff for IsoCon algorithm (default 0.01).')
    parser.add_argument('--gene_length', type=int, help='Estimated gene length for the gene families(s). The longer the estimate the more conservative the coverage value. ')
    parser.add_argument('--abundance_ratio', type=float, help='Sample abundance of the transcript')
    parser.add_argument('--similarity', type=float, default = None, help='Average identity.')
    parser.add_argument('--read_depths', nargs="+", type=int, default = [100,500, 1000, 2000, 3000, 5000, 10000, 20000], help='Different read depths to do calculations over.')
    parser.add_argument('--nr_differences', type=int, default=None, help='Integer. Exact number of mismatches between the gene copies. This will override extimation of gene divergence obtained from the similarity and gene length parameters.')
    parser.add_argument('--replicates', type=int, default=100, help='Number od simulated replicates for calculating abundance.')
    # parser.add_argument('reads_outfile', type=str, help='Generated reads output file')
    # parser.add_argument('ref_outfile', type=str, help='generated ref output file ')
    params = parser.parse_args()
    # print(params.similarity)
    if params.similarity:
        if  not (0.8 < params.similarity < 1):
            print( "Simularity has to be between 0.8 and 1.0")
            assert 0.8 < params.similarity < 1

    # path_, file_prefix = os.path.split(params.reads_outfile)
    # misc_functions.mkdir_p(path_)
    main(params)

