
import os

from modules.input_output import fasta_parser
from modules import edlib_alignment_module

def reverse_complement(string):
    #rev_nuc = {'A':'T', 'C':'G', 'G':'C', 'T':'A', 'N':'N', 'X':'X'}
    # Modified for Abyss output
    rev_nuc = {'A':'T', 'C':'G', 'G':'C', 'T':'A', 'a':'t', 'c':'g', 'g':'c', 't':'a', 'N':'N', 'X':'X', 'n':'n', 'Y':'R', 'R':'Y', 'K':'M', 'M':'K', 'S':'S', 'W':'W', 'B':'V', 'V':'B', 'H':'D', 'D':'H', 'y':'r', 'r':'y', 'k':'m', 'm':'k', 's':'s', 'w':'w', 'b':'v', 'v':'b', 'h':'d', 'd':'h'}

    rev_comp = ''.join([rev_nuc[nucl] for nucl in reversed(string)])
    return(rev_comp)

def find_cutpoint_exact_match(read, barcode_set, barcode_length):
    """
        cut_point = 10 means cutting the first 10 bases 0 to 9, i.e. keeping [10:len(read)] of a string, because strings are 0 indexed.
    """
    match_barcode = None
    cut_point = 0
    if read[:barcode_length] in barcode_set:
        cut_point = barcode_length
        match_barcode = read[:barcode_length]
    else:
        for i in range(50):
            if read[i: barcode_length + i] in barcode_set:
                # print("found start primer further in!", i)
                cut_point = barcode_length + i
                match_barcode = read[i: barcode_length + i]
                break

    if not cut_point:
        barcode_found = False
        for i in range(barcode_length, 8, -1):
            read_beginning = read[0:i]
            for b_code in barcode_set:
                # print(read_beginning,b_code[ barcode_length - len(read_beginning) : ])
                if read_beginning == b_code[ barcode_length - len(read_beginning) : ]:
                    cut_point = i
                    # print("snippet found!", len(read_beginning)) # , b_code[ barcode_length - len(read_beginning) : ], c_seq)
                    match_barcode = b_code
                    barcode_found = True
                    break   
            if barcode_found:
                break

    return cut_point, match_barcode


def find_cutpoint_inexact_match(read, barcode_set, barcode_length):
    match_barcode = None
    cut_point = 0
    best_mismatches = barcode_length
    for barcode in sorted(barcode_set):
        ed, locations, cigar = edlib_alignment_module.edlib_traceback(barcode, read[:30], mode="HW", task="path", k=3)
        # ed_begin, locations_begin, cigar_begin = edlib_alignment_module.edlib_traceback(barcode, read[:30], mode="HW", task="path", k=5)
        # last_global_matches, last_local_matches = 0,0

        # if cigar and cigar[-1] == "=":
        #     try:
        #         last_global_matches = int(cigar[-3:-1])
        #     except ValueError:
        #         last_global_matches = int(cigar[-2:-1])

        # if cigar_begin and cigar_begin[-1] == "=":
        #     try:
        #         # print(ed_begin, locations_begin,cigar)
        #         last_local_matches = int(cigar_begin[-3:-1])
        #     except ValueError:
        #         last_local_matches = int(cigar_begin[-2:-1])
        # print(ed, locations, cigar)
        if 0 <= ed < best_mismatches: # and last_global_matches > 4: # need to have sufficient last matches to get a well defined cuttinng point
            best_mismatches = ed
            cut_point = locations[-1][1]
            match_barcode = barcode
            # if len(locations) > 1:
            #     print(ed, locations, cigar)

        # elif 0 <= ed_begin < best_mismatches and last_local_matches > 4:
        #     best_mismatches = ed
        #     cut_point = locations_begin[-1][1]
        #     match_barcode = barcode
            # print(ed_begin, locations_begin, cigar_begin)

    return cut_point, match_barcode


def remove_barcodes_fcn(read_file, params):
    """
                barcodes:   
                    GGTAGGCGCTCTGTGTGCAGC
                    GGTAGTCATGAGTCGACACTA
                    AGAGTACTACATATGAGATGG
                    CGTGTGCATAGATCGCGATGG
                Require at least 15 matches with indels or a stretch of matches of at least 10bp 
    """

    reads = {acc: seq for (acc, seq) in  fasta_parser.read_fasta(open(read_file, 'r'))}
    barcodes = {acc: seq for (acc, seq) in  fasta_parser.read_fasta(open(params.barcodes, 'r'))}
    barcode_length = 21

    start_barcode_set = set()
    end_barcode_set = set()
    for acc, seq in barcodes.items():
        if acc[0] == "F":
            start_barcode_set.add(seq)
        elif acc[0] == "R":
            seq_r = reverse_complement(seq)
            end_barcode_set.add(seq_r)

    min_seq_len = min([len(seq) for seq in reads.values()])
    reads_before = set([seq for seq in reads.values()])
    print("Nr unique sequences before barcode correction:", len(reads_before))
    print("MIN seq before barcode:", min_seq_len)



    print("total ends to analyze: {0}".format(len(reads)*2))

    reads_filtered = {}
    perfect_start_count = 0
    perfect_end_count = 0
    imperfect_snippet = 0
    imperfect_inexact = 0
    not_found_count = 0

    for c_acc in reads.keys():
        c_seq = reads[c_acc]

        # start of sequence
        start_cut_point, start_barcode_match = find_cutpoint_exact_match(c_seq, start_barcode_set, barcode_length)
        if 21 <= start_cut_point:
            perfect_start_count += 1
        elif 0 < start_cut_point < 21:
            imperfect_snippet += 1
        else:
            start_cut_point, start_barcode_match = find_cutpoint_inexact_match(c_seq, start_barcode_set, barcode_length)

            if 0 < start_cut_point:
                imperfect_inexact += 1
                # print(len(c_seq))

        # end of sequence
        c_seq_rc = reverse_complement(c_seq)
        end_cut_point, end_barcode_match = find_cutpoint_exact_match(c_seq_rc, end_barcode_set, barcode_length)
        if c_acc == "m151209_215352_42146_c100926392550000001823199905121696_s1_p0/137971/0_2460_CCS":
            print(start_cut_point, end_cut_point, len(c_seq))
        if 21 <= end_cut_point:
            perfect_end_count += 1
        elif 0 < end_cut_point < 21:
            imperfect_snippet += 1
        else:
            end_cut_point, end_barcode_match = find_cutpoint_inexact_match(c_seq_rc, end_barcode_set, barcode_length)

            if 0 < end_cut_point:
                imperfect_inexact += 1
                # print(len(c_seq))

        # trim read and add error corrected adapters
        if end_cut_point:
            c_seq_rc = end_barcode_match + c_seq_rc[ end_cut_point : ]
        else:
            not_found_count +=1
        
        cleaned_read = reverse_complement(c_seq_rc)

        if start_cut_point:
            cleaned_read = start_barcode_match + cleaned_read[ start_cut_point : ]
        else:
            not_found_count +=1

        reads_filtered[c_acc] = cleaned_read

    print("Perfect starts:{0}, Perfect ends:{1} imperfect snippets:{2}, imperfect inexact:{3} not found:{4}".format(perfect_start_count, perfect_end_count, imperfect_snippet, imperfect_inexact, not_found_count))
    print("total count: {0}".format(perfect_start_count + perfect_end_count + imperfect_snippet + imperfect_inexact + not_found_count))

    reads_after = set([seq for seq in reads_filtered.values()])
    print("Nr unique sequences after barcode correction:", len(reads_after))

    assert len(reads_filtered) == len(reads)
    corrected_read_file = open(os.path.join(params.outfolder, "reads_barcode_corrected.fa"), "w")
    for acc, seq in list(reads_filtered.items()):
        corrected_read_file.write(">{0}\n{1}\n".format(acc, seq))
    corrected_read_file.close() 
    
    # barcode_corrected_reads_file_name = remove_redundant_after_barcode_removal(read_file_filtered)
    return corrected_read_file.name



# def remove_redundant_after_barcode_removal(read_file_filtered):
#     consensus_file_unique = open(os.path.join(params.outfolder, "querys_unique_after_barcode_removal.fa"), "w")
#     read_file_seq_to_acc = {seq : acc for (acc, seq) in  read_file_filtered.items()}
#     read_file_unique = {acc : seq for (seq, acc) in  read_file_seq_to_acc.items()}
#     min_seq_len = min([len(seq) for seq in read_file_seq_to_acc.keys()])
#     print("MIN seq after barcode:", min_seq_len)
#     for acc, seq in misc_functions.iteritems(read_file_unique):
#         consensus_file_unique.write(">{0}\n{1}\n".format(acc, seq))
#     consensus_file_unique.close() 
#     return consensus_file_unique.name

