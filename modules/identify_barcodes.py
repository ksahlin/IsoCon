
import os

from modules.input_output import fasta_parser
import edlib_alignment_module


def find_inexact_match_cutpoint(read, barcodes):
    best_mismatches = barcode_length
    for barcode in start_barcode_set:
        ed, locations, cigar = edlib_alignment_module.edlib_traceback(barcode, c_seq, mode="HW", task="path", k=2)
        ed_begin, locations_begin, cigar_begin = edlib_alignment_module.edlib_traceback(barcode, c_seq[:30], mode="HW", task="path", k=5)

        if cigar[-1] == "=":
            last_global_matches = int(cigar[-3:-1])

        if cigar_begin[-1] == "=":
            last_local_matches = int(cigar_begin[-3:-1])

        if 0 <= ed < best_mismatches and last_global_matches > 4: # need to have sufficient last matches to get a well defined cuttinng point
            best_mismatches = ed
            # if len(locations) > 1:
            #     print(ed, locations, cigar)

        elif 0 <= ed_begin < best_mismatches and last_local_matches > 4:
            best_mismatches = ed
            print(ed_begin, locations_begin, cigar_begin)                    
        
    if best_mismatches < barcode_length:
        imperfect_count +=1
        start_barcode_found = True

            # print(ed, locations, cigar)    
    return cut_index


def remove_barcodes_fcn(read_file, params):
    """
                barcodes:   
                    GGTAGGCGCTCTGTGTGCAGC
                    GGTAGTCATGAGTCGACACTA
                    AGAGTACTACATATGAGATGG
                    CGTGTGCATAGATCGCGATGG
                Require at least 15 matches with indels or a stretch of matches of at least 10bp 
    """

    read_file = {acc: seq for (acc, seq) in  fasta_parser.read_fasta(open(read_file, 'r'))}
    barcodes = {acc: seq for (acc, seq) in  fasta_parser.read_fasta(open(params.barcodes, 'r'))}
    barcode_length = 21
    perfect_count = 0
    imperfect_count = 0
    not_found_count = 0
    start_barcode_set = set()
    end_barcode_set = set()
    for acc, seq in barcodes.items():
        if acc[0] == "F":
            start_barcode_set.add(seq)
        elif acc[0] == "R":
            end_barcode_set.add(seq)

    # start_barcode_set2 = set(["GGTAGGCGCTCTGTGTGCAGC", "GGTAGTCATGAGTCGACACTA"])
    # end_barcode_set2 = set(["AGAGTACTACATATGAGATGG", "CGTGTGCATAGATCGCGATGG"])
    # assert start_barcode_set == start_barcode_set2
    # assert end_barcode_set == end_barcode_set2

    min_seq_len = min([len(seq) for seq in read_file.values()])
    print("MIN seq before barcode:", min_seq_len)



    print("total ends to analyze: {0}".format(len(read_file)*2))

    read_file_filtered = {}

    for c_acc in read_file.keys():
        c_seq = read_file[c_acc]
        # 1. find if the start and end snippet is part of a barcode with sw (allow edit distance 2 maybe?).
        #    if so, mask. otherwise add whole snippet to insertion
        # initial mask of primers

        # exact match
        start_barcode_found = False
        end_barcode_found = False



        if c_seq[:barcode_length] in start_barcode_set:
            c_seq = c_seq[barcode_length : ]
            # print("found start primer!")
            start_barcode_found = True
            perfect_count += 1
        else:
            for i in range(100):
                if c_seq[i: barcode_length + i] in  start_barcode_set:
                    print("found start primer further in!", i)
                    c_seq = c_seq[barcode_length + i : ]
                    start_barcode_found = True
                    perfect_count += 1
                    break


        if c_seq[-barcode_length :] in end_barcode_set:
            c_seq = c_seq[ : -barcode_length]
            # print("found end primer!")
            end_barcode_found = True
            perfect_count += 1

        else:
            for i in range(100):
                if c_seq[-(barcode_length + i) : -i] in  end_barcode_set:
                    print("found end primer further in!", i)
                    c_seq = c_seq[ : -(barcode_length + i) : ]
                    end_barcode_found = True
                    perfect_count += 1
                    break


        c_seq_start_cutpoint = 0
        c_seq_end_cutpoint = len(c_seq)

        if not start_barcode_found:
            best_mismatches = barcode_length
            for barcode in start_barcode_set:
                ed, locations, cigar = edlib_alignment_module.edlib_traceback(barcode, c_seq, mode="HW", task="path", k=2)
                ed_begin, locations_begin, cigar_begin = edlib_alignment_module.edlib_traceback(barcode, c_seq[:30], mode="HW", task="path", k=5)

                if cigar[-1] == "=":
                    last_global_matches = int(cigar[-3:-1])

                if cigar_begin[-1] == "=":
                    last_local_matches = int(cigar_begin[-3:-1])

                if 0 <= ed < best_mismatches and last_global_matches > 4: # need to have sufficient last matches to get a well defined cuttinng point
                    best_mismatches = ed
                    # if len(locations) > 1:
                    #     print(ed, locations, cigar)

                elif 0 <= ed_begin < best_mismatches and last_local_matches > 4:
                    best_mismatches = ed
                    print(ed_begin, locations_begin, cigar_begin)                    
                
            if best_mismatches < barcode_length:
                imperfect_count +=1
                start_barcode_found = True

                    # print(ed, locations, cigar)


        if not end_barcode_found:
            best_mismatches = barcode_length
            for barcode in end_barcode_set:
                ed, locations, cigar = edlib_alignment_module.edlib_traceback(barcode, c_seq, mode="HW", task="path", k=2)
                ed_end, locations_end, cigar_end = edlib_alignment_module.edlib_traceback(barcode, c_seq[-30:], mode="HW", task="path", k=5)

                if cigar[] == "=":
                    first_global_matches = int(cigar[-3:-1])

                if cigar_begin[-1] == "=":
                    last_local_matches = int(cigar_begin[-3:-1])

                if 0 <= ed < best_mismatches:
                    best_mismatches = ed
                    # if len(locations) > 1:
                    #     print(ed, locations, cigar)

                elif 0 <= ed_end < best_mismatches:
                    best_mismatches = ed
                    print(ed_end, locations_end, cigar_end)     

            if best_mismatches < barcode_length:
                imperfect_count +=1
                end_barcode_found = True


        if not start_barcode_found:
            for i in range(21, 9, -1):
                c_begin = c_seq[0:i]
                for b_code in start_barcode_set:
                    # print(c_begin,b_code[ barcode_length - len(c_begin) : ])
                    if c_begin == b_code[ barcode_length - len(c_begin) : ]:
                        # print("snippet found!", len(c_begin)) # , b_code[ barcode_length - len(c_begin) : ], c_seq)
                        start_barcode_found = True
                        break   
                if start_barcode_found:
                    imperfect_count +=1
                    break


        if not end_barcode_found:
            for i in range(21, 9, -1):
                c_end = c_seq[-i:]
                for b_code in end_barcode_set:
                    if c_end == b_code[ : barcode_length - len(c_end) ]:
                        # print("snippet found!")
                        end_barcode_found = True
                        break    
                if end_barcode_found:
                    imperfect_count +=1
                    break


        if not end_barcode_found:
            not_found_count += 1

        if not start_barcode_found:
            not_found_count += 1


        c_seq = c_seq[ c_seq_start_cutpoint : c_seq_end_cutpoint ]
        read_file_filtered[c_acc] = c_seq
        # print(len(c_seq))

    print("Perfect:{0}, imperfect:{1}, not found:{2}".format(perfect_count, imperfect_count, not_found_count))
    print("total count: {0}".format(perfect_count + imperfect_count + not_found_count))
    assert len(read_file_filtered) == len(read_file)
    consensus_file = open(os.path.join(params.outfolder, "querys_filtered.fa"), "w")
    for acc, seq in misc_functions.iteritems(read_file_filtered):
        consensus_file.write(">{0}\n{1}\n".format(acc, seq))
    consensus_file.close() 
    
    barcode_corrected_reads_file_name = remove_redundant_after_barcode_removal(read_file_filtered)

    return barcode_corrected_reads_file_name

def remove_redundant_after_barcode_removal(read_file_filtered):
    consensus_file_unique = open(os.path.join(params.outfolder, "querys_unique_after_barcode_removal.fa"), "w")
    read_file_seq_to_acc = {seq : acc for (acc, seq) in  read_file_filtered.items()}
    read_file_unique = {acc : seq for (seq, acc) in  read_file_seq_to_acc.items()}
    min_seq_len = min([len(seq) for seq in read_file_seq_to_acc.keys()])
    print("MIN seq after barcode:", min_seq_len)
    for acc, seq in misc_functions.iteritems(read_file_unique):
        consensus_file_unique.write(">{0}\n{1}\n".format(acc, seq))
    consensus_file_unique.close() 
    return consensus_file_unique.name

