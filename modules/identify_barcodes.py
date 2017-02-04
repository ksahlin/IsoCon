

def remove_barcodes_fcn(args):
    """
                barcodes:   
                    GGTAGGCGCTCTGTGTGCAGC
                    GGTAGTCATGAGTCGACACTA
                    AGAGTACTACATATGAGATGG
                    CGTGTGCATAGATCGCGATGG
                Require at least 15 matches with indels or a stretch of matches of at least 10bp 
    """

    # test to see random score
    # for i in range(100):
    #     seq1 = "".join([random.choice("AGCT") for i in range(21)])
    #     seq2 = "".join([random.choice("AGCT") for i in range(21)])
    #     barcode_aligned, start_snippet_aligned, start_snippet_start_snippet, start_snippet_end_snippet, matches, mismatches, deletions, insertions = ssw_alignment( "c_acc", "barcode", seq1, seq2, 2, 2)
    #     print(matches)


    consensus_transcripts = {acc: seq for (acc, seq) in  read_fasta(open(args.consensus_transcripts, 'r'))}

    min_seq_len = min([len(seq) for seq in consensus_transcripts.values()])
    print("MIN seq before barcode:", min_seq_len)

    # remove_redundant_after_barcode_removal(consensus_transcripts)
    # sys.exit()

    barcode_length = 21
    perfect_count = 0
    imperfect_count = 0
    not_found_count = 0
    start_barcode_set = set(["GGTAGGCGCTCTGTGTGCAGC", "GGTAGTCATGAGTCGACACTA"])
    end_barcode_set = set(["AGAGTACTACATATGAGATGG", "CGTGTGCATAGATCGCGATGG"])

    print("total ends to analyze: {0}".format(len(consensus_transcripts)*2))

    consensus_transcripts_filtered = {}

    for c_acc in consensus_transcripts.keys():
        c_seq = consensus_transcripts[c_acc]
        # 1. find if the start and end snippet is part of a barcode with sw (allow edit distance 2 maybe?).
        #    if so, mask. otherwise add whole snippet to insertion
        # initial mask of primers

        # exact match
        start_barcode_found = False
        end_barcode_found = False

        # if "CCATCGCGATCTATGCACACG" in consensus_transcripts[c_acc][:barcode_length] or "CCATCGCGATCTATGCACACG" in consensus_transcripts[c_acc][ : -barcode_length]:
        #     print("OMG")
        # if "CCATCTCATATGTAGTACTCT" in consensus_transcripts[c_acc][:barcode_length] or "CCATCTCATATGTAGTACTCT" in consensus_transcripts[c_acc][ : -barcode_length]:
        #     print("OMG")


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


        potential_start_barcode_found = False
        match_pattern= "[|]+"
        c_seq_start_cutpoint = 0
        c_seq_end_cutpoint = len(c_seq)
        if not start_barcode_found:
            start_snippet = c_seq[:barcode_length]
            best_matches_start = 0
            best_matches_endpoint = 0
            longest_consequtive_stretch = 0
            longest_consequtive_stretch_endpoint = 0
            for barcode in start_barcode_set:
                barcode_aligned, start_snippet_aligned, start_snippet_start_snippet, start_snippet_end_snippet, matches, mismatches, deletions, insertions, match_line, c_seq_begin, c_seq_end, reference_begin, reference_end = ssw_alignment( "c_acc", "barcode", start_snippet, barcode, 2, 2)
                for m in re.finditer(match_pattern, match_line):
                    # print(m.group(0), m.start(0), m.end(0))
                    if m.end(0) - m.start(0) >= longest_consequtive_stretch and m.end(0) - m.start(0) >= 9:
                        longest_consequtive_stretch = m.end(0) - m.start(0)
                        longest_consequtive_stretch_endpoint = c_seq_begin + m.end(0) 
                if matches > best_matches_start and matches > 15:
                    # print(barcode_aligned, start_snippet_aligned, matches)
                    # print(match_line, c_acc, barcode_aligned, start_snippet_aligned, "start", start_snippet_start_snippet, "end", start_snippet_end_snippet, "barcode end snippet", barcode[-len(start_snippet_end_snippet):],  matches)
                    best_matches_start = matches
                    best_matches_endpoint = c_seq_end
                    # if start_snippet_end_snippet == "AGC" and barcode[-len(start_snippet_end_snippet):] == "AGC":
                    #     print("HERE", barcode, c_seq[:barcode_length] )
                    #     sys.exit()
                else:
                    pass
                    # print(barcode_aligned, start_snippet_aligned )
            c_seq_start_cutpoint = max(longest_consequtive_stretch_endpoint, best_matches_endpoint)
            if c_seq_start_cutpoint > 0:
                potential_start_barcode_found = True
                imperfect_count += 1
                # print("start cutpoint", c_seq_start_cutpoint)


        potential_end_barcode_found = False
        if not end_barcode_found:
            end_snippet = c_seq[ -(barcode_length) :]
            best_matches_end = 0
            best_matches_startpoint = barcode_length
            longest_consequtive_stretch = 0
            longest_consequtive_stretch_startpoint = barcode_length

            for barcode in end_barcode_set:
                barcode_aligned, end_snippet_aligned, end_snippet_start_snippet, end_snippet_end_snippet, matches, mismatches, deletions, insertions, match_line, c_seq_begin, c_seq_end, reference_begin, reference_end = ssw_alignment( "c_acc", "barcode", end_snippet, barcode, 2, 2)
                for m in re.finditer(match_pattern, match_line):
                    # print(m.group(0), m.start(0), m.end(0))
                    if m.end(0) - m.start(0) >= longest_consequtive_stretch and m.end(0) - m.start(0) >= 9:
                        longest_consequtive_stretch = m.end(0) - m.start(0)
                        longest_consequtive_stretch_startpoint = c_seq_begin + m.start(0)
                if matches > best_matches_end and matches > 15:
                    # print("THIS", barcode_aligned, end_snippet_aligned, matches)
                    # print(match_line, c_acc, "THIS", barcode_aligned, end_snippet_aligned, "start", end_snippet_start_snippet,  "end", end_snippet_end_snippet,matches)
                    best_matches_end = matches
                    best_matches_startpoint = c_seq_begin
                else:
                    pass


            max_end_offset = barcode_length - min(longest_consequtive_stretch_startpoint, best_matches_startpoint) 
            c_seq_end_cutpoint = len(c_seq) - max_end_offset
            if c_seq_end_cutpoint < len(c_seq):
                potential_end_barcode_found = True
                imperfect_count += 1
                # print("end cutpoint", c_seq_end_cutpoint)

            # if  best_matches_end > 15:
            #     potential_end_barcode_found = True
            #     imperfect_count += 1
            #     print(best_matches_end)

        c_seq = c_seq[ c_seq_start_cutpoint : c_seq_end_cutpoint ]

        if not start_barcode_found and not potential_start_barcode_found:
            not_found_count +=1
            # print("Start not found at all!!", c_acc)
        if not end_barcode_found and not potential_end_barcode_found:
            # print("End not found at all!!", c_acc)
            not_found_count +=1

        consensus_transcripts_filtered[c_acc] = c_seq
        # print(len(c_seq))

    print("Perfect:{0}, imperfect:{1}, not found:{2}".format(perfect_count, imperfect_count, not_found_count))
    print("total count: {0}".format(perfect_count + imperfect_count + not_found_count))
    assert len(consensus_transcripts_filtered) == len(consensus_transcripts)
    consensus_file = open(os.path.join(args.outfolder, "querys_filtered.fa"), "w")
    for acc, seq in misc_functions.iteritems(consensus_transcripts_filtered):
        consensus_file.write(">{0}\n{1}\n".format(acc, seq))
    consensus_file.close() 
    
    remove_redundant_after_barcode_removal(consensus_transcripts_filtered)


def remove_redundant_after_barcode_removal(consensus_transcripts_filtered):
    consensus_file_unique = open(os.path.join(args.outfolder, "querys_unique_after_barcode_removal.fa"), "w")
    consensus_transcripts_seq_to_acc = {seq : acc for (acc, seq) in  consensus_transcripts_filtered.items()}
    consensus_transcripts_unique = {acc : seq for (seq, acc) in  consensus_transcripts_seq_to_acc.items()}
    min_seq_len = min([len(seq) for seq in consensus_transcripts_seq_to_acc.keys()])
    print("MIN seq after barcode:", min_seq_len)
    for acc, seq in misc_functions.iteritems(consensus_transcripts_unique):
        consensus_file_unique.write(">{0}\n{1}\n".format(acc, seq))
    consensus_file_unique.close() 