from collections import defaultdict
import re
import sys
import math

from modules.functions import reverse_complement


class CCS(object):
    """docstring for CCS"""
    def __init__(self, name, seq, qual, np):
        super(CCS, self).__init__()
        self.name = name
        self.seq = seq
        self.qual = qual
        bad_qual_values = [val for val in qual if  val < 0 or val > 93 ]
        if len(bad_qual_values) > 0:
            print(name, qual)
            print(bad_qual_values)
            print("At least one bad quality value in read. Phred quality values are assumed to be larger than 0 and smaller than 94 by protocol.")
            sys.exit()
        self.np = np
        self.subreads = {}
    def positions_with_p_error_higher_than(self, prob):
        # indicator_probs = ""
        indicator_probs = []
        for i, p_error in  enumerate(get_p_error_from_q(self.qual)):
            if p_error > prob:
                # indicator_probs += str(i) + "," + str(p_error)
                indicator_probs.append((i,p_error) )
            # else:
            #     indicator_probs += " "
        return(indicator_probs)


    def read_aln_to_ccs_coord(self, read_aln, pos):
        """
            Finds the correct position in the ccs read to get the base quality prediction from, given a position
            in our alignment from the fasta seq of the read to a candidate.
        """

        fasta_seq = "".join([n for n in read_aln if n != "-"])

        index = self.seq.index(fasta_seq)
        # print("found at index:", index, "length piece:", len(seq_piece), ccs_alignment_vector[pos] == "-")
            
        if (index + pos ) < len(self.seq):
            coord_in_ccs = index + pos  

        elif (index + pos ) == len(self.seq): #deletion occurs after last base pair (corner case and base quality is NA)
            coord_in_ccs = index + pos - 1
            print("Occurred at last position.", index, "length seq_piece:", len(fasta_seq), "length sequence:", len(self.seq))
        else:
            print("Index error:", index, "length seq_piece:", len(fasta_seq), "length sequence:", len(self.seq))
            sys.exit()

        return coord_in_ccs

    def alignment_matrix_pos_to_ccs_coord(self, ccs_alignment_vector, pos):
        """
            Finds the correct position in the ccs read to get the base quality prediction from, given a position
            in our alignment matrix.
        """


        aln_piece = ccs_alignment_vector[:pos+1]        

        # coord_in_ccs = sum([1 for n in aln_piece if n != "-"]) - 1  # go from 1-indexed to 0-indexed
        # if ccs_alignment_vector[pos] == "-": # if deletion, jump one position to the right
        #     coord_in_ccs += 1

        seq_piece = "".join([n for n in aln_piece if n != "-"])
        if len(seq_piece) > 20:
            # print(len(ccs_alignment_vector), pos+1)
            # print("".join([n for n in ccs_alignment_vector if n != "-"])) #print(seq_piece)
            # print(self.seq ) #[ : len(seq_piece)])
            index = self.seq.index(seq_piece)
            # print("found at index:", index, "length piece:", len(seq_piece), ccs_alignment_vector[pos] == "-")
            if ccs_alignment_vector[pos] == "-": # the position is within a deletion in the ccs sequence, uncertainty need to be obtained from the base following the deletion
                
                if (index + len(seq_piece) ) < len(self.seq):
                    coord_in_ccs = index + len(seq_piece)  

                elif (index + len(seq_piece) ) == len(self.seq): #deletion occurs after last base pair (corner case and base quality is NA)
                    coord_in_ccs = index + len(seq_piece) - 1
                else:
                    print("Index error:", index, "length seq_piece:", len(seq_piece), "length sequence:", len(len(self.seq)), ccs_alignment_vector[pos] == "-")
                    sys.exit()

            else:
                coord_in_ccs = index + len(seq_piece) - 1

        else:
            aln_piece_after = ccs_alignment_vector[pos:]        
            seq_piece = "".join([n for n in aln_piece_after if n != "-"])
            index = self.seq.index(seq_piece)
            coord_in_ccs = index
            print("Gah here")

        # below code does not work if substitution within a homopolymenr region..
        # # verify that index is the left most if homopolymer region, otherwise shift left
        # nucl = self.seq[coord_in_ccs]
        # if coord_in_ccs > 0:
        #     if nucl != self.seq[coord_in_ccs - 1]:
        #         return  coord_in_ccs #0-indexed
        #     else:
        #         while coord_in_ccs > 0 and (nucl == self.seq[coord_in_ccs - 1]):
        #             # print("HERRRRE")
        #             coord_in_ccs -= 1
        #             nucl = self.seq[coord_in_ccs]

        #         return coord_in_ccs
        # else:
        #     return  coord_in_ccs #0-indexed


        return coord_in_ccs
        

    def get_p_error_in_base(self, coord):
        q = self.qual[coord]
        p = 10**(-q/10.0)
        return p


def p_error_to_qual(p):
    # print(p)
    q = -10*math.log( p, 10)
    return q

def fix_quality_values(seq, qualities):
    assert len(seq) == len(qualities)
    left_shifted_qualities = []
    homo_pl_region = [ qualities[0] ]
    for i in range(1, len(seq)):
        if seq[i-1] == seq[i]:
            homo_pl_region.append( qualities[i] )
        else:
            left_shift = sorted(homo_pl_region) # homo_pl_region[::-1]
            # if left_shift[0] != min(left_shift):
            #     print(left_shift, seq[max(0,i-6):i+5], qualities[max(0,i-6):i+5], i, len(seq))
                # assert left_shift[0] == min(left_shift)
            left_shifted_qualities.append(left_shift) 
            homo_pl_region =  [ qualities[i] ]

    # last homopolymer region or base
    left_shift = sorted(homo_pl_region) # homo_pl_region[::-1]
    left_shifted_qualities.append(left_shift) 
    qual_values = [nucl_qual for poly_region in left_shifted_qualities for nucl_qual in poly_region]
    return qual_values


def modify_strings_and_acc_fastq(ccs_dict_raw, X_ids, X):
    # print(len(ccs_dict_raw))
    print(len(X_ids), len(X))
    assert len(X_ids) == len(X)
    # print(X.keys())
    # print(sorted(X_ids.keys()))
    # print(sorted(ccs_dict_raw.keys()))

    for q_id in list(ccs_dict_raw.keys()):
        if q_id in X_ids:
            # print(q_id, "in reads!")

            q_acc = X_ids[q_id]
            p = r"strand=-"
            m = re.search(p, q_acc)
            # print(q_acc)
            if m:
                ccs_record = ccs_dict_raw[q_id]
                # seq_rc = reverse_complement(ccs_record.seq)
                # qual_r = ccs_record.qual[::-1]
                qualities = fix_quality_values(ccs_record.seq, ccs_record.qual )
                start_index = ccs_record.seq.index(X[q_acc])
                stop_index = start_index + len(X[q_acc])
                ccs_record.seq = ccs_record.seq[start_index: stop_index]
                ccs_record.qual = qualities[start_index: stop_index]
                # index = ccs_record.seq.find("TCAGCCTCT")
                # if index >= 0:
                    # print("reversed:",  ccs_record.qual[index + 8:index + 14], ccs_record.seq[index + 8:index + 14])
                # print("HETEEE")
                # print(ccs_record.seq)
                # print(ccs_record.qual)
                assert ccs_record.seq == X[q_acc]
                assert len(ccs_record.seq) == len(ccs_record.qual)

            else:
                ccs_record = ccs_dict_raw[q_id]
                start_index = ccs_record.seq.index(X[q_acc])
                stop_index = start_index + len(X[q_acc])
                ccs_record.seq = ccs_record.seq[start_index: stop_index]
                ccs_record.qual = list(ccs_record.qual)[start_index: stop_index]
                # print(ccs_record.seq)
                # print(ccs_record.qual)
                # index = ccs_record.seq.find("TCAGCCTCT")
                # if index >= 0:
                #     print(index, ccs_record.qual[index + 8:index + 14], ccs_record.seq[index + 8:index + 14])
                assert ccs_record.seq == X[q_acc]
                assert len(ccs_record.seq) == len(ccs_record.qual)


            # change bach the name of the ccs_record to match with the flnc reads 
            ccs_obj = ccs_dict_raw[q_id]
            del ccs_dict_raw[q_id]
            new_q_name = X_ids[q_id]
            ccs_obj.name = new_q_name
            ccs_dict_raw[new_q_name] = ccs_obj

        else:
            del ccs_dict_raw[q_id]

    
    print(len(ccs_dict_raw))
    assert len(ccs_dict_raw) == len(X_ids)

    return ccs_dict_raw


def modify_strings_and_acc(ccs_dict_raw, X_ids, X):
    # print(len(ccs_dict_raw))
    print(len(X_ids), len(X))
    assert len(X_ids) == len(X)
    # print(X.keys())
    # print(sorted(X_ids.keys()))
    # print(sorted(ccs_dict_raw.keys()))

    for q_id in list(ccs_dict_raw.keys()):
        if q_id in X_ids:
            # print(q_id, "in reads!")

            q_acc = X_ids[q_id]
            p = r"strand=-"
            m = re.search(p, q_acc)
            if m:
                ccs_record = ccs_dict_raw[q_id]
                seq_rc = reverse_complement(ccs_record.seq)
                qual_r = ccs_record.qual[::-1]
                # print(seq_rc)
                # print(qual_r)
                qualities = fix_quality_values(seq_rc, qual_r )
                start_index = seq_rc.index(X[q_acc])
                stop_index = start_index + len(X[q_acc])
                ccs_record.seq = seq_rc[start_index: stop_index]
                ccs_record.qual = qualities[start_index: stop_index]
                # index = ccs_record.seq.find("TCAGCCTCT")
                # if index >= 0:
                #     print("reversed:",  ccs_record.qual[index + 8:index + 14], ccs_record.seq[index + 8:index + 14])
                assert ccs_record.seq == X[q_acc]
                assert len(ccs_record.seq) == len(ccs_record.qual)

            else:
                ccs_record = ccs_dict_raw[q_id]
                start_index = ccs_record.seq.index(X[q_acc])
                stop_index = start_index + len(X[q_acc])
                ccs_record.seq = ccs_record.seq[start_index: stop_index]
                ccs_record.qual = list(ccs_record.qual)[start_index: stop_index]
                # index = ccs_record.seq.find("TCAGCCTCT")
                # if index >= 0:
                #     print(index, ccs_record.qual[index + 8:index + 14], ccs_record.seq[index + 8:index + 14])
                assert ccs_record.seq == X[q_acc]
                assert len(ccs_record.seq) == len(ccs_record.qual)


            # change bach the name of the ccs_record to match with the flnc reads 
            ccs_obj = ccs_dict_raw[q_id]
            del ccs_dict_raw[q_id]
            new_q_name = X_ids[q_id]
            ccs_obj.name = new_q_name
            ccs_dict_raw[new_q_name] = ccs_obj

        else:
            del ccs_dict_raw[q_id]

    
    print(len(ccs_dict_raw))
    assert len(ccs_dict_raw) == len(X_ids)
    # print("HERE!")

    # lambda_ = 0
    # lambda_2 = 0
    # cntr = 0
    # for q_acc in ccs_dict_raw:
    #     ccs_record = ccs_dict_raw[q_acc]
    #     # print(ccs_record.qual[575:610], ccs_record.seq[575:610])
    #     index = ccs_record.seq.find("TCAGCCTCT")
    #     if index >= 0:
    #         print(ccs_record.qual[index + 9: index + 15], ccs_record.seq[index + 9: index + 15])
    #         p_error = ccs_record.get_p_error_in_base(index + 9)
    #         p_error2 = sum([ccs_record.get_p_error_in_base(index + pos) for pos in range(9,13)])
    #         print(index + 9, p_error)
    #         lambda_ += p_error
    #         lambda_2 += min(1.0, p_error2)
    #         cntr += 1
    # print("tot prob:", lambda_, "tot obs:", cntr)
    # print("tot sum prob homopol:", lambda_2, "tot obs:", cntr)
    
    return ccs_dict_raw


def get_p_error_from_q(qual_list):
    probs = []
    for q in qual_list:
        p = 10**(-q/10.0)
        probs.append(p)
    return probs


def get_p_error_from_char(qual_string):
    probs = []
    for char in qual_string:
        q = ord(char) - 33
        p = 10**(-q/10.0)
        probs.append(p)
    return probs


def get_ccs(ccs_file):  
    ccs_dict = defaultdict(dict)
    # read_id_pattern = r"[\d]+/ccs"
    read_id_pattern = r".+/ccs"
    for read in ccs_file.fetch(until_eof=True):
        m = re.search(read_id_pattern, read.query_name)
        # read_id = m.group(0).split("/")[0]
        read_id = m.group(0)[:-4]
        # print(read_id)
        ccs_read = CCS(read_id, read.query_alignment_sequence, read.query_qualities, read.get_tag("np"))
        ccs_dict[read_id] = ccs_read
        # ccs_dict[read.query_name]["seq"] = read.query_alignment_sequence
        # print(read.query_qualities)
        # sys.exit()
        # ccs_dict[read.query_name]["qual"] = read.query_qualities
        # ccs_dict[read.query_name]["np"] = read.get_tag("np")
        assert len(read.query_alignment_sequence) == len(read.query_qualities)
        
        # if ccs_read.np > 10:
        #     print(ccs_read.np, ccs_read.positions_with_p_error_higher_than(0.01))
    return ccs_dict

