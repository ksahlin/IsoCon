from collections import defaultdict
class CCS(object):
    """docstring for CCS"""
    def __init__(self, name, seq, qual, np):
        super(CCS, self).__init__()
        self.name = name
        self.seq = seq
        self.qual = qual
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

    def alignment_matrix_pos_to_ccs_coord(self, ccs_alignment_vector, pos):
        aln_piece = ccs_alignment_vector[:pos+1]
        coord_in_ccs = sum([1 for n in aln_piece if n != "-"]) - 1  # go from 1-indexed to 0-indexed

        if ccs_alignment_vector[pos] == "-": # if delation, jump one position to the right
            coord_in_ccs += 1
        return  coord #0-indexed

    def get_p_error_in_base(self, coord):
        q = self.qual[coord]
        p = 10**(-q/10.0)
        return p


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

import re

def get_ccs(ccs_file):  
    ccs_dict = defaultdict(dict)
    read_id_pattern = r"[\d]+/ccs"
    for read in ccs_file.fetch(until_eof=True):
        m = re.search(read_id_pattern, read.query_name)
        read_id = m.group(0).split("/")[0]
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

