import argparse
import os 
import sys

import re
import errno
import shutil

import pysam
from collections import defaultdict

def mkdir_p(path):
    try:
        os.makedirs(path)
        print("creating", path)
    except OSError as exc:  # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise


def read_fasta(fasta_file):
    fasta_seqs = {}
    k = 0
    temp = ''
    accession = ''
    for line in fasta_file:
        if line[0] == '>' and k == 0:
            accession = line[1:].strip().replace(" ", "_") #.split()[0]
            fasta_seqs[accession] = ''
            k += 1
        elif line[0] == '>':
            yield accession, temp
            temp = ''
            accession = line[1:].strip().replace(" ", "_") #.split()[0]
        else:
            temp += line.strip()

    if accession != '':
        yield accession, temp


class CCS(object):
    """docstring for CCS"""
    def __init__(self, name, seq, qual, np):
        super(CCS, self).__init__()
        self.name = name
        self.seq = seq
        self.qual = qual
        self.qualstring = "".join([ chr(int(q)+33) for q in self.qual ])

        bad_qual_values = [val for val in qual if  val < 0 or val > 93 ]

        if len(bad_qual_values) > 0:
            print(name, qual)
            print(bad_qual_values)
            print("At least one bad quality value in read. Phred quality values are assumed to be larger than 0 and smaller than 94 by protocol.")
            sys.exit()
        self.np = np
        self.subreads = {}



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

        elif (index + len(seq_piece) ) == len(self.seq): #deletion occurs after last base pair (corner case and base quality is NA)
            coord_in_ccs = index + pos - 1
            print("Occurred at last position.", index, "length seq_piece:", pos, "length sequence:", len(self.seq))
        else:
            print("Index error:", index, "length seq_piece:", len(seq_piece), "length sequence:", len(self.seq))
            sys.exit()

        return coord_in_ccs



def reverse_complement(string):
    #rev_nuc = {'A':'T', 'C':'G', 'G':'C', 'T':'A', 'N':'N', 'X':'X'}
    # Modified for Abyss output
    rev_nuc = {'A':'T', 'C':'G', 'G':'C', 'T':'A', 'a':'t', 'c':'g', 'g':'c', 't':'a', 'N':'N', 'X':'X', 'n':'n', 'Y':'R', 'R':'Y', 'K':'M', 'M':'K', 'S':'S', 'W':'W', 'B':'V', 'V':'B', 'H':'D', 'D':'H', 'y':'r', 'r':'y', 'k':'m', 'm':'k', 's':'s', 'w':'w', 'b':'v', 'v':'b', 'h':'d', 'd':'h'}

    rev_comp = ''.join([rev_nuc[nucl] for nucl in reversed(string)])
    return(rev_comp)

def get_ccs(ccs_file):  
    ccs_dict = defaultdict(dict)
    # read_id_pattern = r"[\d]+/ccs"
    unique_acc = set()
    cntr = 0
    read_id_pattern = r".+/ccs"
    for read in ccs_file.fetch(until_eof=True):
        cntr += 1
        m = re.search(read_id_pattern, read.query_name)
        # read_id = m.group(0).split("/")[0]
        read_id = m.group(0)[:-4]
        # print(read_id)
        ccs_read = CCS(read_id, read.query_alignment_sequence, read.query_qualities, read.get_tag("np"))
        ccs_dict[read_id] = ccs_read
        unique_acc.add(read_id)
        # ccs_dict[read.query_name]["seq"] = read.query_alignment_sequence
        # print(read.query_qualities)
        # sys.exit()
        # ccs_dict[read.query_name]["qual"] = read.query_qualities
        # ccs_dict[read.query_name]["np"] = read.get_tag("np")
        assert len(read.query_alignment_sequence) == len(read.query_qualities)
    
        # if ccs_read.np > 10:
        #     print(ccs_read.np, ccs_read.positions_with_p_error_higher_than(0.01))
    print("UNIQUE ACC:", len(unique_acc), "nr reads:", cntr)
    return ccs_dict

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


def modify_strings_and_acc(ccs_dict_raw, fasta_ids, fasta_dict):
    # print(len(ccs_dict_raw))
    print(len(fasta_ids), len(fasta_dict), len(ccs_dict_raw))
    assert len(fasta_ids) == len(fasta_dict)
    # print(fasta_dict.keys())
    # print(sorted(fasta_ids.keys()))
    # print(sorted(ccs_dict_raw.keys()))

    for q_id in list(ccs_dict_raw.keys()):
        if q_id in fasta_ids:
            # print(q_id, "in reads!")

            q_acc = fasta_ids[q_id]
            p = r"strand=-"
            m = re.search(p, q_acc)
            if m:
                ccs_record = ccs_dict_raw[q_id]
                seq_rc = reverse_complement(ccs_record.seq)
                qual_r = ccs_record.qual[::-1]
                # print(seq_rc)
                # print(qual_r)
                qualities = fix_quality_values(seq_rc, qual_r )
                start_index = seq_rc.index(fasta_dict[q_acc])
                stop_index = start_index + len(fasta_dict[q_acc])
                ccs_record.seq = seq_rc[start_index: stop_index]
                ccs_record.qual = qualities[start_index: stop_index]
                ccs_record.qualstring = "".join( [ chr(int(q)+33) for q in  ccs_record.qual ] )
                # index = ccs_record.seq.find("TCAGCCTCT")
                # if index >= 0:
                #     print("reversed:",  ccs_record.qual[index + 8:index + 14], ccs_record.seq[index + 8:index + 14])
                assert ccs_record.seq == fasta_dict[q_acc]
                assert len(ccs_record.seq) == len(ccs_record.qual)

            else:
                ccs_record = ccs_dict_raw[q_id]
                start_index = ccs_record.seq.index(fasta_dict[q_acc])
                stop_index = start_index + len(fasta_dict[q_acc])
                ccs_record.seq = ccs_record.seq[start_index: stop_index]
                ccs_record.qual = list(ccs_record.qual)[start_index: stop_index]
                ccs_record.qualstring = "".join( [ chr(int(q)+33) for q in ccs_record.qual ] )

                # index = ccs_record.seq.find("TCAGCCTCT")
                # if index >= 0:
                #     print(index, ccs_record.qual[index + 8:index + 14], ccs_record.seq[index + 8:index + 14])
                assert ccs_record.seq == fasta_dict[q_acc]
                assert len(ccs_record.seq) == len(ccs_record.qual)


            # change bach the name of the ccs_record to match with the flnc reads 
            ccs_obj = ccs_dict_raw[q_id]
            del ccs_dict_raw[q_id]
            new_q_name = fasta_ids[q_id]
            ccs_obj.name = new_q_name
            ccs_dict_raw[new_q_name] = ccs_obj

        else:
            del ccs_dict_raw[q_id]

    
    print(len(ccs_dict_raw))
    assert len(ccs_dict_raw) == len(fasta_ids)    
    return ccs_dict_raw



def main(args):
    fasta_dict = {acc : seq for acc, seq in read_fasta(open(args.flnc_fasta, "r"))}
    ccs_file = pysam.AlignmentFile(args.bam, "rb", check_sq=False)
    ccs_dict_raw = get_ccs(ccs_file)
    fasta_ids = { "/".join(acc.split("/")[:2]) : acc for acc in fasta_dict} 
    print("number of fasta reads in fasta file:", len(fasta_ids))
    ccs_dict = modify_strings_and_acc(ccs_dict_raw, fasta_ids, fasta_dict)
    
    fastq_out = open(args.outfile, "w")
    for acc, seq in fasta_dict.items():
        assert fasta_dict[acc] == ccs_dict[acc].seq
        assert len(seq) == len(ccs_dict[acc].qualstring)
        fastq_out.write("@{0}_passes:{1}\n{2}\n+\n{3}\n".format(acc, ccs_dict[acc].np, seq,  ccs_dict[acc].qualstring ))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Generate a fastq file with flnc reads and their base quality values from a fasta file and a ccs bam file.")
    parser.add_argument('flnc_fasta', type=str, help='The fasta file')
    parser.add_argument('bam',  type=str, help='The bam file')
    parser.add_argument('outfile', type=str, help='Output path to fasta file')

    if len(sys.argv) == 1:
        sys.argv.append('-h')

    args = parser.parse_args()
    path_, file_prefix = os.path.split(args.outfile)
    mkdir_p(path_)
    main(args)


