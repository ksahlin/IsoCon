
import os
import errno
import re
import argparse


'''
    Below code taken from https://github.com/lh3/readfq/blob/master/readfq.py
'''

def readfq(fp): # this is a generator function
    last = None # this is a buffer keeping the last unprocessed line
    while True: # mimic closure; is it a bad idea?
        if not last: # the first record or a record following a fastq
            for l in fp: # search for the start of the next record
                if l[0] in '>@': # fasta/q header line
                    last = l[:-1] # save this line
                    break
        if not last: break
        name, seqs, last = last[1:].replace(" ", "_"), [], None
        for l in fp: # read the sequence
            if l[0] in '@+>':
                last = l[:-1]
                break
            seqs.append(l[:-1])
        if not last or last[0] != '+': # this is a fasta record
            yield name, (''.join(seqs), None) # yield a fasta record
            if not last: break
        else: # this is a fastq record
            seq, leng, seqs = ''.join(seqs), 0, []
            for l in fp: # read the quality
                seqs.append(l[:-1])
                leng += len(l) - 1
                if leng >= len(seq): # have read enough quality
                    last = None
                    yield name, (seq, ''.join(seqs)); # yield a fastq record
                    break
            if last: # reach EOF before reading enough quality
                yield name, (seq, None) # yield a fasta record instead
                break


def mkdir_p(path):
    try:
        os.makedirs(path)
        print("creating", path)
    except OSError as exc:  # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise

def main(args):
    from collections import defaultdict
    clusters = defaultdict(list)

    with open(args.clusters) as f:
        for line in f:
            items = line.strip().split()
            acc, cl_id  = items[0], items[1]
            if acc == "from":
                continue
            if "/" in cl_id:
                cl_id = cl_id.split("/")[1]
            clusters[cl_id].append(acc)

    mkdir_p(args.outfolder)
    reads = { acc : (seq, qual) for acc, (seq, qual) in readfq(open(args.fastq, 'r'))}
    
    for cl_id in clusters:
        r = clusters[cl_id]

        if len(r) >= args.N:
            curr_file = open(os.path.join(args.outfolder, str(cl_id) + ".fastq" ), "w")
            for acc in r:
                seq, qual = reads[acc]
                curr_file.write("@{0}\n{1}\n{2}\n{3}\n".format(acc, seq, "+", qual))
            curr_file.close()



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="clusters to fastq", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--version', action='version', version='%(prog)s 0.0.1')

    parser.add_argument('fastq', type=str, help='Path to consensus fastq file(s)')
    parser.add_argument('clusters', type=str, help='The tsv cluster file: from to') 
    parser.add_argument('outfolder', type=str, help='The outfolder file') 
    parser.add_argument('--N', type=int, default = 0, help='Write out clusters with more or equal than N reads')

    args = parser.parse_args()

    main(args)      
