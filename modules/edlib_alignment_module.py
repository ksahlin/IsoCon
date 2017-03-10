import signal
from multiprocessing import Pool
import multiprocessing as mp
import sys

import edlib

def edlib_align_sequences(matches, single_core = False):
    exact_edit_distances = {}
    if single_core:
        for j, s1 in enumerate(matches):
            for i, s2 in enumerate(matches[s1]):
                if s1 in exact_edit_distances:
                    if s2 in exact_edit_distances[s1]:
                        continue
                # print(s1,s2)
                s1, s2, ed = edlib_alignment_helper( ((s1, s2, i, j), {}) )
                if s1 in exact_edit_distances:
                    exact_edit_distances[s1][s2] = ed
                else:
                    exact_edit_distances[s1] = {}
                    exact_edit_distances[s1][s2]  = ed
    else:
        ####### parallelize alignment #########
        # pool = Pool(processes=mp.cpu_count())
        original_sigint_handler = signal.signal(signal.SIGINT, signal.SIG_IGN)
        signal.signal(signal.SIGINT, original_sigint_handler)
        pool = Pool(processes=mp.cpu_count())
        try:
            res = pool.map_async(edlib_alignment_helper, [ ((s1, s2, i,j), {}) for j, s1 in enumerate(matches) for i, s2 in enumerate(matches[s1]) ] )
            alignment_results =res.get(999999999) # Without the timeout this blocking call ignores all signals.
        except KeyboardInterrupt:
            print("Caught KeyboardInterrupt, terminating workers")
            pool.terminate()
            sys.exit()
        else:
            # print("Normal termination")
            pool.close()
        pool.join()
        for s1, s2, ed in alignment_results:
            if s1 in exact_edit_distances:
                exact_edit_distances[s1][s2] = ed
            else:
                exact_edit_distances[s1] = {}
                exact_edit_distances[s1][s2] = ed

    return exact_edit_distances

def edlib_align_sequences_keeping_accession(matches, single_core = False):
    """
        @param: Matches is be a 2D matrix implemented as a dict of dict with the accession of the two sequences as keys. The value is  a tuple (s1,s2) with the sequences.
        return: the same structure but a 3-tuple of seq1,seq2, ed
    """
    exact_matches = {}

    if single_core:
        for j, s1_acc in enumerate(matches):
            for i, s2_acc in enumerate(matches[s1_acc]):
                s1, s2 = matches[s1_acc][s2_acc]
                if s1_acc in exact_matches:
                    if s2_acc in exact_matches[s1_acc]:
                        continue

                s1_acc, s2_acc, (x,y,ed) = edlib_alignment_helper( ((s1, s2, i, j), {"x_acc" : s1_acc, "y_acc" :s2_acc}) )

                if s1_acc in exact_matches:
                    exact_matches[s1_acc][s2_acc] = (s1,s2, ed)
                else:
                    exact_matches[s1_acc] = {}
                    exact_matches[s1_acc][s2_acc]  = (s1,s2, ed)

    else:
        ####### parallelize alignment #########
        # pool = Pool(processes=mp.cpu_count())
        original_sigint_handler = signal.signal(signal.SIGINT, signal.SIG_IGN)
        signal.signal(signal.SIGINT, original_sigint_handler)
        pool = Pool(processes=mp.cpu_count())

        try:
            res = pool.map_async(edlib_alignment_helper, [ ((matches[s1_acc][s2_acc][0], matches[s1_acc][s2_acc][1], i,j), {"x_acc": s1_acc, "y_acc" : s2_acc}) for j, s1_acc in enumerate(matches) for i, s2_acc in enumerate(matches[s1_acc]) ] )
            alignment_results =res.get(999999999) # Without the timeout this blocking call ignores all signals.
        except KeyboardInterrupt:
            print("Caught KeyboardInterrupt, terminating workers")
            pool.terminate()
            sys.exit()
        else:
            # print("Normal termination")
            pool.close()
        pool.join()
        for s1_acc, s2_acc, (x,y,ed) in alignment_results:
            if s1_acc in exact_matches:
                exact_matches[s1_acc][s2_acc] = (x,y,ed)
            else:
                exact_matches[s1_acc] = {}
                exact_matches[s1_acc][s2_acc] = (x,y,ed)

    return exact_matches



def edlib_alignment_helper(arguments):
    args, kwargs = arguments
    return edlib_alignment(*args, **kwargs)

def edlib_alignment(x, y, i,j, x_acc = "", y_acc = ""):
    if i % 10000 == 0 and j % 100 == 0:
        print("processing alignments on y_j with j={0}".format(j+1))

    result = edlib.align(x,y, "NW")
    ed = result["editDistance"]
    # if result["editDistance"] != mismatches + indels:
    #     print(result["editDistance"], mismatches + indels, len(x), len(y), result["cigar"] )
    #     print(x_alignment)
    #     print(y_alignment)
    #     print(x)
    #     print(y)
    if x_acc == y_acc == "":
        return (x,y, ed)
    else:
        return (x_acc, y_acc, (x,y,ed))  

def edlib_traceback(x, y, mode="NW", task="path", k=1):
    result = edlib.align(x, y, mode=mode, task=task, k=k)
    ed = result["editDistance"]
    locations =  result["locations"]
    cigar =  result["cigar"]
    return ed, locations, cigar
