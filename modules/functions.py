import unittest

def position_query_to_alignment(query_aligned, target_aligned, target_alignment_start_position):
    """
        input:      0-indexed target positions
        returns:    target vector is a list of 2*t_len +1 positions to keep insertions between base pairs    
                list of strings of nucleotides for each position, start position in target vector list, end position in target vector list
    """

    query_positioned = [] 
    target_position = target_alignment_start_position # 0-indexed
    temp_ins = ""
    # iterating over alignment positions
    for p in range(len(target_aligned)):
        if target_aligned[p] == "-":
            temp_ins += query_aligned[p]
        else:
            if not temp_ins:
                query_positioned.append("-") #[2*target_position] = "-"
            else:
                query_positioned.append(temp_ins)
                temp_ins = ""

            query_positioned.append(query_aligned[p])

            target_position += 1

    if not temp_ins:
        query_positioned.append("-")
    else:
        query_positioned.append(temp_ins)

    target_vector_start_position = 2*target_alignment_start_position
    target_vector_end_position = 2*(target_position-1) + 2

    return query_positioned, target_vector_start_position, target_vector_end_position


def arrange_query_alignments_to_target(target_accession, target_sequence, query_to_target_alignments):
    """
        input:
            a dictionary query_to_target_alignments of the form
            {query_accession1 : [(target_aligned, query_aligned, target_start, target_end)],
             query_accession2 : [(target_aligned, query_aligned, target_start, target_end)],
             ....
            }
        output: 
            Target vector is a list of 2*len(target_seq) +1 positions
    """

    query_to_target_positioned_dict = {} # query_accession : [] list of length 2l_j +1 to hold insertions
    for query_accession in query_to_target_alignments:
        target_aligned, query_aligned, target_start, target_end = query_to_target_alignments[query_accession]

        query_to_target_positioned, target_vector_start_position, target_vector_end_position = position_query_to_alignment(query_aligned, target_aligned, target_start)
        query_to_target_positioned_dict[query_accession] = (query_to_target_positioned, target_vector_start_position, target_vector_end_position)
        
    return target_accession, query_to_target_positioned_dict


def get_non_overlapping_intervals(ranges):
    """
        example input:     # ranges = [(0,100,'a'),(0,75,'b'),(95,150,'c'),(120,130,'d')]
        output: 
    """
    non_overlapping_parts = []
    # ranges = [(0,100,'a'),(0,75,'b'),(95,150,'c'),(120,130,'d')]
    endpoints = sorted(list(set([r[0] for r in ranges] + [r[1] for r in ranges])))
    start = {}
    end = {}
    for e in endpoints:
        start[e] = set()
        end[e] = set()
    for r in ranges:
        start[r[0]].add(r[2])
        end[r[1]].add(r[2])

    current_ranges = set()
    prev_set_size = 0
    for e1, e2 in zip(endpoints[:-1], endpoints[1:]):
        current_ranges.difference_update(end[e1])
        current_ranges.update(start[e1])

        if prev_set_size > len(current_ranges):
            start_offset = 1
        else:
            start_offset = 0

        next_current_ranges = set(current_ranges)
        next_current_ranges.difference_update(end[e2])
        next_current_ranges.update(start[e2])
        next_set_size = len(next_current_ranges)
        if next_set_size < len(current_ranges):
            stop_offset = 0
        else:
            stop_offset = -1

        if current_ranges:
            # print('%d - %d: %s' % (e1 + start_offset, e2 + stop_offset, ','.join(current_ranges)))
            non_overlapping_parts.append((e1 + start_offset, e2 + stop_offset, tuple(current_ranges)))
        else:
            pass
            # print('%d - %d: %s' % (e1 + start_offset, e2 + stop_offset, ','.join(current_ranges)), "lol")
        
        prev_set_size = len(current_ranges)

    # sys.exit()
    return(non_overlapping_parts)




class TestFunctions(unittest.TestCase):

    def test_positioning(self):
        """
            the start and end of query positions onto the target is always one position 
            before and one position after the start and end of the target coordinates.
        """
        t_seq = "ACGGA"

        q, t, t_start = "ACGGA", "ACGGA", 0
        self.assertEqual(position_query_to_alignment(q, t, t_start), (["-", "A","-", "C", "-","G","-", "G", "-", "A", "-"], 0, 10) )

        q, t, t_start = "TACGGA", "-ACGGA", 0
        self.assertEqual(position_query_to_alignment(q, t, t_start), (["T", "A","-", "C", "-","G","-", "G", "-", "A", "-"], 0, 10) )

        q, t, t_start = "ACGGATTT", "ACGGA---", 0
        self.assertEqual(position_query_to_alignment(q, t, t_start), (["-", "A","-", "C", "-","G","-", "G", "-", "A", "TTT"], 0, 10) )

        q, t, t_start = "ACG", "ACG", 0
        self.assertEqual(position_query_to_alignment(q, t, t_start), (["-", "A","-", "C", "-","G","-"], 0, 6) )

        q, t, t_start = "GGA", "GGA", 2
        self.assertEqual(position_query_to_alignment(q, t, t_start), (["-","G","-", "G", "-", "A", "-"], 4, 10) )

    def test_intervals(self):
        ranges = [(0,100,'a'),(0,75,'b'),(95,150,'c'),(120,130,'d')]
        result = [(0,75,("a", "b")),(76, 94, ("a",)), (95,100, ("a", "c")), (101,119, ("c",)), (120,130,("c","d")), (131,150,("c",)) ]
        self.assertEqual(get_non_overlapping_intervals(ranges), result)



if __name__ == '__main__':
    unittest.main()
