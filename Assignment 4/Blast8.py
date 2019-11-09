import logging

logger = logging.getLogger()

__author__ = "Patrick Schirm, Manuel Glöckler, Finn Mier"

#query, subject, %identical, alignment length, mismatches, gap openings, query start, query end, subject start, subject end, E value, bit score,


class Blast8Reader:

    def __init__(self, subject_lengths = None):
        self.subject_lengths = subject_lengths
        self.reads = {}

    def get_subject_length(self, subject):
        if type(self.subject_lengths) == int:
            return self.subject_lengths
        elif self.subject_lengths is None:
            return None
        else:
            return self.subject_lengths[subject]


    def compute_hit_attributes(self, read, count_occurences = False):
        """ Computes additional useful attributes per found read: the start the match would have to have,
        the direction of the hit, and optionally the number of occurences of this combination of query and subject.

        Parameters:
            read: List of read attributes without query and subject.
            count_occurences (bool, optional): If True adds the amount of occurences of this query/subject combination at the end.

        Returns:
            list of given read attributes and computed attributes.
        """
        query, subject, identical, length, mismatches, gap_openings, query_start, query_end, subject_start, subject_end, E_value, bit_score = read

        direction = '+' if int(subject_start) < int(subject_end) else '-'

        if count_occurences: occurences = self.reads.get((query, subject), [0])[-1] + 1
        #TODO: not sure if the actual_start location is correctly implemented for the negative case
        if direction == "+":
            actual_start = int(subject_start) - int(query_start)
        elif self.get_subject_length(subject):
            actual_start = self.get_subject_length(subject) - int(subject_start) + (int(query_start) - 1) #seems very convoluted
        else:
            actual_start = None

        # actual_start is attribute 11 (index 10)
        if count_occurences: return [*read[2:], actual_start, direction, occurences]
        else: return [*read[2:], actual_start, direction]



    def read_blast8(self, file_path):
        """ Takes the path to a .blat8 file and writes its contents to a dictionary mapping
        a query and subject to all corresponding hits which is a list with the following attributes: 
        [%identical, alignment length, mismatches, gap openings, query start, query end, subject start, 
        subject end, E value, bit_score, hit_start, number hits].
        
        Parameters:
            file_path (str): Path to the .blast8 file.

        Returns:
            {(query, subject) : [[%identical, alignment length, mismatches, gap openings, query start, query end, subject start, subject end, E value, bit, score, hit_start, number hits]]}
        """
        if not (file_path.endswith(".blast8")):
            logger.error("Filenames must end in '.blast8' but was {}".format(file_path))
            raise ValueError
        else:
            with open(file_path) as blast8_file:
                for line_nr, line in enumerate(blast8_file):
                    line = line.rstrip("\n")
                    blast_read = line.split("\t")
                    query, subject = blast_read[0], blast_read[1]
                    if (query, subject) in self.reads:
                        self.reads[(query, subject)].append(self.compute_hit_attributes(blast_read))
                    else:
                        self.reads[(query, subject)] = [self.compute_hit_attributes(blast_read)]
            return self.reads


    def read_probable_blast8(self, file_path, should_replace):
        """ Takes the path to a .blat8 file and writes its contents to a dictionary mapping
        a query and subject to the most probable (according to the given should_replace function)
        corresponding hit which is a list with the following attributes: 
        [%identical, alignment length, mismatches, gap openings, query start, query end, subject start, 
        subject end, E value, bit_score, hit_start, number hits].
        
        Parameters:
            file_path (str): Path to the .blast8 file.
            should_replace ((cur_read, candidate_read) => bool): Function which should return true if current read should be replaced by candidate read.

        Returns:
            {(query, subject) : [%identical, alignment length, mismatches, gap openings, query start, query end, subject start, subject end, E value, bit, score, hit_start, number hits]}
        """
        if not (file_path.endswith(".blast8")):
            logger.error("Filenames must end in '.blast8' but was {}".format(file_path))
            raise ValueError
        else:
            with open(file_path) as blast8_file:
                for line_nr, line in enumerate(blast8_file):
                    line = line.rstrip("\n")
                    blast_read = line.split("\t")
                    query, subject = blast_read[0], blast_read[1]
                    if (query, subject) in self.reads:
                        if should_replace(self.reads[(query, subject)], blast_read):
                            self.reads[(query, subject)] = self.compute_hit_attributes(blast_read, True)
                        else:
                            self.reads[(query, subject)][-1] += 1 #increase occurence counter by 1
                    else:
                        self.reads[(query, subject)] = self.compute_hit_attributes(blast_read, True)
            return self.reads

    def longest_sequence(cur_read, candidate_read):
        return int(cur_read[1]) < int(candidate_read[3])
                

    def e_value(cur_read, candidate_read):
        return float(cur_read[8]) > float(candidate_read[10])

    def identical(cur_read, candidate_read):
        return float(cur_read[0]) < float(candidate_read[2])

    def bit_score(cur_read, candidate_read):
        return float(cur_read[9]) < float(candidate_read[11])

