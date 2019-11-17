
from SAScoringSystem import ScoringSystem
from FastA import read_sequence, write_sequences, write_comparison
import numpy as np

#constants
DIAGONAL = 0
DOWN = 1
RIGHT = 2


class GlobalSequenceLinearGapAligner:
    def __init__(self, scoring_system : ScoringSystem, sequences = None, headers = None, align = False):
        #dict mapping header to sequence
        self.sequences = [] if sequences is None else sequences
        self.headers = [] if headers is None else headers
        self.aligned_sequences = []
        self.scoring_system = scoring_system
        if align: self.align()

    def clear_sequences(self):
        self.aligned_sequences.clear()
        self.sequences.clear()
        self.headers.clear()

    def append_seq_from_file(self, filename):
        headers, sequences = read_sequence(filename)
        self.sequences.extend(sequences)
        self.headers.extend(headers)

    def align(self):
        if len(self.sequences) != 2:
            logger.error("Currently alignment is only implemented for exactly 2 sequences.")
            raise NotImplementedError
        self._initialize_alignment()
        self._compute_matrices()
        return self._traceback()

    def _initialize_alignment(self):
        self.cols = len(self.sequences[0]) + 1
        self.rows = len(self.sequences[1]) + 1

        #matrices are row major
        self.nw_matrix = [[None] * self.cols for i in range(self.rows)]
        self.trace_matrix = [[None] * self.cols for i in range(self.rows)]
        
        #initialize the nw_matrix with 0 and subsequent gap penalties
        self.nw_matrix[0] = [self.scoring_system.gap_penalty(1) * i for i in range(self.cols)]
        self.trace_matrix[0] = [RIGHT] * self.cols
        for i in range(self.rows):
            self.nw_matrix[i][0] = self.scoring_system.gap_penalty(1) * i
            self.trace_matrix[i][0] = DOWN

    def _compute_matrices(self):
        for row_index in range(1, self.rows):
            for column_index in range(1, self.cols):
                top_value = self.nw_matrix[row_index - 1][column_index]
                left_value = self.nw_matrix[row_index][column_index - 1]
                diag_value = self.nw_matrix[row_index - 1][column_index - 1]
                top_value += self.scoring_system.gap_penalty(1) #TODO: currently we assume linear gap penalty
                left_value += self.scoring_system.gap_penalty(1) #TODO: currently we assume linear gap penalty
                diag_value += self.scoring_system.compare(self.sequences[0][column_index - 1], self.sequences[1][row_index - 1])

                # find max_value
                # beware that this has to correspond to the definitions of DIAGONAL, DOWN and RIGHT
                maximum_value = max(diag_value, top_value, left_value)
                if maximum_value == diag_value:
                    self.trace_matrix[row_index][column_index] = DIAGONAL # even if there's multiple maxima we arbitrarily choose the first
                elif maximum_value == top_value:
                    self.trace_matrix[row_index][column_index] = DOWN # even if there's multiple maxima we arbitrarily choose the first
                elif maximum_value == left_value:
                    self.trace_matrix[row_index][column_index] = RIGHT # even if there's multiple maxima we arbitrarily choose the first

                self.nw_matrix[row_index][column_index] = maximum_value
                
    def _traceback(self):
        x_aligned = ""
        y_aligned = ""
        (row_index, column_index) = (self.rows - 1, self.cols - 1)
        while row_index and column_index:
            if self.trace_matrix[row_index][column_index] == DIAGONAL:
                x_aligned += self.sequences[0][column_index - 1]
                y_aligned += self.sequences[1][row_index - 1]
                row_index -= 1
                column_index -= 1
            elif self.trace_matrix[row_index][column_index] == DOWN:
                y_aligned += self.sequences[1][row_index - 1]
                x_aligned += "-"
                row_index -= 1
            elif self.trace_matrix[row_index][column_index] == RIGHT:
                x_aligned += self.sequences[0][column_index - 1]
                y_aligned += "-"
                column_index -= 1
        while column_index: # row index is 0, i.e. we arrived at the top
            y_aligned += "-"
            x_aligned += self.sequences[0][column_index - 1]
            column_index -= 1
        while row_index:
            x_aligned += "-"
            y_aligned += self.sequences[1][row_index - 1]
            row_index -= 1
        self.aligned_sequences.append(x_aligned[::-1])
        self.aligned_sequences.append(y_aligned[::-1])
        return x_aligned[::-1], y_aligned[::-1]

    def get_maximum_score(self):
        return self.nw_matrix[-1][-1]
