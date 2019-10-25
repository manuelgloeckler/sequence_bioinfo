#!/usr/bin/python
import argparse
import numpy as np

"""Basic Needleman-Wunsch"""

__author__ = "Patrick Schirm, Manuel Gloeckler"

"""
The script is created by the above-mentioned two students for the same course 1 year ago
This script implements basic Needleman-Wunsch algorithm.
Gap cost is linear. Command line argument parsing is not included.

	To-dos for those who wants to use this script:
	modify the code so that:
	- you can easily change the scoring scheme
	- the initialization, recursion and traceback steps are clear.
	These will make it way easier to switch to other algorithms
	- it allows commandline options, as follows (set defaults as you see fit):
	#  <input-file> read sequences from named file
	#  <output-file> write aligned sequences to named fasta file
	#  <score>: Set match score.
	#  <score>: Set mismatch score.
	#  <score>: Set gap penalty.
	- it might be able to catch some possible error




"""



class Score:
    """scores are stored here
    Hints: you may not need it
    """

    def __init__(self, match=1, mismatch=-2, gap=3, best_score=None):
        self.match = match
        self.mismatch = mismatch
        self.gap = gap
        self.best_score = best_score


def _tracking_guide(diag, top, left):
    winner = max(diag, top, left)
    # diag as 0, top 1, left 2
    if winner == diag:
        return winner, 0
    elif winner == top:
        return winner, 1
    elif winner == left:
        return winner, 2


def parse_score(path):
    file = open(path, 'r')
    score_dict = {}
    alphabet = []
    i = 0
    for x in file:
        if i == 0:
            alphabet = x.split("\t")
            alphabet = [x.strip() for x in alphabet]
        else:
            row = x.split("\t")
            subst = row[0]
            del row[0]
            for j in range(len(row)):
                score_dict[(alphabet[j], subst)] = row[j].strip()
        i+=1
    return score_dict

def gotho_matrix_builder(x_seq, y_seq, score, d, e):
    ncol = len(x_seq) + 1
    nrow = len(y_seq) + 1
    M = [i[:] for i in [[0] * ncol] * nrow]
    I_x = [i[:] for i in [[0] * ncol] * nrow]
    I_y = [i[:] for i in [[0] * ncol] * nrow]

    # Fill table with values
    for i in range(nrow):
        for j in range(ncol):
            if i == 0 and j == 0:
                M[i][j] = 0
                I_x[i][j] = -np.inf
                I_y[i][j] = -np.inf
            elif i == 0:
                M[i][j] = -np.inf
                I_x[i][j] = -d + (j)
                I_y[i][j] = -np.inf
            elif j == 0:
                matF[i][j] = i * (- score.gap)
                matT[i][j] = 2
            else:
                s = score.match if x_seq[j - 1] == y_seq[i - 1] else score.mismatch
                diag = matF[i - 1][j - 1] + s
                top = matF[i - 1][j] - score.gap
                left = matF[i][j - 1] - score.gap
                matF[i][j], matT[i][j] = _tracking_guide(diag, top, left)

    return True




# compute alignment score
# Initialize F with |y|+1 rows and |x|+1 columns
def matrix_builder(x_seq, y_seq, score):
    ncol = len(x_seq) + 1
    nrow = len(y_seq) + 1
    matF = [i[:] for i in [[0] * ncol] * nrow]
    matT = [i[:] for i in [[0] * ncol] * nrow]

    # Fill table with values
    for i in range(nrow):
        for j in range(ncol):
            if i == 0:
                matF[i][j] = j * (- score.gap)
                matT[i][j] = 1
            elif j == 0:
                matF[i][j] = i * (- score.gap)
                matT[i][j] = 2
            else:
                s = score.match if x_seq[j - 1] == y_seq[i - 1] else score.mismatch
                diag = matF[i - 1][j - 1] + s
                top = matF[i - 1][j] - score.gap
                left = matF[i][j - 1] - score.gap
                matF[i][j], matT[i][j] = _tracking_guide(diag, top, left)
    score.best_score = matF[nrow - 1][ncol - 1]
    print("Best alignment score: {}".format(score.best_score))
    return matT


# perform trace-back
# While doing the traceback, the aligned sequences (i.e. including gaps) are constructed in reverse order
def traceback_to_alignment(matT, x_seq, y_seq):
    x_aligned = ""
    y_aligned = ""
    i, j = len(y_seq), len(x_seq)
    while i and j:
        if matT[i][j] == 0:
            x_aligned += x_seq[j - 1]
            y_aligned += y_seq[i - 1]
            i -= 1
            j -= 1
        elif matT[i][j] == 1:
            x_aligned += "-"
            y_aligned += y_seq[i - 1]
            i -= 1
        elif matT[i][j] == 2:
            x_aligned += x_seq[j - 1]
            y_aligned += "-"
            j -= 1
    while i:
        x_aligned += "-"
        y_aligned += y_seq[i - 1]
        i -= 1
    while j:
        y_aligned += "-"
        x_aligned += x_seq[j - 1]
        j -= 1
    return x_aligned[::-1], y_aligned[::-1]


def run_nw(x_seq, y_seq, score):
    matT = matrix_builder(x_seq, y_seq, score)
    return traceback_to_alignment(matT, x_seq, y_seq)


def seq_file_reader(input_file):
    with open(input_file, "r") as f:
        header = f.readline().strip()
        seq = ""
        line = f.readline()
        while line:
            # if another header
            if line[0] == ">":
                break
            seq += line.strip()
            line = f.readline()
    f.close()
    return header, seq


def init_parser():
    # change 1
    parser = argparse.ArgumentParser(description='Needleman-Wunsch',prog="pqa", formatter_class=argparse.RawDescriptionHelpFormatter)
    optional = parser._action_groups.pop()
    required = parser.add_argument_group('required arguments')
    required.add_argument("-i", "--input", metavar="", type=str, help="File paths to both required files", required=True, nargs=2)
    optional.add_argument("-o", "--output", metavar="", type=str, help="File output path")
    optional.add_argument("-d", "--gap", metavar="", type=check_positive, help="give gap penalty")
    optional.add_argument("-ms", "--match", metavar="", type=check_positive, help="give match score")
    optional.add_argument("-mms", "--mismatch", metavar="", type=check_negative, help="give mismatch score")
    optional.add_argument("-single", metavar="", help="one best output")
    optional.add_argument("-all", metavar="", help="all best alignments as output")
    parser._action_groups.append(optional)
    return parser


def check_negative(value):
    # change 2
    val = int(value)
    if val >=0:
        raise argparse.ArgumentTypeError("%s is an invalid positive in value")
    return val


def check_positive(value):
    val = int(value)
    if val < 0:
        raise argparse.ArgumentTypeError("%s is an invalid negative in value")
    return val


def main():
    # parser for user input
    parser = init_parser()
    args = parser.parse_args()

    print(args.gap)

    print("Program: global aligner with linear gap cost")
    print("Author:", __author__)

    score = Score()
    if args.gap or args.gap == 0: score.gap = args.gap
    if args.match or args.match == 0: score.match = args.match
    if args.mismatch: score.mismatch = args.mismatch

    print("Parameters used in the algorithm:")
    print("\tMatch Score: " + str(score.match))
    print("\tMismatch Score: " + str(score.mismatch))
    print("\tGap Penalty: " + str(score.gap))

    x_header, x_seq = seq_file_reader(args.input[0])
    y_header,y_seq = seq_file_reader(args.input[1])

    # aligned_x, aligned_y = run_nw(x_seq[0:5], y_seq[0:5], score)
    aligned_x, aligned_y = run_nw(x_seq, y_seq, score)

    print("Alignment:")
    print("\n".join([x_header, aligned_x, y_header, aligned_y]))

    if args.output:
        f = open(args.output,"w+")
        f.write("Alignment:" + "\n")
        f.write(x_header + "\n")
        f.write(aligned_x + "\n")
        f.write(y_header + "\n")
        f.write(aligned_y + "\n")
        f.write("Best Score: {}".format(score.best_score))
        f.close()


if __name__ == "__main__":
    main()
