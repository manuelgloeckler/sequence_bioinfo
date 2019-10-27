#!/usr/bin/python
import argparse
import warnings
import numpy as np

"""Basic Needleman-Wunsch"""

__author__ = "Patrick Schirm, Manuel Gloeckler"

"""
The script is created by the above-mentioned two students for the same course 1 year ago
This script implements basic Needleman-Wunsch algorithm.
"""


class Score:
    """scores are stored here
    Hints: you may not need it
    """

    def __init__(self, match=1, mismatch=-1, gap=4, best_score=None):
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
                score_dict[(alphabet[j], subst)] = int(row[j].strip())
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
            elif i == 0 and j != 0:
                M[i][j] = -np.inf
                I_x[i][j] = -np.inf
                I_y[i][j] = -d - (j - 1)*e
            elif j == 0 and i != 0:
                M[i][j] = -np.inf
                I_x[i][j] = -d - (i - 1) *e
                I_y[i][j] = -np.inf
            else:
                if isinstance(score, Score):
                    s = score.match if x_seq[j - 1] == y_seq[i - 1] else score.mismatch
                elif isinstance(score, dict):
                    s = score[(x_seq[j - 1], y_seq[i - 1])]
                else:
                    s = 0
                    warnings.warn("No compatible scoring was given, score = 0")
                M[i][j] = max(M[i - 1][j - 1] + s, I_x[i - 1][j - 1] + s, I_y[i - 1][j - 1] + s)
                I_x[i][j] = max(M[i - 1][j] - d, I_x[i - 1][j] - e)
                I_y[i][j] = max(M[i][j - 1] - d, I_y[i][j - 1] - e)

            # diag = matF[i - 1][j - 1] + s
                # top = matF[i - 1][j] - score.gap
                # left = matF[i][j - 1] - score.gap
                # matF[i][j], matT[i][j] = _tracking_guide(diag, top, left)

    return M, I_x, I_y


def gotho_traceback(M, I_x, I_y, x_seq, y_seq, d, e):
    ncol = len(x_seq) + 1
    nrow = len(y_seq) + 1
    result_x = ""
    result_y = ""

    i = nrow -1
    j = ncol -1

    # Picking max end score and the corresponding matrix where it is reached.
    maxs = [(M[i][j]), (I_x[i][j]), (I_y[i][j] )]
    score = max(maxs)
    matrix = maxs.index(max(maxs))

    while i > 0 and j > 0:
        # Dependent on the current matrix, compute the maximum values of the step before
        if matrix == 0:
            vals = np.array([(M[i-1][j-1]), (I_x[i-1][j-1]), (I_y[i-1][j-1])])
            trace = list(np.where(vals == max(vals))[0])
        elif matrix == 1:
            vals = np.array([M[i-1][j] - d, I_x[i-1][j] -e, - np.inf])
            trace = list(np.where(vals == max(vals))[0])
        elif matrix == 2:
            vals = np.array([M[i][j-1] - d, -np.inf, I_y[i][j-1] - e])
            trace = list(np.where(vals == max(vals))[0])

        # Dependent on the trace of the maximum value and the orignin matrix, add corresponding alignment and
        # change matrix if needed.
        if 0 in trace:
            if matrix == 1:
                result_x += "-"
                result_y += y_seq[i-1]
                i -= 1
            elif matrix == 2:
                result_x += x_seq[j-1]
                result_y += "-"
                j -= 1
            else:
                result_x += x_seq[j-1]
                result_y += y_seq[i-1]
                i -= 1
                j -= 1
            matrix = 0
        elif 1 in trace:
            if matrix == 0:
                result_x += x_seq[j - 1]
                result_y += y_seq[i - 1]
                i -= 1
                j -= 1
            else:
                result_x += "-"
                result_y += y_seq[i-1]
                i -= 1
            matrix = 1
        elif 2 in trace:
            if matrix == 0:
                result_x += x_seq[j - 1]
                result_y += y_seq[i - 1]
                i -= 1
                j -= 1
            else:
                result_x += x_seq[j-1]
                result_y += "-"
                j -= 1
            matrix = 2

    return score, result_x[::-1], result_y[::-1]


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


def run_gotho(x_seq, y_seq, score, d, e):
    # Dummy checks
    if (d < e):
        warnings.warn("You choice d < e does not make sense, rethink your parameters")

    if isinstance(score, Score):
        if not -d - e < score.mismatch:
            warnings.warn("d-e should be smaller than the minimum mismatch score")

    if isinstance(score, dict):
        if not -d - e < min(score.values()):
            warnings.warn("d-e should be smaller than the minimum mismatch score")

    if not (isinstance(score, Score) or isinstance(score, dict)):
        raise Exception("Incompatible scoring function...")

    if not isinstance(x_seq, str) and isinstance(y_seq, str):
        raise Exception("Incompatible sequence formats")

    M, I_x, I_y = gotho_matrix_builder(x_seq, y_seq, score, d, e)
    return gotho_traceback(M, I_x, I_y, x_seq, y_seq, d, e)


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
    parser = argparse.ArgumentParser(description='Needleman-Wunsch',prog="pqa", formatter_class=argparse.RawDescriptionHelpFormatter)
    optional = parser._action_groups.pop()
    required = parser.add_argument_group('required arguments')
    required.add_argument("-i", "--input", metavar="", type=str, help="File paths to both required files", required=True, nargs=2)
    optional.add_argument("-o", "--output", metavar="", type=str, help="File output path")
    # optional.add_argument("-single", metavar="", help="one best output")
    # optional.add_argument("-all", metavar="", help="all best alignments as output")
    optional.add_argument("-d", "--gap", metavar="", type=check_positive, help="give gap penalty")
    optional.add_argument("-ad", "--affine", action="store_true", help="give affine gap penalty, use -open and -extend for penalties as positive integers")
    optional.add_argument("-open", metavar="", type=check_positive, help="If affine penalty, please specify the gap open cost")
    optional.add_argument("-extend", metavar="", type=check_positive, help="If affine penalty, please specify the gap extension cost")
    optional.add_argument("-ms", "--match", metavar="", type=check_positive, help="give match score")
    optional.add_argument("-mms", "--mismatch", metavar="", type=check_negative, help="give mismatch score")
    optional.add_argument('-sm', '--matrix', action="store_true", help="use the blossum62 matrix for (miss) match scoring")
    parser._action_groups.append(optional)
    return parser


def check_negative(value):
    val = int(value)
    if val >=0:
        raise argparse.ArgumentTypeError("%s is an invalid positive in value")
    return val


def check_positive(value):
    val = int(value)
    if val < 0:
        raise argparse.ArgumentTypeError("%s is an invalid negative in value")
    return val


def check_args(args, parser):
    if (not args.affine) and (not args.gap):
        parser.error("--affine or --gap is required")
    if (not args.mismatch or not args.match) and (not args.matrix):
        parser.error("--match and --mismatch together or --matrix is required")
    if args.affine and (not args.open or (not args.extend)):
        parser.error("--affine requires -open and -extend")
    if args.affine and args.gap:
        parser.error("--affine and --gap can not be used together")
    if (args.match or args.mismatch) and (args.matrix):
        parser.error("--match or --mismatch can not be used with --matrix")
    if args.gap and (args.open or args.extend):
        parser.error("when using the -gap parameter, -open and -extend can not be used")
    if (args.extend and args.open) and (args.extend > args.open):
        parser.error("extend cost can not be greater than opening cost")


def main():
    parser = init_parser()
    args = parser.parse_args()
    check_args(args, parser)

    print("Program: global aligner with linear gap cost")
    print("Author:", __author__)

    score = Score()
    if args.gap or args.gap == 0: score.gap = args.gap
    if args.match or args.match == 0: score.match = args.match
    if args.mismatch: score.mismatch = args.mismatch
    if args.matrix:
        score = parse_score('./BLOSUM62.txt')

    x_header, x_seq = seq_file_reader(args.input[0])
    y_header,y_seq = seq_file_reader(args.input[1])

    if args.gap:
        aligned_x, aligned_y = run_nw(x_seq, y_seq, score)
        # print("Alignment with linear cost:")
        # print("\n".join([x_header, aligned_x, y_header, aligned_y]))
        # print(" Best score: " + str(score.best_score))
        if args.output:
            f = open(args.output,"w+")
            f.write("Alignment with linear cost:" + "\n")
            f.write(x_header + "\n")
            f.write(aligned_x + "\n")
            f.write(y_header + "\n")
            f.write(aligned_y + "\n")
            f.write("Best Score: {}".format(score.best_score))
            f.close()
    else:
        M, I_x, I_y = gotho_matrix_builder(x_seq,y_seq,score,args.open,args.extend)
        score, a_x, a_y = gotho_traceback(M, I_x, I_y, x_seq, y_seq,args.open,args.extend)
        # print("Alignment with affine cost:")
        # print(a_x)
        # print(a_y)
        # print("Best score: " + str(score))
        if args.output:
            f = open(args.output,"w+")
            f.write("Alignment with affine cost:" + "\n")
            f.write(x_header + "\n")
            f.write(a_x + "\n")
            f.write(y_header + "\n")
            f.write(a_y + "\n")
            f.write("Best score: " + str(score))
            f.close()

    print("finished")


if __name__ == "__main__":
    main()
