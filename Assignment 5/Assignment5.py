from FastA import *
import numpy as np

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

def hamming_dist(str1, str2):
    assert (len(str1) != len(str2))
    score = 0
    for i in range(len(str1)):
        if str1[i] != str2[i] and (str1[i] != "-" or str1[i] != "-"):
            score +=1
    return score

def kimmura_dist(str1, str2):
    hamm_dist = hamming_dist(str1, str2)
    return -mp.ln(1-hamm_dist - hamm_dist**2/5)

def main():
    blossumMatrix = parse_score('./BLOSUM62.txt')
    headers, sequences = read_sequence('./BB11007_unaligned.fasta')


if __name__ == "__main__":
    main()