from FastA import *
import numpy as np
from Bio import pairwise2

def parse_score(path):
    file = open(path, 'r')
    score_dict = {}
    alphabet = []
    i = 0
    for x in file:
        if i == 0:
            alphabet = x.split()
            alphabet = [x.strip() for x in alphabet]
        else:
            row = x.split()
            subst = row[0]
            del row[0]
            for j in range(len(row)):
                score_dict[(alphabet[j], subst)] = int(row[j].strip())
        i+=1
    return score_dict

def hamming_dist(str1, str2):
    #assert (len(str1) != len(str2))
    score = 0
    for i in range(len(str1)):
        if str1[i] != str2[i] and (str1[i] != "-" or str1[i] != "-"):
            score +=1
    return score/len(str1)

def kimmura_dist(str1, str2):
    hamm_dist = hamming_dist(str1, str2)
    return -np.log(1-hamm_dist - hamm_dist**2/5)

def distance_matrix(strs):
    n = len(strs)
    dist_matrix = np.zeros((n,n))
    for i in range(n):
        for j in range(n):
            dist_matrix[i,j] = kimmura_dist(run_nw(strs[i], strs[j])[:2])
    return dist_matrix

def pair_guided_alignment(align1:list, align2:list, score):
    seq1 = np.random.choice(align1)
    seq2 = np.random.choice(align2)
    new_alignment = pairwise2.align.globaldx(seq1, seq2, score)[0]
    new_seq1 = new_alignment[0]
    new_seq2 = new_alignment[1]

    gaps1 = np.where(np.array(list(seq1)) == "-")[0]
    gaps2 = np.where(np.array(list(seq2)) == "-")[0]
    newgaps1 =  np.where(np.array(list(seq1)) == "-")[0]





    for i in range(len(seq1)):


def main():
    blossumMatrix = parse_score('./blosum62.txt')
    print(blossumMatrix)
    #headers, sequences = read_sequence('./BB11007_unaligned.fasta')
    align = pairwise2.align.globaldx("A-GG", "G--GA", blossumMatrix)
    seq1 = align[0][0]
    seq2 = align[0][1]
    print(seq1,seq2)
    print(hamming_dist(seq1, seq2))
    print(kimmura_dist(seq1,seq2))


if __name__ == "__main__":
    main()