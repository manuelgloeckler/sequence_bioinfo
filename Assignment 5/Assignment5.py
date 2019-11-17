import numpy as np
from Bio import pairwise2
from skbio import DistanceMatrix
from skbio.tree import nj
import logging
from MatrixReader import *
from GlobalSequenceLinearGap import GlobalSequenceLinearGapAligner
from SAScoringSystem import *
import re
import FastA

logger = logging.getLogger()

FORMAT = '%(levelname)s: %(filename)s %(lineno)d, %(funcName)s: %(message)s'
logging.basicConfig(level = logging.INFO, format = FORMAT)

def find_new_gaps(old_sequence, new_sequence):
    old_index = 0
    gap_indices = []
    for new_index in range(len(new_sequence)):
        if old_index < len(old_sequence) and new_sequence[new_index] == old_sequence[old_index]:
            old_index += 1
        else:
            gap_indices.append(new_index)
    return gap_indices




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
    assert (len(str1) == len(str2))
    score = 0
    for i in range(len(str1)):
        if str1[i] != str2[i] and (str1[i] != "-" or str1[i] != "-"):
            score +=1
    return score/len(str1)

def kimmura_dist(str1, str2):
    hamm_dist = hamming_dist(str1, str2)
    exp_val = 1-hamm_dist - ((hamm_dist**2)/5)
    if exp_val > 0:
        return -np.log(exp_val)
    else:
        return 1e10

def distance_matrix(strs, aligner):
    n = len(strs)
    dist_matrix = np.zeros((n,n))
    for i in range(n):
        for j in range(i,n):
            aligner.sequences = [strs[i], strs[j]]
            alignment = aligner.align()

            seq1 = alignment[0]
            seq2 = alignment[1]
            dist_matrix[i,j] = abs(kimmura_dist(seq1, seq2))
            dist_matrix[j,i] = abs(kimmura_dist(seq1, seq2))

    return dist_matrix
icq = []
def pair_guided_alignment(align1:list, align2:list, aligner):
    global icq
    seq1 = np.random.choice(align1)
    seq2 = np.random.choice(align2)
    aligner.sequences = [seq1, seq2]

    new_alignment = aligner.align()
    icq.extend(new_alignment)
    
    new_seq1 = new_alignment[0]
    new_seq2 = new_alignment[1]
    logger.info("old: {},\nnew: {}".format(seq1, new_seq1))

    newgaps1 = find_new_gaps(seq1, new_seq1)
    newgaps2 = find_new_gaps(seq2, new_seq2)

    logger.info("gaps1: {},\ngaps2: {}".format(newgaps1, newgaps2))


    for i in range(len(align1)):
        for col in newgaps1:
            align1[i] = align1[i][:col] + "-" + align1[i][col:]

    for i in range(len(align2)):
        for col in newgaps2:
            align2[i] = align2[i][:col] + "-" + align2[i][col:]

    return align1, align2


def get_joining_list(tree):
    node_names = []
    for node in tree.levelorder():
        node_names.append(node.name)

    node_names.reverse()
    return node_names


def main():
    alphabet, score_matrix = read_scoring_matrix("blosum62matrix.txt", re.compile("\s+"))
    aligner = GlobalSequenceLinearGapAligner(MatrixLinearScoring(alphabet, score_matrix, 4))

    headers, sequences = FastA.read_sequence('./BB11007_unaligned.fasta')
    headers_idx = [str(i) for i in range(len(sequences))]

    data_dm = distance_matrix(sequences, aligner)
    dm = DistanceMatrix(data_dm, headers_idx)
    tree = nj(dm).root_at_midpoint()
    joining_list =  get_joining_list(tree)
    id_sequence_dict = dict(zip(headers_idx,sequences))
    qu = []
    i = 0
    while i < len(joining_list) -1:

        current = joining_list[i]
        next_seq = joining_list[i+1]
        if current in id_sequence_dict.keys() and next_seq in id_sequence_dict.keys():
            seqs1, seqs2 = pair_guided_alignment([id_sequence_dict[current]], [id_sequence_dict[next_seq]], aligner)
            seqs1.extend(seqs2)
            qu.append(seqs1)

            i += 2
        elif current == None and next_seq in id_sequence_dict.keys():
            prof1 = qu[0]
            del qu[0]
            seqs1, seqs2 = pair_guided_alignment(prof1, [id_sequence_dict[next_seq]], aligner)
            seqs1.extend(seqs2)
            qu.append(seqs1)
            i += 2
        elif next_seq == None and current in id_sequence_dict.keys():
            prof1 = qu[0]
            del qu[0]
            seqs1, seqs2 = pair_guided_alignment(prof1, [id_sequence_dict[current]], aligner)
            seqs1.extend(seqs2)
            qu.append(seqs1)
            i += 2
        elif current == None and next_seq == None:
            prof1 = qu[0]
            del qu[0]
            prof2 = qu[0]
            del qu[0]
            seqs1, seqs2 = pair_guided_alignment(prof1, prof2, aligner)
            seqs1.extend(seqs2)
            qu.append(seqs1)
            i += 2
        else:
            i += 1
    FastA.write_comparison("test.fasta", headers, qu[0]) #TODO: sort headers correctly
    FastA.write_comparison("pairwise_alignments.fasta", ["alignment {}".format(i) for i in range(len(icq))], icq) #TODO: sort headers correctly









# print(pair_guided_alignment(subal1, subal2, blossumMatrix))



if __name__ == "__main__":
    main()