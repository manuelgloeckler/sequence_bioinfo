import numpy as np
from Bio import pairwise2
from skbio import DistanceMatrix
from skbio.tree import nj

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

def distance_matrix(strs, score):
    n = len(strs)
    dist_matrix = np.zeros((n,n))
    for i in range(n):
        for j in range(n):
            alignment = pairwise2.align.globaldx(strs[i], strs[j], score)
            seq1 = alignment[0][0]
            seq2 = alignment[0][1]
            dist_matrix[i,j] = abs(kimmura_dist(seq1, seq2))
    return dist_matrix

def pair_guided_alignment(align1:list, align2:list, score):
    seq1 = np.random.choice(align1)
    seq2 = np.random.choice(align2)
    new_alignment = pairwise2.align.globaldx(seq1, seq2, score)[0]
    new_seq1 = new_alignment[0]
    new_seq2 = new_alignment[1]

    newgaps1 =  np.where(np.array(list(new_seq1)) == "-")[0]
    newgaps2 = np.where(np.array(list(new_seq2)) == "-")[0]

    oldgaps1 =  np.where(np.array(list(seq1)) == "-")[0]
    oldgaps2 = np.where(np.array(list(seq2)) == "-")[0]

    newgaps1 = np.delete(newgaps1, np.where(newgaps1 == oldgaps1)[0])
    newgaps2 = np.delete(newgaps2, np.where(newgaps2 == oldgaps2)[0])


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
        print(node.name)

    node_names.reverse()
    return node_names


def main():
    blossumMatrix = parse_score('./blosum62matrix.txt')
    #headers, sequences = read_sequence('./BB11007_unaligned.fasta')
    raw1 = "AGGA"
    raw2 = "GGAC"
    raw3 = "AAAA"
    raw4 = "CCAC"
    raw5 = "GAGA"
    raw6 = "GTAC"
    raw7 = "GTTT"
    raw8 = "GTTG"
    data_dm = distance_matrix([raw1, raw2, raw3,raw4,raw5, raw6,raw7,raw8], blossumMatrix)
    dm = DistanceMatrix(list(data_dm), ["blaa","b","c","d","ganzlangeid","f",'was','wer'])
    tree = nj(dm).root_at_midpoint()
    print(tree.ascii_art())
    joining_list =  get_joining_list(tree)
    id_sequence_dict = dict(zip(["blaa","b","c","d","ganzlangeid","f",'was','wer'],[raw1, raw2, raw3,raw4,raw5, raw6,raw7,raw8]))

    qu = []
    i = 0
    while i < len(joining_list) -1:
        print(i)
        current = joining_list[i]
        next = joining_list[i+1]
        print(current)
        print(next)

        if current in id_sequence_dict.keys() and next in id_sequence_dict.keys():
            seqs1, seqs2 = pair_guided_alignment([id_sequence_dict[current]], [id_sequence_dict[next]], blossumMatrix)
            seqs1.extend(seqs2)
            qu.append(seqs1)
            print("aaaaaaaaaaaaaaa")
            i += 2

        if current == None and next in id_sequence_dict.keys():
            prof1 = qu[0]
            del qu[0]
            seqs1, seqs2 = pair_guided_alignment(prof1, [id_sequence_dict], blossumMatrix)
            seqs1.extend(seqs2)
            qu.append(seqs1)
            i += 2
        if current == None and next == None:
            prof1 = qu[0]
            del qu[0]
            prof2 = qu[0]
            del qu[0]
            seqs1, seqs2 = pair_guided_alignment(prof1, prof2, blossumMatrix)
            seqs1.extend(seqs2)
            qu.append(seqs1)
            i += 2

        i += 1

    print(qu)




# print(pair_guided_alignment(subal1, subal2, blossumMatrix))



if __name__ == "__main__":
    main()