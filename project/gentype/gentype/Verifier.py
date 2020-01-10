import numpy as np
import sys

def get_distances(samples, test_set):
    # compute all occurences for every test sequence
    occurences = {}
    test_set = tuple(test_set)
    for individual in test_set:
        individual = tuple(individual)
        occurences[individual] = occurences.get(individual, 0) + 1

    # assign each sample sequence to one of the test set
    assignments = {}
    assigned_scores = {}
    for sample in samples:
        sample = tuple(sample)
        cur_score = sys.maxsize
        cur_seq = None
        scores = []
        for individual in test_set:
            score = sum([int(x)^int(y) for x, y in zip(sample, individual)])
            scores.append(score)
            if score < cur_score:
                cur_score = score
                cur_seq = individual
                cur_comp = [int(x)^int(y) for x, y in zip(sample, individual)]
        #print("seq1: {}\n seq2: {} \n comparison: {}\n score: {}, scores : {}".format("".join([str(int(x)) for  x in cur_seq]), "".join([str(int(x)) for  x in sample]),"".join([str(int(x)) for  x in cur_comp]) , cur_score, scores))
        assignments[sample] = cur_seq
        assigned_scores[sample] = cur_score
    inverted_assignments = {}
    for key in assignments:
        cur_values = inverted_assignments.get(key, [])
        cur_values.append(key)
        inverted_assignments[tuple(assignments[key])] = cur_values


    # compute similarity of the distributions
    distribution_differences = []
    for individual in occurences:
        test_quota = occurences[individual] / len(test_set)
        generated_quota = len(inverted_assignments.get(individual, [])) / len(samples)
        distribution_differences.append(np.absolute(test_quota - generated_quota))
    distribution_difference = sum(distribution_differences)

    return assigned_scores.values(), distribution_difference