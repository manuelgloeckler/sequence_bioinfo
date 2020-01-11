import numpy as np
import sys

def get_distances(samples, test_set, return_assignments = False):
    # compute all occurences for every test sequence
    occurences = {}
    test_set = tuple(test_set)
    for individual in test_set:
        individual = tuple(individual)
        occurences[individual] = occurences.get(individual, 0) + 1
    sample_occurences = {}
    samples = tuple(samples)
    for sample in samples:
        sample = tuple(sample)
        sample_occurences[sample] = sample_occurences.get(sample, 0) + 1

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
        assignments[sample] = cur_seq
        assigned_scores[sample] = cur_score

    inverted_assignments = {}
    for key in assignments:
        cur_values = inverted_assignments.get(tuple(assignments[key]), [])
        cur_values.append(key)
        inverted_assignments[tuple(assignments[key])] = cur_values


    # compute similarity of the distributions
    distributions = {}
    for individual in occurences:
        test_quota = occurences[individual] / len(test_set)
        assigned_samples = 0
        for sample in inverted_assignments.get(individual, []):
            assigned_samples += sample_occurences[sample]
        generated_quota = assigned_samples / len(samples)
        distributions[individual] = (test_quota, generated_quota)

    if return_assignments:
        return assigned_scores, distributions, assignments
    else:
        return assigned_scores, distributions