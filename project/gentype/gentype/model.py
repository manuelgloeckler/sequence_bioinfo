import numpy as np

class CategoricalMM:
    # Final Model class... TODO implement
    def __init__(self, pi, theta):
        self.pi = np.array(pi)
        self.theta = np.array(theta)

    def p(self, x):
        """ Return probability """
        logpx = x @ np.log(self.theta.T)
        p =  np.sum(self.pi * np.exp(logpx))
        return p

    def simple_sample(self, draws, replace = True):
        k = np.random.choice(range(len(self.pi)), p = self.pi)
        theta = self.theta[k]

        return np.random.choice(range(len(theta)), p = theta, size = draws, replace = replace)

    @staticmethod
    def ML_est_pi(Z):
        K = np.max(Z)+1
        pi = np.array([np.sum(Z == k) for k in range(int(K))])
        pi = pi / pi.sum()

        return pi
    

class AlleleMM(CategoricalMM):
    def __init__(self, Z, theta, inference_matrix, variant_ranges = None, variant_map = None):
        pi = CategoricalMM.ML_est_pi(Z)
        super().__init__(pi, theta)
        self.Z = Z
        self.clusters = int(max(self.Z)) + 1
        self.inference_matrix = inference_matrix
        self._compute_variation_distribution()        
        self.overlaps = {}
        if not (variant_ranges is None or variant_map is None):
            self.compute_variant_overlaps(variant_ranges, variant_map)


    def sample_variations(self, cluster = None):
        if cluster is None:
            k = np.random.choice(range(len(self.pi)), p = self.pi)
        else:
            k = cluster
        theta = np.array(self.theta[k])
        choices, probabilities = self.distributions[k]
        draws = np.random.choice(choices, p = probabilities)
        samples = []
        for i in np.arange(draws):
            sample = np.random.choice(range(len(theta)), p = theta)
            theta[self.overlaps[sample]] = 0
            theta = theta / sum(theta)
            samples.append(sample)
        return samples




    def _compute_variation_distribution(self):
        self.distributions = [None] * self.clusters
        for i in range(self.clusters):
            self.distributions[i] = {}
            variant_count = np.sum(self.inference_matrix[self.Z == i], axis=1)
            for count in variant_count:
                self.distributions[i][count] = self.distributions[i].get(count, 0) + 1 

            choices = []
            probabilities = []
            for choice in self.distributions[i]:
                choices.append(choice)
                probabilities.append(self.distributions[i][choice])
            probabilities = np.array(probabilities) / sum(probabilities)
            self.distributions[i] = (choices, probabilities)

    def compute_variant_overlaps(self, variant_ranges, variant_map):
        while variant_ranges:
            variant_id, start, end = variant_ranges.pop()
            variant_row = variant_map[variant_id]
            self.overlaps[variant_row] = self.overlaps.get(variant_row, [])
            self.overlaps[variant_row].append(variant_row)
            for other_id, other_start, other_end in variant_ranges:
                if start < other_end and other_start < end: # variants are stored with excluded end thus we use < (and not <=)
                    other_row = variant_map[other_id]
                    self.overlaps[variant_row] = self.overlaps.get(variant_row, [])
                    self.overlaps[variant_row].append(other_row)
                    self.overlaps[other_row] = self.overlaps.get(other_row, [])
                    self.overlaps[other_row].append(variant_row)
                

