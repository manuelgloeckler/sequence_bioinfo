import numpy as np

class CategoricalMM:
    """
    General categorical mixture model.

    Args:
        pi (iterable): Iterable containing the probabilities that a cluster is chosen.
        theta (iterable): Iterable containing the probabilities for each of the variants within the clusters.
    """
    def __init__(self, pi, theta):
        self.pi = np.array(pi)
        self.theta = np.array(theta)

    def p(self, x):
        """ Return probability """
        logpx = x @ np.log(self.theta.T)
        p =  np.sum(self.pi * np.exp(logpx))
        return p

    def simple_sample(self, draws, replace = True):
        """
        Draws a given amount of samples from this model.

        Args:
            draws (int): Number of draws.
            replace (bool, optional): If set to False draws will be done without replacement.
                Defaults to True.
        """
        k = np.random.choice(range(len(self.pi)), p = self.pi)
        theta = self.theta[k]

        return np.random.choice(range(len(theta)), p = theta, size = draws, replace = replace)

    @staticmethod
    def ML_est_pi(Z):
        """
        Compute pi (probabilites that a cluster is chosen) from the assignment of individuals to each cluster (Z).

        Args:
            Z (np.ndarray): Assignment for each individual to each respective cluster.
        """
        K = np.max(Z)+1
        pi = np.array([np.sum(Z == k) for k in range(int(K))])
        pi = pi / pi.sum()

        return pi
    

class AlleleMM(CategoricalMM):
    """
    Specialized catigorical mixture model which support drawing sets of variants.

    Args:
        Z (np.ndarray): Assignment for each individual to each respective cluster.
        theta (iterable): Iterable containing the probabilities for each of the variants within the clusters.
        inference_matrix (np.ndarray): Inference matrix which was used to construct the model from which the
            parameters were derived.
        variant_ranges ([(variant_id, variant_start, variant_end)]): List containing tuples assigning each
            variant its start and end position within the reference sequence. If this is not given in the constructor
            self.compute_variant_overlaps must be called before sampling.
        variant_map ({variant_id : idx}): Dictionary mapping each variant id to its index in the inference matrix. 
            If this is not given in the constructor self.compute_variant_overlaps must be called before sampling.
    """
    def __init__(self, Z, theta, inference_matrix, variant_ranges = None, variant_map = None):
        pi = CategoricalMM.ML_est_pi(Z)
        super().__init__(pi, theta)
        self.Z = Z
        self.clusters = int(max(self.Z)) + 1
        self.inference_matrix = inference_matrix
        self._compute_variation_distribution()        
        if not (variant_ranges is None or variant_map is None):
            self.compute_variant_overlaps(variant_ranges, variant_map)


    def sample_variations(self, cluster = None, return_cluster = False):
        """
        Samples a set of variations from this model. Returns them as a bit vector (use the variant_map 
        from the inference_matrix creation to convert them to actual variant ids).

        Args:
            cluster (int, optional): Index of the cluster to sample from. This should only provided if you want to 
                poll from a specific cluster. If not given will sample properly from the entire model.
            return_cluster (bool, optional): If set to true will also return the cluster this sample stems from.

        Returns:
            Bitvector containing a 1 if the variation is expressed and 0 if not (use the variant_map 
                from the inference_matrix creation to convert them to actual variant ids).
                If return_cluster is True this will also return the cluster the sample stems from.
        """
        if cluster is None:
            k = np.random.choice(range(len(self.pi)), p = self.pi)
        else:
            k = cluster
        theta = np.array(self.theta[k])
        choices, probabilities = self.distributions[k]
        draws = np.random.choice(choices, p = probabilities)
        variations = []
        for i in np.arange(draws):
            variation = np.random.choice(range(len(theta)), p = theta)
            theta[self.overlaps[variation]] = 0
            theta = theta / sum(theta)
            variations.append(variation)
        if return_cluster:
            return variations, k
        else:
            return variations


    def _compute_variation_distribution(self):
        """
        Computes the distribution of the amount of variations per individual within each cluster.
        """
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
        """
        Computes the overlaps for each of the variants. This is important so no two overlapping variants will be sampled.

        Args:
            variant_ranges ([(variant_id, variant_start, variant_end)]): List containing tuples assigning each
                variant its start and end position within the reference sequence.
            variant_map ({variant_id : idx}): Dictionary mapping each variant id to its index in the inference matrix. 
        """
        self.overlaps = {}
        variant_ranges = list(variant_ranges) # make copy to not destroy given list
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
                

