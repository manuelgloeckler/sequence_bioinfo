import logging

import numpy as np


logger = logging.getLogger()

def create_hash(alphabet):
    mapping = {a : i for i,a in enumerate(alphabet)}
    kappa = len(alphabet)
    def hash_function(seq):
        k = len(seq) - 1
        return sum([(kappa**(k - i)) * mapping[c] for i,c in enumerate(seq)])
    return hash_function



class HashSearcher:
    """Class which allows searching a query string in a database.

    Parameters:
        database ([str]): list of sequences which will be searched through.
        alphabet (str, optional): List of characters which appear in the database. If it is not given, it is derived from the given database.
        k (int, optional): Length of the k-tuples the database will be stored as. Defaults to 5.
    """
    def __init__(self, database, alphabet = None, k = 5):
        self.database = database
        if alphabet is None: self._get_alphabet_from_db()
        else: self.alphabet = alphabet
        self.k = k
        self.hash_function = create_hash(self.alphabet)
        self.build_hashtable()
        

    def build_hashtable(self):
        """Function to build a hashtable from this instances database.
        The result is a dictionary mapping hash values to a list containing pairs describing the sequence as 
        well as the position the k-tuple was found at.
        """
        self.hash_table = {}
        for sequence_index, sequence in enumerate(self.database):
            for starting_position in range(0, len(sequence), self.k):
                k_tuple = sequence[starting_position:starting_position + self.k]
                hash_value = self.hash_function(k_tuple)
                if hash_value not in self.hash_table: self.hash_table[hash_value] = []
                hash_bucket = self.hash_table[hash_value]
                hash_bucket.append((sequence_index, starting_position))


    def search_sequence(self, query):
        # create all overlapping k-tuples from the query
        k_tuples = [query[i:i+self.k] for i in range(len(query) - self.k)]
        hits = []

        # collect all hits from the hash table
        for starting_index, k_tuple in enumerate(k_tuples):
            hash_value = self.hash_function(k_tuple)
            hash_bucket = self.hash_table.get(hash_value, [])
            hits.extend([(seq_id, pos - starting_index, pos) for (seq_id, pos) in hash_bucket])

        # combine consecutive hits
        hits.sort()
        greatest_hits = []
        cur_hit = None
        for hit in hits:
            if cur_hit and hit[0] == cur_hit[0] and hit[1] == cur_hit[1] and hit[2] == cur_hit[2] + self.k: # consecutive
                cur_hit[-1] += self.k
            else:
                if cur_hit: greatest_hits.append(tuple(cur_hit))
                cur_hit = [*hit, self.k]
        if cur_hit: greatest_hits.append(tuple(cur_hit))

        return greatest_hits
                
    
    def _get_alphabet_from_db(self):
        chars = set()
        for sequence in self.database:
            chars.update(sequence)
        self.alphabet = list(chars)
        logger.info("No alphabet was given. Recovered {} as alphabet.".format(self.alphabet))




if __name__ == "__main__":
    hash_function = create_hash("ACGT")
    print(hash_function("AGCA"))