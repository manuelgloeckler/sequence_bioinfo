import logging

import numpy as np


logger = logging.getLogger()



class HashSearcher:
    """Class which allows searching a query string in a database.

    Parameters:
        database ([str]): list of sequences which will be searched through.
        alphabet (str, optional): List of characters which appear in the database. If it is not given, it is derived from the given database.
        k (int, optional): Length of the k-tuples the database will be stored as. Defaults to 5.
    """
    def __init__(self, database, alphabet = None, k = 5):
        self.database = database
        db_alphabet = self._get_alphabet_from_db()
        if alphabet is None: 
            logger.info("No alphabet was given. Recovered {} as alphabet.".format(self.alphabet))
            alphabet = db_alphabet
        else:
            if not set(alphabet).issubset(db_alphabet):
                logger.error("Given alphabet: {} did not correspond to alphabet of database: {}".format(alphabet, db_alphabet))
                raise ValueError
            else: self.alphabet = alphabet

        self.mapping = {a : i for i,a in enumerate(self.alphabet)}
        self.k = k
        self.build_hashtable()

    def hash_function(self, sequence): 
        kappa = len(self.alphabet)
        k = len(sequence) - 1
        return sum([(kappa**(k - i)) * self.mapping[c] for i,c in enumerate(sequence)])

    def near_hashes(self, hash_value, sequence):
        near_hashes = []
        kappa = len(self.alphabet)
        k = len(sequence) - 1
        for i in range(len(sequence)):
            base_hash = hash_value - (self.mapping[sequence[i]] * (kappa**(k-i)))
            near_hashes.extend([base_hash + (j * (kappa**(k-i))) for j in range(kappa)])
        near_hashes = set(near_hashes); near_hashes = list(near_hashes) # remove duplicates
        return near_hashes


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


    def search_sequence(self, query, allow_mismatch = False):
        #sanity check given query
        query_alphabet = set()
        query_alphabet.update(query)
        if not set(self.alphabet).issubset(query_alphabet):
            logger.error("Query alphabet: {} did not correspond to alphabet of database: {}".format(query_alphabet, self.alphabet))
            raise ValueError

        # create all overlapping k-tuples from the query
        k_tuples = [query[i:i+self.k] for i in range(len(query) - self.k)]
        hits = []

        # collect all hits from the hash table
        for starting_index, k_tuple in enumerate(k_tuples):
            hash_value = self.hash_function(k_tuple)
            if allow_mismatch:
                near_hashes = self.near_hashes(hash_value, k_tuple)
                for hash_value in near_hashes:
                    hash_bucket = self.hash_table.get(hash_value, [])
                    hits.extend([(seq_id, pos - starting_index, pos) for (seq_id, pos) in hash_bucket])
            else:
                hash_bucket = self.hash_table.get(hash_value, [])
                hits.extend([(seq_id, pos - starting_index, pos) for (seq_id, pos) in hash_bucket])

        # combine consecutive hits
        hits.sort()
        greatest_hits = []
        cur_hit = None
        for hit in hits:
            if cur_hit and hit[0] == cur_hit[0] and hit[1] == cur_hit[1] and hit[2] == cur_hit[2] + cur_hit[-1]: # consecutive
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
        return chars




if __name__ == "__main__":
    hash_function = create_hash("ACGT")
    print(hash_function("AGCA"))