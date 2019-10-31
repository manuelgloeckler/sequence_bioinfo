import logging

logger = logging.getLogger()

def hash_seq(alphabet, seq):
    mapping = {a : i for i,a in enumerate(alphabet)}
    kappa = len(alphabet)
    k = len(seq) - 1
    return sum([(kappa**(k - i)) * mapping[c] for i,c in enumerate(seq)])



class HashSearcher:
    """
    Parameters:
        database ([str]): list of sequences which will be searched through
    """
    def __init__(self, database, k = 5):
        self.database = database
        self.k = k
        self.build_hashtable()
        

    def build_hashtable(self):
        for sequence in self.database:
            pass




if __name__ == "__main__":
    print(hash_seq("ACGT", "AGCA"))