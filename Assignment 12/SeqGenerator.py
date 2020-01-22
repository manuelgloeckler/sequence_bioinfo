alphabet = ["A", "C", "G", "T"]

import numpy as np
class SeqGenerator:
    """
    Class to randomly generate sequences. Based on Jukes-Cantor model.
    Sequences will generally be represented as vectors of ints in [0, length of alphanet].
    In order to transform them to nucleotide sequences call vector_to_seq.

    Args:
        mutation_rate(float): Rate of mutation for Jukes-Cantor model.
    """
    def __init__(self, mutation_rate):
        self.mutation_rate = mutation_rate


    def generate_sequence(self, sequence_length):
        """
        Function to randomly generate a sequence from scratch. Each nucleotide is assumend to be equally likely.

        Args:
            sequence_length(int): Length of the sequence that is to be generated.        
  
        Returns:
            Randomly generated sequence.
        """
        return np.random.randint(0,len(alphabet),sequence_length)

    def mutate_sequence(self, reference, passed_time):
        """
        Function to mutate a given sequence based on the Jukes-Cantor model.

        Args:
            reference(np.ndarray): Reference sequence represented as vector of ints in [0, length of alphabet].
            time(float): Time passed used to Jukes-Cantor model.

        Returns:
            Mutated sequence based on given time and reference and Jukes-Cantor model.

        """
        changes = np.random.rand(len(reference)) < 3/4 * (1 - np.exp(-4/3 * self.mutation_rate * passed_time))
        new_bases = np.random.randint(1,len(alphabet), np.sum(changes))
        new_sequence = reference.copy()
        new_sequence[changes] += new_bases
        new_sequence[changes] %= len(alphabet)
        return new_sequence

    def vector_to_seq(self, v):
        """
        Function to turn vector of ints into nucleotides.

        Args:
            v(np.ndarray): Vector of ints to be turned into nucleotide sequence.

        Returns:
            string representing the nucleotide sequence represented by v.

        """
        return "".join([alphabet[i] for i in v])

if __name__ == "__main__":
    seq_gen = SeqGenerator(0.1)
    seq = seq_gen.generateSequence(100)
    print(seq_gen.vector_to_seq(seq))
    mut = seq_gen.mutate_sequence(seq, 0.5)
    print(seq_gen.vector_to_seq(mut))
