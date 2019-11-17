import logging
import re
from typing import Pattern
logger = logging.getLogger()

__author__ = "Yanpeng Li, Josua Stadelmaier, Finn Mier"

def read_scoring_matrix(filename, delimiter = "\t"):
    """Reads a scoring_matrix from a file.
    Independent of the ordering of the file, this scoring matrix will be symmetric.

    Args:
        filename (str): Path to the file from which to read the matrix. The first row and column of the matrix
            should be the corresponding scored letters.
        delimiter (str, optional): Delimiter for the values in the matrix. Defaults to tab.

    Returns:
        A tuple containing a dictionary mapping scored letters to indices and a list of list containing the scores for each pair of letters. 
    """

    with open(filename) as scoring_file:
        if isinstance(delimiter, Pattern):
            alphabet = re.split(delimiter, scoring_file.readline().rstrip("\n"))[1:]
        else:
            alphabet = scoring_file.readline().rstrip("\n").split(delimiter)
        abc_to_idx = {alphabet[i] : i for i in range(len(alphabet))}
        scoring_matrix = [None] * len(alphabet)
        for line_nr, line in enumerate(scoring_file):
            if isinstance(delimiter, Pattern):
                separated = re.split(delimiter, line.rstrip("\n"))
            else:
                separated = line.rstrip("\n").split(delimiter)
            scoring_matrix[abc_to_idx[separated[0]]] = [int(x) for x in separated[1:]]
        return abc_to_idx, scoring_matrix
    

if __name__ == "__main__":
    print(read_scoring_matrix("Material-A2/BLOSUM62.txt"))