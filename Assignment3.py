#!/usr/bin/python

"""Basic Needleman-Wunsch"""

__author__ = "Yanpeng Li, Finn Mier"



import logging
from FastA import *
import argparse
from HashSearcher import HashSearcher
import Alphabets

logger = logging.getLogger()


def parse_args():
    parser = argparse.ArgumentParser(description="Basic aligner using hashing")
    parser.add_argument('-i', '--input', nargs = 2, type=str, required=True, help="Databasefile and query file")
    parser.add_argument('-o', '--output', nargs = 1, default=["out.out"], type=str, help="Name of the Output file. Defaults to alignment.fasta. If the file already exists an error will be thrown.")
    parser.add_argument('-l', '--loglevel', nargs = 1, default="warning", type=str, choices=['debug', 'info', 'warning', 'error', 'critical'], help="Sets logging level. Defaults to warning.")
    parser.add_argument('-k', '--tuple_length', nargs = 1, default=5, type=int, help="Length of tuples used to find seeds.")
    parser.add_argument('-a', '--alphabet', nargs = 1, type=str, choices=["amino", "DNA"], help="Sets the alphabet. If it is not given the alphabet will be recovered from the database.")
    parser.add_argument('-n', '--near_perfect',  action="store_true", help="When enabled will also search for hits with one mismatch.")
    parser.add_argument('-?', action="help", help="Shows this help message and exits.")

    args = parser.parse_args()
    return args


def logging_setup(loglevel : str): 
    # setting logging stuff
    loglevel = loglevel.lower()
    logging_map = {"debug" : logging.DEBUG, "info" : logging.INFO, "warning" : logging.WARNING, "error" : logging.ERROR, "critical" : logging.CRITICAL}
    level = logging_map.get(loglevel, logging.INFO)
    FORMAT = '%(levelname)s: %(filename)s %(lineno)d, %(funcName)s: %(message)s'
    logging.basicConfig(level = level, format = FORMAT)


def main():
    #parse input arguments
    args = parse_args()
    logging_setup(args.loglevel[0])

    # print some information about the script
    print("Program: global aligner with affine gap cost")
    print("Author: {}".format(__author__))

    headers, sequences = read_sequence(args.input[0])
    if args.alphabet and args.alphabet[0] == "DNA": alphabet = Alphabets.DNA
    elif args.alphabet and args.alphabet[0] == "amino": alphabet = Alphabets.amino_acid
    else: alphabet = None
    hash_searcher = HashSearcher(sequences, alphabet, args.tuple_length[0])

    headers, queries = read_sequence(args.input[1])
    print("Writing results to {}".format(args.output[0]))
    with open(args.output[0], "w") as file:
        for query_index, query in enumerate(queries):
        #seq_id, pos - starting_index, pos, length
            for seq_id, offset, pos, length in hash_searcher.search_sequence(query, args.near_perfect):
                print("{}\t{}\t{}\t{}\t{}".format(query_index, seq_id, -(offset - pos), pos, length), file = file)
        


if __name__ == "__main__":
    main()
