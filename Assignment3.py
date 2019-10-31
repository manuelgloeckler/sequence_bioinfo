#!/usr/bin/python

"""Basic Needleman-Wunsch"""

__author__ = "Yanpeng Li, Finn Mier"



import logging
from FastA import *
import argparse
from HashSearcher import HashSearcher

logger = logging.getLogger()


def parse_args():
    parser = argparse.ArgumentParser(description="Basic aligner using hashing")
    parser.add_argument('-i', '--input', nargs = 2, type=str, required=True, help="Databasefile and query file")
    parser.add_argument('-o', '--output', nargs = 1, default=["out.out"], type=str, help="Name of the Output file. Defaults to alignment.fasta. If the file already exists an error will be thrown.")
    parser.add_argument('-l', '--loglevel', nargs = 1, default="warning", type=str, choices=['debug', 'info', 'warning', 'error', 'critical'], help="Sets logging level. Defaults to warning.")
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
    hash_searcher = HashSearcher(sequences)

    headers, queries = read_sequence(args.input[1])
    print("Writing results to {}".format(args.output))
    with open(args.output[0], "w") as file:
        for query_index, query in enumerate(queries):
            #seq_id, pos - starting_index, pos, length
            for seq_id, offset, pos, length in hash_searcher.search_sequence(query):
                print("{}\t{}\t{}\t{}\t{}".format(query_index, seq_id, -(offset - pos), pos, length), file = file)


if __name__ == "__main__":
    main()
