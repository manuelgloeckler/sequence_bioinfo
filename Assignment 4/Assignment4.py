#!/usr/bin/python

"""Analysis tool for blat output files"""

__author__ = "Patrick Schirm, Manuel Glöckler, Finn Mier"

import argparse
from Blast8 import *
from TruePositive import read_TP
from FastA import read_sequence

import logging 
logger = logging.getLogger()


def parse_args():
    parser = argparse.ArgumentParser(description="Analysis tool for blat output files")
    parser.add_argument('-b', '--blat', type=str, required=True, help="blat output file")
    parser.add_argument('-c', '--criterion', type=str, default=None, choices=['e_value', 'length', 'bit_score', 'identical'], help="how to evaluate the 'best' hit.")
    parser.add_argument('-t', '--true_positive', type=str, required=True, help="true positive file")
    parser.add_argument('-d', '--data_base', type=str, help="fasta database file. Either this or data_base_length is required to compute the opposite alignments")
    parser.add_argument('-n', '--data_base_length', type=str, help="database length. Either this or data_base is required to compute the opposite alignments")
    parser.add_argument('-o', '--output', default="analysis.txt", type=str, help="Name of the Output file. Defaults to analysis.txt. If the file already exists an error will be thrown.")
    parser.add_argument('-l', '--loglevel', default="warning", type=str, choices=['debug', 'info', 'warning', 'error', 'critical'], help="Sets logging level. Defaults to warning.")
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
    logging_setup(args.loglevel)

    # print some information about the script
    print("Author: {}".format(__author__))

    true_positives = read_TP(args.true_positive)


    if args.data_base:
        headers, sequences = read_sequence(args.data_base)
        blast_reader = Blast8Reader({headers[i].lstrip(">"): len(sequences[i]) for i in range(len(sequences))})
    elif args.data_base_length:
        blast_reader = Blast8Reader(int(args.data_base_length))
    else:
        logger.warning("No data_base or data_base_length was passed. Thus alignments in negative directions can not be checked.")
        blast_reader = Blast8Reader()

#TODO: First Assignment task doesn't ask for perfect ones but for any!
    if args.criterion is None:
        blast8_reads = blast_reader.read_blast8(args.blat)
    if args.criterion == "e_value":
        blast8_reads = blast_reader.read_probable_blast8(args.blat, Blast8Reader.e_value)
    elif args.criterion == "length":
        blast8_reads = blast_reader.read_probable_blast8(args.blat, Blast8Reader.longest_sequence)
    elif args.criterion == "identical":
        blast8_reads = blast_reader.read_probable_blast8(args.blat, Blast8Reader.identical)
    elif args.criterion == "bit_score":
        blast8_reads = blast_reader.read_probable_blast8(args.blat, Blast8Reader.bit_score)

    found_positives = 0
    false_positives = 0
    found_keys = set()
    for tp_key in true_positives:
        if tp_key in blast8_reads:
            if type(blast8_reads[tp_key][0]) == list:
                found = False
                for read in blast8_reads[tp_key]:
                    if read[10] == true_positives[tp_key][0] and read[11] == true_positives[tp_key][1]:
                        found_keys.add(tp_key)
                        found_positives += 1
                        found = True
                        break
                if found == False: false_positives += 1
            else:
                if blast8_reads[tp_key][10] == true_positives[tp_key][0] and blast8_reads[tp_key][11] == true_positives[tp_key][1]:
                    found_keys.add(tp_key)
                    found_positives += 1
                else:
                    false_positives += 1


    print("Found positives: {}/{} with percentage: {}.".format(found_positives, len(true_positives), found_positives/len(true_positives)))

    uniquely_mapped = 0
    correctly_uniquely_mapped = 0
    for read_key in blast8_reads:
        if (type(blast8_reads[read_key][-1]) == list and len(blast8_reads[read_key]) == 1) or blast8_reads[read_key][-1] == 1:
            uniquely_mapped += 1
            if read_key in found_keys:
                correctly_uniquely_mapped += 1

    print("Uniquely mapped reads: {}, of which {} were true positives.".format(uniquely_mapped, correctly_uniquely_mapped))

    # number of true local alignments found relative to the total number of alignments reported
    sensitivity = found_positives / len(blast8_reads)
    
    # %of true local alignment which are correctly identified
    specificity = found_positives / len(true_positives)

    true_negatives = found_positives * (len(true_positives) - 1) #maybe theres a better way to define the true negatives
    # TN / (TN + FP)
    other_specificity = true_negatives / (true_negatives + false_positives)
    print("Sensitivity: {}, Specificity (by Definition): {}, Specificity (by TN): {}.".format(sensitivity, specificity, other_specificity))


if __name__ == "__main__":
    main()
