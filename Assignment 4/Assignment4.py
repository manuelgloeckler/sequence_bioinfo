#!/usr/bin/python

"""Analysis tool for blat output files"""

__author__ = "Patrick Schirm, Manuel Glöckler, Finn Mier"

import argparse
from Blast8 import read_blast8

import logging 
logger = logging.getLogger()


def parse_args():
    parser = argparse.ArgumentParser(description="Analysis tool for blat output files")
    parser.add_argument('-i', '--input', nargs = 2, type=str, required=True, help="pass your blat output and true positive file for analysis.")
    parser.add_argument('-o', '--output', nargs = 1, default=["analysis.txt"], type=str, help="Name of the Output file. Defaults to analysis.txt. If the file already exists an error will be thrown.")
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
    print("Author: {}".format(__author__))

    blast8_reads = read_blast8(args.input[0])

    print("Writing results to {}".format(args.output[0]))
       


if __name__ == "__main__":
    main()
