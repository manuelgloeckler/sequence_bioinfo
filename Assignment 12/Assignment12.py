from JC_Tree import *
from FastA import *
from SeqGenerator import *
import argparse


def parse_args():
    parser = argparse.ArgumentParser(description="Artificial phylogenetic tree construction based on Jukes Cantor.")
    parser.add_argument('-n', '--newick', type=str, required=True, help="Name of the file containing newick representation of the tree.")
    parser.add_argument('-m', '--mutation_rate', type=float, required=True, help="Rate of mutation for the sequences.")
    parser.add_argument('-l', '--sequence_length', type=int, required=True, help="Length to the generated sequences.")
    parser.add_argument('-o', '--output', default="assignment12.fasta", type=str, help="Name of the Output file. Defaults to assignment12.fasta.")
    parser.add_argument('-?', action="help", help="Shows this help message and exits.")

    args = parser.parse_args()
    return args

def run(newick_string, sequence_length, mutation_rate, output):
    # parse newick string
    parser = Newick_Parser()
    tree_root = parser.parse_newick(newick_string)

    # generate root sequence
    seq_generator = SeqGenerator(mutation_rate)
    tree_root.sequence = seq_generator.generate_sequence(sequence_length)

    # generated sequences for the rest of the tree    
    worklist = [tree_root]
    leafs = []
    while worklist:
        node = worklist.pop()
        if node.sequence is None: node.sequence = seq_generator.mutate_sequence(node.parent.sequence, node.distance)
        worklist.extend(node)
        if len(node.children) == 0: leafs.append(node)
    
    # write generated sequences to fasta
    headers = ["> {}".format(leaf.name) for leaf in leafs]
    sequences = [seq_generator.vector_to_seq(leaf.sequence) for leaf in leafs]
    write_sequences(output, headers, sequences)
    print("Output written to {}.".format(output))


if __name__ == "__main__":
    args = parse_args()
    with open(args.newick) as newick_file:
        newick_string = newick_file.read()
    run(newick_string, args.sequence_length, args.mutation_rate, args.output)