import logging

logger = logging.getLogger()

__author__ = "Yanpeng Li, Josua Stadelmaier, Finn Mier"

def read_sequence(sequence_path):
    if not (sequence_path.endswith(".fasta") or sequence_path.endswith(".fna")):
        logger.error("Filenames must end in '.fasta' or '.fna' but was {}".format(sequence_path))
        raise ValueError
    else:
        sequences = []
        headers = []
        # list containing each line read of the current sequence. Will be joined at the end.
        current_sequence = []
        # header describing the current sequence
        current_header = None
        with open(sequence_path) as seq_file:
            for line_nr, line in enumerate(seq_file):
                if line[0] == ">": # header
                    if not current_header is None:
                        sequences.append("".join(current_sequence))
                        headers.append(current_header)
                        current_sequence = []
                    current_header = line.rstrip("\n")
                elif line[0] == ";": # comment
                    pass
                else: # part of sequence
                    if current_header is None:
                        logger.warning("Found sequence in line {} in file {} before header. This sequence will be ignored.".format(line_nr, sequence_path))
                    else:
                        current_sequence.append(line.rstrip("\n"))
            # add final sequence
            if not current_header is None:
                sequences.append("".join(current_sequence))
                headers.append(current_header)
        return headers, sequences
            

def write_sequences(file_path, headers, sequences):
    if len(headers) != len(sequences):
        logger.error("Tried to write sequences but length of given headers({}) and sequences({}) was not equal.".format(len(headers), len(sequences)))
        raise ValueError
    with open(file_path, 'w') as target:
        for i, sequence in enumerate(sequences):
            header = headers[i]
            seperated_header = [header[i:i + 80] for i in range(0, len(header), 80)] #fasta requires line length < 80
            for line in seperated_header:
                print(line, file=target)
            seperated_sequence = [sequence[i:i + 80] for i in range(0, len(sequence), 80)] #fasta requires line length < 80
            for line in seperated_sequence:
                print(line, file=target)

def write_comparison(file_path, headers, sequences):
    if len(headers) != len(sequences):
        logger.error("Tried to write sequences but length of given headers({}) and sequences({}) was not equal.".format(len(headers), len(sequences)))
        raise ValueError
    with open(file_path, 'w') as target:
        for i, sequence in enumerate(sequences):
            header = headers[i]
            print(header, file=target)
            print(sequence, file=target)

