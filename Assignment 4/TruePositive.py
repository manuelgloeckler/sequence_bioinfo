import logging

logger = logging.getLogger()

__author__ = "Patrick Schirm, Manuel Glöckler, Finn Mier"

#subject, query, true_position_start, direction


def read_TP(file_path):
    """ Takes the path to a .TP file and writes its contents to a dictionary mapping
    a query and subject to a list containing all true locations.
    
    Parameters:
        file_path (str): Path to the .TP file.

    Returns:
        {(query, subject) : (true_position_start, direction)}
    """
    true_locations = {}
    with open(file_path) as tp_file:
        for line_nr, line in enumerate(tp_file):
            if line_nr == 0: continue #skip first line
            line = line.rstrip("\n")
            line = line.lstrip(">")
            tp_info = line.split("\t")
            subject, query = tp_info[0], tp_info[1]
            true_locations[(query, subject)] = (int(tp_info[2]), tp_info[3])
    return true_locations
            