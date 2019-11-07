import logging

logger = logging.getLogger()

__author__ = "Patrick Schirm, Manuel Glöckler, Finn Mier"

def read_blast8(file_path):
    if not (file_path.endswith(".blast8")):
        logger.error("Filenames must end in '.blast8' but was {}".format(file_path))
        raise ValueError
    else:
        blast_reads = []
        with open(file_path) as blast8_file:
            for line_nr, line in enumerate(blast8_file):
                line = line.rstrip("\n")
                blast_read = line.split("\t")
                logger.debug(blast_read)
                
        return blast_reads
            