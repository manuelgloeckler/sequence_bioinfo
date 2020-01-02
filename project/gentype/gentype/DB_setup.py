import sqlite3
import logging
logger = logging.getLogger()

def ensure_tables(cursor):
    """
    Esures that all neccesary tables exist in the sqlite database.
    If they don't they are created.

    Args:
        cursor (sqlite3.Cursor): Cursor to the database which is to be checked.
    """
    ensure_individuals_table(cursor)
    ensure_populations_table(cursor)
    ensure_individuals_populations_table(cursor)
    ensure_referenceSets_table(cursor)
    ensure_reference_sequences_table(cursor)
    ensure_variants_table(cursor)
    ensure_individuals_variants_table(cursor)

    
def ensure_individuals_table(cursor):
    """
    Esures that the individuals table exist in the sqlite database.
    It has the entries (name, gender).
    If it doesn't exist it will be created.

    Args:
        cursor (sqlite3.Cursor): Cursor to the database which is to be checked.
    """
    try:
        cursor.execute("""
        CREATE TABLE IF NOT EXISTS individuals (
            name VARCHAR(64) PRIMARY KEY,
            gender VARCHAR(6)
        )
        """)
    except sqlite3.OperationalError as e:
        logger.warning("Creating individuals table failed with: '{}'. Trying to continue anyways.".format(e))

def ensure_populations_table(cursor):
    """
    Esures that the populations table exist in the sqlite database.
    It has the entries (name, size, description).
    If it doesn't exist it will be created.

    Args:
        cursor (sqlite3.Cursor): Cursor to the database which is to be checked.
    """
    try:
        cursor.execute("""
        CREATE TABLE IF NOT EXISTS populations (
            name VARCHAR(64) PRIMARY KEY,
            size INT,
            description VARCHAR(64)
        )
        """)
    except sqlite3.OperationalError as e:
        logger.warning("Creating populations table failed with: '{}'. Trying to continue anyways.".format(e))

def ensure_individuals_populations_table(cursor):
    """
    Esures that the individuals_populations table exist in the sqlite database.
    It has the entries (individual -> individuals, population -> populations).
    If it doesn't exist it will be created.

    Args:
        cursor (sqlite3.Cursor): Cursor to the database which is to be checked.
    """
    try:
        cursor.execute("""
        CREATE TABLE IF NOT EXISTS individuals_populations (
            individual VARCHAR(64),
            population VARCHAR(64),
            UNIQUE(individual, population),
            FOREIGN KEY (individual) REFERENCES individuals(name),
            FOREIGN KEY (population) REFERENCES populations(name)
        )
        """)
    except sqlite3.OperationalError as e:
        logger.warning("Creating individuals_populations table failed with: '{}'. Trying to continue anyways.".format(e))

def ensure_referenceSets_table(cursor):
    """
    Esures that the referenceSets table exist in the sqlite database.
    It has the entries (id, description, sourceURI, name, ncbiTaxonId).
    If it doesn't exist it will be created.

    Args:
        cursor (sqlite3.Cursor): Cursor to the database which is to be checked.
    """
    try:
        cursor.execute("""
        CREATE TABLE IF NOT EXISTS referenceSets (
            id VARCHAR(16) PRIMARY KEY,
            description VARCHAR(64),
            sourceURI VARCHAR(256),
            name VARCHAR(64),
            ncbiTaxonId INT
        )
        """)
    except sqlite3.OperationalError as e:
        logger.warning("Creating referenceSets table failed with: '{}'. Trying to continue anyways.".format(e))

def ensure_reference_sequences_table(cursor):
    """
    Esures that the reference_sequences table exist in the sqlite database.
    It has the entries (id, isPrimary, sourceURI, name, isDerived, length, ncbiTaxonId).
    If it doesn't exist it will be created.

    Args:
        cursor (sqlite3.Cursor): Cursor to the database which is to be checked.
    """
    try:
        cursor.execute("""
        CREATE TABLE IF NOT EXISTS reference_sequences (
            id VARCHAR(64) PRIMARY KEY,
            isPrimary BIT,
            sourceURI VARCHAR(256),
            name VARCHAR(64),
            isDerived BIT,
            length INT,
            ncbiTaxonId INT
        )
        """)
    except sqlite3.OperationalError as e:
        logger.warning("Creating reference_sequences table failed with: '{}'. Trying to continue anyways.".format(e))

def ensure_variants_table(cursor):
    """
    Esures that the variants table exist in the sqlite database.
    It has the entries (id, updated, created, start, end, referenceName, referenceBases, variantSetId, name, alternateBases).
    If it doesn't exist it will be created.

    Args:
        cursor (sqlite3.Cursor): Cursor to the database which is to be checked.
    """
    try:
        cursor.execute("""
        CREATE TABLE IF NOT EXISTS variants (
            id VARCHAR(64) PRIMARY KEY,
            updated BIGINT,
            created BIGINT,
            start INT,
            end INT,
            referenceName VARCHAR(64),
            referenceBases VARCHAR(64),
            variantSetId INT,
            name VARCHAR(64),
            alternateBases VARCHAR(64)
        )
        """)
    except sqlite3.OperationalError as e:
        logger.warning("Creating variants table failed with: '{}'. Trying to continue anyways.".format(e))

def ensure_individuals_variants_table(cursor):
    """
    Esures that the individuals_variants table exist in the sqlite database.
    It has the entries (variant, individual, expression1, expression2).
    
    If it doesn't exist it will be created.

    Args:
        cursor (sqlite3.Cursor): Cursor to the database which is to be checked.
    """
    try:
        cursor.execute("""
        CREATE TABLE IF NOT EXISTS individuals_variants (
            variant VARCHAR(64),
            individual VARCHAR(64),
            expression1 BIT,
            expression2 BIT,
            UNIQUE(individual, variant),
            FOREIGN KEY (individual) REFERENCES individuals(name),
            FOREIGN KEY (variant) REFERENCES variants(id)
        )
        """)
    except sqlite3.OperationalError as e:
        logger.warning("Creating individuals_variants table failed with: '{}'. Trying to continue anyways.".format(e))



if __name__ == "__main__":
    connector = sqlite3.connect("potato2.db")
    cursor = connector.cursor()
    ensure_samples_tables(cursor)
    connector.commit()
    connector.close()