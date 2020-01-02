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

    try:
        cursor.execute("""
        CREATE TABLE IF NOT EXISTS individuals (
            name VARCHAR(64) PRIMARY KEY,
            gender VARCHAR(6)
        )
        """)
    except sqlite3.OperationalError as e:
        logger.warning("Creating individuals table for {} failed with: '{}'. Trying to continue anyways.".format(database, e))

    try:
        cursor.execute("""
        CREATE TABLE IF NOT EXISTS populations (
            name VARCHAR(64) PRIMARY KEY,
            size INT,
            description VARCHAR(64)
        )
        """)
    except sqlite3.OperationalError as e:
        logger.warning("Creating populations table for {} failed with: '{}'. Trying to continue anyways.".format(database, e))

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
        logger.warning("Creating populations table for {} failed with: '{}'. Trying to continue anyways.".format(database, e))

if __name__ == "__main__":
    connector = sqlite3.connect("potato2.db")
    cursor = connector.cursor()
    ensure_samples_tables(cursor)
    connector.commit()
    connector.close()