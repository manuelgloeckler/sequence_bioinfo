import pickle
import sqlite3
import sys
from EnsemblClient import EnsemblClient
import DB_setup
import logging
FORMAT = '%(levelname)s: %(filename)s %(lineno)d, %(funcName)s: %(message)s'
logging.basicConfig(level=logging.WARNING, format=FORMAT)
logger = logging.getLogger()


class Model:
    def __init__(self, client, local_database):
        """
        Model to manage all relevant information and access to it.

        Args:
            client (EnsemblClient): Object with which to make request to the server from which to fetch the data
            local_database (str): Name of the local database to be used. Currently we use sqlite.
        """
        self.client = client
        
        #connect to local database
        self.db_connector = sqlite3.connect(local_database)
        self.db_cursor = self.db_connector.cursor()
        DB_setup.ensure_tables(self.db_cursor)
        self.db_connector.commit()
        self.individuals = {}
        pass


    def fetch_individuals(self, pop = "ALL", project = "1000GENOMES:phase_3"):
        """
        Send a request to the ensembl database, requesting all individuals for the given population and project. 
        Updates the database and returns the stored information.

        Args:
            population (str): Name of the population the individuals should be returned for.
            project (str): Name of the project to use (default will be correct in almost all cases).

        Returns:
            List of all individuals (by population). Will also be stored in the database.
        """

        ext = "/info/variation/populations/human/"+ project + ":" + pop
        populations = self.client.perform_rest_action(ext)
        # I think the way we currently query ensembl it will be always len(populations) = 1, still a for loop seems more semantically correct
        for population in populations:
            # make sure that the population is present in the local database
            try:
                population_info = (population['name'], population['size'], population['description'])
                self.db_cursor.execute("INSERT INTO populations(name, size, description) VALUES (?, ?, ?)", population_info)
            except sqlite3.IntegrityError as e:
                logger.info("Tried to insert {} into populations table, but got error: {}.".format(population_info, e))

            #insert individuals into local db
            individuals = population["individuals"]
            def individual_generator():
                for individual in individuals:
                    yield (individual['name'], individual['gender'])
            if sqlite3.sqlite_version > "3.24":
                try:
                    self.db_cursor.executemany("INSERT INTO individuals(name, gender) VALUES (?, ?) ON CONFLICT DO NOTHING", individual_generator())
                except sqlite3.OperationalError as e:
                    logger.warning("Updating individuals table failed with error: '{}'. Maybe try updating your SQLite Version.".format(e))

            else:
                for individual_info in individual_generator():
                    try:
                        self.db_cursor.execute("INSERT INTO individuals(name, gender) VALUES (?, ?)", individual_info)
                    except sqlite3.IntegrityError as e:
                        logger.info("Tried to insert {} into individuals table, but got error: {}.".format(individual_info, e))

            # connect individuals with their population
            def individual_population_generator():
                for individual in individuals:
                    yield (individual['name'], population['name'])
            if sqlite3.sqlite_version > "3.24":
                try:
                    self.db_cursor.executemany("INSERT INTO individuals_populations(individual, population) VALUES (?, ?) ON CONFLICT DO NOTHING", individual_population_generator())
                except sqlite3.OperationalError as e:
                    logger.warning("Updating individuals_populations table failed with error: '{}'. Maybe try updating your SQLite Version.".format(e))

            else:
                for individual_population_info in individual_population_generator():
                    try:
                        self.db_cursor.execute("INSERT INTO individuals_populations(individual, population) VALUES (?, ?)", individual_population_info)
                    except sqlite3.IntegrityError as e:
                        logger.info("Tried to insert {} into individuals_populations table, but got error: {}.".format(individual_population_info, e))

            self.db_connector.commit()
        return populations


    def fetch_populations(self, pop_filter = "LD", species = "homo_sapiens"):
        """
        Sends a request for the populations specified by given species and filter to the ensembl database.
        The result will be written to the local database and returned.

        Args:
            pop_filter (str): TODO: not sure what this does.
            species (str): Species the populations should belong to.

        Returns:
            List of populations belonging to given species adhering to given filter.
        """
        ext = "/info/variation/populations/" + species
        params = {'filter' : pop_filter} if pop_filter else None
        populations = self.client.perform_rest_action(ext, params=params)

        #insert populations into local db
        def population_generator():
            for population in populations:
                yield (population['name'], population['size'], population['description'])
        if sqlite3.sqlite_version > "3.24":
            try:
                self.db_cursor.executemany("INSERT INTO populations(name, size, description) VALUES (?, ?, ?) ON CONFLICT DO NOTHING", population_generator())
            except sqlite3.OperationalError as e:
                logger.warning("Updating populations table failed with error: '{}'. Maybe try updating your SQLite Version.".format(e))
        else:
            for population_info in population_generator():
                try:
                    self.db_cursor.execute("INSERT INTO populations(name, size, description) VALUES (?, ?, ?)", population_info)
                except sqlite3.IntegrityError as e:
                    logger.info("Tried to insert {} into populations table, but got error: {}.".format(population_info, e))
            
        self.db_connector.commit()
        
        return populations


if __name__ == '__main__':
    client = EnsemblClient()
    model = Model(client, "test_db1.db")
    #print(model.fetch_populations(pop_filter=None))
    print(model.fetch_individuals())


