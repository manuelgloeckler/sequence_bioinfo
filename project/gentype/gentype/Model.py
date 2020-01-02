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
        Updates the database.

        Args:
            population (str): Name of the population the individuals should be returned for.
            project (str): Name of the project to use (default will be correct in almost all cases).

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
                    yield (individual['name'].split(":")[2], individual['gender'])
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
                    yield (individual['name'].split(":")[2], population['name'])
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


    def fetch_populations(self, pop_filter = "LD", species = "homo_sapiens"):
        """
        Sends a request for the populations specified by given species and filter to the ensembl database.
        The result will be written to the local database.

        Args:
            pop_filter (str): TODO: not sure what this does.
            species (str): Species the populations should belong to.

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

    def fetch_reference_set(self):
        """
        Sends a request for all reference sets in Ensembl. These sets are associated with reference sequences and are needed to access them.
        The result will be written to the local database.

        """
        ext = "/ga4gh/referencesets/search"
        headers={ "Content-Type" : "application/json", "Accept" : "application/json"}
        referenceSets = self.client.perform_rest_action(ext, headers, data = {})
        # referenceSets generally will only contain one entry, but this seems more robust
        for referenceSet in referenceSets['referenceSets']:
            try:
                referenceSet_info = (referenceSet['id'], referenceSet['description'], referenceSet['sourceURI'], referenceSet['name'], referenceSet['ncbiTaxonId'])
                self.db_cursor.execute("INSERT INTO referenceSet(id, description, sourceURI, name, ncbiTaxonId) VALUES (?, ?, ?, ?, ?)", referenceSet_info)
            except sqlite3.IntegrityError as e:
                logger.info("Tried to insert {} into populations table, but got error: {}.".format(population_info, e))
                
        self.db_connector.commit()

    def fetch_reference_sequences(self, set_id):
        """
        Sends a request for all reference sequences in Ensembl. These sequences serve as basis for variants.
        The result will be written to the local database.

        Args:
            set_id (str): Id of the reference set for which to request all associated reference sequences.

        """
        ext = "/ga4gh/references/search"
        headers={ "Content-Type" : "application/json", "Accept" : "application/json"}
        next_page = 0
        data = {"referenceSetId" : set_id, "pageSize" : 10}
        while not next_page is None:
            data['pageToken'] = next_page
            answer = self.client.perform_rest_action(ext, headers, data = data)
            next_page = answer['nextPageToken']
            def sequence_generator():
                for sequence in answer['references']:
                    yield (sequence['id'], sequence['isPrimary'] == 'true', sequence['sourceURI'], sequence['name'], sequence['isDerived'] == 'true', sequence['length'], sequence['ncbiTaxonId'])
            if sqlite3.sqlite_version > "3.24":
                try:
                    self.db_cursor.executemany("INSERT INTO reference_sequences(id, isPrimary, sourceURI, name, isDerived, length, ncbiTaxonId) VALUES (?, ?, ?, ?, ?, ?, ?) ON CONFLICT DO NOTHING", population_generator())
                except sqlite3.OperationalError as e:
                    logger.warning("Updating reference_sequences table failed with error: '{}'. Maybe try updating your SQLite Version.".format(e))
            else:
                for reference_sequence in sequence_generator():
                    try:
                        self.db_cursor.execute("INSERT INTO reference_sequences(id, isPrimary, sourceURI, name, isDerived, length, ncbiTaxonId) VALUES (?, ?, ?, ?, ?, ?, ?)", reference_sequence)
                    except sqlite3.IntegrityError as e:
                        logger.info("Tried to insert {} into reference_sequences table, but got error: {}.".format(reference_sequence, e))
                
        self.db_connector.commit()

    def fetch_variants(self, start, end, reference_name, variant_set_id = 3):
        """
        Sends a request for all variants contained in the reference sequence (identified by reference_name) from start to end 
        and part of the variant set (identified by variantSetId, always using 1 seems sufficient) in Ensembl. 
        The result will be written to the local database.

        Args:
            start (int): Starting position, 0 indexed, included.
            end (int): End position, 0 indexed, excluded.
            reference_name (str): Name of the reference sequence for which to get the variant.
            variantSetId (int, optional): Id of the variant set from which to get the variant. 

        """
        ext = "/ga4gh/variants/search"
        headers={ "Content-Type" : "application/json", "Accept" : "application/json"}
        next_page = 0
        data = {"end" : end, "referenceName" : reference_name, "start" : start, "variantSetId" : variant_set_id, "pageSize" : 10}
        while True:
            logger.info(data.get('pageToken', ""))
            answer = self.client.perform_rest_action(ext, headers, data = data)
            data['pageToken'] = answer['nextPageToken']
            for variant in answer['variants']:
                variant_info = (variant['id'], variant['updated'], variant['created'], variant['start'], variant['end'], variant['referenceName'], str(variant['referenceBases']), variant['variantSetId'], variant['names'][0], "".join(variant["alternateBases"]))
                try:
                    self.db_cursor.execute("INSERT INTO variants(id, updated, created, start, end, referenceName, referenceBases, variantSetId, name, alternateBases) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)", variant_info)
                except sqlite3.IntegrityError as e:
                    logger.info("Tried to insert {} into variants table, but got error: {}.".format(variant_info, e))
                try:
                    def individual_variant_generator():
                        for call in variant["calls"]:
                            yield (variant['id'], call['callSetName'], sum(call['genotype']))
                    self.db_cursor.executemany("INSERT INTO individuals_variants(variant, individual, expression) VALUES (?, ?, ?)", individual_variant_generator())
                except sqlite3.IntegrityError as e:
                    logger.info("Tried to insert {} into individuals_variants table, but got error: {}. Make sure to fetch individuals first.".format(variant["calls"], e))
                    
                
            self.db_connector.commit()
            if data['pageToken'] is None: break

if __name__ == '__main__':
    client = EnsemblClient()
    model = Model(client, "test_db1.db")
    model.fetch_populations(pop_filter=None)
    model.fetch_individuals()
    model.fetch_variants(17190024, 17190824, "22")


