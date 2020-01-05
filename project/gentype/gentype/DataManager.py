import pickle
import sqlite3
import sys
try:
    from .EnsemblClient import EnsemblClient
except:
    from EnsemblClient import EnsemblClient
try:
    from . import DB_setup
except:
    import DB_setup
import logging
import numpy as np
FORMAT = '%(levelname)s: %(filename)s %(lineno)d, %(funcName)s: %(message)s'
logging.basicConfig(level=logging.WARNING, format=FORMAT)
logger = logging.getLogger()


class DataManager:
    def __init__(self, client, local_database):
        """
        Class to manage all relevant information and access to it.

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
                    try:
                        yield (individual['name'].split(":")[2], individual['gender'])
                    except IndexError:
                        logger.warning("Name '{}' did not adhere to naming convention. Was skipped.".format(individual['name']))

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
                    try:
                        yield (individual['name'].split(":")[2], population['name'])
                    except IndexError:
                        logger.warning("Name '{}' did not adhere to naming convention. Was skipped.".format(individual['name']))
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
                self.db_cursor.execute("INSERT INTO referenceSets(id, description, sourceURI, name, ncbiTaxonId) VALUES (?, ?, ?, ?, ?)", referenceSet_info)
            except sqlite3.IntegrityError as e:
                logger.info("Tried to insert {} into referenceSets table, but got error: {}.".format(referenceSet_info, e))
                
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

    def fetch_variants(self, start, end, reference_name, variant_set_id = 1):
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
                            yield (variant['id'], call['callSetName'], call['genotype'][0], call['genotype'][1])
                    self.db_cursor.executemany("INSERT INTO individuals_variants(variant, individual, expression1, expression2) VALUES (?, ?, ?, ?)", individual_variant_generator())
                except sqlite3.IntegrityError as e:
                    logger.info("Tried to insert {} into individuals_variants table, but got error: {}. Make sure to fetch individuals first.".format(variant["calls"], e))
                    
                
            self.db_connector.commit()
            if data['pageToken'] is None: break

    def generate_inference_matrix(self, start = 0, end = None, population = "ALL", project = "1000GENOMES:phase_3", sum_allels = False):
        """
        Generates the inference matrix for the section specified by start and end,
        based on the database and the given project and population.
        Make sure to have fetched the variants and individuals beforehand.
        Each row in the matrix represents an individual, each column a variant.

        Args:
            start (int, optional): Begining (inclusive) of the section for which to consider variants.
                Defaults to 0.
            end (int, optional): End (exclusive) of the section for which to consider variants.
                Defaults to sys.maxsize.
            population (str, optional): Name of the population for which to
                generate the matrix. Defaults to ALL.
            project (str, optional): Name of the project for which to generate the matrix.
                Defaults to 1000GENOMES:phase_3. Our algorithms does not support all 
                naming conventions.
            sum_allels (bool, optional): When True, expression of a variant will be collected per individual,
                summing the expression per strand (if expressed on both -> 2, on one -> 1, on neither -> 0).
                Defaults to False.

        Returns:
            inference matrix. Each row in the matrix represents an individual,
                each column a variant.
        """
        if end is None: end = sys.maxsize
        population = "{}:{}".format(project, population)
        variants_individuals = self.db_cursor.execute("""
        SELECT variant, IV.individual, expression1, expression2 
        FROM individuals_variants IV 
        JOIN individuals_populations IP ON IV.individual = IP.individual
        JOIN variants V ON IV.variant = V.id
        WHERE population = ? AND start >= ? AND end < ?;
        """, (population, start, end))

        # setup map to numerical values
        variants_individuals = list(variants_individuals) #make the iterable permanent
        number_variants = 0
        number_individuals = 0
        variants_map = {}
        individuals_map = {}
        for entry in variants_individuals:
            variant, individual, expr1, expr2 = entry
            if not variant in variants_map:
                variants_map[variant] = number_variants
                number_variants += 1
            if not individual in individuals_map:
                individuals_map[individual] = number_individuals
                number_individuals += 1
        
        # create matrix
        if sum_allels:
            inference_matrix = np.ndarray((number_individuals, number_variants))
            for entry in variants_individuals:
                variant, individual, expr1, expr2 = entry
                inference_matrix[individuals_map[individual]][variants_map[variant]] = expr1 + expr2
        else:
            inference_matrix = np.ndarray((2 * number_individuals, number_variants)) 
            for entry in variants_individuals:
                variant, individual, expr1, expr2 = entry
                inference_matrix[2 * individuals_map[individual]][variants_map[variant]] = expr1
                inference_matrix[2 * individuals_map[individual] + 1][variants_map[variant]] = expr2

        
        return inference_matrix

    def get_variation_distribution(self, start = 0, end = None, population = "ALL", project = "1000GENOMES:phase_3"):
        """
        Returns a dictionary mapping number of variants within the specified region 
        to the number of strands with that number of variants.

        Args:
            start (int, optional): Begining (inclusive) of the section for which to consider variants.
                Defaults to 0.
            end (int, optional): End (exclusive) of the section for which to consider variants.
                Defaults to sys.maxsize.
            population (str, optional): Name of the population for which to
                generate the distribution. Defaults to ALL.
            project (str, optional): Name of the project for which to generate the distribution.
                Defaults to 1000GENOMES:phase_3. Our algorithms does not support all 
                naming conventions.

        Returns:
            Distribution via dictionary. Maps number of variants within the specified region 
                to the number of strands with that number of variants. I.e. {n : #strands with n variations}.
        """
        if end is None: end = sys.maxsize
        population = "{}:{}".format(project, population)
        variants_individuals = self.db_cursor.execute("""
        SELECT IV.individual, SUM(expression1), SUM(expression2) 
        FROM individuals_variants IV 
        JOIN individuals_populations IP ON IV.individual = IP.individual
        JOIN variants V ON IV.variant = V.id
        WHERE population = ? AND start >= ? AND end < ?
        GROUP BY IV.individual;
        """, (population, start, end))

        distribution = {}
        for entry in variants_individuals:
            individual, sum_expr1, sum_expr2 = entry
            distribution[sum_expr1] = distribution.get(sum_expr1, 0) + 1
            distribution[sum_expr2] = distribution.get(sum_expr2, 0) + 1

        return distribution

if __name__ == '__main__':
    client = EnsemblClient()
    model = DataManager(client, "test_db1.db")
    #model.fetch_populations(pop_filter=None)
    #model.fetch_individuals()
    #model.fetch_individuals("CHB", "1000GENOMES:phase_3")
    #matrix = model.generate_inference_matrix("CHB")
    #np.set_printoptions(edgeitems = max(matrix.shape))
    #print(matrix)
    #model.fetch_reference_set()
    #model.fetch_reference_sequences("GRCh37.p13")
    #model.fetch_variants(17671934, 17675934, "22")
    model.get_variation_distribution(17671934, 17675934, "CHB")

