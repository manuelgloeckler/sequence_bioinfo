import pickle
import sqlite3
import ast
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


    def fetch_individuals(self, pop = "ALL", project = "1000GENOMES:phase_3", report_progress = False):
        """
        Send a request to the ensembl database, requesting all individuals for the given population and project. 
        Updates the database.

        Args:
            population (str): Name of the population the individuals should be returned for.
            project (str): Name of the project to use (default will be correct in almost all cases).
            report_progress (bool): If set to true fetching progress will be printed out.

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
            if report_progress: print("Finished fetching individuals for population: {}.".format(population.get("description", "unknown")))


    def fetch_populations(self, pop_filter = "LD", species = "homo_sapiens", report_progress = False):
        """
        Sends a request for the populations specified by given species and filter to the ensembl database.
        The result will be written to the local database.

        Args:
            pop_filter (str): TODO: not sure what this does.
            species (str): Species the populations should belong to.
            report_progress (bool): If set to true fetching progress will be printed out.

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
        if report_progress: print("Finished fetching populations for species: {}.".format(species))

    def fetch_reference_set(self, report_progress = False):
        """
        Sends a request for all reference sets in Ensembl. These sets are associated with reference sequences and are needed to access them.
        The result will be written to the local database.

        Args:
            report_progress (bool): If set to true fetching progress will be printed out.

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
        if report_progress: print("Finished fetching reference set.")

    def fetch_reference_sequences(self, set_id, report_progress = False):
        """
        Sends a request for all reference sequences in Ensembl. These sequences serve as basis for variants.
        The result will be written to the local database.

        Args:
            set_id (str): Id of the reference set for which to request all associated reference sequences.
            report_progress (bool): If set to true fetching progress will be printed out.

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
                    self.db_cursor.executemany("INSERT INTO reference_sequences(id, isPrimary, sourceURI, name, isDerived, length, ncbiTaxonId) VALUES (?, ?, ?, ?, ?, ?, ?) ON CONFLICT DO NOTHING", sequence_generator())
                except sqlite3.OperationalError as e:
                    logger.warning("Updating reference_sequences table failed with error: '{}'. Maybe try updating your SQLite Version.".format(e))
            else:
                for reference_sequence in sequence_generator():
                    try:
                        self.db_cursor.execute("INSERT INTO reference_sequences(id, isPrimary, sourceURI, name, isDerived, length, ncbiTaxonId) VALUES (?, ?, ?, ?, ?, ?, ?)", reference_sequence)
                    except sqlite3.IntegrityError as e:
                        logger.info("Tried to insert {} into reference_sequences table, but got error: {}.".format(reference_sequence, e))
                
        self.db_connector.commit()
        if report_progress: print("Finished fetching reference sequences for set: {}.".format(set_id))

    def fetch_variants(self, start, end, reference_name, variant_set_id = 1, report_progress = False):
        """
        Sends a request for all variants contained in the reference sequence (identified by reference_name) from start to end 
        and part of the variant set (identified by variantSetId, always using 1 seems sufficient) in Ensembl. 
        The result will be written to the local database.

        Args:
            start (int): Starting position, 0 indexed, included.
            end (int): End position, 0 indexed, excluded.
            reference_name (str): Name of the reference sequence for which to get the variant.
            variantSetId (int, optional): Id of the variant set from which to get the variant. 
            report_progress (bool): If set to true fetching progress will be printed out.

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
                variant_info = (variant['id'], variant['updated'], variant['created'], variant['start'], variant['end'], variant['referenceName'], str(variant['referenceBases']), variant['variantSetId'], variant['names'][0], str(variant["alternateBases"]))
                try:
                    self.db_cursor.execute("INSERT INTO variants(id, updated, created, start, end, referenceName, referenceBases, variantSetId, name, alternateBases) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)", variant_info)
                except sqlite3.IntegrityError as e:
                    logger.info("Tried to insert {} into variants table, but got error: {}.".format(variant_info, e))
                try:
                    def individual_variant_generator():
                        for call in variant["calls"]:
                            if len(call['genotype']) != 2: call['genotype'].append(None)
                            else: call['genotype'][1] = 0 != call['genotype'][1]
                            yield (variant['id'], call['callSetName'], 0 != call['genotype'][0], call['genotype'][1])
                    self.db_cursor.executemany("INSERT INTO individuals_variants(variant, individual, expression1, expression2) VALUES (?, ?, ?, ?)", individual_variant_generator())
                except sqlite3.IntegrityError as e:
                    logger.info("Tried to insert {} into individuals_variants table, but got error: {}. Make sure to fetch individuals first.".format(variant["calls"], e))
                    
                
            self.db_connector.commit()
            if report_progress: print("Finished fetching 10 variants next is {}.".format(data['pageToken']))
            if data['pageToken'] is None: break


    def fetch_phenotypes(self, start, end, reference_name, species = "homo_sapiens", feature_type = "Variation", report_progress = False):
        """
        Fetches all phenotype information in the specified region (via start & end) 
        on the given reference sequence (chromosome) for the given species (currently only homo_sapiens is sensible).
        Args:
            start (int): Start of the region for which to fetch phenotype information.
            end (int) : End of the region for which to fetch phenotype information.
            reference_name (str): Name of the reference sequence for which to fetch the phenotype information.
            species (str, optional): Name of the species for which to fetch the phenotype information.
                Defaults to 'homo_sapiens' which currently is the only sensible choice.
            feature_type (str, optional): Options are: Variation, StructuralVariation, Gene, QTL.
            report_progress (bool): If set to true fetching progress will be printed out.
        """
        region = "{}:{}-{}".format(reference_name, start, end)
        ext = "/phenotype/region/{}/{}".format(species, region)
        params = {"feature_type" : feature_type}
        phenotypes = self.client.perform_rest_action(ext, params = params)
        for phenotype in phenotypes:
            ensembl_id = phenotype["id"]
            for association in phenotype["phenotype_associations"]:
                location_split = association["location"].split(":")[1].split("-")
                start = location_split[0]
                end = location_split[1]
                attributes = association.get('attributes', {})
                p_value = attributes.get("p_value", None)
                associated_gene = attributes.get("associated_gene", None)
                risk_allele = attributes.get("risk_allele", None)
                external_reference = attributes.get("external_reference", None)
                phenotype_info = (ensembl_id, start, end, reference_name, association["description"], association["source"], p_value, associated_gene, risk_allele, external_reference)
                try:
                    self.db_cursor.execute("INSERT INTO phenotypes(ensembl_id, start, end, referenceName, description, source, p_value, associated_gene, risk_allele, external_reference) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)", phenotype_info)
                except sqlite3.IntegrityError as e:
                    logger.info("Tried to insert {} into variants table, but got error: {}.".format(phenotype_info, e))
               
        if report_progress: print("Finished fetching phenotypes from {} to {} for sequence {}.".format(start, end, reference_name))
        self.db_connector.commit()

    def fetch_region(self, start, end, reference_name, species = "homo_sapiens"):
        """
        Fetches the region of the reference sequence specified by reference name, start and end.

        Args:
            start (int): Start of the region which to fetch.
            end (int) : End of the region which to fetch.
            reference_name (str): Name of the reference sequence which to fetch.
            species (str, optional): Name of the species for which to fetch the sequence information.
                Defaults to 'homo_sapiens' which currently is the only sensible choice.
        """
        region = "{}:{}-{}".format(reference_name, start, end)
        ext = "/sequence/region/{}/{}".format(species, region)
        region = self.client.perform_rest_action(ext)
        #TODO: store to DB
        return region

    def fetch_all(self, start, end, reference_sequence, project="1000GENOMES:phase_3", **kwargs):
        """
        Fetches all information required to construct the model for the given section specified by
        chromosome, start and end.

        Args:
            start (int): Starting position of the section for which to fetch the information.
            end (int): Ending position of the section for which to fetch the information.
            reference_sequence (str): Name of the reference sequence (usuallly name of the chromosome)
                for which to fetch the information.
            kwargs: Any additional keyword arguments will be passed onto the specific fetch operations.
        """
        self.fetch_reference_set(**kwargs)
        self.fetch_reference_sequences("GRCh38", **kwargs)
        self.fetch_populations(pop_filter = None, **kwargs)

        for pop in self.get_populations(project=project):
            self.fetch_individuals(pop, project, **kwargs) 

        self.fetch_variants(start, end, reference_sequence, **kwargs)
        self.fetch_phenotypes(start, end, reference_sequence, **kwargs)


    def get_associated_phenotypes(self, start, end):
        return list(self.db_cursor.execute("""
        SELECT *
        FROM phenotypes
        WHERE start >= ? and end <= ?;
        """, (start, end)))

    def get_individuals(self):
        individuals = self.db_cursor.execute("""
        SELECT name FROM 'individuals'
        """)

        return sorted(list(map(lambda x: x[0], individuals)))

    def get_populations(self, project="1000GENOMES:phase_3"):
        populations = self.db_cursor.execute("""
        SELECT name FROM 'populations'
        """)
        pop_name = []
        for pop in populations:
            name = pop[0]
            if project in name:
                pop_name.append(name.split(":")[-1])

        return sorted(pop_name)

    def generate_individual_population_map(self):
        ind_pop = self.db_cursor.execute("""
        SELECT individual, population FROM 'individuals_populations'
        """)
        ind_pop_map = dict()
        for ind, pop in ind_pop:
            if ind not in ind_pop_map:
                ind_pop_map[ind] = [pop.split(":")[-1]]
            else:
                ind_pop_map[ind] += [pop.split(":")[-1]]
        return ind_pop_map



    def generate_inference_matrix(self, start = 0, end = None, reference_name = None, population = "ALL", project = "1000GENOMES:phase_3", sum_allels = False, shuffle = True):
        """
        Generates the inference matrix for the section specified by start and end,
        based on the database and the given reference sequence (chromosome), project and population.
        Make sure to have fetched the variants and individuals beforehand.
        Each row in the matrix represents an individual, each column a variant.
        Args:
            start (int, optional): Beginning (inclusive) of the section for which to consider variants.
                Defaults to 0.
            end (int, optional): End (exclusive) of the section for which to consider variants.
                Defaults to sys.maxsize.
            reference_name (str, optional): Name of the reference sequence the variants should belong to.
                If None is given it is assumed all variants in the DB belong to the same reference sequence.
            population (str, optional): Name of the population for which to
                generate the matrix. Defaults to ALL.
            project (str, optional): Name of the project for which to generate the matrix.
                Defaults to 1000GENOMES:phase_3. Our algorithms does not support all 
                naming conventions.
            sum_allels (bool, optional): When True, expression of a variant will be collected per individual,
                summing the expression per strand (if expressed on both -> 2, on one -> 1, on neither -> 0).
                Defaults to False.
            shuffle (bool, optional): When True, shuffles data fetched from the database before creating the matrix.
        Returns:
            inference matrix. Each row in the matrix represents an individual,
                each column a variant.
            individuals map. Maps each individual name to the corresponding row 
                in the inference matrix.
            variants map. Maps each variant to the corresponding column
                in the inference matrix.
        """
        if end is None: end = sys.maxsize
        population = "{}:{}".format(project, population)
        if reference_name is None:
            reference_condition = ""
            sql_args = (population, start, end)
        else:
            reference_condition = " AND referenceName = ?"
            sql_args = (population, start, end, reference_name)
        variants_individuals = self.db_cursor.execute("""
        SELECT variant, IV.individual, expression1, expression2 
        FROM individuals_variants IV 
        JOIN individuals_populations IP ON IV.individual = IP.individual
        JOIN variants V ON IV.variant = V.id
        WHERE population = ? AND start >= ? AND end < ?{};
        """.format(reference_condition), sql_args)

        # setup map to numerical values
        variants_individuals = list(variants_individuals) #make the iterable permanent
        number_variants = 0
        number_individuals = 0
        variants_map = {}
        individuals_map = {}
        if shuffle: np.random.shuffle(variants_individuals)
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

        
        return inference_matrix, individuals_map, variants_map

    def get_variation_distribution(self, start = 0, end = None, reference_name = None, population = "ALL", project = "1000GENOMES:phase_3"):
        """
        Returns a dictionary mapping number of variants within the specified region 
        to the number of strands with that number of variants.

        Args:
            start (int, optional): Begining (inclusive) of the section for which to consider variants.
                Defaults to 0.
            end (int, optional): End (exclusive) of the section for which to consider variants.
                Defaults to sys.maxsize.
            reference_name (str, optional): Name of the reference sequence the variants should belong to.
                If None is given it is assumed all variants in the DB belong to the same reference sequence.
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
        if reference_name is None:
            reference_condition = ""
            sql_args = (population, start, end)
        else:
            reference_condition = " AND referenceName = ?"
            sql_args = (population, start, end, reference_name)
        variants_individuals = self.db_cursor.execute("""
        SELECT IV.individual, SUM(expression1), SUM(expression2) 
        FROM individuals_variants IV 
        JOIN individuals_populations IP ON IV.individual = IP.individual
        JOIN variants V ON IV.variant = V.id
        WHERE population = ? AND start >= ? AND end < ?{}
        GROUP BY IV.individual;
        """.format(reference_condition), sql_args)

        distribution = {}
        for entry in variants_individuals:
            individual, sum_expr1, sum_expr2 = entry
            distribution[sum_expr1] = distribution.get(sum_expr1, 0) + 1
            distribution[sum_expr2] = distribution.get(sum_expr2, 0) + 1

        return distribution


    def get_variation_range(self, start = 0, end = None, reference_name = None, population = "ALL", project = "1000GENOMES:phase_3"):
        """
        Returns a dictionary mapping number of variants within the specified region 
        to the number of strands with that number of variants.

        Args:
            start (int, optional): Begining (inclusive) of the section for which to consider variants.
                Defaults to 0.
            end (int, optional): End (exclusive) of the section for which to consider variants.
                Defaults to sys.maxsize.
            reference_name (str, optional): Name of the reference sequence the variants should belong to.
                If None is given it is assumed all variants in the DB belong to the same reference sequence.
            population (str, optional): Name of the population for which to
                generate the distribution. Defaults to ALL.
            project (str, optional): Name of the project for which to generate the distribution.
                Defaults to 1000GENOMES:phase_3. Our algorithms does not support all 
                naming conventions.

        Returns:
            Set containing tuples of form (variant_id, variant_start, variant_end).
        """
        if end is None: end = sys.maxsize
        population = "{}:{}".format(project, population)
        if reference_name is None:
            reference_condition = ""
            sql_args = (population, start, end)
        else:
            reference_condition = " AND referenceName = ?"
            sql_args = (population, start, end, reference_name)
            
        #TODO: Could probably optimize the sql query
        variant_ranges = self.db_cursor.execute("""
        SELECT V.id, V.start, V.end
        FROM individuals_variants IV 
        JOIN individuals_populations IP ON IV.individual = IP.individual
        JOIN variants V ON IV.variant = V.id
        WHERE population = ? AND start >= ? AND end < ?{};
        """.format(reference_condition), sql_args)
        variant_ranges = set(variant_ranges)
        return variant_ranges

    def get_variation_alternate(self, start = 0, end = None, reference_name = None, population = "ALL", project = "1000GENOMES:phase_3"):
        """
        Returns a dictionary mapping variation ids to their alternate base, their start, and their end.

        Args:
            start (int, optional): Begining (inclusive) of the section for which to consider variants.
                Defaults to 0.
            end (int, optional): End (exclusive) of the section for which to consider variants.
                Defaults to sys.maxsize.
            reference_name (str, optional): Name of the reference sequence the variants should belong to.
                If None is given it is assumed all variants in the DB belong to the same reference sequence.
            population (str, optional): Name of the population for which to
                generate the distribution. Defaults to ALL.
            project (str, optional): Name of the project for which to generate the distribution.
                Defaults to 1000GENOMES:phase_3. Our algorithms does not support all 
                naming conventions.

        Returns:
            Dictionary mapping variant id to (alternate_base, start, end) of that variant.
        """
        if end is None: end = sys.maxsize
        population = "{}:{}".format(project, population)
        if reference_name is None:
            reference_condition = ""
            sql_args = (population, start, end)
        else:
            reference_condition = " AND referenceName = ?"
            sql_args = (population, start, end, reference_name)
            
        #TODO: Could probably optimize the sql query
        variant_bases = self.db_cursor.execute("""
        SELECT V.id, V.alternateBases, V.start, V.end
        FROM individuals_variants IV 
        JOIN individuals_populations IP ON IV.individual = IP.individual
        JOIN variants V ON IV.variant = V.id
        WHERE population = ? AND start >= ? AND end < ?{};
        """.format(reference_condition), sql_args)
        variant_bases = set(variant_bases)
        alternates_map = {}
        for variant_base in variant_bases:
            alternates_map[variant_base[0]] = (ast.literal_eval(variant_base[1]), variant_base[2], variant_base[3])
        return alternates_map




            

if __name__ == "__main__":
    Database_Name = "Gentype_DB.db"
    client = EnsemblClient()
    data_manager = DataManager(client, Database_Name)
    data_manager.get_variation_range(start = 29941260, end = 29945884, population = "ALL")