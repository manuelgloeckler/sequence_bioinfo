import warnings
import pickle
from .fetch_ensembl import *


class Ensembl():
    def __init__(self):
        # Save fetched information for reuse
        self.samples = list()
        self.populations = list()
        self.population_info = dict()
        self.samples_info = dict()
        self.map_samples_to_pop = dict()
        self.map_pop_to_samples = dict()
        self.fetched_variants = dict()
        self.sequence_by_id = dict()

        # New information

    """ Save state of the collector """
    def save(self, path):
        with  open(path, 'w') as file1: 
            pickle.dump(self, file1)

    """ Static Function to load a state """
    @staticmethod
    def load(path):
        with open(path, 'rb') as file1:
            new_datacolecter = pickle.load(file1)
            return new_datacolecter


    """ This method will return all populations of ensembl with the filter set to LD -> 1000 Genome phase 3"""
    def get_populations(self):
        if len(self.populations) == 0:
            populations = fetch_populations()
            for entry in populations:
                key = entry["name"].split(":")[-1]
                self.populations.append(key)
                self.population_info[key] = entry

        return self.populations

    """ This method will return all samples from 1000GENOMES:phase_3 for all populations """
    def get_samples(self):
        if len(self.populations) == 0:
            
            samples = fetch_samples()[0]["individuals"]
            for entry in samples:
                key = entry["name"].split(":")[-1]
                self.samples.append(key)
                self.samples_info[key] = entry 
        return self.samples

    def _create_pop_sample_mappings(self):
        if len(self.samples) == 0:
            self.get_samples()

        if len(self.populations) == 0:
            self.get_populations()

        if len(self.map_pop_to_samples) == 0:
            for pop in self.populations:
                samples_pop = fetch_samples(pop = pop)[0]["individuals"]
                
                for entry in samples_pop:
                    key = entry["name"].split(":")[-1]
                    self.map_samples_to_pop[key] = pop
                    item = self.map_pop_to_samples.get(pop, [])
                    item.append(key)
                    self.map_pop_to_samples[pop] = item
    
    """ Get samples by population """
    def get_samples_by_pop(self, population):
        self._create_pop_sample_mappings()

        if population in self.map_pop_to_samples.keys():
            return self.map_pop_to_samples[population]
        else:
            warnings.warn("This population does not exist, returned all samples")
            return self.samples

    """ This method will return the population of a sample """
    def get_pop_by_sample(self, sample):
        self._create_pop_sample_mappings()

        if sample in self.map_samples_to_pop.keys():
            return self.map_samples_to_pop[sample]
        else:
            warnings.warn("This sample does not exist, returned ALL")
            return "ALL"

    """ Returns Variants from a genomic region of definied samples """
    def get_variants(self, chromosome, start, end, samples):
        key = (chromosome, start, end, samples)
        if not key in self.fetched_variants:
            self.fetched_variants[key] = list(map(lambda x: Variant(x),fetch_variants(chromosome, start, end, samples)))
        
        return self.fetched_variants[key]
        
    """ Returns a sequence by id"""
    def get_sequence_by_id(self, id, format = "plain"):
        key = (id, format)
        if not key in self.sequence_by_id:
            self.sequence_by_id[key] = fetch_sequence_by_id(id, format=format)
        
        return self.sequence_by_id[key]

    # TODO add sequence by region and so one...


class Individual:
    def __init__(id, population, variants = []):
        self.id = id
        self.population = population
        self.variants = variants
        self.allele = list()
        self.genotype = list()

    def generate_allele(variants):
        pass

    def generate_genotype(variants):
        pass




class Variant():
    def __init__(self, kwargs):
        self.INFO = dict()
        for key, value in kwargs.items():
            if key == "start":
                self.start = value
            elif key == "end":
                self.end = value
            elif key == "id":
                self.id = value
            elif key == "referenceBases":
                self.referenceBases = value
            elif key == "alternateBases":
                self.alternateBases = value
            elif key == "calls":
                self.calls = value
            else:
                self.INFO[key] = value

    def __str__(self):
        var = self.id[0] + ":" + str(self.start) + "..." + str(self.end) + ":" + self.referenceBases + "|"
        for alternate in self.alternateBases:
            var += alternate +","
    
        return var



