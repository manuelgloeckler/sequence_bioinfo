import requests, sys
import warnings
import pickle


class Ensembl():
    def __init__(self):
        self.samples = list()
        self.populations = list()
        self.population_info = dict()
        self.samples_info = dict()
        self.map_samples_to_pop = dict()
        self.map_pop_to_samples = dict()

    """ Save state of the collector """
    def save(self, path):
        with  open(path, 'w') as file1: 
            pickle.dump(self, file1)

    """ Static Function to load a state """
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
    
    """ Get samples by population """
    def get_samples_by_pop(self, population):
        if len(self.samples) == 0:
            self.get_samples()

        if len(self.populations) == 0:
            self.get_populations()

        if len(self.map_pop_to_samples) == 0:
            for pop in self.populations:
                samples_pop = fetch_samples(pop = "1000GENOMES:phase_3:" + population)[0]["individuals"]
                
                for entry in samples_pop:
                    key = entry["name"].split(":")[-1]
                    self.map_samples_to_pop[key] = pop
                    item = self.map_pop_to_samples.get(pop, [])
                    item.append(key)
                    self.map_pop_to_samples[pop] = item

        if population in self.map_pop_to_samples.keys():
            return self.map_pop_to_samples[population]
        else:
            warnings.warn("This population does not exist, returned all samples")
            return self.samples

    """ This method will return the population of a sample """
    def get_pop_by_sample(self, sample):
        if len(self.map_samples_to_pop) == 0:
            this.get_samples_by_pop()

        if sample in map_samples_to_pop.keys():
            return map_samples_to_pop[sample]
        else:
            warnings.warn("This sample does not exist, returned ALL")
            return "ALL"

    # TODO make classes of sample, with corresponding variation and so one.
    def get_variants(self, chromosome, start, end, samples):
        return fetch_variants(chromosome, start, end, samples)

    
""" Fetches all samples from ENSEMBL and reprot dictionary """
def fetch_samples(pop = "1000GENOMES:phase_3:ALL"):
    server = "https://grch37.rest.ensembl.org"
    ext = "/info/variation/populations/human/" + pop + "?"

    r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})

    if not r.ok:
        r.raise_for_status()
        sys.exit()

    return r.json()


""" Fetches all variants for a defined start and endpoint in defined samples """            
def fetch_variants(chromosome,
                start,
                end,
                samples,
                pageSize = 10,
                variantSetId = 3, 
                referenceName = 22):
    decoded = _fetch_variants_worker(chromosome, start, end, samples, pageSize, variantSetId, referenceName, None)
    variants = decoded["variants"]
    nextPage = decoded["nextPageToken"]

    while(nextPage != None):
        decoded = _fetch_variants_worker(chromosome, start, end, samples, pageSize, variantSetId, referenceName, nextPage) 
        variants.append(decoded["variants"])
        nextPage = decoded["nextPageToken"]

    return variants

def _fetch_variants_worker(chromosome,
                start,
                end,
                samples,
                pageSize,
                variantSetId, 
                referenceName,
                nextPageToken):
    server = "https://grch37.rest.ensembl.org"
    ext = "/ga4gh/variants/search"
    headers={ "Content-Type" : "application/json", "Accept" : "application/json"}

    callSetIds = str([str(chromosome) + ':' + sample  for sample in samples]).replace("'", '"')
    if nextPageToken == None:
        nextPageToken = '"None"'
    else:
        nextPageToken = str(nextPageToken)

    r = requests.post(server+ext,
        headers=headers,
        data='{"pageSize": ' + str(pageSize) +',"variantSetId": ' + str(variantSetId) + ', "callSetIds" :' + callSetIds +\
                ', "referenceName": ' + str(referenceName) + ', "start": '+ str(start) +' , "end": ' + str(end) +\
                ', "pageToken": '+ str(nextPageToken) + ' }')
    
    if not r.ok:
        r.raise_for_status()
        sys.exit()
    
    return r.json()
            
""" Fetches all populations from ENSEMBL """
def fetch_populations(filter = "filter=LD"):
    server = "https://grch37.rest.ensembl.org"
    ext = "/info/variation/populations/homo_sapiens?" + filter
    r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})

    if not r.ok:
        r.raise_for_status()
        sys.exit()

    return r.json()




