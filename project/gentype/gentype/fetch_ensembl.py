import sys, requests


""" Fetches all samples from ENSEMBL and reprot dictionary """
def fetch_samples(pop = "ALL", project = "1000GENOMES:phase_3"):
    server = "https://rest.ensembl.org"
    ext = "/info/variation/populations/human/"+ project + ":" + pop + "?"

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

    variants = []
    nextPage = None
    i = 0
    while True:
        decoded = _fetch_variants_worker(chromosome, start, end, samples, pageSize, variantSetId, referenceName, nextPage) 
        variants.extend(decoded["variants"])
        nextPage = decoded["nextPageToken"]
        i+=pageSize
        print("Fetched " + str(i) + " Variants")

        if nextPage == None:
            break


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
def fetch_populations(filter = "filter=LD", species = "homo_sapiens"):
    server = "https://grch37.rest.ensembl.org"
    ext = "/info/variation/populations/" + species + "?" + filter
    r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})

    if not r.ok:
        r.raise_for_status()
        sys.exit()

    return r.json()


def ping():
    server = "https://grch37.rest.ensembl.org"
    ext = "/info/ping?"
    
    r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})
    
    if not r.ok:
        r.raise_for_status()
        sys.exit()
    
    return bool(r.json()["ping"])

def fetch_sequence_by_id(id, format = 'plain'):
    server = "https://rest.ensembl.org"
    ext = "/sequence/id/" + id + "?"
    
    r = requests.get(server+ext, headers={ "Content-Type" : "text/" + format})
    
    if not r.ok:
        r.raise_for_status()
        sys.exit()
    
    
    return r.text

def fetch_sequence_by_region(chromosome, start, end, direction, species = "human"):
    server = "https://rest.ensembl.org"
    ext = "/sequence/region/" + species + "/" + str(chromosome) + ":" + str(start) + ".." + str(end) + ":" + str(direction)+ "?"
 
    r = requests.get(server+ext, headers={ "Content-Type" : "text/plain"})
 
    if not r.ok:
        r.raise_for_status()
        sys.exit()
 
    return r.text
