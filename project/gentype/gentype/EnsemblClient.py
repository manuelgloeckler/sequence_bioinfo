#!/usr/bin/env python

# code is slightly altered version from https://github.com/Ensembl/ensembl-rest/wiki/Example-Python-Client

import sys
import json
import time
import logging

logger = logging.getLogger()
# Python 2/3 adaptability
try:
    from urllib.parse import urlparse, urlencode
    from urllib.request import urlopen, Request
    from urllib.error import HTTPError
except ImportError:
    from urlparse import urlparse
    from urllib import urlencode
    from urllib2 import urlopen, Request, HTTPError

class EnsemblClient(object):
    def __init__(self, server='https://grch37.rest.ensembl.org', reqs_per_sec=15):
        self.server = server
        self.reqs_per_sec = reqs_per_sec
        self.req_count = 0
        self.last_req = 0

    def perform_rest_action(self, endpoint, hdrs=None, params=None, data=None):
        if hdrs is None:
            hdrs = {}

        if 'Content-Type' not in hdrs:
            hdrs['Content-Type'] = 'application/json'

        if params:
            endpoint += '?' + urlencode(params)

        # check if we need to rate limit ourselves
        if self.req_count >= self.reqs_per_sec:
            delta = time.time() - self.last_req
            if delta < 1:
                time.sleep(1 - delta)
            self.last_req = time.time()
            self.req_count = 0
        
        try:
            if not data is None: data = str(json.dumps(data)).encode()
            request = Request(self.server + endpoint, headers=hdrs, data=data)
            response = urlopen(request)
            content = response.read()
            if content:
                data = json.loads(content)
            self.req_count += 1

        except HTTPError as e:
            # check if we are being rate limited by the server
            if e.code == 429:
                if 'Retry-After' in e.headers:
                    retry = e.headers['Retry-After']
                    logger.warning("Got time restricted by Ensembl-Server. Will retry in {}.".format(retry))
                    time.sleep(float(retry))
                    self.perform_rest_action(endpoint, hdrs, params = {})
            else:
                sys.stderr.write('Request failed for {0}: Status code: {1.code} Reason: {1.reason}\n'.format(endpoint, e))
           
        return data

    def get_variants(self, species, symbol):
        genes = self.perform_rest_action(
            endpoint='/xrefs/symbol/{0}/{1}'.format(species, symbol), 
            params={'object_type': 'gene'}
        )
        if genes:
            stable_id = genes[0]['id']
            variants = self.perform_rest_action(
                '/overlap/id/{0}'.format(stable_id),
                params={'feature': 'variation'}
            )
            return variants
        return None


def run(species, symbol):
    client = EnsemblRestClient()
    #data = client.perform_rest_action('/variation/human/rs56116432', params = {'genotypes' : 1})
    #print(data)
    variants = client.get_variants(species, symbol)
    if variants:
        for v in variants: #currently we only take the first 10 to save some time
            if (v['id'] == 'rs56116432'): print("OI")
            data = client.perform_rest_action('/variation/{}/{}'.format(species, v['id']), params = {'genotypes' : 1})
            if data:
                print(data['genotypes'])

if __name__ == '__main__':
    if len(sys.argv) == 3:
        species, symbol = sys.argv[1:]
    else:
        species, symbol = 'human', 'ABO'

    run(species, symbol)