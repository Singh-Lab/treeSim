#Functions for retrieving and manipulating orthogroups from OrthoDB

from urllib2 import urlopen
import argparse
import json
from time import sleep
import sys
import numpy as np
import matplotlib.pyplot as plt
import pickle
from random import randint
from ConfigParser import ConfigParser as CP

baseUrl = 'http://www.orthodb.org/'
searchUrl = baseUrl + 'search?'
fastaUrl = baseUrl + 'fasta?'

#Finds orthogroup IDs related to the query gene within the orthogroup scope
def findOrthogroupID(gene, level='9347', species='9347'):
    """
    Finds Orthogroup IDs related to the query gene within the orthogroup scope.
    
    Args:
        gene    (str ): The gene name whose orthogroup we'd like to find 
                        (e.g. ENSG00000196670)
        level   (str ): NCBI Taxonomy browser level ID
        species (str ): The NCBI Taxonomy browser ID for the species group,
                        usually the same as the level ID

    Output:
        ID (str ): The orthogroup ID ('XXXXatYYYY') the gene is a part of
    """
    url = searchUrl + 'query=' + gene + '&level=' + level + '&species=' + species
    response = json.loads(urlopen(url).read())
    if len(response['data']) == 0:
        return ''
    return str(response['data'][0])

def getFasta(oID, level='9347', species='9347'):
    """
    Retrieves raw fasta of an orthogroup given its ID. 
    Useful Orthogroup levels:
    Human: Species ID = 9606
    Primates: level = species = 9443
    Eutheria: level = species = 9347

    Args:
        oID (str ): The target orthogroup's ID
        level   (str ): NCBI Taxonomy browser level ID
        species (str ): The NCBI Taxonomy browser ID for the species group,
                        usually the same as the level ID
    
    Output:
        A list of strings, one for each line of the fasta file (headers + seqs)
    """
    url = fastaUrl + 'id=' + oID + '&level=' + level + '&species=' + species
    sleep(1)
    return urlopen(url).read().split('\n')

def getGroup(oid, src):
    prefix = CP.DATAPATH + src + '/' #pylint: disable=no-member
    return list(open(prefix + oid))

def fastaToSeqs(fasta):
    """
    Extracts sequences, sequence names, and species IDs per sequence from a 
    fasta file in list format

    Args:
        fasta (list): A list of strings in the format
                      <HEADER_1>
                      <SEQUENCE_1>
                      ...
                      <HEADER_N>
                      <SEQUENCE_N>

    Output:
        names      (list): ith entry is the PubGene ID of the ith sequence in the file.
        sequences  (list): the ith sequence in the file
        speciesIDs (list): ith entry is the NCBI Taxonomy browser ID of 
                           the ith sequence
    """
    headers = []
    sequences = []

    skipped = False
    for i in range(len(fasta)):
        if skipped:
            skipped = False
            continue
        if fasta[i] == '' or ('pub_gene_id' not in fasta[i] and i % 2 == 0):
            skipped = True
            continue
        if i % 2 == 0:
            headers.append(fasta[i])
        else:
            sequences.append(fasta[i])

    speciesIDs = [line.split(":")[0][1:] for line in headers]
    names = [line.split('"pub_gene_id":')[1].split(',')[0].strip('"') for line in headers]
    names = [thing.split(';')[1] if ';' in thing else thing for thing in names]

    return (names, sequences, speciesIDs)