#Contains functions for manipulating and analyzing orthogroup files
import os
import numpy as np
from TreeUtils import findDomains
from ConfigParser import ConfigParser
from RetrieveOrthogroup import fastaToSeqs
from Similarity import domainSim, sequenceSim
from matplotlib import pyplot as plt
import seaborn as sns
from Fasta import FastaFile

def groupDomains(names, sequences, hmmfile):
    """
    Takes a list of input sequences and returns a list of domain strings for each.
    Leaves an empty string at position i of the jth list if the jth sequence does 
    not have a copy of domain i. This aligns all existing domains and makes it clear
    which domains are present in which sequence

    Example (domains marked as xxx):

    sequences = ["AAxxxAAxxxAAAAAAA",
                 "AAAAAAAxxxAAxxxAA",
                 "AAxxxAAAAAAAxxxAA"]


    grouped = [[dom1, dom2, ''  ],
              [''  , dom2, dom3],
              [dom1, ''  , dom3]]

    Args:
        sequences (list): A list of sequences 
        hmmfile   (str ): The name of the hmmfile containing the desired domain model

    Returns:
        grouped   (list): A list of lists of all domain sequences from each domain
        domNames  (list): A list of lists of domain names for each domain in each sequence
    """

    domStarts = [findDomains(i, hmmfile)[0] for i in sequences]
    domNames = []
    allStarts = sorted(list(set.union(*[set(i) for i in domStarts])))
    grouped = []
    for i in range(len(domStarts)):
        domains = ['' for _ in range(len(allStarts))]
        dnames = ['' for _ in range(len(allStarts))]
        for start in domStarts[i]:
            domains[allStarts.index(start)] = sequences[i][start: start+23]
            dnames[allStarts.index(start)] = names[i] + "_" + str(start)
        grouped.append(domains)
        domNames.append(dnames)
    return grouped, domNames

def filter(ff, hmmfile):
    """
    Filters an orthogroup to remove every entry with no domain instances.
    
    Args:
        ff (FastaFile): A fastafile object representation of the orthogroup
        hmmfile (str ): The file path to the hmm representation of the domain
    """

    for i in range(ff.length()):
        if len(ff.getDomains(i)) == 0:
            ff.delete(i)

def removeDuplicates(names, sequences, speciesIDs):
    """
    Removes duplicate entries, where bot the name and speciesID match

    Args:
        names: gene names
        sequences: raw gene sequences
        speciesIDs: species IDs

    """

    used = set()
    for i in range(len(names)-1, -1, -1):
        if (names[i], speciesIDs[i]) in used:
            names.pop(i)
            sequences.pop(i)
            speciesIDs.pop(i)
        else:
            used.add((names[i], speciesIDs[i]))
            names[i] += "_" + speciesIDs[i]  

def oneToOneSequences(ff):
    """
    Given a fasta file of an orthogroup, checks whether it is one to one.
    This means that the same number of sequences exist per species (this
    could mean two human and two chimp sequences, but not two human and 
    three chimp sequences. Returns True/False if the set is/isn't 1-1
    """
    sequences = ff.getAllSequences()
    speciesIDs = ff.getAllSpeciesIDs()

    seqsBySpecies = {}
    for i in range(len(sequences)):
        if speciesIDs[i] in seqsBySpecies:
            seqsBySpecies[speciesIDs[i]].append(sequences[i])
        else:
            seqsBySpecies[speciesIDs[i]] = [sequences[i]]

    
    firstLen = len(seqsBySpecies[seqsBySpecies.keys()[0]])
    for key in seqsBySpecies.keys():
        if len(seqsBySpecies[key]) != firstLen:
            return False

    return True

def oneToOneDomains(ff):
    domains = [ff.getDomains(i) for i in ff.length()]
    firstlen = len(domains[0])

    for line in domains:
        if len(line) != firstlen:
            return False
    
    return True

def selfSimilarity(name, sequence, hmmfile, heatmap=False):
    """
    Given a single sequence, checks the level of self similarity between 
    its constituent domains. Optionally creates a heatmap of this similarity

    Args:
        sequence (str ): An amino acid string representing a protein
        hmmfile  (str ): File path of the hmm used to find domains
        heatmap  (bool): (optional, default False) If true, displays a heatmap
                         of self similarity between domains on the sequence

    Output:
        simMatrix (list): A 
    """
    domains = findDomains(sequence, hmmfile)[2]
    numDomains = len(domains)
    simMatrix = np.zeros((numDomains, numDomains))

    for i in range(numDomains):
        for j in range(i, numDomains):
            simMatrix[i][j] = domainSim(domains[i], domains[j])
            simMatrix[j][i] = simMatrix[i][j]

    if heatmap:
        sns.heatmap(simMatrix, cmap='viridis')
        plt.savefig('tmp/' + name + '.pdf')
        plt.close()

    return simMatrix

def columnSimilarity(sequences, hmmfile):
    """
    Groups domains into columns based on msa and computes pairwise similarity
    within each column
    
    Args:
        sequences (list): raw sequences from an orthogroup to find domains in
        hmmfile   (str ): path to hmm representation of domain

    Output:
        A list of avg similarities per column
    """
    names = [str(i) for i in range(len(sequences))]
    g = groupDomains(names, sequences, hmmfile)[0]
    columns = [0. for i in range(len(g[0]))]
    for dom in range(len(g[0])):
        col = [o[dom] for o in g if o[dom] != '']
        count = 0
        for i in col:
            for j in col:
                columns[dom] += domainSim(i, j)
                count += 1

        columns[dom] /= count

    return columns

def averageSimilarity(ff):
    domains = ff.getAllDomains()


if __name__ == "__main__":
    data = ConfigParser.ORTHOGROUP_PATH #pylint: disable=no-member
    group = data + os.listdir(data)[0]
    hmmfile = ConfigParser.HMM_PATH #pylint: disable=no-member
    names, sequences, speciesIDs = fastaToSeqs(list(open(group)))
    grouped, domNames = groupDomains(names, sequences, hmmfile)
