#Contains functions for manipulating and analyzing orthogroup files
import os
from TreeUtils import findDomains
from ConfigParser import ConfigParser
from RetrieveOrthogroup import fastaToSeqs

#
def isValid(domain):
    """Checks if the input string is a valid zf-C2H2 domain"""
    valid = len(domain) == 23 and domain[2] == "C" and domain[5] == "C"
    valid &= domain[18] == "H" and domain[22] == "H"
    return valid

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

def filter(names, sequences, speciesIDs, hmmfile):
    """
    Filters an orthogroup to remove every entry with no domain instances.
    
    Args:
        names: The name of each sequence
        sequences: The sequences to be filtered
        speciesIDs:
        
    """

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

def oneToOneSequences(sequences, speciesIDs):
    """
    Given a fasta file of an orthogroup, checks whether it is one to one.
    This means that the same number of sequences exist per species (this
    could mean two human and two chimp sequences, but not two human and 
    three chimp sequences. Returns True/False if the set is/isn't 1-1
    """
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

def oneToOneDomains(sequences, hmmfile):
    #finds all domains occurences in each sequence
    domains = []
    for seq in sequences:
        domains.append(findDomains(seq, hmmfile)[2])

    firstlen = len(domains[0])


if __name__ == "__main__":
    data = ConfigParser.ORTHOGROUP_PATH #pylint: disable=no-member
    group = data + os.listdir(data)[0]
    hmmfile = ConfigParser.HMM_PATH #pylint: disable=no-member
    names, sequences, speciesIDs = fastaToSeqs(list(open(group)))
    grouped, domNames = groupDomains(names, sequences, hmmfile)