#Performs basic sequence evolution on a single domain.
#Only performs substitutions, ignoring indels

import numpy as np
from TreeUtils import findDomains, isValid
from stats import exp, drawFromDiscrete
from math import log
from pyvolve import Model, Partition, Evolver, read_tree

#####################################
#                                   #
#        Atomic Sequence Ops        #
#                                   #
#####################################

def duplicate(sequence, hmmfile, domNumber):
    """
    Given a sequence and a domain number (ith domain in sequence), duplicates this 
    domain in the sequence.

    Args:
        sequence: The full gene sequence (including both domain and linker regions)
        hmmfile: hmm file containing hmm of domain to be duplicated
        domNumber: The position of the domain to be duplicated w.r.t. the other domains

    Returns:
        sequence (str ): The sequence after the specified duplication.	
    """
    BASELINKER = 'TGEVK'
    ends, domSeqs = findDomains(sequence, hmmfile)[1:]
    sequence = sequence[:ends[domNumber]+1] + BASELINKER + \
                    domSeqs[domNumber] + sequence[ends[domNumber]+1:]
    return sequence

def remove(sequence, hmmfile, domNumber):
    """
    Deletes a specified domain in the input sequence. 

    Args:
        sequence: The full gene sequence (including both domain and linker regions)
        hmmfile: hmm file containing hmm of domain to be removed
        domNumber: The position of the domain to be duplicated w.r.t. the other domains
    Returns sequence with specified duplication
    """
    starts, ends = findDomains(sequence, hmmfile)[:2]

    #Removes one of the linkers if necessary
    if domNumber > 0:
        sequence = sequence[:starts[domNumber] - 5] + sequence[ends[domNumber]+1:]
    elif domNumber < len(starts) - 1:
        sequence = sequence[:starts[domNumber]] + sequence[ends[domNumber]+1+5:]
    else:
        sequence = sequence[:starts[domNumber]] + sequence[ends[domNumber]+1:]
    return sequence

#####################################
#                                   #
#        Sequence Evolution         #
#                                   #
#####################################

alphabet = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 
                'T', 'V', 'W', 'Y']

#Performs evolutionary process on a domain according to its hmm. Currently assumes zf_shilpa
def evolveDomain(sequence, rate, branchLength, emissionProbabilities):

    invalid = True
    
    #Returns the number of mutations that occur on a branch with time t
    def numMutations(t):
        count = 0
        t -= exp(rate)
        while t > 0:
            count += 1
            t -= exp(rate)
        return count

    def ic(line):
        return [-1 * p * log(p) for p in line]

    ics = [sum(ic(i)) for i in emissionProbabilities]
    ics = ics / sum(ics)
        
    nMuts = numMutations(branchLength)
    while invalid:
        seqCopy = sequence
        for i in range(nMuts):
            position = drawFromDiscrete(ics)
            character = drawFromDiscrete(emissionProbabilities[position])
            seqCopy = seqCopy[:position] + alphabet[character] + seqCopy[position+1:]
        invalid = not isValid(seqCopy)

    return seqCopy

def evolveLinker(sequence, branchLength):
    """
    Evolves non-domain sequence a specified distance using pyvolve. Simulates substitutions only
    (no indels). branchLength * sequence is the expected fraction of positions to mutate (with
    replacement). Returns sequence post modification.
    """
    m = Model("JTT")
    p = Partition(models = m, root_sequence = sequence)
    t = read_tree(tree = "(A:" + str(branchLength) + ",B:" + str(branchLength) + ");") #(A:BL,b:BL)
    e = Evolver(partitions=p, tree=t)
    e()
    return e.get_sequences()["A"]

def evolveEmptyLinker(branchLength):
    """Inserts and evolves a template linker when no linker exists"""
    BASELINKER = 'TGEVK'
    return evolveLinker(BASELINKER, branchLength)

#Simulates evolution of a full sequence including both domains and linker regions
#Splits sequence into domain and linker regions, deals with each independently
def evolveSequence(sequence, rate, branchLength, emissionProbabilities, hmmfile):
    """
    TODO: Fill in this comment.
    """

    #Find domains, check if sequence begins and/or ends with a domain
    domains = findDomains(sequence, hmmfile)[2]

    #split on all domains
    for seq in domains:
        sequence = sequence.replace(seq, "xxx")
    sequences = sequence.split("xxx")

    #Evolve sequence fragments individually
    for i in range(len(domains)):
        domains[i] = evolveDomain(domains[i], rate, branchLength, emissionProbabilities)

    for i in range(len(sequences)):
        if sequences[i] == '':
            sequences[i] = evolveEmptyLinker(branchLength)
        else:
            sequences[i] = evolveLinker(sequences[i], branchLength)

    #Reassemble full sequence post evolution
    sequence = ''
    for i in range(len(domains)):
        sequence += sequences[i] + domains[i]

    sequence += sequences[-1] if len(sequences) > len(domains) else domains[-1]

    return sequence