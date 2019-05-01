#Performs basic sequence evolution on a single domain.
#Only performs substitutions, ignoring indels

from ete3 import Tree
import numpy as np
from TreeUtils import findDomains, isValid, printDomSeq
from stats import exp, drawFromDiscrete
from math import log
from pyvolve import Model, Partition, Evolver, read_tree
from OrthoAnalysis import selfSimilarity

#####################################
#                                   #
#        Atomic Sequence Ops        #
#                                   #
#####################################

def duplicate(sequence, hmmfile, domNumber, length):
    """
    Given a sequence and a domain number (ith domain in sequence), duplicates this 
    domain in the sequence.

    Args:
        sequence: The full gene sequence (including both domain and linker regions)
        hmmfile: hmm file containing hmm of domain to be duplicated
        domNumber: The position of the domain to be duplicated w.r.t. the other domains
		length: Length (#domains) involved in duplication

    Returns:
        sequence (str ): The sequence after the specified duplication.	
    """
    BASELINKER = 'TGEVK'
    starts, ends = findDomains(sequence, hmmfile)[:2]

    sequence = sequence[:ends[domNumber + length - 1] + 1] + BASELINKER + \
                    sequence[starts[domNumber]:]

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
def evolveDomain(sequence, rate, branchLength, emissionProbs):
    """
    Performs evolutionary process on a domain according to its hmm. Currently assumes zf_shilpa

    Args:
        sequence (str ):       The domain sequence to be evolved
        rate (float):          rate at which mutations occur
        branchLength (float):  distance to evolve sequence
        emissionsProbs (list): matrix with dimensions (n x 20) where n is the length of 
                               the domain. Each row contains the probability of each 
                               aa appearing at that position (in pfam hmm order)
    """

    invalid = True
    
    #Returns the number of mutations that occur on a branch with time t
    def numMutations(t):
        count = 0
        t -= exp(rate)
        while t > 0:
            count += 1
            t -= exp(rate)
        return count

    #information content of a probability distribution
    def ic(line):
        return [-1 * p * log(p) for p in line]

    ics = [sum(ic(i)) for i in emissionProbs]
    ics = ics / sum(ics)
        
    nMuts = numMutations(branchLength)
    while invalid:
        seqCopy = sequence
        for i in range(nMuts):
            position = drawFromDiscrete(ics)
            character = drawFromDiscrete(emissionProbs[position])
            seqCopy = seqCopy[:position] + alphabet[character] + seqCopy[position+1:]
        invalid = not isValid(seqCopy)
        if invalid:
            print 'yikes'

    #assert(isValid(seqCopy))
    return seqCopy

def evolveLinker(sequence, branchLength):
    """
    Evolves non-domain sequence a specified distance using pyvolve. Simulates substitutions only
    (no indels). branchLength * sequence is the expected fraction of positions to mutate (with
    replacement). Returns sequence post modification.
    """
    m = Model("JTT")
    p = Partition(models = m, root_sequence = sequence)
    #t = (A:BL,b:BL)
    t = read_tree(tree = "(A:" + str(branchLength) + ",B:" + str(branchLength) + ");") 
    e = Evolver(partitions=p, tree=t)
    e()
    return e.get_sequences()["A"]

def evolveEmptyLinker(branchLength):
    """Inserts and evolves a template linker when no linker exists"""
    BASELINKER = 'TGEVK'
    return evolveLinker(BASELINKER, branchLength)

#Simulates evolution of a full sequence including both domains and linker regions
#Splits sequence into domain and linker regions, deals with each independently
def evolveSequence(sequence, rate, branchLength, emissionProbs, hmmfile):
    """
    Putting the previous steps together, simulates evolutoin of a full sequence including 
    both domains and non-domain sequence. 

    Args:
        sequence (str): The full sequence to be evolved
        emissionsProbs: matrix with dimensions (n x 20) where n is the length of 
                        the domain. Each row contains the probability of each 
                        aa appearing at that position (in pfam hmm order) 
    """

    #Find domains, check if sequence begins and/or ends with a domain
    domains = findDomains(sequence, hmmfile)[2]

    #split on all domains
    for seq in domains:
        sequence = sequence.replace(seq, "xxx")
    sequences = sequence.split("xxx")

    #Evolve sequence fragments individually
    for i in range(len(domains)):
        domains[i] = evolveDomain(domains[i], rate, branchLength, emissionProbs)

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

#TODO: Need a better name for this function
def domainOrder(sequence, rate, hmmfile, emissionProbs, tree, hnodename):
    """
    Evolves sequence along domain subtree within a single host node. Assumes that every
    node in the domain tree has an event assigned. Requires that every root has a position
    and the root sequence has the correct number of domains (matches the tree)

    Args:
        
    """
    starts, ends, seqs = findDomains(sequence, hmmfile)

    jobs = []
    workingSet = set()

    #Assumes that the root node given is a dummy node
    #Also assumes that all children of the fake root have a 'position' field
    for child in tree.children:
        #child.domain = seqs[child.position]
        jobs.append([child.dist, child])
        workingSet.add(child)

    #process these roots one at a time
    while len(jobs) != 0:
        starts, ends, seqs = findDomains(sequence, hmmfile)
        jobs.sort() #sorts jobs by distance from the root
        dist, node = jobs.pop(0)
        workingSet.remove(node)
        if dist > 0:
            print 'evolving for:', dist, '\n'
            sequence = evolveSequence(sequence, rate, dist, emissionProbs, hmmfile)

        #node.domain = seqs[node.position]

        for job in jobs:
            job[0] -= dist

        #Only need to deal with duplication and loss events explicitly;
        #speciation and leaf nodes require no work
        if node.event == 'DUPLICATION':

            print (node.name, node.dupNumber), [(thing[1].name, thing[1].dupNumber) for thing in jobs]

            #Find all other nodes marked in the same tandem duplication
            td = [(node.position, node)]
            for i in range(len(jobs) - 1, -1, -1):
                otherNode = jobs[i][1]
                if otherNode.event == 'DUPLICATION' and otherNode.dupNumber == node.dupNumber:
                    td.append((otherNode.position, jobs.pop(i)[1]))

            td.sort()
            size = len(td)

            sequence = duplicate(sequence, hmmfile, td[0][1].position, size)

            #increment position of all the succeeding domains
            for otherNode in workingSet:
                if otherNode.position > td[-1][1].position:
                    #should already have a position feature, if not there is a problem
                    otherNode.position += size 

            #add position to children
            for otherNode in [i[1] for i in td]:
                a,b = otherNode.children
                a.add_feature('position', otherNode.position)
                b.add_feature('position', otherNode.position + size)
                workingSet.add(a)
                workingSet.add(b)
                jobs.append([a.dist, a])
                jobs.append([b.dist, b])

            print 'DUPLICATION, size ' + str(size) + 'at ' + str(td[0][1].position) + \
                    ' ,now ' + str(len(starts) + size) + ' domains'
            printDomSeq(sequence, hmmfile)
            print tree.get_ascii(attributes=['position'])

        if node.event == 'LOSS':
            sequence = remove(sequence, hmmfile, node.position)
            for otherNode in workingSet:
                if otherNode.position >= node.position: #Why does this have to be >= instead of > ?
                    otherNode.position -= 1
                    
            print 'LOSS at ' + str(node.position) + ' ,now ' + str(len(starts) - 1) + ' domains'
            printDomSeq(sequence, hmmfile)
            print tree.get_ascii(attributes=['position'])

        #selfSimilarity(hnodename, sequence, hmmfile, True)

    return sequence

def evolveAlongTree(host, guest, reverseMap, rootSequence, hmmfile, emissionProbs):
    """
    Evolves a root sequence along an entire host tree, taking into account the domain level 
    events present in the guest tree (duplication, loss, speciation)

    Args:
        host (Tree)       : The host tree (ete3 format) inside which the guest tree evolved
        guest (Tree)      : The guest tree over which to evolve a sequence 
        reverseMap (dict) : mapping from nodes in the host node -> guest nodes 
        rootSequence (str): Initial sequence to evolve. Should contain sequence with ONE domain
        hmmfile (str )    : path to hmmfile used to identify domains
        emissionProbs     : matrix with dimensions (n x 20) where n is the length of 
                            the domain. Each row contains the probability of each 
                            aa appearing at that position (in pfam hmm order) 
    """

    for node in host.traverse():
        node.add_feature('sequence', "")

    for hostNode in host.traverse():
        tempSequence = rootSequence if hostNode == host else hostNode.up.sequence

        #No events occured at this node
        if hostNode not in reverseMap:
            hostNode.sequence = evolveSequence(tempSequence, 0.1, hostNode.dist, \
                                    emissionProbs, hmmfile)
            continue

        allGuestNodes = reverseMap[hostNode]
        allGuestNodesSet = set(allGuestNodes)
        upAncestors, leafChildren = {}, {}

        for guestNode in allGuestNodes:
            if guestNode.up not in allGuestNodesSet:
                upAncestors[guestNode] = guestNode.up
                #pass positional information on from the previous species
                if guestNode.up != None:
                    guestNode.add_feature('position', guestNode.up.position)
                guestNode.up = None
            if guestNode.children != [] and guestNode.children[0] not in allGuestNodesSet:
                leafChildren[guestNode] = guestNode.children
                guestNode.children = []

        if hostNode != host:
            t = Tree()
            t.dist = 0
            t.children = upAncestors.keys()
            for guestNode in upAncestors.keys():
                guestNode.up = t

        else:
            t = guest

        #Actually do the work
        tempSequence = domainOrder(tempSequence, 0.1, hmmfile, emissionProbs, t, hostNode.name)
        hostNode.sequence = tempSequence

        #Reconnect all root and leaf nodes to the rest of the guest tree
        for node in upAncestors:
            node.up = upAncestors[node]
        for node in leafChildren:
            node.children = leafChildren[node]