#Performs basic sequence evolution on a single domain.
#Only performs substitutions, ignoring indels

from ete3 import Tree
import numpy as np
from TreeUtils import isValid, printDomSeq
from stats import exp, drawFromDiscrete
from math import log
from pyvolve import Model, Partition, Evolver, read_tree
from OrthoAnalysis import selfSimilarity
from scipy.linalg import expm
from ConfigParser import ConfigParser as CP

#####################################
#                                   #
#        Atomic Sequence Ops        #
#                                   #
#####################################

def duplicate(sequence, starts, ends, domNumber, length):
    """
    Given a sequence and a domain number (ith domain in sequence), duplicates this 
    domain in the sequence.

    Args:
        sequence: The full gene sequence (including both domain and linker regions)
        starts: Starting position of every domain in the sequence
        ends: ending position of every domain in the sequence
        domNumber: The position of the domain to be duplicated w.r.t. the other domains
		length: Length (#domains) involved in duplication

    Returns:
        sequence (str ): The sequence after the specified duplication.	
    """
    BASELINKER = 'TGEVK'
    amount_added = ends[domNumber + length - 1] - starts[domNumber] + 1 + 5

    sequence = sequence[:ends[domNumber + length - 1] + 1] + BASELINKER + \
                    sequence[starts[domNumber]:]

    bstarts, bends = starts[:domNumber + length], ends[:domNumber + length]
    estarts = [i + amount_added for i in starts[domNumber:]]
    eends = [i + amount_added for i in ends[domNumber:]]

    starts, ends = bstarts + estarts, bends + eends

    return sequence, starts, ends

def remove(sequence, starts, ends, domNumber):
    """
    Deletes a specified domain in the input sequence. 

    Args:
        sequence: The full gene sequence (including both domain and linker regions)
        starts: Starting position of every domain in the sequence
        ends: ending position of every domain in the sequence
        domNumber: The position of the domain to be duplicated w.r.t. the other domains
    Returns sequence with specified duplication
    """
    #Removes one of the linkers if necessary

    amount_removed = ends[domNumber] - starts[domNumber] + 1

    if domNumber > 0:
        sequence = sequence[:starts[domNumber] - 5] + sequence[ends[domNumber]+1:]
        amount_removed += 5
    elif domNumber < len(starts) - 1:
        sequence = sequence[:starts[domNumber]] + sequence[ends[domNumber]+1+5:]
        amount_removed += 5
    else:
        sequence = sequence[:starts[domNumber]] + sequence[ends[domNumber]+1:]

    #Remove domain from list, adjust positions of the others
    starts.pop(domNumber), ends.pop(domNumber)
    for i in range(domNumber, len(starts)):
        starts[i] -= amount_removed
        ends[i] -= amount_removed

    return sequence, starts, ends

#####################################
#                                   #
#        Sequence Evolution         #
#                                   #
#####################################

alphabet = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 
                'T', 'V', 'W', 'Y']

def genTransitionMatrix(emissionProbs, transmat, branchLength):
    """
    Uses emission probabilities to re-weight transition matrix.
    Generates one transition matrix per position in the list of emission probabilities

    Args:
        emissionProbs (list): matrix with dimensions (n x 20) where n is the length of 
                              the domain. Each row contains the probability of each 
                              aa appearing at that position (in pfam hmm order)
        transmat (list):      aa transition matrix with dimensions (20 x 20)
    """

    #Turn into transition matrix
    transmat = np.asarray(transmat)
    transmat *= branchLength
    transmat = expm(transmat)
    
    #Zero out probability of staying in the same place
    for i in range(len(transmat)):
        transmat[i][i] = 0
    transmat = [i / sum(i) for i in transmat]

    return [transmat] * 100
    
    """
    out = []
    for position in emissionProbs:
        temp = [np.multiply(position, i) for i in transmat]
        temp = [i / sum(i) for i in temp]
        out.append(temp)

    return out
    """

def evolveDomain(sequence, rate, branchLength, emissionProbs, transmat, hmmfile):
    """
    Performs evolutionary process on a domain according to its hmm. Currently assumes zf_shilpa

    Args:
        sequence (str ):       The domain sequence to be evolved
        rate (float):          rate at which mutations occur
        branchLength (float):  distance to evolve sequence
        emissionsProbs (list): matrix with dimensions (n x 20) where n is the length of 
                               the domain. Each row contains the probability of each 
                               aa appearing at that position (in pfam hmm order)
        transmat (list):       aa transition matrix with dimensions (20 x 20) 
    """

    invalid = True
    
    #Returns the number of mutations that occur on a branch with time t
    def numMutations(t):
        count = 0
        t -= exp((1./rate) * 1./len(sequence))
        while t > 0:
            count += 1
            t -= exp(rate * 1./len(sequence))
        return count

    #Entropy of a probability distribution
    def entropy(line):
        return [-1 * p * log(p, 2) for p in line]

    entropies = [sum(entropy(i)) for i in emissionProbs]
    entropies = [i / sum(entropies) for i in entropies]

    #TODO: change so we only generate matrices for the required positions?
    transitions = genTransitionMatrix(emissionProbs, transmat, branchLength)
        
    nMuts = numMutations(branchLength)
    invalidCounter = 0
    while invalid:
        seqCopy = sequence
        for i in range(nMuts):
            #position = drawFromDiscrete(entropies)
            #character = drawFromDiscrete(emissionProbs[position])
            
            position = np.random.choice(range(len(seqCopy)))
            oldChar = alphabet.index(sequence[position])
            row = transitions[position][oldChar]
            character = drawFromDiscrete(row)
            
            seqCopy = seqCopy[:position] + alphabet[character] + seqCopy[position+1:]
        
        #invalid = not isValid(seqCopy, hmmfile)
        invalid = False
        invalidCounter += 1
        if invalidCounter >= 25:
            raise ValueError

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
def evolveSequence(sequence, starts, ends, hmmfile, rate, branchLength, emissionProbs, transmat):
    """
    Putting the previous steps together, simulates evolutoin of a full sequence including 
    both domains and non-domain sequence. 

    Args:
        sequence (str):  The full sequence to be evolved
        starts: Starting position of every domain in the sequence
        ends: ending position of every domain in the sequence (+1)
        emissionsProbs:  matrix with dimensions (n x 20) where n is the length of 
                         the domain. Each row contains the probability of each 
                         aa appearing at that position (in pfam hmm order) 
        transmat (list): aa transition matrix with dimensions (20 x 20)
    """
    #FOR TESTING
    original_sequence = sequence
    #END TESTING BLOCK

    domains = [sequence[starts[i]:ends[i]+1] for i in range(len(starts))]

    #split on all domains
    for seq in domains:
        sequence = sequence.replace(seq, "xxx")
    sequences = sequence.split("xxx")

    #Evolve sequence fragments individually
    for i in range(len(domains)):
        domains[i] = evolveDomain(domains[i], rate, branchLength, emissionProbs, transmat, hmmfile)

    for i in range(len(sequences)):
        if sequences[i] == '':
            sequences[i] = evolveEmptyLinker(branchLength)
        else:
            sequences[i] = evolveLinker(sequences[i], branchLength)

    #Reassemble full sequence post evolution
    sequence = ''
    try:
        for i in range(len(domains)):
            sequence += sequences[i] + domains[i]
    except:
        print "EVOLVE SEQUENCE ERROR"
        print original_sequence
        print domains
        print sequences
        printDomSeq(original_sequence, hmmfile)
        raise Exception

    sequence += sequences[-1] if len(sequences) > len(domains) else domains[-1]

    return sequence

#TODO: Need a better name for this function
def domainOrder(sequence, starts, ends, rate, emissionProbs, hmmfile, tree, hnodename, transmat):
    """
    Evolves sequence along domain subtree within a single host node. Assumes that every
    node in the domain tree has an event assigned. Requires that every root has a position
    and the root sequence has the correct number of domains (matches the tree)

    Args:
        
    """
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
        jobs.sort() #sorts jobs by distance from the root
        dist, node = jobs.pop(0)
        workingSet.remove(node)
        if dist > 0:
            #print 'evolving for:', dist, '\n'
            sequence = evolveSequence(sequence, starts, ends, hmmfile, rate, dist, emissionProbs, transmat)

        for job in jobs:
            job[0] -= dist

        #Only need to deal with duplication and loss events explicitly;
        #speciation and leaf nodes require no work
        if node.event == 'DUPLICATION':

            #Find all other nodes marked in the same tandem duplication
            td = [(node.position, node)]
            for i in range(len(jobs) - 1, -1, -1):
                otherNode = jobs[i][1]
                if otherNode.event == 'DUPLICATION' and otherNode.dupNumber == node.dupNumber:
                    td.append((otherNode.position, jobs.pop(i)[1]))

            td.sort()
            size = len(td)

            sequence, starts, ends = duplicate(sequence, starts, ends, td[0][1].position, size)

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

        if node.event == 'LOSS':
            sequence, starts, ends = remove(sequence, starts, ends, node.position)
            for otherNode in workingSet:
                if otherNode.position >= node.position: #Why does this have to be >= instead of > ?
                    otherNode.position -= 1

    return sequence, starts, ends

def evolveAlongTree(host, guest, reverseMap, rootSequence, rootStarts, rootEnds, 
                        hmmfile, emissionProbs, transmat):
    """
    Evolves a root sequence along an entire host tree, taking into account the domain level 
    events present in the guest tree (duplication, loss, speciation)

    Args:
        host (Tree)       : The host tree (ete3 format) inside which the guest tree evolved
        guest (Tree)      : The guest tree over which to evolve a sequence 
        reverseMap (dict) : mapping from nodes in the host node -> guest nodes 
        rootSequence (str): Initial sequence to evolve. Should contain sequence with ONE domain
        rootStarts (list) :
        rootEnds (list)   :
        hmmfile (str )    : path to hmmfile used to identify domains
        emissionProbs     : matrix with dimensions (n x 20) where n is the length of 
                            the domain. Each row contains the probability of each 
                            aa appearing at that position (in pfam hmm order) 
    """

    for node in host.traverse():
        node.add_feature('sequence', "")
        node.add_feature('starts', [])
        node.add_feature('ends', [])

    for hostNode in host.traverse():
        tempSequence = rootSequence if hostNode == host else hostNode.up.sequence
        starts = rootStarts if hostNode == host else hostNode.up.starts
        ends = rootEnds if hostNode == host else hostNode.up.ends

        #No events occured at this node
        if hostNode not in reverseMap:
            hostNode.sequence = evolveSequence(tempSequence, starts, ends, hmmfile, 1, \
                                    hostNode.dist, emissionProbs, transmat)
            hostNode.starts = starts
            hostNode.ends = ends
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
        tempSequence, starts, ends = domainOrder(tempSequence, starts, ends, \
                                                .75, emissionProbs, hmmfile, \
                                                t, hostNode.name, transmat)
        hostNode.sequence = tempSequence
        hostNode.starts = starts
        hostNode.ends = ends

        #Reconnect all root and leaf nodes to the rest of the guest tree
        for node in upAncestors:
            node.up = upAncestors[node]
        for node in leafChildren:
            node.children = leafChildren[node]

if __name__ == '__main__':

    from GuestTreeGen import expon, buildGuestTree
    from HostTreeGen import createRandomTopology
    from stats import gaussNoise
    from rootSequence import genRandomSequence2 as grs
    import pickle
    from ConfigParser import ConfigParser as CP

    eppath = CP.EP_PATH #pylint: disable=no-member
    emissionProbs = pickle.load(open(eppath))
    hmmfile = CP.HMM_PATH #pylint: disable=no-member

    def s2(x):
        denom = 1 + exp(10-x) if x < 10 else 1
        return 1 - .6 / denom

    def expfunc(minimum=1, maximum=3):
        """Exponential distribution with lambda=0.75 and min/max parameters"""
        return min(maximum, max(minimum, int(expon(0.75).rvs())))

    def selfSim(seqs):
        return selfSimilarity('asdf', seqs, hmmfile, False)

    sd = 5 #startingDomains

    hostTree = createRandomTopology(1, 1, lambda x: x)
    guestTree, nodeMap = buildGuestTree(hostTree, s2, expfunc, .2, gaussNoise, sd)

    rootSequence = grs(sd)
    transmat = pickle.load(CP.TRANSMAT) #pylint: disable=no-member
    #evolveAlongTree(hostTree, guestTree, nodeMap, rootSequence, hmmfile, emissionProbs, transmat)
    selfSimilarity('asdf', hostTree.sequence, hmmfile, True) #pylint: disable=no-member
