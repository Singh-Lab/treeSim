import sys
from os import system
#sys.path.append('treeSim/')
from ete3 import Tree
from HostTreeGen import createRandomTopology
from OrthoAnalysis import selfSimilarity, groupDomains, columnSimilarity
from OrthoNJ import prune, mlTree, mlTree2
from GuestTreeGen import buildGuestTree, exp, expon
from CreateSequences import evolveAlongTree, evolveSequence
from stats import gaussNoise
from rootSequence import genRandomSequence2 as grs
from TreeUtils import writeTree, writeFasta, findDomains
from random import randint
import pickle
import numpy as np

eppath = '/home/caluru/Documents/shilpa/treeSimulation/simulation/zf_shilpa_probabilities.pickle'
emissionProbs = pickle.load(open(eppath))
hmmfile = '/home/caluru/Data/hmmfiles/zf_shilpa_232.hmm'

def s2(x):
    denom = 1 + exp(7-x)
    return 1 - .6 / denom

def expfunc(minimum=1, maximum=3):
    """Exponential distribution with lambda=0.75 and min/max parameters"""
    return min(maximum, max(minimum, int(expon(0.75).rvs())))

def selfSim(seqs):
    return selfSimilarity('asdf', seqs, hmmfile, False)

def generateIQTree():

    sd = 1 #startingDomains

    hostTree = createRandomTopology(1, 1, lambda x: x)
    guestTree, nodeMap = buildGuestTree(hostTree, s2, expfunc, .2, gaussNoise, sd)

    print hostTree.get_ascii()

    rootSequence = grs(sd)
    evolveAlongTree(hostTree, guestTree, nodeMap, rootSequence, hmmfile, emissionProbs)

    seqs = findDomains(hostTree.sequence, hmmfile)[2] #pylint: disable=no-member
    names = [(leaf.position, leaf.name) for leaf in guestTree if leaf.event != 'LOSS']
    names.sort()
    names = [name[1] for name in names]
    guestTree = prune(guestTree, names)
    outgroup = Tree()
    outgroup.up = guestTree
    guestTree.children.append(outgroup)
    outgroup.name = 'Outgroup'
    outseq = evolveSequence(rootSequence, .1, 2, emissionProbs, hmmfile)
    outseq = findDomains(outseq, hmmfile)[2][0]
    outgroup.add_feature('sequence', outseq)
    seqs.insert(0, outseq)
    names.insert(0, 'Outgroup')

    guestTree.write(outfile = 'testtree.nwk')

    writeFasta(names, seqs, 'testfasta.fa', False)
    mlTree('testfasta.fa', 'testtree.nwk', True)
    iqtree = Tree('testfasta.fa.treefile')
    iqtree.set_outgroup(iqtree&('Outgroup'))

    writeFasta(names, seqs, 'testfasta2.fa', False)
    mlTree2('testfasta2.fa', True)
    iqtree2 = Tree('testfasta2.fa.treefile')
    iqtree2.set_outgroup(iqtree2&('Outgroup'))

    print guestTree
    print iqtree
    print iqtree2

    return hostTree, guestTree, iqtree, iqtree2

def fivenumberstatistic(sim):
    flat = []
    for i in range(len(sim)):
        for j in range(i+1, len(sim)):
            flat.append(sim[i][j])

    flat.sort()

    fq = int(len(flat) / 4)
    median = int(len(flat) / 2)
    tq = int(3 * len(flat) / 4)

    return flat[0], flat[fq], flat[median], flat[tq], flat[-1]

failcounter = 0
f = open('output.txt','w')
f.write('RF\tMaxRF\tMin\tQ1\tMedian\tQ3\tMax\n')
i = 0
while i < 200:
    try:
        hostTree, guestTree, iqtree, iqtree2 = generateIQTree()
        system('rm testfasta*')
        rf, maxrf = guestTree.robinson_foulds(iqtree)[:2]
        sim = selfSim(hostTree.sequence) #pylint: disable=no-member
        out = str(rf) + '\t' + str(maxrf) + '\t'
        fns = fivenumberstatistic(sim)
        for stat in fns:
            out += str(stat) + '\t'
        f.write(out + '\n')
        i += 1
    except:
        failcounter += 1
        continue

print '\nFAILCOUNTER WAS', failcounter
f.close()

"""
selfSim(hostTree.sequence) #pylint: disable=no-member

g = groupDomains(names, seqs, hmmfile)
val = 0 #randint(0,9)

selfSimilarity(names[val], seqs[val], hmmfile, True)
print columnSimilarity(seqs, hmmfile)


writeFasta(names, seqs, 'testfasta.fa')
guestTree.write(outfile='testtree.nwk', format=1)
"""
