import sys
from os import system
#sys.path.append('treeSim/')
from ete3 import Tree
from HostTreeGen import createRandomTopology
from OrthoAnalysis import selfSimilarity, groupDomains, columnSimilarity
from OrthoNJ import prune, mlTree, mlTree2, recScore
from GuestTreeGen import buildGuestTree, exp, expon
from CreateSequences import evolveAlongTree, evolveSequence
from stats import gaussNoise
from rootSequence import genRandomSequence2 as grs
from TreeUtils import writeTree, writeFasta, findDomains
from random import randint
import pickle
import numpy as np
from Utils import printProgressBar
from datetime import datetime

eppath = '/home/caluru/Documents/shilpa/treeSimulation/simulation/zf_shilpa_probabilities.pickle'
emissionProbs = pickle.load(open(eppath))
hmmfile = '/Users/chaitanya/Data/hmmfiles/zf_shilpa_232.hmm'

def s2(x):
    denom = 1 + exp(10-x) if x < 10 else 1
    return 1 - .5 / denom

def expfunc(minimum=1, maximum=3):
    """Exponential distribution with lambda=0.75 and min/max parameters"""
    return min(maximum, max(minimum, int(expon(0.75).rvs())))

def selfSim(seqs):
    return selfSimilarity('asdf', seqs, hmmfile, False)

def treeskew(tree):
    #outputs height of tree, trying to measure imbalance in the tree
    farthest = tree.get_farthest_leaf(topology_only=True)
    return str(farthest[1])

def findLeaves(nodes):
    leaves = []
    for node in nodes:
        if node.children == []:
            leaves.append(node)
    return leaves

def generateIQTree():

    sd = 1 #startingDomains

    hostTree = createRandomTopology(1, 1, lambda x: x)
    guestTree, nodeMap = buildGuestTree(hostTree, s2, expfunc, .2, gaussNoise, sd)

    rootSequence = grs(sd)
    evolveAlongTree(hostTree, guestTree, nodeMap, rootSequence, hmmfile, emissionProbs)

    names, seqs = [], []

    for node in hostTree:
        seqs += findDomains(node.sequence, hmmfile)[2]
        gnodes = findLeaves(nodeMap[node])
        n = [(leaf.position, leaf.name) for leaf in gnodes if leaf.event != 'LOSS']
        n.sort()
        names += [name[1] for name in n]

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
    hostTree.write(outfile='hosttree.nwk')
    addRandomTree('testtree.nwk')

    writeFasta(names, seqs, 'testfasta.fa', False)
    mlTree('testfasta.fa', 'testtree.nwk', True)
    iqtree = Tree('testfasta.fa.treefile')
    iqtree.set_outgroup(iqtree&('Outgroup'))

    return hostTree, guestTree, iqtree

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

def noHost(samples):
    failcounter = 0
    f = open('output.txt','w')
    f.write('Filename\tRF\tMaxRF\tPct RF\tMin\tQ1\tMedian\tQ3\tMax\tReal Tree Height \
                \tIQTree Height\tIQTree Score\tReal Score\tRandom Score \t')
                #\tIQTree Rec Score\tReal Rec Score\n')
    #write out actual and reconstructed tree as well as host sequence
    #Filenames: 1.tree, 1.iqtree, 1.seq
    treefilepath = 'output_trees/'
    i = 0
    while i < samples:
        printProgressBar(i+1, samples, suffix=str(failcounter))

        try:
            #Generate Data
            hostTree, guestTree, iqtree = generateIQTree()
            guestTree.write(outfile = treefilepath + str(i) + '.tree')
            iqtree.write(outfile= treefilepath + str(i) + '.iqtree')

            #Collect Stats
            iqscore, realscore, randscore = parseIQOutput('thingy.txt')
            system('rm testfasta*')
            rf, maxrf = guestTree.robinson_foulds(iqtree)[:2]
            thing = [node for node in hostTree][0]
            sim = selfSim(thing.sequence) #pylint: disable=no-member
            out = str(i) + '.txt\t' + str(rf) + '\t' + str(maxrf) + '\t'
            out += str(round(float(rf)/maxrf,2)) + '\t'
            fns = fivenumberstatistic(sim)
            for stat in fns:
                out += str(stat) + '\t'
            out += treeskew(guestTree) + '\t'
            out += treeskew(iqtree) + '\t'
            out += iqscore + '\t' + realscore + '\t' + randscore + '\t'
            #out += recScore('hosttree.nwk', 'output_trees/' + str(i) + '.iqtree') + '\t'
            #out += recScore('hosttree.nwk', 'output_trees/' + str(i) + '.tree')
            f.write(out + '\n')
            
            #Write out sequence file
            seqfile = treefilepath + str(i) + '.seq'
            g = open(seqfile, 'w')
            g.write(hostTree.sequence) #pylint: disable=no-member
            g.close()

            i += 1

        except:
            failcounter += 1

    f.close()

def withHost():
    sd = 1 #startingDomains

    hostTree = createRandomTopology(1, .5, lambda x: x)
    guestTree, nodeMap = buildGuestTree(hostTree, s2, expfunc, .1, gaussNoise, sd)

    hostTree.write(outfile='host.nwk')
    guestTree.children[0].write(outfile='guest.nwk')
    print guestTree

    rootSequence = grs(sd)
    evolveAlongTree(hostTree, guestTree, nodeMap, rootSequence, hmmfile, emissionProbs)
    names = [(leaf.position, leaf.name) for leaf in guestTree if leaf.event != 'LOSS']
    names.sort()
    names = [i[1] for i in names]
    seqs = findDomains(hostTree.sequence, hmmfile)[2] #pylint: disable=no-member
    writeFasta(names, seqs, 'sequences.fa')

def parseIQOutput(filename):
    f = list(open(filename))
    oll = 0 #index of optimal log likelihood
    test = 0 #index of real and random log likelihood

    for i in range(len(f)):
        if 'FINALIZING TREE SEARCH' in f[i]:
            f = f[i:]
            break
    
    for i in range(len(f)):
        if 'Optimal log-likelihood' in f[i]:
            oll = i
        #Assumes that the real score and the score of one random tree were evaluated
        if 'Reading trees in testtree.nwk' in f[i]:
            test = i+2

    optscore = f[oll].split()[-1]
    realscore = f[test].split()[-1]
    randscore = f[test+1].split()[-1]

    return optscore, realscore, randscore

def addRandomTree(treefile):
    #Adds a second tree to the given treefile with the same leaf names but a random topology
    t = Tree(treefile)
    names = set([leaf.name for leaf in t])
    names.discard('Outgroup')
    g = Tree()
    g.populate(len(names))
    for leaf in g:
        leaf.name = names.pop()
    out = Tree()
    out.populate(2)
    out.children[0].name = 'Outgroup'
    out.children[1] = g
    g.up = out
    outtree = out.write(format=9)
    f = open(treefile, 'a')
    f.write(outtree)
    f.close()

if __name__ == "__main__":
    print datetime.now().strftime('%Y-%m-%d %H:%M:%S'), '\n'
    withHost()
    print '\n', datetime.now().strftime('%Y-%m-%d %H:%M:%S')
