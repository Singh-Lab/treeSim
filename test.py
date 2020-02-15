import sys
from os import system, listdir
from ete3 import Tree
from HostTreeGen import createRandomTopology
from OrthoAnalysis import selfSimilarity
from OrthoNJ import prune, mlTree, mlTree2, recScore
from GuestTreeGen import buildGuestTree, exp, expon
from CreateSequences import evolveAlongTree, evolveSequence
from stats import gaussNoise
from rootSequence import genRandomSequence2 as grs
from TreeUtils import findDomains, raxml, raxml_score, raxml_score_from_file
from TreeUtils import writeMapping, writeFasta, writeTree, genMap
from random import randint
import pickle
import numpy as np
from Utils import printProgressBar
from datetime import datetime
from ConfigParser import ConfigParser as CP

#eppath = '/home/caluru/Documents/shilpa/treeSimulation/simulation/zf_shilpa_probabilities.pickle'
eppath = CP.EP_PATH #pylint: disable=no-member
transmat = pickle.load(open(CP.TRANSMAT)) #pylint: disable=no-member
emissionProbs = pickle.load(open(eppath))
hmmfile = CP.HMM_PATH #pylint: disable=no-member

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
    evolveAlongTree(hostTree, guestTree, nodeMap, rootSequence, hmmfile, emissionProbs, transmat)

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
    outseq = evolveSequence(rootSequence, .1, 2, emissionProbs, hmmfile, transmat)
    outseq = findDomains(outseq, hmmfile)[2][0]
    outgroup.add_feature('sequence', outseq)
    seqs.insert(0, outseq)
    names.insert(0, 'Outgroup')

    guestTree.write(outfile = 'testtree.nwk')
    hostTree.write(outfile='hosttree.nwk')
    addRandomTrees('testtree.nwk')

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
            iqtree.write(outfile = treefilepath + str(i) + '.iqtree')

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

def withHost(numLeaves = 4, bl = .5, hostTree = None):
    sd = 1 #startingDomains
    extralen = .05

    if hostTree is None:
        hostTree = createRandomTopology(numLeaves, bl, lambda x: x)
        for leaf in hostTree:
            leaf.dist += extralen

    guestTree, nodeMap = buildGuestTree(hostTree, s2, expfunc, .1, gaussNoise, sd)
    
    for leaf in guestTree:
        leaf.dist += extralen

    rootSequence = grs(sd)
    evolveAlongTree(hostTree, guestTree, nodeMap, rootSequence, hmmfile, emissionProbs, transmat)
    names = [(leaf.position, leaf.name) for leaf in guestTree if leaf.event != 'LOSS']
    names.sort()
    names = [i[1] for i in names]
    names.sort()

    seqs = []
    hnodes = sorted([i.name for i in hostTree]) 
    for node in hnodes:
        seqs += findDomains((hostTree&node).sequence, hmmfile)[2]

    for node in hostTree.traverse():
        node.del_feature('leaves')

    guestTree = guestTree.children[0]
    guestTree.up = None
    writeTree(hostTree, 'host.nwk')
    writeTree(guestTree, 'guest.nwk')
    writeFasta(names, seqs, 'sequences.fa')

    return hostTree, guestTree, names, seqs

def generateTestCases(n=500):

    hostCases = 0
    while hostCases < n / 10:

        try:
            host = withHost(8, .3)[0]
        except:
            continue

        guestCases = 0
        while guestCases < 10:

            printProgressBar(hostCases * 10 + guestCases, n)

            try:
                guest = withHost(8, .3, host)[1]
            except:
                 continue

            writeMapping(genMap(host, guest), 'guest.map')
            folder_name = 'examples/' + str(hostCases*10 + guestCases) + '/'
            system('mkdir ' + folder_name)
            system('mv host.nwk guest.nwk sequences.fa guest.map ' + folder_name)

            guestCases += 1

        hostCases += 1

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

def addRandomTrees(treefile, n=1):
    #Adds n trees to the given treefile with the same leaf names but a random topology
    t = Tree(treefile)
    outtrees = []
    #names.discard('Outgroup')
    for _ in range(n):
        names = set([leaf.name for leaf in t])
        g = Tree()
        g.populate(len(names))
        for leaf in g:
            leaf.name = names.pop()
        """
        out = Tree()
        out.populate(2)
        out.children[0].name = 'Outgroup'
        out.children[1] = g
        g.up = out
        """
        #outtree = out.write(format=9)
        outtrees.append(g.write(format=9))
    f = open(treefile, 'a')
    for tree in outtrees:
        f.write(tree + '\n')
    f.close()

def createRandomTrees(tree, n=1, outfile=None):
    """
    Takes a tree as input and generates n trees with the same leaf names but random topologies
    """
    outtrees = []

    def name(tree):
        i = 100
        for node in tree.traverse():
            if node.name == '':
                node.name = str(i)
                i += 1

    for _ in range(n):
        names = set([leaf.name for leaf in tree])
        g = Tree()
        g.populate(len(names))
        for leaf in g:
            leaf.name = names.pop()
        name(g)
        outtrees.append(g)
    
    if outfile != None:
        f = open(outfile,'w')
        for tree in outtrees:
            f.write(tree.write(format=1) + '\n')
        f.close()
    
    return outtrees

def testLikelihood(checkpoint=False):
    """
    Parameter sweep over number of leaves in the host tree vs length of the tree 
    to see which combinations produce domain trees that are not significantly worse
    than RAxML_bestTree
    """

    from TreeSearch import reroot
    from matplotlib import pyplot as plt
    import seaborn as sns
    sns.set()

    def genMap(host, guest):
        #{guest -> host}
        nodemap = {}
        for leaf in guest:
            hname = 'h' + leaf.name.split("_")[0][1:]
            nodemap[leaf] = host&hname
        return nodemap

    def name(tree):
        i = 100
        for node in tree.traverse():
            if node.name == '':
                node.name = str(i)
                i += 1

    nLeaves = []
    bLengths = []
    colors = [] 
    available_colors = ["#440154FF", "#404788FF", "#287D8EFF", "#29AF7FFF", "#95D840FF", "#FDE725FF"]

    if not checkpoint:
        for nl in [4,6,8,10,12]:
            for bl in np.arange(.05, .55, .05):

                #it = ((nl - 4) + 1) / 2 * 10 + bl / .5 * 5
                #printProgressBar(it, 10)

                successes = 0

                for _ in range(5):

                    host = None
                    while host is None:
                        try:
                            host, guest, names, seqs = withHost(nl, bl)
                        except:
                            continue

                    guest = guest.children[0]
                    guest.up = None
                    name(guest)
                    name(host)

                    f = open('sequences.fa','w')
                    for i in range(len(names)):
                        f.write(">" + names[i] + '\n')
                        f.write(seqs[i] + '\n')
                    f.close()

                    f = open('host.nwk','w')
                    f.write(host.write(format=1) + '\n')
                    f.close()

                    f = open('guest.nwk','w')
                    f.write(guest.write(format=1) + '\n')
                    f.close()

                    #Run RAxML
                    if "RAxML_bestTree.nwk" in listdir('.'):
                        system('rm RAxML*')
                    raxml('sequences.fa', 'nwk')
                    raxml_tree = Tree('RAxML_bestTree.nwk')
                    name(raxml_tree)
                    raxml_tree = reroot(host, raxml_tree, genMap(host, raxml_tree))

                    result = raxml_score_from_file('RAxml_bestTree.nwk', 'guest.nwk', 'sequences.fa')
                    score = result[1][0]
                    if score == 0:
                        successes += 1

                print (nl, bl), successes

                colors.append(available_colors[successes])
                nLeaves.append(nl)
                bLengths.append(bl)
        
        pickle.dump((nLeaves, bLengths, colors), open('plot_data.pickle','w'))

    else:
        nLeaves, bLengths, colors = pickle.load(open('plot_data.pickle'))

    plt.scatter(nLeaves, bLengths, color=colors)
    plt.title('Acceptable BL Ranges for Each Host Tree Size')
    plt.xlabel('# of Host Leaves')
    plt.ylabel('Branch Length')
    plt.show()

def treeSearchTest(n=100):
    """
    Tests the effectiveness of the tree search algo in TreeSearch. Produces a plot
    of rec score x raxml score for the following trees:
    1) The optimal tree found by raxml
    2) The actual tree generated by GuestTreeGen
    3) n randomly generated trees with the same leaves
    4) n trees one spr move away from the raxml optimum tree
    """

    from TreeSearch import reroot, pick_sprs, reconcileDL
    from matplotlib import pyplot as plt
    import seaborn as sns
    sns.set()

    def genMap(host, guest):
        #{guest -> host}
        nodemap = {}
        for leaf in guest:
            hname = 'h' + leaf.name.split("_")[0][1:]
            nodemap[leaf] = host&hname
        return nodemap

    def name(tree):
        i = 100
        for node in tree.traverse():
            if node.name == '':
                node.name = str(i)
                i += 1

    #Generate host and guest tree
    print "Generating Trees"
    host, guest, names, seqs = withHost(8, .3) #Experimentally determined from testLikelihood()
    guest = guest.children[0]
    guest.up = None
    name(guest)
    name(host)

    f = open('sequences.fa','w')
    for i in range(len(names)):
        f.write(">" + names[i] + '\n')
        f.write(seqs[i] + '\n')
    f.close()

    f = open('host.nwk','w')
    f.write(host.write(format=1) + '\n')
    f.close()

    f = open('guest.nwk','w')
    f.write(guest.write(format=1) + '\n')
    f.close()

    #Run RAxML
    print "Running RAxML"
    if "RAxML_bestTree.nwk" in listdir('.'):
        system('rm RAxML*')
    raxml('sequences.fa', 'nwk')
    raxml_tree = Tree('RAxML_bestTree.nwk')
    name(raxml_tree)
    raxml_tree = reroot(host, raxml_tree, genMap(host, raxml_tree))

    #Get ML score of RAxML tree
    f = list(open('RAxML_info.nwk'))
    raxml_mlscore = 0
    for line in f:
        if "Final GAMMA-based Score of best tree" in line:
            raxml_mlscore = float(line.strip().split()[-1])

    #Check if guest tree is significantly worse than raxml tree
    print "Evaluating Guest Tree"
    result = raxml_score_from_file('RAxml_bestTree.nwk', 'guest.nwk', 'sequences.fa')
    guest_score = result[0][0]
    score = result[1][0]
    if score == 0:
        print "Guest Tree is Not Significantly Worse than RAxML Tree"
    else:
        print "Guest Tree is Significantly Worse than RAxML Tree"

    #Generate/Score n random trees
    print "Random Trees"
    random = createRandomTrees(guest, n, 'badTrees.nwk')
    random_scores = raxml_score(raxml_tree, random, 'sequences.fa')[0]
    random = [reroot(host, g, genMap(host, g)) for g in random]
    random_rscores = [reconcileDL(host, g, genMap(host, g))[0] for g in random]

    #Generate/Score n 1spr moves
    print "1SPR Trees"
    one_spr = pick_sprs(raxml_tree, n)
    ospr_scores = raxml_score(raxml_tree, one_spr, 'sequences.fa')[0]
    one_spr = [reroot(host, i, genMap(host, i)) for i in one_spr]
    ospr_rscores = [reconcileDL(host, g, genMap(host, g))[0] for g in one_spr]

    #Score RAxML Tree/Guest Tree
    print "Finishing Up 12 (???)"
    raxml_rscore = reconcileDL(host, raxml_tree, genMap(host, raxml_tree))[0]
    guest_rscore = reconcileDL(host, guest, genMap(host, guest))[0]
    
    #Plot
    plt.figure()
    a = plt.scatter(random_rscores, random_scores, color='blue')
    b = plt.scatter([guest_rscore], [guest_score], color='green')
    c = plt.scatter(ospr_rscores, ospr_scores, color='purple')
    d = plt.scatter([raxml_rscore], [raxml_mlscore], color='red')
    plt.xlabel('Reconciliation Score')
    plt.ylabel('RAxML Score')
    plt.legend((a,b,c,d),('random','real','1SPR','RAxML best'))
    plt.show()

def emMatTest(bl=1):
    from CreateSequences import genTransitionMatrix
    from matplotlib import pyplot as plt
    from scipy.linalg import expm
    import seaborn as sns
    sns.set()
    
    eps = pickle.load(open(CP.EP_PATH)) #pylint:disable=no-member
    transmat = pickle.load(open(CP.TRANSMAT)) #pylint:disable=no-member

    outs = genTransitionMatrix(eps, transmat, bl)
    transmat = expm(transmat * bl)

    for i in range(20):
        transmat[i][i] /= 10
    transmat = [i / sum(i) for i in transmat]

    i = 0
    for (ep, out) in zip(eps, outs):
        plt.figure(figsize=(16,4))
        plt.subplot(1, 3, 1)
        sns.heatmap(transmat, cmap='viridis')
        plt.title('Original Transition Matrix')

        plt.subplot(1, 3, 2)
        sns.heatmap([ep], cmap='viridis')
        plt.title('Mask (Position ' + str(i) + ")")

        plt.subplot(1, 3, 3)
        sns.heatmap(out, cmap='viridis')
        plt.title('Output')

        plt.show()
        i += 1

def seqGenTest(n=500, bl=1):

    failcount = 0

    for i in range(n):
        seq = grs(5)
        try:
            evolveSequence(seq, .1, bl, emissionProbs, hmmfile, transmat)
        except ValueError:
            failcount += 1
        suff = str(failcount) + ' / ' + str(i+1) + ' failures'
        printProgressBar(i+1, n, suffix=suff)

def seqDiff(n=10, bl=1):

    RED = '\033[91m'
    NORMAL = '\033[0m'

    seq = grs(1)
    dom = findDomains(seq, hmmfile)[2][0]
    print dom

    iterations = 0
    while iterations < n:
        try:
            temp = evolveSequence(seq, .05, bl, emissionProbs, hmmfile, transmat)
            tempdom = findDomains(temp, hmmfile)[2][0]
            out = ""
            nMuts = 0
            for i in range(len(dom)):
                if tempdom[i] == dom[i]:
                    out += tempdom[i]
                else:
                    out += RED + tempdom[i] + NORMAL
                    nMuts += 1
            totalMuts = 0
            for i in range(len(temp)):
                if temp[i] != seq[i]:
                    totalMuts += 1
            print out, nMuts, totalMuts, round(totalMuts / float(nMuts) / (len(temp) / 23.))
        except:
            continue
        iterations += 1


if __name__ == "__main__":
    print datetime.now().strftime('%Y-%m-%d %H:%M:%S'), '\n'
    host, guest, names, seqs = withHost()
    #treeSearchTest()
    #emMatTest()
    #for bl in [.1, .25, .5, .75, 1]:
        #seqGenTest(100, bl)
    #seqDiff(bl=.5)
    #testLikelihood()
    #generateTestCases(10)
    print '\n', datetime.now().strftime('%Y-%m-%d %H:%M:%S')
