#Builds a neighbor joining-like tree where starting subtrees are taken from columns
#and a full tree is constructed from this. Intended as input for iqtree.

from ete3 import Tree
from OrthoAnalysis import groupDomains
from TreeUtils import readTree
from ConfigParser import ConfigParser as CP
import numpy as np
import os

#TODO: Untested
def columnOccupancy(column):
    #Checks percentage of entries with a domain in a column
    count = 0.

    for entry in column:
        if entry != '':
            count += 1

    return count / len(column)

#TODO: Untested
def prune(tree, domNames):
    """
    Prunes away leaves not included in domNames, as well as any other nodes with one
    child in order to maintain the complete binary tree property.
    """
    #Mark
    delete = []
    for leaf in tree:
        if leaf.name not in domNames:
            delete.append(leaf)

    #Sweep
    while delete != []:
        leaf = delete.pop(0)
        if leaf.children == []:
            parent = leaf.up
            if parent not in delete:
                delete.append(parent)
            parent.children.remove(leaf)
            leaf.up = None
        else: #only one child
            if leaf == tree:
                tree = leaf.children[0]
                tree.up = None
                continue
            leaf.up.children.remove(leaf)
            leaf.up.children.append(leaf.children[0])
            leaf.children[0].up = leaf.up
            leaf.up = None
            leaf.children = []

    return tree

#TODO: Untested
def fullColumn(hostcopy, domNames, mapping=None):
    """
    Creates the guest subtree topology for a column of domains with high occupancy. 
    Takes a copy of the host topology, renames the leaves to match the domain names, 
    and returns this subtree

    Args:
        hostcopy (Tree): A copy of the host tree topology which retains branch lengths
                         and host node names. 
        domNames (list): contains the names of the domains in this column 
        mapping  (dict): A mapping from domain name -> name of host it came from. If 
                         no mapping is specified, assumes that domain with name G<n>_XXX 
                         comes from host with name H<n>.

    Returns:
        hostcopy (Tree): The modified copy of the host tree.
    """

    #Step 1: Generate mapping (if necessary)
    domNames = [name for name in domNames if name != '']

    if mapping is None:
        mapping = {}
        for name in domNames:
            #Assumes G<n>_XXX format
            hostName = "H" + name[1:].split("_")[0]
            mapping[name] = hostName

    #Step 2: Rename leaves appropriately
    for name in domNames:
        hostName = mapping[name]
        hostNode = hostcopy&hostName
        hostNode.name = name

    #Step 3: Prune away unecessary branches/nodes
    return prune(hostcopy, domNames)

#TODO: Untested
def emptyColumn():
    """
    For columns with low occupancy (<20%). Tries to group them into the minimum number
    of high occupancy clades within the host tree.
    """
    pass

def splitByClade(hostCopy, domNames, threshold=0.8, mapping=None):
    """
    Works on all cases regardless of column occupancy. Takes a host tree and maps the 
    domains in a column onto the tree. After, returns the smallest number of clades 
    such that each clade has an occupancy at least as high as the threshold.
    """
    #Step 1: Generate mapping (if necessary)
    domNames = [name for name in domNames if name != '']

    if mapping is None:
        mapping = {}
        for name in domNames:
            #Assumes G_<n>_XXX format
            hostName = "H_" + name.split("_")[1]
            mapping[name] = hostName

    #Step 2: Rename leaves appropriately
    for name in domNames:
        hostName = mapping[name]
        hostNode = hostCopy&hostName
        hostNode.name = name

    #Step 3: Count occupancy by clade
    for leaf in hostCopy:
        leaf.add_feature('clade_size', 1)
        if leaf.name in domNames:
            leaf.add_feature('clade_occupancy', 1)
        else:
            leaf.add_feature('clade_occupancy', 0)

    jobs = [node for node in hostCopy.traverse()][::-1]
    for node in jobs:
        if node.children == []:
            continue
        a, b = node.children
        node.add_feature('clade_size', a.clade_size + b.clade_size)
        node.add_feature('clade_occupancy', a.clade_occupancy + b.clade_occupancy)

    #Step 4: Break off subtrees
    subtrees = []
    jobs = [hostCopy]
    while jobs != []:
        node = jobs.pop()
        if float(node.clade_occupancy) / node.clade_size >= threshold:
            node.up = None
            subtrees.append(node)
        else:
            jobs += node.children

    return [prune(subtree, domNames) for subtree in subtrees]

#TODO: Untested
def createOrthoTree(hostTree, names, sequences, filename):
    """
    Given a host tree and a guest tree, creates a first pass orthogroup - based
    reconstruction of a likely domain tree.

    Args:
        hostTree (str ): The file name of the host tree to construct the guest inside of
        names    (list): The names of each sequence
        sequences(list): Sequences from a orthogroup to build domain tree from
        filename (str ): The name of  the file to be written out
    """

    host = readTree(hostTree)
    hmmfile = CP.HMM_FILE #pylint: disable=no-member
    
    grouped, domNames = groupDomains(names, sequences, hmmfile)
    grouped = np.asarray(grouped).T
    domNames = np.asarray(domNames).T

    subtrees = []

    for i in range(len(grouped)):
        #column = grouped[i]
        names = domNames[i]
        subtrees += splitByClade(host.copy('newick'), domNames)

        #TODO: Finish constructing tree from clades

def mlTree(msafile, treefile, outgroup=False):
    """Runs IQTree on the input file given ("""
    cmd = "iqtree -s " + msafile + " -z " + treefile #TODO: Change -z back to -t!
    if outgroup: cmd += " -o Outgroup"
    os.system(cmd + " > thingy.txt")

def mlTree2(msafile, outgroup=False):
    """Runs IQTree on the input file given ("""
    cmd = "iqtree -s " + msafile + " -pre " + msafile
    if outgroup: cmd += " -o Outgroup"
    os.system(cmd + " > /dev/null 2>&1")

def recScore(hostfile, guestfile):
    """Runs Notung and returns the reconciliation score with D = 2, L = 1"""
    def createRecTree(treefile):
        #Takes in a gene tree and converts it into notung format. Also removes the outgroup
        t = Tree(treefile)
        #Remove outgroup
        if t.children[0].name == 'Outgroup':
            t = t.children[1]
        elif t.children[1].name == 'Outgroup':
            t = t.children[0]
        #Correct naming scheme
        for leaf in t:
            name = leaf.name.split("_")
            leaf.name = name[1] + "_h" + name[0][1:]
        return t
    t = createRecTree(guestfile)
    t.write(outfile='guest.nwk')
    notung = '/home/caluru/Downloads/Notung-2.9/Notung-2.9.jar --reconcile --parsable'
    cmd = 'java -jar ' + notung + ' -s ' + hostfile + ' -g guest.nwk'
    os.system(cmd + "> notungoutput.txt")
    f = list(open('guest.nwk.reconciled.parsable.txt'))[0]
    score = f.split()[0]
    os.system('rm guest.nwk.reconciled guest.nwk.reconciled.parsable.txt guest.nwk')
    return score

def logLikelihood(msafile, treefile):
    #Takes in an msa and a tree topology and uses IQTree to generate a log likelihood
    cmd = "iqtree -s " + msafile + " -z " + treefile
    os.system(cmd)

if __name__ == '__main__':
    #Test Pruning
    if False:
        t = Tree()
        t.populate(10)

        print t.get_ascii()
        
        keep = [leaf.name for leaf in t]
        np.random.shuffle(keep)
        keep = keep[:7]

        #keep = ['A','B','C','D','G','H','J']
        
        i = 1
        for node in t.traverse():
            if node.name == '':
                node.name = str(i)
                i += 1
        
        print keep
        print prune(t, keep)

    #Test split by clade
    if False:
        t = Tree('test.nwk')

        i = 1
        for node in t.traverse():
            if node.name == '':
                node.name = "H_" + str(i)
                i += 1
            else:
                node.name = "H_" + node.name

        print t.get_ascii()
        keep = ['A','B','C','D','G','H','I','N','O','P']
        keep = ["G_" + i for i in keep]
        print keep
        stuff = splitByClade(t, keep, 10/16.)
        print [thing.name for thing in stuff]

    if True:
        mlTree('testfasta.fa', 'testtree.nwk')
