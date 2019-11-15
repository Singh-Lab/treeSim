#Provides several utility functions useful for reading/writing trees, etc
from ete3 import Tree
from math import exp
from Utils import sortBy
import numpy as np
import os

#Tree I/O 
def readTree(filename):
    """
    Reads a tree in newick or NHX format to an ete3 object. To read trees not created by 
    this package, see readHostTree in HostTreeGen.py
    """
    return Tree(filename, format=1)

def writeTree(tree, filename):
    """Writes an ete3 tree to the given filename in NHX format"""
    output = tree.write(format=1, features=[])[:-1]

    #For whatever reason ete3 doesn't include root name + features, so this adds it on
    output += tree.name + "[&&NHX"
    for feature in tree.features:
        output += ":" + str(feature) + '=' + str(getattr(tree, feature))
    output += '];'

    outputHandle = open(filename, 'w')
    outputHandle.write(output)
    outputHandle.close()

def readMapping(host, guest, mapfile):
    """
    Requires that the host and guest tree that the mapping refers to have already 
    been read into memory. See writeMapping for mapfile format

    Args:
        host    (Tree): The host tree
        guest   (Tree): The guest tree
        mapfile (str ): Name of the file containing mapping between host and guest nodes

    Output:
        nodeMap (dict): A mapping of host -> [guest] nodes
    """

    nodemap = {}
    for node in host.traverse():
        nodemap[node] = []

    nodemapFile = list(open(mapfile))

    for line in nodemapFile:
        line = line.split()
        hostNode = host&line[0]
        mapped = []
        for guestName in line[1:]:
            guestNode = guest&guestName
            mapped.append(guestNode)
        nodemap[hostNode] = mapped

    return nodemap

def writeMapping(nodemap, filename):
    """
    Writes out a mapping between host nodes and guest nodes. Each line of the output
    consists of a host node name followed by the names of each of the guest nodes 
    that maps to it (all separated by spaces).
    """
    outputHandle = open(filename, 'w')

    for node in nodemap:
        out = node.name + ' ' + ' '.join([i.name for i in nodemap[node]])
        outputHandle.write(out + '\n')

    outputHandle.close()

def writeFasta(names, sequences, filename, outgroup=False):
    """
    Writes out a fasta sequence to a file named <filename>. If outgroup is true,
    creates a fake outgroup (length 23 "AAAAAAAAAAAAAAAAAAAAAAA") and adds as the 
    last sequence in the fasta file
    """
    f = open(filename, 'w')

    #TODO: Make sure that this is actually a valid outgroup in all cases (should be)
    if outgroup:
        alphabet = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', \
                        'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
        fakeseq = ''.join(np.random.choice(alphabet, 23))
        f.write(">Outgroup\n" + fakeseq + '\n')

    for i in range(len(names)):
        f.write(">" + names[i] + '\n')
        f.write(sequences[i] + '\n')

    f.close()

#Duplication rate functions for use with guest tree generation
def s(x):
    """Sigmoid function designed to quickly reduce losses as domains are lost"""
    denom = 1 + exp(8-x)
    return .7 - .3 / denom

def s2(x):
    denom = 1 + exp(10-x)
    return .7 - .6 / denom

def s3(x):
    if x >= 13:
        return 0
    if x <= 7:
        return 1
    return s2(x)

#Utilities for finding domains

#Finds all domains in infile according to hmm in hmmfile using hmmsearch232
def findDomainsFile(infile, hmmfile):
    #domName = os.popen("grep 'NAME' " + hmmfile).read().split()[1]
    seqName = list(open(infile))[0][1:].strip()
    os.system("hmmsearch232 --domE 0.00001 " + hmmfile + " " + infile 
                + " | grep '^ [ ]*" + seqName + "' > tmp/grepoutput.txt")
    hits = list(open('tmp/grepoutput.txt'))
    starts, ends, seqs = [], [], []

    for line in hits:
        line = line.split()
        if not line[1].isdigit():
            continue
        starts.append(int(line[1]) - 1) #HMMER is 1-indexed :(
        ends.append(int(line[3]) - 1)
        seqs.append(line[2].upper())

    #HMMER doesn't sort by position in the sequence. This fixes that
    seqs = sortBy(seqs, starts)
    starts.sort()
    ends.sort()
        
    return starts, ends, seqs

#Finds all domains in the input sequence according to hmm in hmmfile using hmmsearch232	
def findDomains(sequence, hmmfile):
    g = open('tmp/tmp.fa','w')
    g.write('>seq\n' + sequence)
    g.close()
    
    return findDomainsFile('tmp/tmp.fa', hmmfile)

def printDomSeq(sequence, hmmfile, minimal_rep = False):
    """
    prints the sequence with domains highlighted in red 
    (first character highlighted in green)
    """

    #Escape sequences used to colorize terminal output
    RED = '\033[91m'
    GREEN = '\033[92m'
    NORMAL = '\033[0m'

    #Find domains, check if sequence begins and/or ends with a domain
    domains = findDomains(sequence, hmmfile)[2]

    #split on all domains
    for domain in domains:
        sequence = sequence.replace(domain, "xxx")
    sequences = sequence.split("xxx")

    if minimal_rep:
        out = ''
        for i in range(len(domains)):
            out += '---' + RED + "XXX" + NORMAL
        out += '---'

        print out
        return
    
    #Reassemble full sequence post evolution
    out = ''
    for i in range(len(domains)):
        out += sequences[i] + GREEN + domains[i][0] + RED + domains[i][1:] + NORMAL
    out += sequences[-1] #if len(sequences) > len(domains) else RED + domains[-1] + NORMAL

    print out

def isValid(domain):
    """Checks if the input string is a valid zf-C2H2 domain"""
    valid = len(domain) == 23 and domain[2] == "C" and domain[5] == "C"
    valid &= domain[18] == "H" and domain[22] == "H"
    return valid
