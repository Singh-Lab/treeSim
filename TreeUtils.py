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
    """Tighter sigmoid than s(x)"""
    denom = 1 + exp(10-x)
    return .7 - .6 / denom

def s3(x):
    """s2(x) but with strict cutoffs at 7 and 13 domains"""
    if x >= 13:
        return 0
    if x <= 7:
        return 1
    return s2(x)

#Utilities for finding domains

def findDomainsFile(infile, hmmfile):
    """
    Finds all domains in infile according to hmm in hmmfile using hmmsearch232.
    Assumes only one sequence is in the infile

    Args:
    infile  (str): The path to the fasta file to search for domains in
    hmmfile (str): Path to file containing hmm to search for

    Output:
    starts (list): List of starting positions of domains found
    ends   (list): List of ending positions of domains found
    seqs   (list): List of domain sequences found 
    """
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
    """
    Finds domains in the given sequence

    Args:
    sequence (str): sequence to search for domains in
    hmmfile (str): Path to file containing hmm to search for

    Output:
    starts (list): List of starting positions of domains found
    ends   (list): List of ending positions of domains found
    seqs   (list): List of domain sequences found 
    """
    g = open('tmp/tmp.fa','w')
    g.write('>seq\n' + sequence)
    g.close()
    
    return findDomainsFile('tmp/tmp.fa', hmmfile)

def printDomSeq(sequence, hmmfile, minimal_rep = False):
    """
    prints the sequence with domains highlighted in red 
    (first character highlighted in green)

    Args:
    sequence (str): Protein sequence to print
    hmmfile  (str): Path to hmm file of domain to highlight
    mimimal_rep (bool): If true, prints a string of dashes and X's (nondomain and domain
                        sequences) rather than the full highlighted sequences
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

def raxml(infile, outext):
    """
    Runs raxml with set parameters on the input 

    Args:
    infile (str): fasta file to build tree from
    outext (str): output file extension to use. Output tree will be at
                  RAxML_bestTree.<outext>
    """
    command = '/home/caluru/Downloads/standard-RAxML-master/raxmlHPC-PTHREADS-AVX2 -s '
    command += infile + ' -n ' + outext + ' -m PROTGAMMAJTT -T 8 -p ' + str(np.random.randint(2000))
    command += ' > raxml_log.txt'
    os.system(command)

def raxml_score(benchfile, testfile, seqfile):
    """
    Uses RAxML to compute the SH score between the benchmark tree constructed by RAxML and
    a set of other trees. Outputs a binary vector where the ith entry is 1 if tree i is 
    significantly worse than the benchmark and 0 otherwise

    Args:
    benchfile (str): path to newick file containing the best tree found by RAxML
    testfile  (str): path to newick file containing all trees to test, one per line
    seqfile   (str): path to fasta file containing leaf sequences

    Output:
    scores (list): A list of 0/1 entries specifying whether or not each tree is
                   significantly worse than the benchmark or not.
    """
    #Run RAxML to find if tree in benchfile is significantly better than those in testfile
    command = '/home/caluru/Downloads/standard-RAxML-master/raxmlHPC-PTHREADS-AVX2 '
    #Switch to -f h if this takes too long
    command += '-f H -t' + benchfile + ' -z ' + testfile + ' -s ' + seqfile + ' -m PROTGAMMAJTT -T 8 -n sco'
    command += ' > raxml_log.txt' 
    #TODO: Read results and select tree
    os.system(command)

    #Parse resulting logfile
    f = list(open('RAxML_info.sco'))
    start = 0
    for i in range(len(f)):
        if ' trees in File ' in f[i]:
            start = i

    f = f[start+3:]
    scores = []
    for line in f:
        answer = line.split('Significantly Worse: ')[1].split()[0]
        scores.append(1 if answer == 'Yes' else 0)

    return scores
