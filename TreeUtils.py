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

def writeReconciliation(host, guest, mapping):
    """Writes out a basic output version of the input trees and mapping without any
    information from the trees in nwk format, and writes a guest -> host mapping file.
    Automatically names all nodes that do not have a name.
    """
    os.system('mkdir output')

    i=20
    for node in host.traverse():
        if node.name == '':
            node.name = str(i)
            i += 1

    i=20
    for node in guest.traverse():
        if node.name == '':
            node.name = str(i)
            i += 1

    f = open('output/host.nwk','w')
    f.write(host.write(format=1)[:-1] + host.name + ';')
    f.close()

    g = open('output/guest.nwk','w')
    g.write(guest.write(format=1)[:-1] + guest.name + ';')
    g.close()

    mapfile = open('output/mapping.txt','w')
    for key in mapping:
        out = key.name + '\t' + mapping[key].name
        mapfile.write(out + '\n')
    mapfile.close()

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

def readMapping(host, guest, mapfile=None):
    """
    Requires that the host and guest tree that the mapping refers to have already 
    been read into memory. See writeMapping for mapfile format

    Args:
        host    (Tree): The host tree
        guest   (Tree): The guest tree
        mapfile (str ): Name of the file containing mapping between host and guest nodes
                        If None, mapping is inferred from host and guest node names

    Output:
        nodeMap (dict): A mapping of host -> [guest] nodes
    """

    if mapfile is not None:
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

    else:
        pass

def genMap(host, guest, names=False):
    """
    Works only in the case that node gx_y maps to node hx
    """
    #{guest -> host}
    nodemap = {}
    for leaf in guest:
        gname = leaf.name
        hname = 'h' + gname.split("_")[0][1:]
        if names:
            nodemap[gname] = hname
        else:
            nodemap[leaf] = host&hname
    return nodemap

def writeMapping(nodemap, filename):
    """
    Writes out a mapping between host nodes and guest nodes. Each line of the output
    consists of a host node name followed by the names of each of the guest nodes 
    that maps to it (all separated by spaces).
    """
    outputHandle = open(filename, 'w')

    for node in nodemap:
        out = node.name + '\t' + nodemap[node].name + '\n'
        #out = node.name + ' ' + ' '.join([i.name for i in nodemap[node]]) + '\n'
        outputHandle.write(out)

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

def findDomainsFile(infile, hmmfile, version=2):
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
    if version == 2:
        os.system("hmmsearch232 --domE 0.00001 " + hmmfile + " " + infile 
                    + " | grep '^ [ ]*" + seqName + "' > tmp/grepoutput.txt")
    else:
        os.system("hmmsearch --domE 0.00001 " + hmmfile + " " + infile 
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

def findMotifsFile(infile, mfile):
    """
    Finds all motif matches in infile according to motif in mfile using mast.
    Assumes only one sequence is in the infile

    Args:
    infile  (str): The path to the fasta file to search for domains in
    mfile (str): Path to file containing motif to search for

    Output:
    starts (list): List of starting positions of motifs found
    ends   (list): List of ending positions of motifs found
    seqs   (list): List of motif sequences found 
    """
    starts, ends, seqs = [], [], []
    sequence = list(open(infile))[1].strip()

    os.system("mast -hit_list " + mfile + " " + infile + " > tmp/mast_output.txt")
    f = list(open('tmp/mast_output.txt'))

    for line in f:
        if line[0] != "#":
            temp = line.split()
            starts.append(int(temp[4]) - 1) #MAST is 1-indexed :(
            ends.append(int(temp[5]) - 1)
            seqs.append(sequence[int(temp[4]) : int(temp[5])])

    return starts, ends, seqs

def findDomains(sequence, hmmfile, version=2):
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
    
    return findDomainsFile('tmp/tmp.fa', hmmfile, version)

def findMotifs(sequence, mfile):
    """
    Finds motifs in the given sequence

    Args:
    sequence (str): sequence to search for motifs in
    hmmfile (str): Path to file containing hmm to search for

    Output:
    starts (list): List of starting positions of motifs found
    ends   (list): List of ending positions of motifs found
    seqs   (list): List of motif sequences found 
    """
    g = open('tmp/tmp.fa','w')
    g.write('>seq\n' + sequence)
    g.close()
    
    return findMotifsFile('tmp/tmp.fa', mfile)

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

def isValid(domain, hmmfile):
    """Checks if the input string is a valid zf-C2H2 domain"""

    #Can HMMER find it?
    seqs = findDomains(domain, hmmfile)[2]
    if len(seqs) > 0 and domain != seqs[0]:
        return False

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

    #Remove previous run if it exists (RAxML will not clobber existing results)
    if "RAxML_bestTree.nwk" in os.listdir('.'):
        os.system('rm RAxML_bestTree.nwk')

    command = 'raxml -s '
    command += infile + ' -n ' + outext + ' -m PROTGAMMAJTT -T 8 -p ' + str(np.random.randint(2000) + 1)
    command += ' > raxml_log.txt'
    os.system(command)

def raxml_score_from_file(benchfile, testfile, seqfile):
    """
    Uses RAxML to compute the SH score between the benchmark tree constructed by RAxML and
    a set of other trees. Outputs a binary vector where the ith entry is 1 if tree i is 
    significantly worse than the benchmark and 0 otherwise

    Args:
    benchfile (str): path to newick file containing the best tree found by RAxML
    testfile  (str): path to newick file containing all trees to test, one per line
    seqfile   (str): path to fasta file containing leaf sequences

    Output:
    scores (list): The list of ML scores given by raxml for each input tree
    worse  (list): A list of 0/1 entries specifying whether or not each tree is
                   significantly worse than the benchmark or not. 
    """
    #Run RAxML to find if tree in benchfile is significantly better than those in testfile
    #command = '/home/caluru/Downloads/standard-RAxML-master/raxmlHPC-PTHREADS-AVX2 '
    command = 'raxml '
    #Switch to -f h if this takes too long
    command += '-f H -t ' + benchfile + ' -z ' + testfile + ' -s ' + seqfile + ' -m PROTGAMMAJTT -T 8 -n sco'
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
    worse = []
    for line in f:
        if 'Tree: ' not in line:
            continue
        score = float(line.split()[3])
        scores.append(score)
        answer = line.split('Significantly Worse: ')[1].split()[0]
        worse.append(1 if answer == 'Yes' else 0)

    os.system('rm *.sco')
    return scores, worse

def raxml_score(benchTree, testTrees, seqfile):
    """
    Uses RAxML to compute the SH score between the benchmark tree constructed by RAxML and
    a set of other trees. Outputs a binary vector where the ith entry is 1 if tree i is 
    significantly worse than the benchmark and 0 otherwise. Takes in ete3 objects, writes
    them to the appropriate files and calls raxml_score_file

    Args:
    benchTree (Tree): The best tree found by raxml
    testTree  (list): List of trees to test against the best tree
    seqfile   (str): path to fasta file containing leaf sequences

    Output:
    scores (list): The list of ML scores given by raxml for each input tree
    worse  (list): A list of 0/1 entries specifying whether or not each tree is
                   significantly worse than the benchmark or not. 
    """
    g = open('bestTree.nwk', 'w')
    g.write(benchTree.write(format = 9) + '\n')
    g.close()

    g = open('otherTrees.nwk', 'w')
    for tree in testTrees:
        g.write(tree.write(format=9) + '\n')
    g.close()

    scores, worse = raxml_score_from_file('bestTree.nwk', 'otherTrees.nwk', seqfile)
    os.system('rm bestTree.nwk')
    os.system('rm otherTrees.nwk')

    return scores, worse

#TODO: only tested with Tree and dictionary inputs, not file inputs
def run_treefix(host, guest, lmap, sequences):
    os.system('mkdir -p tfix/config/; mkdir -p tfix/data/0/')

    if type(host) == str:
        os.system('cp ' + host + ' tfix/config/host.stree')
    else:
        writeTree(host, 'tfix/config/host.stree')

    #guest is a path to a tree file
    if type(guest) == str:
        os.system('cp ' + guest + ' tfix/data/0/0.nt.raxml.tree')
    #guest is a tree object
    else:
        writeTree(guest, 'tfix/data/0/0.nt.raxml.tree')

    #lmap is a file
    if type(lmap) == str:
        os.system('cp ' + lmap + 'tfix/config/guest.smap')
    #lmap is a dictionary of guest -> host nodes
    else:
        f = open('tfix/config/guest.smap', 'w')
        for key in lmap:
            f.write(key.name + '\t' + lmap[key].name + '\n')
        f.close()

    #sequences is always a path to a fasta file
    os.system('cp ' + sequences + ' tfix/data/0/0.nt.align')

    #Run TreeFix
    cmd = 'treefix -s tfix/config/host.stree'
    cmd += ' -S tfix/config/guest.smap' 
    cmd += ' -A nt.align '
    cmd += ' -o nt.raxml.tree' 
    cmd += ' -n nt.raxml.treefix.tree'
    cmd += ' -V 0'
    cmd += ' -l data/0/0.nt.raxml.treefix.log'
    cmd += ' -e " -m PROTGAMMAJTT"'
    cmd += ' tfix/data/0/0.nt.raxml.tree'

    print cmd
    os.system(cmd)

    out = Tree('tfix/data/0/0.nt.raxml.treefix.tree')
    os.system('rm -r tfix')
    return out
