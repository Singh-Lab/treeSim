#Provides several utility functions useful for reading/writing trees, etc
from ete3 import Tree
from math import exp
from Utils import sortBy
import numpy as np
import os

#Tree I/O TODO: switch to NHX format to store node annotations
def read(filename):
    """Reads a tree in newick format to an ete3 object"""
    return Tree(filename, format=1)

def write(tree, filename):
    """Writes an ete3 tree to the given filename in NHX format"""
    tree.write(outfile=filename, format=1)

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
	os.system("hmmsearch232 -E 0.00001 " + hmmfile + " " + infile 
				+ " | grep '^ [ ]*" + seqName + "' > tmp/grepoutput.txt")
	hits = list(open('tmp/grepoutput.txt'))
	starts, ends, seqs = [], [], []

	for line in hits:
		line = line.split()
		if not line[1].isdigit():
			continue
		starts.append(int(line[1]) - 1) #HMMER is 1-indexed :/
		ends.append(int(line[3]) - 1)
		seqs.append(line[2])

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

#TODO: Change so that it takes one domain per column rather than one domain per column
def genRandomSequence(numDoms, datapath, hmmfile):
    """
    Generates a random zf-C2H2 protein sequence with the given number of domains.
    Creates sequences using a set of orthogroups (each in fasta format, multiple 
    sequence alignments allowed)found in the datapath folder.

    Args:
        numDoms  (int ): The number of domains in the output sequence
        datapath (str ): the path to the folder containing orthogroups to emulate
        hmmfile  (str ): path to the hmmfile, in HMMER format
    """
    files = os.listdir(datapath)
    f = list(open(datapath + np.random.choice(files)))[1::2]
    sequence = np.random.choice(f).strip()
    
    starts, ends, seqs = findDomains(sequence, hmmfile)
    if len(starts) < numDoms:
        return genRandomSequence(numDoms, datapath, hmmfile)
    prefix = sequence[:starts[0]]
    suffix = sequence[ends[-1]:]
    if prefix == '' or suffix == '':
        return genRandomSequence(numDoms, datapath, hmmfile)
    linkers = []
    for i in range(len(starts)-1):
        linkers.append(sequence[ends[i]+1:starts[i+1]])
    
	middle = ''
    for _ in range(numDoms - 1):
        middle += np.random.choice(seqs) + np.random.choice(linkers)
    middle += np.random.choice(seqs)

    newSeq = prefix + middle + suffix
    return ''.join(newSeq.split('-'))
