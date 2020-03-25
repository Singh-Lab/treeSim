#Draws from an orthogroup, creates a root sequence with desired properties

from TreeUtils import findDomains
from TreeUtils import printDomSeq as pds
from os import listdir as ls
from numpy.random import choice
from ConfigParser import ConfigParser as CP
from OrthoAnalysis import selfSimilarity
import string

DATAPATH = CP.ORTHOGROUP_PATH #pylint: disable=no-member
hmmfile = CP.HMM_PATH #pylint: disable=no-member

def fixzf(dom):
    #Adds in the C's and H's if they aren't there
    dom = list(dom)
    dom[2] = "C"
    dom[5] = "C"
    dom[18] = "H"
    dom[22] = "H"
    return ''.join(dom)

#TODO: Bug causing more domains than required occasionally (grs2 works normally).
def genRandomSequence(numDoms):
    """
    Generates a random zf-C2H2 protein sequence with the given number of domains.
    Includes sequence before and after zf domain.
    """
    files = ls(DATAPATH)
    f = list(open(DATAPATH + choice(files)))[1::2]
    sequence = choice(f).strip()
    sequence.translate(None, '-')
    
    starts, ends, seqs = findDomains(sequence, hmmfile)
    if len(starts) < numDoms:
        return genRandomSequence(numDoms)
    prefix = sequence[:starts[0]]
    suffix = sequence[ends[-1]:]
    if prefix == '' or suffix == '':
        return genRandomSequence(numDoms)
    linkers = []
    for i in range(len(starts)-1):
        linkers.append(sequence[ends[i]+1:starts[i+1]])
    
    middle = ''
    for _ in range(numDoms - 1):
        middle += choice(seqs) + choice(linkers)
    middle += choice(seqs)

    newSeq = prefix + middle + suffix
    newSeq = ''.join(newSeq.split('-'))

    #Deletes all lowercase letters
    newSeq = newSeq.translate(None, string.ascii_lowercase)
    #Deletes all illegal AA characters
    newSeq = newSeq.translate(None, 'BJOUXZ')

    return newSeq

def genRandomSequence2(numDoms):
    """
    Generates a zf-C2H2 protein sequence with the given number of domains.
    Picks 1 domain per sequence, one sequence per orthogroup
    """
    files = ls(DATAPATH)
    for i in range(len(files))[::-1]:
        if '.fa' not in files[i]:
            files.pop(i)
    pool = []
    i = 0
    while i < numDoms:
        f = list(open(DATAPATH + choice(files)))[1::2]
        f = choice(f).strip()
        if len(findDomains(f, hmmfile)[0]) > 1:
            pool.append(string.translate(f, None, '-'))
            i += 1

    starts, ends = findDomains(pool[0], hmmfile)[:2]
    prefix = pool[0][:starts[0]]
    suffix = pool[0][ends[-1]:]
    if prefix == '' or suffix == '':
        return genRandomSequence2(numDoms)

    i,j = starts[0], starts[1]
    middle = fixzf(pool[0][i:j])

    for sequence in pool[1:]:
        starts = findDomains(sequence, hmmfile)[0]
        i, j = starts[0], starts[1]
        middle += fixzf(sequence[i:j])

    newSeq = prefix + middle + suffix
    newSeq = newSeq.translate(None, '-') #''.join(newSeq.split('-'))
    newSeq = newSeq.translate(None, string.ascii_lowercase)
    newSeq = newSeq.translate(None, 'BJOUXZ')

    if len(findDomains(newSeq, hmmfile)[0]) != numDoms:
        return genRandomSequence2(numDoms)

    return newSeq

if __name__ == '__main__':
    crapcounter = 0
    for i in range(1000):
        s = genRandomSequence2(1)
        #pds(s, hmmfile)
        x = selfSimilarity('d', s, hmmfile, False)
        if len(x) != 1:
            print 'fuck'
            crapcounter += 1

    print crapcounter
