#Draws from an orthogroup, creates a root sequence with desired properties

from CreateSequences import findDomains
from os import listdir as ls
from numpy.random import choice
from TreeUtils import printDomSeq as pds

DATAPATH = '/home/caluru/Data/eutheria/'
hmmfile = '/home/caluru/Data/hmmfiles/zf_shilpa_232.hmm'

def genRandomSequence(numDoms):
    """
    Generates a random zf-C2H2 protein sequence with the given number of domains.
    Includes sequence before and after zf domain.
    """
    files = ls(DATAPATH)
    f = list(open(DATAPATH + choice(files)))[1::2]
    sequence = choice(f).strip()
    
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
    return ''.join(newSeq.split('-'))

def genRandomSequence2(numDoms):
    """
    Generates a zf-C2H2 protein sequence with the given number of domains.
    Picks 1 domain per sequence, one sequence per orthogroup
    """
    files = ls(DATAPATH)
    pool = []
    for _ in range(numDoms):
        f = list(open(DATAPATH + choice(files)))[1::2]
        pool.append(choice(f).strip())

    starts, ends = findDomains(pool[0], hmmfile)[:2]
    prefix = pool[0][:starts[0]]
    suffix = pool[0][ends[-1]:]
    if prefix == '' or suffix == '':
        return genRandomSequence(numDoms)

    i,j = starts[0], starts[1]
    middle = pool[0][i:j]

    for sequence in pool[1:]:
        starts = findDomains(sequence, hmmfile)[0]
        i, j = starts[0], starts[1]
        middle += sequence[i:j]

    newSeq = prefix + middle + suffix
    return ''.join(newSeq.split('-'))

if __name__ == '__main__':
    print
    pds(genRandomSequence(5), hmmfile)
