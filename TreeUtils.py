#Provides several utility functions useful for reading/writing trees, etc
from ete3 import Tree
from scipy.stats import expon
from math import exp
import numpy as np

def read(filename):
    """Reads a tree in newick format to an ete3 object"""
    return Tree(filename, format=1)

def write(tree, filename):
    """Writes an ete3 tree to the given filename in newick format"""
    tree.write(outfile=filename, format=1)

#Useful statistical distributions for varying branch lengths
def expfunc(minimum=1, maximum=3):
    """Exponential distribution with lambda=0.75 and min/max parameters"""
    return min(maximum, max(minimum, int(expon(0.75).rvs())))

def gaussNoise(n):
    """"Returns a float drawn from N(n, .1*n)"""
    return np.random.normal(n, .1*n)

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
