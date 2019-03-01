#Implements various models of branch relaxation
from ete3 import Tree
import stats

def uncorrelatedRelaxation(tree, method='exponential'):
    """
    Relaxes each branch length individually without an explicit time process.
    Multiplies each branch length by a rate parameter drawn from the specified
    distribution. Available distributions are 'exponential', 'gamma', or 'gaussian'

    Args:
        tree   (Tree): The ete3 Tree with branch lengths to relax
        method (str ): The distribution to draw rates from. Options are 'exponential',
                       'gamma', and 'gaussian'. Default is exponential, others will
                       default to exponential
    """

    rateFunc = stats.exp
    #if method == 'gamma': rateFunc = stats.gamma
    if method == 'gaussian' : rateFunc = stats.gaussian
    
    for node in tree.traverse():
        rate = rateFunc(node.dist)
        node.add_feature('rate', rate)
        node.dist = rate
        
def autocorrelatedRelaxation(tree):
    """
    Relaxes branches using an autocorrelated process. In this process, each node is
    assigned a rate parameter R_i, and the rate multiplier on a branch of length l between nodes
    i and j is l * (R_i + R_j) / 2. The root node's rate is set to 1.
    """

    tree.add_feature('rate', 1)

    for node in tree.traverse():
        if node == tree:
            continue
        node.add_feature('rate', stats.logN(node.up.rate, node.dist))
        node.dist *= (node.up.rate + node.rate) / 2.

    