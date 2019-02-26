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

    rateFunc = stats.exp(1)
    if method == 'gamma': rateFunc = stats.gamma
    if method == 'gaussian' : rateFunc = stats.gaussian()