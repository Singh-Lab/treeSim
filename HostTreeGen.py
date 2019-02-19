from ete3 import Tree
import numpy as np

#####################################
#                                   #
#       Host Tree Generation        #
#                                   #
#####################################

def assignBranchLengths(tree, treeHeight, heightDistribution):
    """
    Assigns branch lengths to all nodes of the tree such that the sum of branch lengths
    on any root to leaf path is equal to hostLength. Only used on a whole tree (host
    trees here, for the full version see shilpa/scripts/simulation.py). 
    """

    leaves = [leaf for leaf in tree]
    treeHeight = float(treeHeight)

    for node in tree.traverse():
        node.dist = 0
        lset = []
        for leaf in leaves:
            if node.get_common_ancestor(leaf) == node:
                lset.append(leaf)
        height = np.min([leaf.get_distance(node, topology_only=True) for leaf in lset])
        node.add_feature('height', height)

    for node in tree.traverse():
        remainingDist = treeHeight - (node.get_distance(tree) + tree.dist)
        if node.height == 0:
            node.dist = remainingDist
        else:
            node.dist = heightDistribution(remainingDist / (node.height + 1.))

def genHostTree(numLeaves, treeHeight, heightDistribution):
    """
    Generates a host tree with the specified number of leaves and approximate overall height.
    Args:
        numLeaves (int ): The number of leaves desired in the host tree
        treeHeight (float): The average overall length of a root to leaf path
        heightDistribution (func): 
    Output:
        host (Tree): The host tree in ete3 Tree format
    """
    host = Tree()
    host.populate(numLeaves)
    nameCounter = 0

    for node in host.traverse():
        node.name = "H" + str(nameCounter)
        nameCounter += 1
    
    assignBranchLengths(host, treeHeight, heightDistribution)

    return host

def readHostTree(treeFile, treeHeight = -1):
    """ Reads input file and ensures that every node has a label (or labels it otherwise) and
        optionally rescales the branch lengths so that the average height is the desired one

        Args:
            treeFile (str ): name of file to read
            treeHeight (float ): Rescales branch lengths such that the mean root to leaf path
                                 matches the desired tree height. Ignored if treeHeight <= 0
        Output:
            host (Tree): The host tree in ete3 Tree format
    """
    host = Tree(treeFile)

    #Ensure all nodes are named
    nameCounter = 0
    for node in host.traverse():
        if node.name == '':
            node.name = "H" + str(nameCounter)
            nameCounter += 1

    #If rescaling is desired
    if treeHeight > 0:
        height = 0.0
        for leaf in host:
            height += leaf.get_distance(host)
        height /= len([leaf for leaf in host])
        height += host.dist #the "stem" before the root isn't included in a root to leaf path
    
        for node in host.traverse():
            node.dist /= height

    return host

#Test cases
if __name__ == "__main__":
    pass