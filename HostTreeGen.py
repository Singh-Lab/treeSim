from ete3 import Tree
import numpy as np
import stats

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

def createRandomTopology(numLeaves, treeHeight, heightDistribution):
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

def birthDeathTree(birthRate, deathRate, treeHeight):
    """
    Generates a tree topology according to the birth-death model.
    
    Args:
        birthRate (float): birth rate 
        deathRate (float): death rate
        treeHeight (float): The average overall length of a root to leaf path
        numLeaves (int): The number of leaves desired at the end of the run. If 
                         the input is <= 0, this parameter is ignored and   
    """
    birthRate = float(birthRate)
    deathRate = float(deathRate)

    host = Tree()
    host.dist = 0
    lineages = [(host, treeHeight)]

    while lineages != []:
        #waiting time is exp(b + d), P(b) = b/(b+d), P(d) = 1 - P(b)
        node, height = lineages.pop(0)
        eventTime = stats.exp(1./(1./birthRate + 1./deathRate))

        #event occurs
        if eventTime <= height:
            #duplication
            if np.random.random() < birthRate / (birthRate + deathRate):
                left = node.add_child(dist=eventTime)
                right = node.add_child(dist=eventTime)
                lineages.append((left, height - eventTime))
                lineages.append((right, height - eventTime))
            #loss: Remove from queue, delete node later (cleanup process)
            else:
                node.dist *= -1
        #If no event occurs, credit remaining branch length to this node
        else:
            node.dist += height
                
    if host.children == []:
        host.name = "H0"
        return host

    #remove lost nodes
    for node in host.traverse():
        if node.dist < 0:
            node.up.remove_child(node)

    #remove nodes with only one child (ensure full binary tree)
    for node in [a for a in host.traverse()]:
        if len(node.children) == 1:
            #This is the root node
            if node.up == None:
                host.children[0].dist += host.dist
                host = host.children[0]
                host.up = None
            else:
                parent = node.up
                child = node.children[0]
                child.dist += node.dist
                child.up = parent
                parent.remove_child(node)
                parent.children.append(child)

    nameCounter = 0
    for node in host.traverse():
        node.name = "H" + str(nameCounter)
        nameCounter += 1

    return host

def BDMinMax(birthRate, deathRate, treeHeight, minLeaves = 0, maxLeaves = float('inf')):
    """
    Creates a birth death tree within the desired range of leaves. Produces a warning if
    the input birth/death rates are unlikely to produce a tree with the desired number of 
    leaves.
    """
    #TODO: ADD WARNING IF PARAMETERS ARE MISMATCHED -> sum of exps is gaussian (goddammit)
    """
    mean = 1./(1./birthRate - 1./deathRate)

    if minLeaves > 2 * mean or maxLeaves < 0:
        print """

    incorrect = True

    while(incorrect):
        host = birthDeathTree(birthRate, deathRate, treeHeight)
        numLeaves = len([leaf for leaf in host])
        incorrect = numLeaves <= minLeaves or numLeaves >= maxLeaves

    return host
    

def readHostTree(treeFile, treeHeight = -1):
    """
    Reads input file and ensures that every node has a label (or labels it otherwise) and
    optionally rescales the branch lengths so that the average height is the desired one. 
    Unlike readTree in TreeUtils, this function is intended to read trees not created by
    this package.

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
    t = birthDeathTree(0.3,0.1,1)
    print t.get_ascii()
    print t.get_ascii(attributes=['dist'])