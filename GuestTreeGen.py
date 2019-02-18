from ete3 import Tree
import numpy as np

LOSS_CODE = -1
DUPLICATION_FACTOR = 3

#####################################
#                                   #
#         Helper Functions          #
#                                   #
#####################################

def gaussNoise(n):
    """"Returns a float drawn from N(n, .1*n)"""
    return np.random.normal(n, .1*n)

def createNode():
    """Creates a domain node with required fields precreated"""
    node = Tree()
    node.name = 'placeholder'
    node.add_feature('pos', 0)
    node.add_feature('event', 'SPECIATION')
    node.dist = 0
    return node

def createTree(n):
    """Creates a starter domain tree with n leaves"""
    tree = Tree()
    tree.populate(n)

    i = 0
    for node in tree.traverse():
        node.add_feature('event', 'DUPLICATION')
        node.name = "G0_" + str(i)
        i += 1

    i = 0
    for leaf in tree:
        leaf.add_feature('pos', i)
        i += 1
        leaf.add_feature('event', 'SPECIATION')
        leaf.dist = 0
    return tree

def addcherry(node, hostName = '', domCounter = 0):
    """Takes a node with no children and adds two children"""
    left = createNode()
    right = createNode()

    if hostName != '':
        left.name = "G" + hostName[1:] + "_" + str(domCounter)
        right.name = "G" + hostName[1:] + "_" + str(domCounter + 1)
        domCounter += 2

    left.add_feature('pos', 0)
    right.add_feature('pos', 1)

    node.add_child(left)
    node.add_child(right)

    return domCounter

def clean(host, guest):
    """
    Takes in host and guest trees and removes all attributes other than name, dist, and event
    """
    allowedFeatures = ['support','dist','name','event']

    for node in host:
        features = node.features
        for feature in features:
            if feature not in allowedFeatures:
                node.del_feature(feature)

    for node in guest:
        features = node.features
        for feature in features:
            if feature not in allowedFeatures:
                node.del_feature(feature)        

#####################################
#                                   #
#       Guest Tree Generation       #
#                                   #
#####################################

def buildGuestNode(startingTree, dupRateFunc, dupfunc, hostName = '', 
                    branchLength = 1, branchFunc = gaussNoise, eventDist = .05):
    """
    Creates topology of guest tree within a single node of the host tree.
    Starting tree nodes must have a pos feature (position of domain in sequence).
    Number of events at this node is determined by branchLength and eventDistribution

    Args:
        startingTree (Tree): The symbolic tree to evolve. For example, a tree with 
                             three leaves will be treated as a sequence with three 
                             domains
        dupRateFunc (func): function with # domains as input and a duplication rate output
        dupfunc (func): A function which returns the size of a tandem duplication.
                        Must be of the form dupfunc(x,y), where x and y are the 
                        min and max duplication sizes respectively
        hostName (str ): If given, hostName (of the form H<HOST_ID>) will be used to label 
                         guest nodes in the form G<HOST_ID>_<GENE_ID>
        branchLength (float): used to determine the number of events that occur here
        eventDist (float): Average distance between events
        branchFunc (func): takes eventDist as input, returns actual distance between two events

    Returns:
        branchLength (float): 
    """
    domCounter = 0 #if hostName is given, names of dom nodes will be 'hostName_domCounter'
    leaves = [leaf for leaf in startingTree]
    for leaf in leaves:
        leaf.add_feature('bl', branchLength)
    if hostName != '':
        for leaf in leaves:
            leaf.name = 'G' + hostName[1:] + "_" + str(domCounter)
            domCounter += 1
    leaves.sort(key=lambda x: x.pos)
    dist = branchFunc(eventDist)

    while(branchLength > dist and len(leaves) > 0):
        event = np.random.random()
    
        #Duplication
        if event < dupRateFunc(len(leaves)):
            branchLength -= dist
            numDomains = len(leaves)
            size = dupfunc(1, len(leaves))
            start = np.random.randint(numDomains - size + 1)
            

            for i in range(start, start+size):
                node = leaves[i]
                node.event = "DUPLICATION"
                domCounter = addcherry(node, hostName, domCounter)
                node.children[0].dist = dist
                node.children[1].dist = DUPLICATION_FACTOR * dist
                node.children[0].add_feature('bl', branchLength)                
                node.children[1].add_feature('bl', branchLength)

                leaves[i] = node.children[0]
                node.children[0].pos = node.pos
                node.children[1].pos = node.pos + size
                leaves.insert(i + size, node.children[1])

            for i in range(start + 2*size, len(leaves)):
                leaves[i].pos += size
                
        #Loss
        else:
            position = np.random.randint(len(leaves))
            leaves[position].pos = LOSS_CODE
            leaves[position].event = 'LOSS'
            leaves.pop(position)
            for i in range(position, len(leaves)):
                leaves[i].pos -= 1

        dist = branchFunc()

    for node in leaves:
        if 'bl' in node.features:
            node.dist += node.bl

    return branchLength

def buildGuestTree(host, dupRateFunc, dupfunc, eventDist, branchFunc, startSize):
    """
    Creates a guest tree topology inside of the host tree given parameters

    Args:
        host (Tree): The host tree (ete3 format) to evolve inside
        dupRateFunc (func): A function that takes the current number of domains as input
                            and returns the rate at which duplications occur
        dupfunc (func): A function which returns the size of a tandem duplication.
                        Must be of the form dupfunc(x,y), where x and y are the 
                        min and max duplication sizes respectively
        eventDist (func): function that takes no arguments and returns the distance to 
                          the next evolutionary event
        startSize (int ): The number of leaves in the initial guest, prior to any
                          evolutionary event.

    Returns:
        guest (Tree):   The complete guest tree topology
        nodemap (dict): host -> guest mapping of nodes
    """
    guest = createTree(startSize)
    guest.name = "G" + host.name[1:] + "_0"
    guest.add_feature('pos', 0)
    nodemap = {}
    branchLength = 0

    for node in host:
        node.add_feature('bl', node.dist)
        nodemap[node] = []

    for hostNode in host.traverse():
        if hostNode == host: #Root node
            opstack, branchLength = buildGuestNode(guest, dupRateFunc, dupfunc, hostNode.name, hostNode.dist)
            hostNode.add_feature('ops', opstack)
            hostNode.bl = branchLength
            nodemap[hostNode] = [node for node in guest.traverse()]
            hostNode.add_feature('leaves', [leaf for leaf in guest if leaf.pos != LOSS_CODE])
            
        else:
            #Find last occurrence of domain nodes above me
            gleaves = hostNode.up.leaves
            hostNode.add_feature('leaves', [])
            if gleaves == []:
                continue
            newTrees = []
            treepositions = []


            #create a subtree for each gleaf
            for i in range(len(gleaves)):
                if gleaves[i].pos != LOSS_CODE:
                    t = createNode()
                    t.pos = gleaves[i].pos
                    newTrees.append(t)
                    treepositions.append(i)
            if len(treepositions) == 0:
                continue
            supernode = Tree()
            supernode.children = newTrees

            opstack, branchLength = buildGuestNode(supernode, dupRateFunc, dupfunc, 
                                    hostNode.name, hostNode.dist, branchFunc, eventDist)
            hostNode.bl = branchLength

            hostNode.leaves += [leaf for leaf in supernode if leaf != LOSS_CODE]
            nodemap[hostNode] = [n for n in supernode.traverse()]
            nodemap[hostNode].pop(0)

            for i in range(len(treepositions)):
                gleaves[treepositions[i]].add_child(newTrees[i])
    
    clean(host, guest)
    return guest, nodemap