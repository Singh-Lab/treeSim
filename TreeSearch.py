#Given a starting tree, perform search around tree to find tree with high likelihood and low reconciliation score

from ete3 import Tree
import numpy as np
import os
from ilp_generator import createTreeRepresentation, createDistMatrix, createEqns, write
from gurobipy import *

def spr(tree, subtree, new_sibling):
    """
    Performs a subtree prune and regraft operation moving a subtree to 
    a new location in the tree

    Arguments:
    tree (Tree): The full tree in which to perform an spr move
    subtree (Tree): The subtree that will be pruned and grafted
    new_sibling (Tree): The new sibling of the subtree being moved
    """
    
    if subtree == tree or subtree.up == tree:
        return

    #Add node between new_sibling and its parent
    temp = Tree()
    temp.up = new_sibling.up
    temp.children = [subtree, new_sibling]
    new_sibling.up = temp
    if temp.up.children[0] == new_sibling:
        temp.up.children[0] = temp
    else:
        temp.up.children[1] = temp

    #Remove subtree from its current location
    old_parent = subtree.up
    subtree.up = temp

    #Remove old parent
    ancestor = old_parent.up
    if old_parent.children[0] == subtree:
        other_child = old_parent.children[1]
    else:
        other_child = old_parent.children[0]

    other_child.up = ancestor
    if ancestor.children[0] == old_parent:
        ancestor.children[0] = other_child
    else:
        ancestor.children[1] = other_child

def pick_spr(tree):
    """
    Picks (and performs) a random spr move by picking a subtree to move 
    and a location to move it to

    Arguments:
    tree: The tree to perform a random spr for
    """
    nodes = [i for i in tree.traverse()][3:]
    subtree = np.random.choice(nodes[:3])
    invalid = set([i for i in subtree.traverse()])
    remaining = set(nodes) - invalid
    ns = np.random.choice(list(remaining))
    print 'subtree', subtree.name, 'ns', ns.name
    spr(tree, subtree, ns)

def generate_rootings(tree):
    """
    Takes a rooted tree as input and generates all possible rerootings 
    of that tree

    Arguments:
    tree (Tree): The rooted tree to produce alternate rootings of
    """
    trees = []
    copy = tree.copy()
    copy.unroot()

    #Every node needs a unique name for this to work
    i = 0
    for node in copy.traverse():
        if node.name == '':
            node.name = str(i)
            i += 1

    print copy.get_ascii()
    for node in copy.traverse():
        if node == copy:
            continue
        temp = copy.copy()
        temp.set_outgroup(temp&(node.name))
        print temp.get_ascii()
        trees.append(temp)
    return trees

def reconcile(host, guest, leafmap):
    h = createTreeRepresentation(host)
    g = createTreeRepresentation(guest)
    d = createDistMatrix(host)
    eqnames, rhs, coldict = createEqns(host, guest, h, g, leafmap, d)
    write('treesolve.mps', eqnames, rhs, coldict)

    m = read('treesolve.mps')
    m.optimize()
    cost = m.getObjective.getValue()

    os.system('rm treesolve.mps')
    return cost

def reroot(host, guest, leafmap):
    trees = generate_rootings(guest)
    costs = []
    
    for tree in trees:
        #Generate new mapping
        newmap = {}
        for node in leafmap:
            newguest = tree&(node.name)
            newmap[newguest] = leafmap[node]
            costs.append(reconcile(host, tree, newmap))

    index = np.argmin(costs)
    return trees[index]

"""
t = Tree()
t.populate(200)

i = 0
for node in t.traverse():
    node.name = str(i)
    i += 1

for i in range(300):
    #print t.get_ascii()
    pick_spr(t)
    print i
    #print t.get_ascii()
"""

t = Tree()
t.populate(5)

i = 0
for node in t.traverse():
    node.name = str(i)
    i += 1

a = generate_rootings(t)
print t.get_ascii()
