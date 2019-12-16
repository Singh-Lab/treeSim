#Given a starting tree, perform search around tree to find tree with high likelihood and low reconciliation score

from ete3 import Tree
import numpy as np
import os
import logging
from ilp_generator import *
from gurobipy import *
from TreeUtils import raxml, raxml_score

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
    subtree = np.random.choice(nodes)
    invalid = set([i for i in subtree.traverse()])
    remaining = set(nodes) - invalid
    ns = np.random.choice(list(remaining))
    print 'subtree', subtree.name, 'ns', ns.name
    spr(tree, subtree, ns)
    return (subtree.label, ns.label)

def pick_sprs(tree, n):
    """
    Tries to pick n unique sprs of a tree and return a list of those sprs.
    If that many are difficult to find, may give up early and return fewer

    Arguments:
    tree: The tree to find random sprs for
    n: The number of unique trees to attempt to find

    Output:
    the list of sprs
    """
    i = 0
    for node in tree.traverse():
        node.add_feature('label', str(i))
        i += 1

    sprs = []
    seen = set()
    i = 0
    failcount = 0

    while i < n:

        newTree = tree.copy()
        used = pick_spr(newTree)

        if used not in seen:
            seen.add(used)
            sprs.append(newTree)
            i += 1
            failcount = 0

        else:
            failcount += 1
            if failcount >= 100:
                return sprs
                
    logging.debug('Found ' + str(len(sprs)) + ' sprs out of ' + str(n))
    return sprs

def generate_rootings(tree):
    """
    Takes a rooted tree as input and generates all possible rerootings 
    of that tree. Does not change the input tree

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

    for node in copy.traverse():
        if node == copy:
            continue
        temp = copy.copy()
        temp.set_outgroup(temp&(node.name))
        print temp.get_ascii()
        trees.append(temp)

    return trees

#TODO: Test this function
def reconcile(host, guest, leafmap):
    """
    Performs TDL reconciliation on the guest tree and returns both the cost and the full
    mapping of guest -> host.

    Args:
    host (Tree): The host tree in ete3 format
    guest (Tree): The guest tree in ete3 format
    leafmapping (dict): guest -> host mapping containing all leaf nodes in the guest tree

    Output:
    cost (float): The cost of the min cost reconciliation with L = 1 and 
                  D(n) = 2 + .5(n-1)
    fullmap (dict): guest -> mapping containing all nodes in the guest tree
    """
    h = createTreeRepresentation(host)
    g = createTreeRepresentation(guest)
    d = createDistMatrix(host)
    mapping = createMapping(leafmap)

    eqnames, rhs, coldict = createEqns(host, guest, h, g, mapping, d)
    write('treesolve.mps', eqnames, rhs, coldict)

    m = read('treesolve.mps') # pylint: disable=undefined-variable
    m.Params.outputflag = 0
    m.optimize()
    cost = m.getObjective().getValue()

    fullmap = extractMapping(m, host, guest)

    os.system('rm treesolve.mps')
    return cost, fullmap

#TODO: Test this function
def reconcileDL(host, guest, leafmap):
    """
    Performs DL reconciliation on the guest tree and returns both the cost and the full
    mapping of guest -> host.

    Args:
    host (Tree): The host tree in ete3 format
    guest (Tree): The guest tree in ete3 format
    leafmapping (dict): guest -> host mapping containing all leaf nodes in the guest tree

    Output:
    cost (float): The cost of the min cost reconciliation with L = 1 and D = 2
    fullmap (dict): guest -> mapping containing all nodes in the guest tree
    """
    fullmap = {}
    cost = 0
    seen = set()
    jobs = [leaf for leaf in guest]

    #Perform Reconciliation
    while jobs != []:
        job = jobs.pop(0)
        if job in leafmap:
            fullmap[job] = leafmap[job]
        else:
            a, b = job.children
            hostmap = fullmap[a].get_common_ancestor(fullmap[b])
            fullmap[job] = hostmap
        seen.add(job)
        guestParent = job.up
        if job.up == None:
            continue
        a, b = guestParent.children
        if a in seen and b in seen:
            jobs.append(guestParent)

    #Compute Cost
    for node in guest.traverse():
        #CASE 1: Leaf, no cost 
        if node.children == []:
            continue
        #CASE 2: Both children mapped to same node (DUPLICATION) easy duplication
        myhost = fullmap[node]
        lchild = fullmap[node.children[0]]
        rchild = fullmap[node.children[1]]

        if myhost == lchild and myhost == rchild:
            cost += 2
        #CASE 3: either lchild or rchild is mapped to the same but not the other (DUPLOSS)
        elif myhost == lchild:
            cost += 2 + myhost.get_distance(rchild)
        elif myhost == rchild:
            cost += 2 + myhost.get_distance(lchild)
        #CASE 4: Neither child is mapped to myhost (SPECIATION + LOSS)
        else:
            cost += myhost.get_distance(rchild) + myhost.get_distance(lchild) - 2

    return cost, fullmap

def reroot(host, guest, leafmap, recModule=reconcile):
    """
    Roots the input guest tree by checking all possible rootings for the lowest 
    reconciliation score. If the midpoint rooting has the lowest possible score,
    it will be used. #TODO: make the second sentence true.

    Args:
    host (Tree): The host tree to reconcile against
    guest (Tree): The guest tree to reroot by reconciliation
    leafmap (dict): guest -> host mapping including all guest leaves
    recModule (func): Function to reconcile guest with host. Must take as input
                      the host, guest, and leafmap arguments passed here and 
                      output (cost, fullmap) where cost is the float cost of the 
                      reconciliation and fullmap is a guest -> host mapping 
                      including all guest nodes.

    Output:
    The rooting of the guest tree with the lowest reconciliation cost.
    """
    trees = generate_rootings(guest)
    print len(trees)
    costs = []

    copy = guest.copy()
    copy.set_outgroup(copy.get_midpoint_outgroup())
    trees.append(copy)
    
    for tree in trees:
        #Generate new mapping
        newmap = {}
        for node in leafmap:
            newguest = tree&(node.name)
            newmap[newguest] = leafmap[node]
        costs.append(recModule(host, tree, newmap)[0])

    index = np.argmin(costs)
    return trees[index]

def perform_search(sequences, host, guest, leafmap, num_iter=100):
    """
    Performs a search in tree space surrounding the highest scoring guest 
    tree according to raxml. Tries to find a guest with a similar likelihood 
    value and a better reconciliation score.

    Steps:
    1) Determine reconciliation score of guest w.r.t. host
    2) Generate 100 SPR moves, test all for reconciliation and raxml scores
    3) Pick a better tree, move to it
    4) Repeat 2-3 100 times
    """

    bestTree = guest
    bestScore = reconcileDL(host, guest, leafmap)

    #Base nodemap on names, not actual nodes
    nodemap = {}
    for node in leafmap:
        nodemap[node.name] = leafmap[node].name

    for iteration in range(num_iter):
        logging.debug('Iteration number ' + str(iteration))
        sprs = pick_sprs(guest, 100)
        scores = raxml_score(bestTree, sprs, sequences)
        good_trees = [sprs[i] for i in range(len(sprs)) if scores[i] == 0]
        logging.debug('Found ' + str(len(good_trees)) + ' close trees')

        rec_scores = []
        for i in range(len(good_trees)):
            tree = good_trees[i]
            lmap = {}
            for node in tree:
                lmap[node] = host&(nodemap[node.name])
            good_trees[i] = reroot(host, tree, lmap, reconcileDL) #Maybe use reconcile() instead? (Test performance)

            tree = good_trees[i]
            lmap = {}
            for node in tree:
                lmap[node] = host&(nodemap[node.name])
            rec_scores.append(reconcileDL(host, tree, leafmap)[0]) #Replace this with reconcile() after testing

        index = np.argmax(rec_scores)
        newScore = rec_scores[index]

        if newScore <= bestScore:
            logging.debug('Did not find a better tree')
        else:
            logging.debug('Found better tree, new: ' + str(newScore) + ' old: ' + str(bestScore))
            bestTree = good_trees[index]
            bestScore = newScore

    return bestTree
        
if __name__ == '__main__':
    logging.basicConfig(level=logging.DEBUG, format="%(asctime)15s %(message)s")
    logging.info('This is a logging test')
    t = Tree('RAxML_bestTree.out')
    t.set_outgroup(t.get_midpoint_outgroup())
    seqs = 'ENSG00000196724.fa'
    a = pick_sprs(t, 100)
    logging.info('found %s sprs', str(len(a)))
    logging.info(str(raxml_score(t, a, seqs)))

    #perform_search(seqs, a, t, )
