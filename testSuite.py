import unittest
from TreeSearch import reconcile, reconcileDL, generate_rootings
from ete3 import Tree

class TestTreeSearch(unittest.TestCase):
    
    def setUp(self):
        self.host1 = Tree('((5,4)1,((1,0)1,2)1);')
        self.guest1 = Tree('(((11,10)1,((7,(3,2)1)1,(5,(1,0)1)1)1)1,(13,12)1);')
        tempmap = {11: 0, 10: 4, 13: 1, 12: 1, 1: 1, 0: 4, 3: 4, 2: 2, 5: 4, 7: 4}
        self.map1 = {}
        for key in tempmap:
            self.map1[self.guest1&str(key)] = self.host1&(str(tempmap[key]))

        self.host2 = Tree('(A,B)C;', format=1)
        self.guest2 = Tree('(((1,2)7,(3,4)8)10,(5,6)9)11;', format=1)
        tempmap = {'1':"A",'2':"B",'3':"A",'4':"B",'5':"A",'6':"A"}
        self.map2 = {}
        for key in tempmap:
            self.map2[self.guest2&str(key)] = self.host2&(str(tempmap[key]))

    def tearDown(self):
        pass

    def testReconcileDL(self):
        cost, mapping = reconcileDL(self.host2, self.guest2, self.map2)
        namemap = {}
        for key in mapping:
            namemap[key.name] = mapping[key].name

        self.assertEqual(cost, 7)
        self.assertEqual(namemap, {'1':"A",'2':"B",'3':"A",'4':"B",'5':"A",'6':"A",'7':"C",'8':"C",'9':"A",'10':"C",'11':"C"})

    def testReconcile(self):
        cost, mapping = reconcile(self.host2, self.guest2, self.map2)
        namemap = {}
        for key in mapping:
            namemap[key.name] = mapping[key].name

        self.assertEqual(cost, 6.5)
        self.assertEqual(namemap, {'1':"A",'2':"B",'3':"A",'4':"B",'5':"A",'6':"A",'7':"C",'8':"C",'9':"C",'10':"C",'11':"C"})

    def testSPR(self):
        """
        Three cases:
            1) Child of internal node to child of internal node
            2) Child of root to child of internal node
            3) Child of internal node to child of root

        Exceptions (ValueError):
            1) Root to any
            2) Node to self
            3) Node to parent
            4) Node to sibling
            5) Node to any node in its subtree
        """

        from TreeSearch import spr

        t1 = Tree('(((A,B)F,C)H,(D,E)G)R;', format=1)
        t2 = Tree('((A,(D,E)G)R,(B,C)F)H;', format=1)
        t3 = Tree('(((A,C)H,(D,E)G)R,B)F;', format=1)

        t = Tree('((A,(B,C)F)H,(D,E)G)R;', format=1)
        t = spr(t, t&"B", t&"A")
        assert(t1.robinson_foulds(t)[0] == 0)

        t = Tree('((A,(B,C)F)H,(D,E)G)R;', format=1)
        t = spr(t, t&"G", t&"A")
        assert(t2.robinson_foulds(t)[0] == 0)

        t = Tree('((A,(B,C)F)H,(D,E)G)R;', format=1)
        t = spr(t, t&"B", t&"R")
        assert(t3.robinson_foulds(t)[0] == 0)

        t = Tree('((A,(B,C)F)H,(D,E)G)R;', format=1)
        self.assertRaises(ValueError, spr, t, t&"R", t&"F")
        self.assertRaises(ValueError, spr, t, t&"F", t&"F")
        self.assertRaises(ValueError, spr, t, t&"A", t&"F")
        self.assertRaises(ValueError, spr, t, t&"H", t&"B")

if __name__ == '__main__':
    unittest.main()

    """
    host1 = Tree('((5,4)1,((1,0)1,2)1);')
    guest1 = Tree('(((11,10)1,((7,(3,2)1)1,(5,(1,0)1)1)1)1,(13,12)1);')
    tempmap = {11: 0, 10: 4, 13: 1, 12: 1, 1: 1, 0: 4, 3: 4, 2: 2, 5: 4, 7: 4}
    map1 = {}
    for key in tempmap:
        map1[guest1&str(key)] = host1&(str(tempmap[key]))
    cost, map1 = reconcileDL(host1, guest1, map1)
    print cost
    from TreeUtils import writeReconciliation
    writeReconciliation(host1, guest1, map1)
    """
"""
import HostTreeGen, GuestTreeGen
from Utils import printProgressBar
from TreeUtils import s3
from stats import expfunc, gaussNoise
import numpy as np
"""

"""
RUN_ALL = False

#Basic test to ensure errors aren't encountered
if(False or RUN_ALL):
    print "\nRandom Topology + BuildGuestTree"
    for i in range(100):
        treeHeight = np.random.random()
        host = HostTreeGen.createRandomTopology(i, treeHeight, gaussNoise)
        guest, nodemap = GuestTreeGen.buildGuestTree(host, s3, expfunc, treeHeight / 5, gaussNoise, 10)
        printProgressBar(i+1, 100)

if(True or RUN_ALL):
    print "\nBirthDeath Host + BuildGuestTree"
    sizes = []
    for i in range(100):
        treeHeight = np.random.random()
        host = HostTreeGen.birthDeathTree(.3, .1, treeHeight)
        sizes.append(len([leaf for leaf in host]))
        guest, nodemap = GuestTreeGen.buildGuestTree(host, s3, expfunc, treeHeight / 5, gaussNoise, 10)
        printProgressBar(i+1, 100)
    print sizes
    print np.mean(sizes)
"""