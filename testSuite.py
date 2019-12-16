#Test suite for HostTreeGen and GuestTreeGen. Change each if statement to true to run that test
import unittest
from TreeSearch import reconcile, reconcileDL, generate_rootings
from ete3 import Tree

class TestTreeSearch(unittest.TestCase):
    
    def setUp(self):
        self.host1 = Tree('((5:1,4:1)1:1,((1:1,0:1)1:1,2:1)1:1);')
        self.guest1 = Tree('(((11:1,10:1)1:1,((7:1,(3:1,2:1)1:1)1:1,(5:1,(1:1,0:1)1:1)1:1)1:1)1:1,(13:1,12:1)1:1);')
        tempmap = {11: 0, 10: 4, 13: 1, 12: 1, 1: 1, 0: 4, 3: 4, 2: 2, 5: 4, 7: 4}
        self.map1 = {}
        for key in tempmap:
            self.map1[self.guest1&str(key)] = self.host1&(str(tempmap[key]))

        self.host2 = Tree('(A,B)C;', format=1)
        self.guest2 = Tree('(((1,2)7,(3,4)8)10,(5,6)9)11;', format=1)
        tempmap = {'1':"A",'2':"C",'3':"A",'4':"C",'5':"B",'6':"B"}
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

        self.assertEqual(cost, 13)
        self.assertEqual(namemap, {'1':"A",'2':"C",'3':"A",'4':"C",'5':"B",'6':"B",'7':"C",'8':"C",'9':"B",'10':"C",'11':"C"})

    def testReconcile(self):
        cost, mapping = reconcile(self.host2, self.guest2, self.map2)
        namemap = {}
        for key in mapping:
            namemap[key.name] = mapping[key].name

        self.assertEqual(cost, 11)
        self.assertEqual(namemap, {'1':"A",'2':"C",'3':"A",'4':"C",'5':"B",'6':"B",'7':"C",'8':"C",'9':"C",'10':"C",'11':"C"})

if __name__ == '__main__':
    unittest.main()
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