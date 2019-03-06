#Test suite for HostTreeGen and GuestTreeGen. Change each if statement to true to run that test
import HostTreeGen, GuestTreeGen
from Utils import printProgressBar
from TreeUtils import s3
from stats import expfunc, gaussNoise
import numpy as np

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
        treeHeight = 5 * np.random.random()
        host = HostTreeGen.birthDeathTree(.3, .1, treeHeight)
        sizes.append(len([leaf for leaf in host]))
        guest, nodemap = GuestTreeGen.buildGuestTree(host, s3, expfunc, treeHeight / 5, gaussNoise, 10)
        printProgressBar(i+1, 100)
    print sizes
    print np.mean(sizes)