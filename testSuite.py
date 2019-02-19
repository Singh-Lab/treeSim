#Test suite for HostTreeGen and GuestTreeGen
from HostTreeGen import genHostTree
from GuestTreeGen import buildGuestTree
from Utils import printProgressBar
from TreeUtils import s3, expfunc, gaussNoise
import numpy as np

#Basic test to ensure errors aren't encountered
for i in range(100):
    treeHeight = np.random.random()
    host = genHostTree(i, treeHeight, gaussNoise)
    guest = buildGuestTree(host, s3, expfunc, treeHeight / 5, gaussNoise, 10)
    printProgressBar(i+1, 100)