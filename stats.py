#Holds useful statistical models for generating trees
import numpy as np
from scipy import stats
import math

def expfunc(minimum=1, maximum=3):
    """Exponential distribution with lambda=0.75 and min/max parameters"""
    return min(maximum, max(minimum, int(stats.expon(0.75).rvs())))

def gaussNoise(n):
    """"Returns a float drawn from N(n, .1*n)"""
    return np.random.normal(n, .1*n)

def exp(m):
    #Draws from an exponential distribution with mean m
    return stats.expon(scale=m).rvs()

"""
#TODO: Implement this model maybe?
def gamma(m):
    return stats.gamma(m).rvs()

"""

def gaussian(m):
    return np.random.normal(m, m)

def logN(m, l):
    """
    Calculates the log of the logN distribution with parameters such that the log of 
    the mean of this distribution is equal to m

    Args:
        m (float): The desired mean of the log of the logN distribution
        l (float): The actual variance of the logN distribution
    """

    m = m - (l ** 2 / 2.)
    l = math.sqrt(l)

    return stats.lognorm(s=l, scale=m).rvs()

if __name__ == '__main__':

    passed = True
    print "\n*** EXP TESTS ***"
    for i in range(10):
        dist = [exp(i) for _ in range(1000)]
        mean, var = np.mean(dist), np.std(dist)
        print "Expect " + str(i), mean, var
        passed &= abs(i - mean) < (.05 * i) and abs(i - var) < (.05 * i)
    if passed:
        print "TESTS PASSED"
    else:
        print "TESTS FAILED"

    passed = True
    print "\n*** GAMMA TESTS ***"
    for i in range(10):
        dist = [exp(i) for _ in range(1000)]
        print "Expect " + str(i), np.mean(dist), np.std(dist)

    passed = True
    print "\n*** GAUSSIAN TESTS ***"
    for i in range(10):
        dist = [gaussian(i) for _ in range(1000)]
        print "Expect " + str(i), np.mean(dist), np.std(dist)