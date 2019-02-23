#Holds useful statistical models for generating trees
import numpy as np
from scipy.stats import expon

def expfunc(minimum=1, maximum=3):
    """Exponential distribution with lambda=0.75 and min/max parameters"""
    return min(maximum, max(minimum, int(expon(0.75).rvs())))

def gaussNoise(n):
    """"Returns a float drawn from N(n, .1*n)"""
    return np.random.normal(n, .1*n)

def exp(x):
    return expon(-1 * x)