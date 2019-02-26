#Holds useful statistical models for generating trees
import numpy as np
from scipy import stats

def expfunc(minimum=1, maximum=3):
    """Exponential distribution with lambda=0.75 and min/max parameters"""
    return min(maximum, max(minimum, int(stats.expon(0.75).rvs())))

def gaussNoise(n):
    """"Returns a float drawn from N(n, .1*n)"""
    return np.random.normal(n, .1*n)

def exp(l):
    return lambda x: l * stats.expon(-1 * l * x).rvs()

def gamma(x):
    return stats.gamma(x, x**2).rvs()

def gaussian():
    return lambda x: np.random.normal(x, x)