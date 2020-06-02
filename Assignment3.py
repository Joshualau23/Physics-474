# -*- coding: utf-8 -*-
"""
Created on Wed Oct 16 22:26:05 2019

@author: joshu
"""

from scipy import optimize as op
import numpy as np
import matplotlib.pyplot as plt
import datetime
from astropy.io import ascii
import random
from scipy.stats import linregress
from astropy import constants
from math import e,log


def f(x):
    func = log(x+1)/x - 1/(1+x)
    func_deri = ((x-(x+1)*np.log(x+1))/(x*x*(x+1))) - ((2*x+1)/(x*x*(x+1)*(x+1))) + (1/(x*x))
    vp = 2*func**(0.5)*func_deri
    return vp

def velocity(x):
    func = log(x+1)/x - 1/(1+x)
    v = func**(0.5)
    return v


r = op.bisect(f,1,3)
vmax = velocity(r)

print(vmax)






