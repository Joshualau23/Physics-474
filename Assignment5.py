# -*- coding: utf-8 -*-
"""
Created on Sat Nov 23 22:47:14 2019

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
from math import e,log,sin,cos


def f(x):
    func = (sin(x) * (x - sin(x))) / ((1.0 - cos(x))**2) + 2.182
    return func

n = op.bisect(f,3,5)

print(n)
