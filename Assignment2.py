# -*- coding: utf-8 -*-
"""
Created on Sun Oct  6 11:12:22 2019

@author: joshu
"""

import numpy as np
import matplotlib.pyplot as plt
import datetime
from astropy.io import ascii
import random
from scipy.stats import linregress
from astropy import constants
from math import e,sqrt
from scipy import optimize


### Q2 b)


data=ascii.read("phy474_ass2_SNdata.txt")

cz_list = data['log10cz']

mmax_list = data['mmax']

'''
def lineofbest(m,x,b):
    return m*x + b
best_fit = optimize.curve_fit(lineofbest,cz_list, mmax_list)

slope_fit = best_fit[0]
f = np.poly1d(slope_fit)
error = best_fit[1]
slope_err = error[0][0]
intercept_err = error[1][1]

y_fit = f(cz_list)


plt.figure()
plt.scatter(cz_list, mmax_list,color = "deepskyblue")
plt.plot(cz_list, y_fit,color = "Gold")
plt.xlabel("log(cz)")
plt.ylabel('m_max')
plt.title("Maximum Magnitude versus Log(cz)" )
plt.grid()
plt.savefig('A2_2b.png')
plt.show()

print slope_fit , sqrt(slope_err), sqrt(intercept_err)

'''
### Q2 c)

'''
def lineofbest2(x,b):
    return 5*x + b

best_fit = optimize.curve_fit(lineofbest2,cz_list, mmax_list)

intercept_fit = best_fit[0][0]
f2 = np.poly1d([5,intercept_fit])
error = best_fit[1][0][0]


y_fit2 = f2(cz_list)


plt.figure()
plt.scatter(cz_list, mmax_list,color = "deepskyblue")
plt.plot(cz_list, y_fit2,color = "Black")
plt.xlabel("log(cz)")
plt.ylabel('m_max')
plt.title("Maximum Magnitude versus Log(cz)" )
plt.grid()
plt.savefig('A2_2c.png')
plt.show()

print intercept_fit , sqrt(error)

stddev_list =[]
for i in range(0,12):
    ymu = mmax_list[i]-y_fit2[i]
    stddev_list.append(ymu**2)

print(stddev_list)

sigma = np.sqrt((1/12.0)*(sum(stddev_list)))
print(sigma)
'''

### Q3 a)

from scipy import optimize as op

def surface_brightness(r):
    return ((np.exp(-r))*(1 - r)) - 0.5

r1 = 0.0
r2 = 10.0

root = op.bisect(surface_brightness,r1,r2)

print(root)

### Q3 b)
