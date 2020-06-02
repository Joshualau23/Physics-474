# -*- coding: utf-8 -*-
"""
Created on Fri Nov  8 13:27:43 2019

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
from math import e,log,pi,sqrt


### Part a)

Msun = 1.99e30
ro = 4e6*Msun 
a = 1.0
G = 6.67e-11*(1/3.086e16)*(1/3.086e16)*(1/3.086e16)

def velocity_dispersion_a(r):
    x = r/a
    vel = sqrt((2.0/9.0)*pi*G*ro*(a**2)*(1.0/sqrt(1+(x**2))))
    vel_kms = vel*3.086e13
    return vel_kms




r = np.arange(0.1,10,0.1)
velocity_list = []

for i in range(0,len(r)):
    v = velocity_dispersion_a(r[i])
    velocity_list.append(v)
    
plt.figure()
plt.plot(r,velocity_list,color = "deepskyblue")
plt.title('Velocity Dispersion Part (a)')
plt.xlabel('Distance r (pc)')
plt.ylabel('Velocity Dispersion (km/s)')
plt.grid()
plt.savefig('2c1.png')

### Part b)

Mbh = 4e6*Msun

def velocity_stars(r):
    x = r/a
    vel = sqrt((2.0/9.0)*pi*G*ro*(a**2)*(1.0/sqrt(1+(x**2))))
    vel_kms = vel*3.086e13
    return vel_kms

def velocity_bh(r):
    vel = ((G*Mbh)/(3.0*r))*(1.0 / ((1.0+r**2))**(5.0/2.0))*((-8.0/3.0)+((8*r**4 + 12*r**2 + 3.0)/(3.0*r*((1+r**2))**(3/2.0))))
    vel_kms = sqrt(vel)*3.086e13
    return vel_kms
    #return vel



def Total_velocity(r):
    return velocity_stars(r) + velocity_bh(r)

r = np.arange(0.1,10,0.1)
velocity_bhlist = []
velocity_starslist = []
velocity_total = []

for i in range(0,len(r)):
    v_stars = velocity_stars(r[i])
    v_bh = velocity_bh(r[i])
    v_total = Total_velocity(r[i])
    velocity_starslist.append(v_stars)
    velocity_bhlist.append(v_bh)
    velocity_total.append(v_total)

radius = []

for i in range(0,len(r)):
    if velocity_bhlist[i] - velocity_starslist[i] >= 0:
        radius.append(r[i])
        
print (radius[-1])

plt.figure()
plt.plot(r,velocity_total,color = "deepskyblue")
plt.plot(r,velocity_bhlist,color = "lawngreen")
plt.plot(r,velocity_starslist,color = "gold")
plt.title('Velocity Dispersion Part (b)')
plt.xlabel('Distance r (pc)')
plt.ylabel('Velocity Dispersion (km/s)')
plt.grid()
plt.savefig('2c2.png')




