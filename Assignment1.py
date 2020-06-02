# -*- coding: utf-8 -*-
"""
Created on Sun Sep 15 19:42:15 2019

@author: joshu
"""


### Q1 a)

import numpy as np
import matplotlib.pyplot as plt
import datetime
from astropy.io import ascii
import random
from scipy.stats import linregress
from astropy import constants
from math import e

'''

Mv_list = []
R_list = []

Mv_sun = 4.83
sigma = 0.3


for i in range(0,2000):
    abs_mag =  np.random.normal(Mv_sun, sigma)
    R = random.randrange(70,130,1)
    Mv_list.append(abs_mag)
    R_list.append(R)
    
plt.figure()
plt.scatter(R_list, Mv_list,color = "gold")
plt.xlabel("r[pc]")
plt.ylabel('Mv')
plt.title("Mv versus Radius" )
plt.grid()
plt.savefig('A1a.png')
plt.show()

### Q1 b)



Mv_StarsWeSee = []
R_StarsWeSee = []


def apparentmag_formula(Mv,d):
    return Mv + 5*np.log10((d/10))

for i in range(0,2000):
    if apparentmag_formula(Mv_list[i],R_list[i]) >= 10 :
        continue
    else:
        Mv_StarsWeSee.append(Mv_list[i])
        R_StarsWeSee.append(R_list[i])
    
plt.figure()
plt.scatter(R_list, Mv_list,color = "gold")
plt.scatter(R_StarsWeSee, Mv_StarsWeSee,color = "deepskyblue")
plt.xlabel("r[pc]")
plt.ylabel('Mv')
plt.title("For Stars we can see: Mv versus Radius" )
plt.grid()
plt.savefig('A1b.png')
plt.show() 


print sum(Mv_StarsWeSee) / float(len(Mv_StarsWeSee))
    
### Q1 c)

#average distance of the stars we can see
print sum(R_StarsWeSee) / float(len(R_StarsWeSee))

#average distance of the sample of stars
print sum(R_list) / float(len(R_list))

'''

### Q2


'''

from scipy.integrate import quad
from math import e

def star_formation(t,tau):
    return e**(-t/tau)


I = quad(star_formation, Tgal - Tms,Tgal,args=(tau))

tau_list = np.linspace(2,10,41)


Tgal = 10

# 2M stars

Tms = 1.1
fraction_list = []

for i in range(0,41):
    tau = tau_list[i]
    I2 = quad(star_formation, Tgal - Tms,Tgal,args=(tau))
    I1 = quad(star_formation, Tgal,0,args=(tau))
    fraction = I2[0]*-1 / I1[0]
    fraction_list.append(fraction)
    print tau
    print fraction



plt.figure()
plt.scatter(tau_list, fraction_list,color = "gold")
plt.xlabel("Tau")
plt.ylabel('Fraction')
plt.title("Fraction of Stars That are Still Visible Today (2M)" )
plt.grid()
plt.savefig('A1_2_1.png')
plt.show() 




# 3M stars


Tms = 0.35
fraction_list = []

for i in range(0,41):
    tau = tau_list[i]
    I2 = quad(star_formation, Tgal - Tms,Tgal,args=(tau))
    I1 = quad(star_formation, Tgal,0,args=(tau))
    fraction = I2[0]*-1 / I1[0]
    fraction_list.append(fraction)
    #print tau
    #print fraction



plt.figure()
plt.scatter(tau_list, fraction_list,color = "gold")
plt.xlabel("Tau")
plt.ylabel('Fraction')
plt.title("Fraction of Stars That are Still Visible Today (3M)" )
plt.grid()
plt.savefig('A1_2_2.png')
plt.show() 


# 1M stars


Tms = 9.8
fraction_list = []

for i in range(0,41):
    tau = tau_list[i]
    I2 = quad(star_formation, Tgal - Tms,Tgal,args=(tau))
    I1 = quad(star_formation, Tgal,0,args=(tau))
    fraction = I2[0]*-1 / I1[0]
    fraction_list.append(fraction)
    #print tau
    #print fraction



plt.figure()
plt.scatter(tau_list, fraction_list,color = "gold")
plt.xlabel("Tau")
plt.ylabel('Fraction')
plt.title("Fraction of Stars That are Still Visible Today (1M)" )
plt.grid()
plt.savefig('A1_2_3.png')
plt.show()


# 5M stars


Tms = 0.094
fraction_list = []

for i in range(0,41):
    tau = tau_list[i]
    I2 = quad(star_formation, Tgal - Tms,Tgal,args=(tau))
    I1 = quad(star_formation, Tgal,0,args=(tau))
    fraction = I2[0]*-1 / I1[0]
    fraction_list.append(fraction)
    #print tau
    #print fraction



plt.figure()
plt.scatter(tau_list, fraction_list,color = "gold")
plt.xlabel("Tau")
plt.ylabel('Fraction')
plt.title("Fraction of Stars That are Still Visible Today (5M)" )
plt.grid()
plt.savefig('A1_2_4.png')
plt.show()




# 9M stars


Tms = 0.026
fraction_list = []

for i in range(0,41):
    tau = tau_list[i]
    I2 = quad(star_formation, Tgal - Tms,Tgal,args=(tau))
    I1 = quad(star_formation, Tgal,0,args=(tau))
    fraction = I2[0]*-1 / I1[0]
    fraction_list.append(fraction)
    #print tau
    #print fraction



plt.figure()
plt.scatter(tau_list, fraction_list,color = "gold")
plt.xlabel("Tau")
plt.ylabel('Fraction')
plt.title("Fraction of Stars That are Still Visible Today (9M)" )
plt.grid()
plt.savefig('A1_2_5.png')
plt.show()




# 15M stars


Tms = 0.012
fraction_list = []

for i in range(0,41):
    tau = tau_list[i]
    I2 = quad(star_formation, Tgal - Tms,Tgal,args=(tau))
    I1 = quad(star_formation, Tgal,0,args=(tau))
    fraction = I2[0]*-1 / I1[0]
    fraction_list.append(fraction)
    #print tau
    #print fraction



plt.figure()
plt.scatter(tau_list, fraction_list,color = "gold")
plt.xlabel("Tau")
plt.ylabel('Fraction')
plt.title("Fraction of Stars That are Still Visible Today (15M)" )
plt.grid()
plt.savefig('A1_2_6.png')
plt.show()

'''

### Q3 b)

from math import log

Z_sun = 0.0122

p = (Z_sun*(1-0.15))/ log(50/13)


### Q3 c)

massfraction = 1 - e**(-Z_sun*(0.25-0.15)/p)


print massfraction







