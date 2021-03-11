#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  5 16:09:37 2021

@author: joelbjervig
"""


import numpy as np
import math
import matplotlib.pyplot as plt


import sys
sys.path.insert(1, '/Users/joelbjervig/documents/universitet/kurser/pagaende/computational physics/computational-physics/homeExc')

from library import bode
from library import bissection

# def real(b):
#     retrun 2*(np.arcsin(b/rmax) -np.pi/2-np.arccos(b/rmax*np.sqrt(E/(E-V))))



# Constants
N=100001
V=1
a=3
# r>b = 1otherwise theta1 will be imaginary
# r>b*1/srqt(1-V/E) =  1*otherwise theta2 will be imaginary
# 1/sqrt(1-V/E) must be real, so E>V
b=1 
rmax=3*a

# Energy is varying
E = V*np.arange(1.1,10.1,0.1)

# Leonard Jones potential
pot = lambda r: 4*V*((a/r)**12-(a/r)**6)

# first integrand term for b<r<r_max
theta1 = lambda r: 2*b*(r**(-2)*1/(np.sqrt(1-b**2/r**2)))
# second integrand term for r_min<r<r_max
theta2 = lambda r, E: 2*b*(r**(-2)*1/(np.sqrt(1-b**2/r**2-pot(r)/E)))

# Now when the potential is dependent on r, so is rmin
rminfunc = lambda r,E: 1-b**2/r**2-pot(r)/E
rmin = np.zeros(len(E))

# Find the positive root of the function
for i in range(0,len(E)):
    rminfuncE = lambda r: 1-b**2/r**2-pot(r)/E[i]
    rmin[i] = bissection(0.1, rmax, rminfuncE, 100)
    
# numberical solution
def theta(b,E,rmax,rmin):
    return np.subtract(bode(theta1,b,rmax,N), bode(theta2,rmin,rmax,N))

sol = np.zeros(len(E))

print(theta1(2))
#print(theta2(2))
"""
for i in range(0,len(E)):

    print("b = ", b, " | E[i] = ",E[i], " | rmax = ", rmax," | rmin[i] = ",rmin[i])
    print("THETA SHAPE:", np.shape(theta(b, E[i],rmax,rmin[i])))

    sol[i]= theta(b, E[i], rmax, rmin[i])
 """   
plt.figure()
plt.plot(E,sol,label="Numerical solution")
plt.xlabel("Energy E")
plt.ylabel('Deflection angle [radians]')
plt.legend(loc='best')
plt.show()
