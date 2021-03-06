#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  3 16:48:07 2021

@author: joelbjervig
"""

import time
import numpy as np
from math import *
import matplotlib.pyplot as plt

import sys
sys.path.insert(1, '/Users/joelbjervig/documents/universitet/kurser/pagaende/computational physics/computational-physics/homeExc')

from library import bode


# constants
V = 3
E = 4
a = 10
r_max = 3*a
b = np.arange(0,r_max/2-1,0.1) # so arcsin wont let b go up to r_max, only b<r_max/2
r_min = b*np.sqrt(1-V/E)

print(r_min)

# PROBLEMS: some values are not accepted in sqrt and in "true_divide". need to figure out range?
# first integrand term for b<r<r_max


theta1 = lambda r: (1/r**(2)*1/np.sqrt((1-b**2/r**2)))
# second integrand term for r_min<r<r_max
theta2 = lambda r: (1/r**(2)*1/np.sqrt((1-b**2/r**2)))

#theta2 = lambda r: 1/r**(-2)*1/((1-(b**2/r**2)-V/E))**(0.5)

N = 401

def theta(b,rmax,rmin):
    return 2*b*np.subtract(bode(theta1,b+1e-5,rmax,N), bode(theta2,rmin+1e-5,rmax,N))
    #return 2*b*np.subtract(bode(theta1,b+1e-5,rmax,N), bode(theta2,rmin+1e-5,rmax,N))

print(theta(b,r_max,r_min))
#r = np.arange(0.01,10,0.01)



"""
# analytical solution integral
anSol = lambda b: 2*( (np.arcsin(b/r_max*1/((1-V/E)**0.5)) - pi/2) - ( np.arcsin(b/r_max)-pi/2 ) )






# plot analytical solution
plt.figure(1)
plt.xlabel("r")
plt.ylabel("Theta")
plt.plot(r,theta1(r)+theta2(r))

plt.legend()
plt.show()
"""

"""
# plot analytical solution
plt.figure(1)
plt.title("Analytical solution")
plt.xlabel("impact paramenter b")
plt.ylabel("Theta")
plt.plot(b,anSol(b))

plt.legend()
plt.show()
"""