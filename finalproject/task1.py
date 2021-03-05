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
E = 400
a = 10
r_max = 3*a
b = np.arange(0.01,r_max/2-1,0.1) # so arcsin wont let b go up to r_max, only b<r_max/2
r_min = b*(1-V/E)


# PROBLEMS: some values are not accepted in sqrt and in "true_divide". need to figure out range?
# first integrand term for b<r<r_max
theta1 = lambda r: 2*b*(r**(-2)*1/((1-b**2/r**2)**(0.5)))
# second integrand term for r_min<r<r_max
theta2 = lambda r: 2*b*(r**(-2)*1/((1-b**2/r**2-V/E)**(0.5)))


# analytical solution integral
anSol = lambda b: 2*( (np.arcsin(b/r_max*1/((1-V/E)**0.5)) - pi/2) - ( np.arcsin(b/r_max)-pi/2 ) )


N = 401

I1 = bode(theta1,b,r_max,N)
I2 = bode(theta2,r_min,r_max,N)

# plot analytical solution
plt.figure(1)
plt.title("Analytical solution")
plt.xlabel("impact paramenter b")
plt.ylabel("Theta")
plt.plot(b,180/pi*anSol(b))

plt.legend()
plt.show()