#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  4 15:50:03 2021

@author: Joel Bjervig, Anton Palm, Thomas Herard
"""

import matplotlib.pyplot as plt
import numpy as np
from math import *

import sys
sys.path.insert(1, '/Users/joelbjervig/documents/universitet/kurser/pagaende/computational physics/computational-physics/homeExc')

from library import numerov


# function defenitions
#denisty distribution function
rho = lambda r:1/(8*pi)*exp(-r)
# k^2 term
k = lambda x:0
# source term
S = lambda r:-1/2*r*exp(r)
# analytical solution
phi_exact = lambda r: 1-1/2*(r+2)*exp(-r)

# grid
a = 0;
b = 30;
N = 1000
x = np.linspace(a,b,N)
h = (b-a)/N

# IC
y0 = 0
y1 = phi_exact(h)

# create empty vector
y = np.zeros(x.shape)

# call numerovvv
y = numerov(k,S,y0,y1,x)

plt.figure(1)
plt.title("Solution to Second Order ODE Using Numerov Method")
plt.xlabel("Radius r")
plt.ylabel("Electric potential phi")
#plt.plot(t,p2,label="Analytical solution")
plt.plot(x,y,label="Numerov method")
plt.legend()
plt.show()