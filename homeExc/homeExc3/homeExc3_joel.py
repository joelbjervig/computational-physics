#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 17 14:29:12 2021

@author: joelbjervig
"""

import numpy as np
from math import *
import matplotlib.pyplot as plt

import sys
sys.path.insert(1, '/Users/joelbjervig/documents/universitet/kurser/pagaende/computational physics/computational-physics/homeExc/homeExc1')

from library import P
from library import multipleroots
from library import bissection
from library import bode
from library import simpsons2

a = -1
b = 1

N = 14 # Legendre polynomial of order 10 has 9 roots

f = lambda x: np.sqrt(1-x**2)

# legendre polynomials of order N and N+1 (for weights)
p = lambda x: P(x,N)
p1 = lambda x: P(x,N+1)

# weights
w = lambda x, p: 2*(1-x**2)/((N+1)**2*(p1(x))**2)


integGL = 0
for L in range(1,N,4):

    # gauss legendre 
    P_roots = multipleroots(a,b,p,N

    for i in range(N):
        integGL = integGL + w(P_roots[i],p(P_roots[i]))*f(P_roots[i])
        

    integBode = bode(f,a,b,L)
    integSimp = simpsons2(f,a,b,L)

print("Integral using Gauss-Legendre quadreture: ", integGL)
print("Integral using Bode method = \t ", integBode)
print("Integral using Simpsons method = \t ", integSimp)
