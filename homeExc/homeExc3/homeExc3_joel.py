#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 17 14:29:12 2021

@author: joelbjervig
"""

import time
import numpy as np
from math import *
import matplotlib.pyplot as plt

import sys
sys.path.insert(1, '/Users/joelbjervig/documents/universitet/kurser/pagaende/computational physics/computational-physics/homeExc')

from library import P
from library import multipleroots
from library import bissection
from library import bode
from library import simpsons2
from library import gausslegendre


# function of interest
f = lambda x: np.sqrt(1-x**2)
a = -1
b = 1

Nlp = 14 # Legendre polynomial of order N has N-1 roots
# legendre polynomials of order N and N+1 (for weights)
p = lambda x: P(x,N)
p1 = lambda x: P(x,N+1)

# weights
w = lambda x, p: 2*(1-x**2)/((N+1)**2*(p1(x))**2)

print()
start = time.process_time()     # meassures time taken to loop
integGL = gausslegendre(f, Nlp)
print("Runtime for Gauss Legendre quadrature with N = 14: T = ",time.process_time() - start)
print("Gauss Legendre accuracy = ", abs(integGL-pi/2))

"""
integGL = 0
for L in range(1,Nlp,4):
    print("N = ", L)
    integGL = gausslegendre(f, L)
    integBode = bode(f,a,b,L)
    integSimp = simpsons2(f,a,b,L)

print("Integral using Gauss-Legendre quadreture: ", integGL)
print("Integral using Bode method = \t ", integBode)
print("Integral using Simpsons method = \t ", integSimp)
"""
eps = 1.7*10**(-6)
N=np.arange(5,1000+2,4) 
integBodeVec = np.zeros(N.shape)
integGLVec = np.full(N.shape,integGL)


start = time.process_time()     # meassures time taken to loop
integBode = bode(f,a,b,489)
print("Runtime for Bodes rule with N = 489: T = ", time.process_time() - start)
print("Bode accuracy = ", abs(integBode-pi/2))

"""
for i in range(len(N)):
    integBodeVec[i] = bode(f,a,b,N[i])

plt.figure(1)
plt.title("Error bode")
plt.xlabel("number of points")
plt.ylabel("Error")
plt.semilogy(N,abs(integBodeVec-np.pi/2),label="Error Bode rule")
plt.semilogy(N,abs(integGLVec-np.pi/2),label=("Error Gauss Legendre quadrature for N = ", Nlp))

plt.legend()
plt.show()
"""
