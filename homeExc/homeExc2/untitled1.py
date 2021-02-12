#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 12 09:37:12 2021

@author: joelbjervig
"""

import numpy as np
from math import *
import matplotlib.pyplot as plt

import sys
sys.path.insert(1, '/Users/joelbjervig/documents/universitet/kurser/pagaende/computational physics/computational-physics/homeExc/homeExc1')

from library import bode

# const def
a = 4;
r_max = 30

n = 10
N = 4*n

# function def
phi_GT  = lambda r: -1/sqrt(2*pi)*(np.exp(-a*r))
phi_LT  = lambda r: 1/sqrt(2*pi)*(np.exp(a*r)-np.exp(-a*r))
S       = lambda r: -1/2*r*np.exp(-r)

r = np.linspace(0,r_max,N)
phi = np.zeros(r.shape)

for x in r:
    phi[x] = phi_GT(x)*bode(phi_LT(r)*S(r),0,r,N) + phi_LT(x)*bode(phi_GT(r)*S(r),r,r_max,N) # N or n? and are the bounds 0 ro r?

plt.figure(2)
plt.title("Title")
plt.xlabel("r")
plt.ylabel("Amplitude")
plt.plot(r,phi,label="solution")
plt.legend()
plt.show()