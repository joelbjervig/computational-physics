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

n = 100
N = 4*n

# function def
phi_GT  = lambda r: -1/sqrt(2*a)*(np.exp(-a*r))
phi_LT  = lambda r: 1/sqrt(2*a)*(np.exp(a*r)-np.exp(-a*r))
S       = lambda r: -1/2*r*np.exp(-r)

phi_GTS = lambda r: -1/sqrt(2*a)*(np.exp(-a*r)) * (-1/2*r*np.exp(-r))
phi_LTS = lambda r: 1/sqrt(2*a)*(np.exp(a*r)-np.exp(-a*r)) * (-1/2*r*np.exp(-r))

phi_true = lambda r : (1/(1-16)**2)*(np.exp(-4*r)-np.exp(-r)*(1+0.5*(1-16)*r))

r = np.linspace(0,r_max,N)
phi = np.zeros(r.shape)
phi_t = phi_true(r)

for i in range(len(r)):
    phi[i] = phi_GT(r[i])*bode(0,r[i],N,phi_LTS) + phi_LT(r[i])*bode(r[i],r_max,N,phi_GTS) # N or n? and are the bounds 0 ro r?
    

plt.figure(1)
plt.title("error")
plt.xlabel("Amplitude")
plt.ylabel("r")
plt.plot(r,abs(phi-phi_t),label="error ye")


plt.legend()
plt.show()