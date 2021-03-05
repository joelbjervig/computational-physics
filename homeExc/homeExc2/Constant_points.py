#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 12 09:37:12 2021

@author: Finn Joel Bjervig, Anton Palm Ekspong, Thomas Herard
"""
import time
import numpy as np
from math import *
import matplotlib.pyplot as plt

import sys
sys.path.insert(1, '/Users/joelbjervig/documents/universitet/kurser/pagaende/computational physics/computational-physics/homeExc')

from library import bode

# const def
a = 4;
r_max = 30

n = 100
N = 4*n

# function def

#homogenous solutions
phi_GT  = lambda r: -1/sqrt(2*a)*(np.exp(-a*r))
phi_LT  = lambda r: 1/sqrt(2*a)*(np.exp(a*r)-np.exp(-a*r))

#source terms
S       = lambda r: -1/2*r*np.exp(-r)

#hom. sol multiplied with the source term for defining the 
phi_GTS = lambda r: -1/sqrt(2*a)*(np.exp(-a*r)) * (-1/2*r*np.exp(-r))
phi_LTS = lambda r: 1/sqrt(2*a)*(np.exp(a*r)-np.exp(-a*r)) * (-1/2*r*np.exp(-r))

#analytical solution of the differential equation
phi_true = lambda r : (1/(1-16)**2)*(np.exp(-4*r)-np.exp(-r)*(1+0.5*(1-16)*r))

# create vectors for the numerical solution
r = np.linspace(0,r_max,N)
phi = np.zeros(r.shape)

# boundary values
phi_0 = 0
phi_inf = 0
phi[0] = phi_0
phi[-1] = phi_inf


#loops over r_m

start = time.process_time()     # meassures time taken to loop
for i in range(1,len(r)-1):
    phi[i] = phi_GT(r[i])*bode(0,r[i],N,phi_LTS) + phi_LT(r[i])*bode(r[i],r_max,N,phi_GTS) # N or n? and are the bounds 0 ro r?
print(time.process_time() - start)

# plot of error
plt.figure(1)
plt.title("Deviation from analytical solution")
plt.xlabel("r")
plt.ylabel("Amplitude")
plt.plot(r,abs(phi-phi_true(r)))

plt.legend()
plt.show()


# plot of homogenous solutions
print(r)
plt.figure(2)
plt.title("Solutions to the homogenous problem")
plt.xlabel("r")
plt.ylabel("Amplitude")
plt.plot(r, phi_GT(r), label="PHI_>")
plt.plot(r, phi_LT(r), label="PHI_<")

plt.legend()
plt.show()
