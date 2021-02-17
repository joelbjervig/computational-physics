#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 17 11:33:43 2021

@author: joelbjervig
"""

import matplotlib.pyplot as plt
import numpy as np
from math import *
import random

#number of points
N = 20

# function of circle for calculating pi
circle = lambda x: np.sqrt(y**2-1)

# weight
w0 = lambda x: 1
w1 = lambda x: (4-2*x)/3

# uniformily distributed vector. no need to sort since MC method doesnt require it
y = np.random.uniform(0,1,(1,N))

#change of vaiables for each respective  weight w0 and w1
x0 = w0(y)
x1 = w1(y)

# function we want to integrate using MC - monte carlo
f = lambda x: 1/(1+x**2)


def monte_carlo(N, f, w, x):
    I = 0;
    #sum the fractions 
    for (i in range(N)):
        I = I + f(x[i])/w(x[i])
    
    # take average
    I = I/N
    
    return I

def varianceSqrd(N,f,y):
    S1 = 0
    S2 = 0
    S = 0
    
    #sum the fractions 
    for (i in range(N)):
        S1 = S1 + (f(y[i])**2)/N    # first summation
        S2 = S2 + ((f(y[i]))/N)**2  # second summation
        S = S1-S2                   # take the differance
        
    return S/N # = Variance squared
    
    