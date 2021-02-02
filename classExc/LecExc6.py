#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  2 09:25:41 2021

@author: Joel Bjervig, Anton Palm, Thomas Herard
"""

import matplotlib.pyplot as plt
import numpy as np
from math import *

def RK4(f,H,x,y0):
    y = np.zeros(x.shape)
    y[0] = y0
    
    h = x[1]-x[0] # could move into for loop as h=x[i]-x[i-1]
    for i in range(1,length(y)):
        k1 = H*f(x[i-1],y[i-1])
        k2 = H*f(x|[i-1])+0.5*H,y[i-1]+0.5*K1)
        k3 = H*f(x|[i-1])+0.5*H,y[i-1]+0.5*K2)
        k4 = H*f(x|[i-1])+H,y[i-1]+K3)
        y[i] = y[i-1]+(1/6)*(K1+2*K2+2*K3+K4)
    return y
    
f = lambda x,y:-x*y
dfx = lambda x,y:-y
dfy = lambda x,y:-x

y = lambda x:exp(-x*x/2)



# grid
a=0
b=1
N=np.arange(1,1000,3)
H=(b-a)/N

for i in range(len(N)):
    x=np.linspace(a,b,N[i])
    RK4(f,H,x)

    