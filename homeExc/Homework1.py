#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  3 10:01:20 2021

@author: Joel Bjervig, Thomas Herrad, Anton Palm Ekspong
"""
from library import RK2
from library import RK3
from library import RK4
import numpy as np
from math import *

f   = lambda p,y: p         
dfx = lambda x,y:-y             #df/dx
dfy = lambda x,y:-x             #df/dy
y_true = lambda x:exp(-x*x/2)   # analytical solution
y0=1    # initial condition

# grid
a=0
b=1
N=np.arange(1,1000,10)  # num internals. Num points is N+1
H=(b-a)/N

# initialize vectors

rk4 = np.zeros(N.shape)

# iterate through number of intervals. Q: how does error change with changing N (H)
for i in range(len(N)):
    x = np.linspace(a,b,N[i])
    rk4[i] = RK4(f,x,y0)[-1]
