#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  3 10:01:20 2021

@author: Joel Bjervig, Thomas Herrad, Anton Palm Ekspong
"""
import numpy as np
from math import *
import matplotlib.pyplot as plt
from library import RK2
from library import RK3
from library import RK4
from library import RK4_SYS









# exacto solution of problem
p_exact = lambda t: c1*cos(2*pi*t)- 2*pi*c2*sin(2*pi*t)
y_exact = lambda t: c2*cos(2*pi*t) + c1/(2*pi)*sin(2*pi*t)

# ODE
f = lambda t,y,p:p
g = lambda t,y,p: -4*pi**2*y
y0=1    # initial condition just trying omething
p0=1

a=0
b=1.5
N=np.arange(1,1000,10)  # num internals. Num points is N+1


# initialize vectors
t = np.linspace(a,b,N[99])
y=RK4_SYS(f,g,y0,p0,t)[0]
# y = np.zeros(N.shape)
# p = np.zeros(N.shape)

# # iterate through number of intervals. Q: how does error change with changing N (H)
# for i in range(len(N)):
#     t = np.linspace(a,b,N[i])
    

#     y[i] = RK4_SYS(f,g,y0,p0,t)[0]


