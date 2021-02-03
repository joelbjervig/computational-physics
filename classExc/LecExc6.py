#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  2 09:25:41 2021

@authors: Joel Bjervig, Anton Palm, Thomas Herard
"""

import matplotlib.pyplot as plt
import numpy as np
from math import *

# Runge-Kutta 2nd order
def RK2(f, x, y0):
    y = np.zeros(x.shape)
    y[0] = y0
    
    for i in range(len(x) - 1):
        h = x[i+1] - x[i]
        
        k = h*f(y[i], x[i])
        
        y[i+1] = y[i] + h*f(x[i]+h/2, y[i]+k/2)
        
    return y

# Runge-Kutta 3rd order
def RK3(f, x, y0):
    y = np.zeros(x.shape)
    y[0] = y0
    
    for i in range(1,len(y)):
        h = x[i]-x[i-1]
        
        k1 = h*f(x[i-1],    y[i-1])
        k2 = h*f(x[i-1]+h/2,y[i-1]+k1/2)
        k3 = h*f(x[i-1]+h,  y[i-1]-k1+2*k2)
        
        y[i]=y[i-1]+(k1+4*k2+k3)/6
    return y

# Runge-Kutta 4th order
def RK4(f, x, y0):
    y = np.zeros(x.shape)
    y[0] = y0
    
    for i in range(len(x)-1):
        h = x[i+1]-x[i]
        
        k1 = h*f(x[i],y[i])
        k2 = h*f(x[i]+0.5*h  ,y[i]+0.5*k1)
        k3 = h*f(x[i]+0.5*h  ,y[i]+0.5*k2)
        k4 = h*f(x[i]+    h  ,y[i]    +k3)
    
        y[i+1] = y[i]+(1/6)*(k1+2*k2+2*k3+k4)
        
    return y

# functions definition
f   = lambda x,y:-x*y           
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
rk2 = np.zeros(N.shape)
rk3 = np.zeros(N.shape)
rk4 = np.zeros(N.shape)

# iterate through number of intervals. Q: how does error change with changing N (H)
for i in range(len(N)):
    x = np.linspace(a,b,N[i])
    
    rk2[i] = RK2(f,x,y0)[-1] # why the [-1]? /joel
    rk3[i] = RK3(f,x,y0)[-1]
    rk4[i] = RK4(f,x,y0)[-1]


# plot abs of errors
plt.figure(1)
plt.xlabel("Size of interval h")
plt.ylabel("Error")
plt.loglog(H,abs(rk2-y_true(b)),label="RK2")
plt.loglog(H,abs(rk3-y_true(b)),label="RK3")
plt.loglog(H,abs(rk4-y_true(b)),label="RK4")
plt.legend()
plt.show()

# larger grid
b=3

# re-initialize vectors
rk2 = np.zeros(N.shape)
rk3 = np.zeros(N.shape)
rk4 = np.zeros(N.shape)

for i in range(len(N)):
    x=np.linspace(a,b,N[i])
    
    rk2[i] = RK2(f,x,y0)[-1]
    rk3[i] = RK3(f,x,y0)[-1]
    rk4[i] = RK4(f,x,y0)[-1] # why the [-1]? /joel

