#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  3 10:01:20 2021

@author: Joel Bjervig, Thomas Herrad, Anton Palm Ekspong
"""
import numpy as np
from math import *
import matplotlib.pyplot as plt
from library import RK2_SYS
from library import RK3_SYS
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
t = np.linspace(a,b,N[50])

p2 = RK2_SYS(f,g,y0,p0,t)[0]
y2 = RK2_SYS(f,g,y0,p0,t)[1]

p3 = RK3_SYS(f,g,y0,p0,t)[0]
y3 = RK3_SYS(f,g,y0,p0,t)[1]

p4 = RK4_SYS(f,g,y0,p0,t)[0]
y4 = RK4_SYS(f,g,y0,p0,t)[1]



plt.figure(1)
plt.title("Solution to ODE using Runge-Kutta of order 2")
plt.xlabel("Time t")
plt.ylabel("Amplitude")
plt.plot(t,p2,label="Momentum p")
plt.plot(t,y2,label="Displacement y")
plt.legend()
plt.show()

plt.figure(2)
plt.title("Solution to ODE using Runge-Kutta of order 3")
plt.xlabel("Time t")
plt.ylabel("Amplitude")
plt.plot(t,p3,label="Momentum p")
plt.plot(t,y3,label="Displacement y")
plt.legend()
plt.show()

plt.figure(3)
plt.title("Solution to ODE using Runge-Kutta of order 4")
plt.xlabel("Time t")
plt.ylabel("Amplitude")
plt.plot(t,p4,label="Momentum p")
plt.plot(t,y4,label="Displacement y")
plt.legend()
plt.show()

# how large is the difference between rk3 and rk4?
plt.figure(4)
plt.title("Abosulute difference between Runge-Kutta methods of order four and three")
plt.xlabel("Time t")
plt.ylabel("Amplitude")
plt.plot(t,abs(p4-p3),label="abs(rk4-rk3) momentum")
plt.plot(t,abs(y4-y3),label="abs(rk4-rk3) displacement")
plt.legend()
plt.show()

# # iterate through number of intervals. Q: how does error change with changing N (H)
# for i in range(len(N)):
#     t = np.linspace(a,b,N[i])
    

#     y[i] = RK4_SYS(f,g,y0,p0,t)[0]


