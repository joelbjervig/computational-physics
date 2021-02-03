#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  3 10:01:20 2021

@author: Joel Bjervig, Thomas Herrad, Anton Palm
"""
<<<<<<< Updated upstream

from math import *
from library import bode

# ODE
dydt = lambda p,p
dpdt = lambda y, -4*pi*y


# exacto solution of problem
p_exact = lambda t, c1*cos(2*pi*t)- 2*pi*c2*sin(2*pi*t)
y_exact = lambda t, c2*cos(2*pi*t) + c1/(2*pi)*sin(2*pi*t)

=======

from math import *
from library import bode
from library import RK2
from library import RK3
from library import RK4





=======
# ODE
dydt = lambda p,p
dpdt = lambda y, -4*pi*y
y0=1    # initial condition just trying omething


# exacto solution of problem
p_exact = lambda t, c1*cos(2*pi*t)- 2*pi*c2*sin(2*pi*t)
y_exact = lambda t, c2*cos(2*pi*t) + c1/(2*pi)*sin(2*pi*t)

a=0
b=1
N=np.arange(1,1000,10)  # num internals. Num points is N+1


# initialize vectors
rk2 = np.zeros(N.shape)
rk3 = np.zeros(N.shape)
rk4 = np.zeros(N.shape)

# iterate through number of intervals. Q: how does error change with changing N (H)
for i in range(len(N)):
    t = np.linspace(a,b,N[i])
    
    rk2[i] = RK2(f,t,y0)[-1] # why the [-1]? /joel
    rk3[i] = RK3(f,t,y0)[-1]
    rk4[i] = RK4(f,t,y0)[-1]


>>>>>>> Stashed changes
