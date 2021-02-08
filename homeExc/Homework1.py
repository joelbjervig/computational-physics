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

fig = plt.figure()
ax = fig.gca(projection='3d')


# exacto solution of problem
c1 = 0; # corresponds to 
c2 = 1/(2*pi); # y0 and p0

p_exact = lambda t: c2*2*pi*np.cos(2*pi*t)- 2*pi*c1*np.sin(2*pi*t)
y_exact = lambda t: c1*np.cos(2*pi*t) + c2*np.sin(2*pi*t)

# ODE
f = lambda t,y,p:p
g = lambda t,y,p: -4*pi**2*y
y0=0    # initial condition just trying omething
p0=1

a=0
b=2*pi

l = 1       # min num of intervals
u = 10000    # max num of intervals
p = 100      # number of partitions
N=np.arange(l,u,p)  # num internals. Num points is N+1
H=(b-a)/N

# initialize vectors
size = len(N)-1
t = np.linspace(a,b,N[size])

y2,p2 = RK2_SYS(f,g,y0,p0,t)
y3,p3 = RK3_SYS(f,g,y0,p0,t)
y4,p4 = RK4_SYS(f,g,y0,p0,t)



y_true = y_exact(t);
p_true = p_exact(t);


print("l = %1.3f" %l)
print("u = %1.3f" %u)
print("p = %1.3f" %p)
print("length of N = %1.3f" %len(N))
print("N array:")
print(N)


# plot analytical solution
plt.figure(2)
plt.title("Analytical solution, steplength h = %1.3f" %H[size])
plt.xlabel("Time t")
plt.ylabel("Amplitude")
plt.plot(t,y_true,label="Displacement")
plt.plot(t,p_true,label="Momentum")
plt.legend()
plt.show()

#plot RK 2,3,4 approximative solutions
plt.figure(2)
plt.title("Solutions by Runge Kutta Methods of order 2,3,4. Steplength h = %1.3f" %H[size])
plt.xlabel("Time t")
plt.ylabel("Amplitude")
plt.plot(t,y2,label="RK2: Displacement")
plt.plot(t,p2,label="RK2: Momentum")
plt.plot(t,y3,label="RK3: Displacement")
plt.plot(t,p3,label="RK3: Momentum")
plt.plot(t,y4,label="RK4: Displacement")
plt.plot(t,p4,label="RK4: Momentum")
plt.legend()
plt.show()

y4h = np.zeros((len(N), len(t)))
p4h = np.zeros((len(N), len(t)))

# iterate through number of intervals. Q: how does error change with changing N (H)
for i in range(len(N)):
    t = np.linspace(a,b,N[i])

    y4h[i] = RK4_SYS(f,g,y0,p0,t)[0]
    p4h[i] = RK4_SYS(f,g,y0,p0,t)[1]
    
    
"""
# plot abs of errors
plt.figure(1)
plt.title("Error of RK4 method: limits a = {a}, b = {b}π,. Number of intervals: Min: I_min = {l}, max: I_max = {u}, num partitions p = {p}".format(a= a, b = int(b/pi), l = l, u = u, p = p))
plt.xlabel("Size of interval h")
plt.ylabel("Error")
plt.loglog(H,abs(y4h-y_exact(b)),label="Displacement")
plt.loglog(H,abs(p4h-p_exact(b)),label="Momentum")
plt.legend()
plt.show()
"""

H, t = np.meshgrid(H, t)
surf = ax.plot_surface(H,t,abs(y4h-y_exact(t)))#, cmap=cm.coolwarm, linewidth=0, antialiased=False)





