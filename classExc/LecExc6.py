#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  2 09:25:41 2021

@author: Joel Bjervig, Anton Palm, Thomas Herard
"""
import matplotlib.pyplot as plt
import numpy as np
from math import *

def rk2(f, y0, x):
    
    n = len(x)
    y=np.zeros(x.shape)
    y[0]=y0
    for i in range( n - 1 ):
        h = x[i+1] - x[i]
        k= h*f(y[i], x[i]) 
        y[i+1] = y[i] + h*f(y[i]+k/2, x[i]+h/2)

    return y



a=0
b=1
f = lambda x,y:-x*y
y_true = lambda x:exp(-x*x/2)
y0=1


N=np.arange(1,1000,1)
rk=np.zeros(N.shape)
for i in range(len(N)):
    x=np.linspace(a,b,N[i])
    rk[i]=rk2(f,y0,x)[-1]
    

H=(b-a)/N

plt.figure(1)
plt.xlabel("h")
plt.ylabel("error")
plt.plot(H,abs(rk-y_true(b)),label="RK2")


b=3

rk=np.zeros(N.shape)
for i in range(len(N)):
    x=np.linspace(a,b,N[i])
    rk[i]=rk2(f,y0,x)[-1]
