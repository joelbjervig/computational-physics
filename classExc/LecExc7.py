#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  4 15:50:03 2021

@author: Joel Bjervig, Anton Palm, Thomas Herard
"""

import matplotlib.pyplot as plt
import numpy as np
from math import *

def numerov(x,y0,y1,kfun,sfun):
    y=np.zeros(x.shape)
    y[0]=y0
    y[1]=y1
    
    for i in range(2,len(x)):
        h=x[i]-x[i-1]
        
        a=1+h**2*kfun(x[i])**2/12
        b=1-5*h**2*kfun(x[i-1])**2/12
        c=1+h**2*kfun(x[i-2])**2/12
        d=h**2*(sfun(x[i])+10*sfun(x[i-1])+sfun(x[i-2]))/12
        
        y[i]=2*b*y[i-1]/a - c*y[i-2]/a + d/a
        
        
    return y

  sfun = lambda x:-0.5*x*np.exp(-x)
kfun= lambda x:0

def phitrue(x):
    return 1-0.5*(x+2)*np.exp(-x)

r=np.linspace(0,30,10000)

phi1=phitrue(r[1]-r[0])
phi1=1.05*phi1

phi=numerov(r,0,phi1,kfun,sfun)


plt.figure(1)
plt.plot(r,phi,label="numerical")
plt.plot(r,phitrue(r),label="exact")
plt.legend()
