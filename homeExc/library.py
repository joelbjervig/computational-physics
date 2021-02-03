#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  3 10:02:52 2021

@author: Joel Bjervig, Thomas Herrad, Anton Palm Ekspong
"""
import numpy as np
from math import *
###############################
######### ODE SOLVERS #########
###############################

# Runge-Kutta 2nd order
def RK2(f, x, y0):
    y = np.zeros(x.shape)
    y[0] = y0
    
    for i in range(len(x) - 1):
        h = x[i+1] - x[i]
        
        k = h*f(y[i], x[i])
        
        y[i+1] = y[i] + h*f(x[i]+H/2, y[i]+k/2)
        
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