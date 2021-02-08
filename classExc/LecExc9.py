#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb  8 14:16:41 2021

@author: joelbjervig Thomas Herard
"""
import numpy as np
from math import *
import matplotlib.pyplot as plt

# guess trial eigenvalue as something small
epsilon1 = 10**(-8)
epsilon2 = 10**(-6)
#k_init = 0.001 
x0 = 0.1
dx = 0.1
y0 = 0
y1 = y0+epsilon

x = np.zeros(100)

f = lambda x:sin(pi*x)
k = lambda x:pi
S = lambda x:0


while( (abs(search(f,x0,dx,epsilon2)[0]-numerov(k,S,y0,y1,x)[-1]>epsilon2)) or (dx>epsilon1) )
    x0,dx = search(f,x0,dx,epsilon2)
    y = numerov(k,S,y0,y1,x)
