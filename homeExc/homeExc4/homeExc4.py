#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 23 14:17:24 2021

@author: joelbjervig
"""
import numpy as np
import pylab
import random
from math import *

# number of steps to take
n = 100
N = 10
# length of step ie lattice resolution



def distance(x1,y1,x2,y2):
    return sqrt((x2-x1)**2+(y2-y1)**2)

def avrg(num):
    sum = 0
    for i in range(0,len(num)):
        sum += num[i]
    return sum/len(num)

def walk(n):
    for j in range(0,n):
        
        # x and y coorditnates for each step taken
        x = np.zeros(n)
        y = np.zeros(n)
        
        # filling coordinates with random steps
        # moving one step in the direction of either north south east or west
        
        for i in range(1,n):
            
            rand = random.randint(1,4)
            
            if rand == 1:
                x[i] = x[i-1] + 1
                y[i] = y[i-1]
            elif rand == 2:
                x[i] = x[i-1] - 1
                y[i] = y[i-1]
            elif rand == 3:
                x[i] = x[i-1]
                y[i] = y[i-1] + 1
            else:
                x[i] = x[i-1]
                y[i] = y[i-1] - 1
        
    return distance(x[0],y[0],x[-1],y[-1])

r = np.zeros(N)
r_avrg = np.zeros(n)

for i in range(1,n): # loop over number of steps taken
    for j in range(1,N): # loop N times to average out
        r[j] = walk(i)
    r_avrg[i] = avrg(r)

pylab.title("how it scales")
pylab.plot(r_avrg)
pylab.show()
"""
pylab.title("Random Walk ($n = " + str(n) + "$ steps)") 
pylab.plot(x,y)
pylab.show()
"""

