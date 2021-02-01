#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 29 08:45:27 2021

@author: joelbjervig
"""

from matplotlib import pyplot as plt
import math as m
import numpy as np

def f(x):
    return m.exp(x)

def trapz(f,a,b,n):
    h = (b-a)/n
    k = 0.0
    
    x = a+h
    for i in range(1,n-1):
        k+=f(x)
        x+=h
        
    return h*(0.5*(f(a)+f(b))+k)

def simpsons(f,a,b,N):
    
    h=(b-a)/N
    k=0.0
    
    x=a+h
    for i in range(1,int(N/2+1)):
        x += 2*h
        k += 4*f(x)
        
    x = a+2*h
    for i in range(2,int((N/2)-1)):
        x += 2*h
        k += 2*f(x)

    return (h/3)*(f(a)+f(b)+k)


def bode(f,a,b,N):
    h = (b-a)/N

    vector = np.arange(a,b-h,4*h)
    # print(vector)
    sum = 0
    
    for x in vector:
        sum = sum + 2*h/45*(7*f(x) + 32*f(x+h) + 12*f(x+2*h)+ 32*f(x+3*h) + 7*f(x+4*h))
    return sum



a = 0
b = 1

N = np.array([4, 40, 400, 4000, 40000])
length = np.size(N)

trapz_e = np.empty(length)
simps_e = np.empty(length)
bode_e = np.empty(length)

for i in range(0,length):
    trapz_e[i] =   abs((f(1)-1)-trapz(f,a,b,N[i]))
    simps_e[i] =   abs((f(1)-1)-simpsons(f,a,b,N[i]))
    bode_e[i] =    abs((f(1)-1)-bode(f,a,b,N[i]))
    
#plot the result
plt.loglog(N,trapz_e, label = 'Trapezodial rule');
plt.loglog(N,simps_e, label = 'Simpsons rule');
plt.loglog(N,bode_e, label = 'Bode rule');

plt.legend()
plt.show()

for i in range(1,4):
    print(i)