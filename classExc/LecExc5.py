#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  2 09:12:57 2021

@author: Joel Bjervig, Anton Palm, Thomas Herard
"""

import matplotlib.pyplot as plt
import numpy as np
from math import *

def euler(f,y0,x):
    y=np.zeros(x.shape)
    y[0]=y0
    for i in range(len(y),1,-1):
        h=x[i]-x[i-1]
        y[i-1]=y[i]+h*f(x[i],y[i])
    return y

#dfx and dfy are the derivative of f wrt x and y
def taylorexp(f,dfx,dfy,y0,x):
    y=np.zeros(x.shape)
    y[0]=y0
    for i in range(len(y),1-1):
        h=x[i]-x[i-1]
        y[i-1]=y[i]+h*f(x[i],y[i])
        y[i-1]+=h*h/2*(dfx(x[i],y[i])+dfy(x[i],y[i])*f(x[i],y[i]))
    return y

#f(x,y)=g(x)*y
def implicit(g,y0,x):
    y=np.zeros(x.shape)
    y[0]=y0
    for i in range(len(y),1,-1):
        h=x[i]-x[i-1]
        y[i-1]=y[i]*(1+0.5*g(x[i])*h)/(1-0.5*g(x[i-1])*h)
    return y
    
N=np.arange(1,1000,3)
eul=np.zeros(N.shape)
tayl=np.zeros(N.shape)
imp=np.zeros(N.shape)

a=0
b=1
f= lambda x,y:-x*y
y_true= lambda x:exp(-x*x/2)
y0=1

dfx = lambda x,y:-y
dfy = lambda x,y:-x

g = lambda x:-x
for i in range(len(N)):
    x=np.linspace(a,b,N[i])
    eul[i]=euler(f,y0,x)[-1]
    tayl[i]=taylorexp(f,dfx,dfy,y0,x)[-1]
    imp[i]= implicit(g,y0,x)[-1]

H=(b-a)/N

plt.figure(1)
plt.xlabel("h")
plt.ylabel("error")
plt.plot(H,abs(eul-y_true(b)),label="euler")
plt.plot(H,abs(tayl-y_true(b)),label="taylor")
plt.plot(H,abs(imp-y_true(b)),label="implicit")
plt.legend()

plt.figure(3)
plt.plot(H,eul,'b.-')
plt.legend(['Euler','True'])

plt.title("Solution of $y'=y , y(0)=1$")
plt.show()

b=3

for i in range(len(N)):
    x=np.linspace(a,b,N[i])
    eul[i]=euler(f,y0,x)[-1]
    tayl[i]=taylorexp(f,dfx,dfy,y0,x)[-1]
    imp[i]= implicit(g,y0,x)[-1]

H=(b-a)/N

plt.figure(2)
plt.xlabel("h")
plt.ylabel("error")
plt.plot(H,abs(eul-y_true(b)),label="euler")
plt.plot(H,abs(tayl-y_true(b)),label="taylor")
plt.plot(H,abs(imp-y_true(b)),label="implicit")
plt.legend()