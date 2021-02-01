#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 29 15:31:28 2021

@author: joelbjervig
"""
import math as m

def f(x):
    return x*x-5
def df(x):
    return 2*x

def search(f,x0,dx,stop):

    if(f(x0+dx)>=0):
        x0 = x0/2
    else:
        x0 = x0+dx
    return x0

def newtonraphson(f,Df,x0,root,acc,max_iter):
    xn = x0
    i = 0
    while abs(f(x)>acc) and i<max_iter:
        x = 
        
        
    for n in range(0,max_iter):
        fxn = f(xn)
        if abs(fxn) < epsilon:
            print('Found solution after',n,'iterations.')
            return xn
        Dfxn = Df(xn)
        if Dfxn == 0:
            Dfxn = rpsilon;
        print('Error Value:', root-x0)
        xn = xn - fxn/Dfxn
    print('Exceeded maximum iterations. No solution found.')
    return None

def secant(f,x0):
    return

x0 = 0.2            # initial guess
dx = 0.1            # steplength
eps = 10^(-10)      # precission
max_iter = 100      # max iteration
root = m.sqrt(5)    # positive root

newton(f,df,x0,root,eps,max_iter)