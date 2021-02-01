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

def search(f,x0,dx,accu):   # guess trial value smaller than solution
    n = 0   #number of iterations

    while(abs(f(x0))>accu):     # while value is larger or equal to than accuracy 
        print("x value: ", x0, "function value: f(x) = ", f(x0))
        
        x0 = x0+dx          # increase x by smal value.
        if(f(x0+dx)>=0):    # if f(x) has changes sign,
            x0 = x0-dx      # take a step back
            
            dx = 0.5*dx     
            x0 = x0+dx      # and take a half step forward
        n+=1            # increace iteration counter
        
    print('Found solution after',n,'iterations.')
    return x0

def newtonraphson(f,df,x0,root,accu,max_iter):
    xn = x0 
    for n in range(0,max_iter):

        fxn = f(xn)
        Dfxn = df(xn)
        
        print('Error Value:', root-xn)
        if abs(abs(fxn) < accu):
            print('Found solution after',n,'iterations.')
            return x0
        
        if Dfxn == 0:
            Dfxn = accu;
            
        xn = xn - fxn/Dfxn
        
    print('Exceeded maximum iterations. No solution found.')
    return None

def secant(f,x0,x1,root,accu,max_iter):
    
    for n in range(0,max_iter):

        
        if f(x0) == f(x1):
            print('Divide by zero error!')
            break
        
        # same as newton but uses approximate derivative
        x2 = x0 - f(x0)*(x0-x1)/(f(x0)-f(x1)) 
        
        print('Error Value:', root-x2)
        if abs(f(x2)) < accu: #originally fxn instead of root-xn
            print('Found solution after',n,'iterations.')
            return x2
        
        # update points
        x0 = x1
        x1 = x2
        
    print('Exceeded maximum iterations. No solution found.')
    return None

x0 = 0.2            # initial guess
x1 = 0.1            # second initial guess (for secant method)
dx = 0.01           # steplength
eps = 10**-10       # precission
max_iter = 100      # max iteration
root = m.sqrt(5)    # positive root

newtonraphson(f,df,x0,root,eps,max_iter)
search(f,x0,dx,eps)
secant(f,x0,x1,root,eps,max_iter)