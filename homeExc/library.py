#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  3 10:02:52 2021

@author: Joel Bjervig, Thomas Herrad, Anton Palm Ekspong
"""
import numpy as np
from math import *



###############################
######### RUNGE-KUTTA #########
###############################

# Runge-Kutta 2nd order
def RK2(f, x, y0):
    y = np.zeros(x.shape)
    y[0] = y0
    
    for i in range(len(x) - 1):
        h = x[i+1] - x[i]
        
        k = h*f(y[i], x[i])
        
        y[i+1] = y[i] + h*f(x[i]+h/2, y[i]+k/2)
        
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



###############################
###### OTHER ODE SOLVERS ######
###############################

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

###############################
######### ROOT FINDING ########
###############################


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

###############################
######### ROOT FINDING ########
###############################



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