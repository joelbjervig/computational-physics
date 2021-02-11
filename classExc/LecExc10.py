# -*- coding: utf-8 -*-
"""
Created on Mon Feb  8 14:41:07 2021

@author: anton
"""
import numpy as np
from math import *
import matplotlib.pyplot as plt


# Legendre polynomial 
def P(x,n):  
    if(n == 0): 
        return 1 # P0 = 1 
    elif(n == 1): 
        return x # P1 = x 
    else: 
        return ((2*n-1)*x*P(x, n-1)-(n-1)*P(x, n-2))/(n)

# Derivative of Legendre Polynomial
def P_prime(x,n):
    if(n==0):
        return 0
    elif(n==1):
        return 1
    else:
        return (-n*x*P(x,n)+n*P(x,n-1))/(1-x**2)
    
#bissection function finding one root in one specific interval
def bissection(a,b,f,maxit):
    if f(a)*f(b) >= 0:
        print("Doesn't change sign on interval")
        return None
    c=(a+b)/2
    for i in range(maxit):
        c=(a+b)/2
        if f(c)*f(a)<0:
            b=c
        elif f(c)*f(b)<0:
            a=c
        elif f(c)==0.0:
            print("Exact solution found.")
            return c
    return (a+b)/2

def multipleroots(a,b,f,numroots):
    for scale in range(10,1000,10):
        X=np.linspace(a,b,numroots*scale)
        x0=np.zeros((numroots,))
        j=0
        for i in range(1,len(X)):
            if f(X[i-1])*f(X[i])<0:
                x0[j]=bissection(X[i-1],X[i],f,1000)
                j+=1
            if f(X[i-1])*f(X[i]) ==0:
                if f(X[i-1]) == 0:
                    x0[j]=X[i-1]
                    j+=1
                else:
                    x0[j]=X[i]
                    j+=1
        if j==numroots:
            print(scale)
            return x0
    print("Number of roots not consistent with function provided.")
    return None

x=np.arange(-1,1,0.001)
plt.figure(1)
plt.title("Legendre Polynomials")
plt.xlabel("x")
plt.ylabel("y")
plt.plot(x,P(x,3),label="P_3")
plt.plot(x,P(x,4),label="P_4")
plt.plot(x,P(x,5),label="P_5")
plt.legend()
plt.show()

plt.figure(2)
plt.title("Legendre Polynomials Derivatives")
plt.xlabel("x")
plt.ylabel("y")
plt.plot(x,P_prime(x,3),label="P_3_prime")
plt.plot(x,P_prime(x,4),label="P_4_prime")
plt.plot(x,P_prime(x,5),label="P_5_prime")
plt.legend()
plt.show()


a=-1
b=1
numroots=3

P3 = lambda x:P(x,3)
roots=multipleroots(a,b,P3,3)
print(roots)

