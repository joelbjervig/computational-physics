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




