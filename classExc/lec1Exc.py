#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 28 16:30:02 2021

@author: joelbjervig
"""
import numpy as np
from matplotlib import pyplot as plt

# spacing
h = np.array([0.1, 0.01, 0.001, 0.0001, 0.00001])


# create x datapoints
xc = np.array([1, 1, 1, 1, 1])
x = np.array([ xc-2*h, xc-1*h, xc, xc+1*h, xc+2*h ])
# print matrix
print("THIS IS X")
print(x);

# define sine function
f = np.sin(x)
fprime = np.cos(x)

# print matrices
#print(f)
print("THIS IS F' at x=1")
print(fprime[2])

def threepoint(A):
    return (A[3]-A[1])/(2*h)

def forward(A):
    return (A[3]-A[2])/(h)

def backward(A):
    return (A[2]-A[1])/(h)

def fivepoint(A):
    return (A[0]-8*A[1]+8*A[3]-A[4])/(12*h)


tp = threepoint(f)
fw = forward(f)
bw = backward(f)
fp = fivepoint(f)

print("THESE ARE THE APPROXIMATIONS")
print(tp)
print(fw)
print(bw)
print(fp)

e_tp = abs(tp-fprime[2])
e_fw = abs(fw-fprime[2])
e_bw = abs(bw-fprime[2])
e_fp = abs(fp-fprime[2])

print("THESE ARE THE ERRORS")
print(e_tp)
print(e_fw)
print(e_bw)
print(e_fp)

#plot the result

plt.loglog(h,e_tp, label = 'Three point');
plt.loglog(h,e_fw, label = 'Forward');
plt.loglog(h,e_bw, label = 'Backward');
plt.loglog(h,e_fp, label = 'Five point');
plt.legend()

plt.show()