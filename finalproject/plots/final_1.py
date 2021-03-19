# -*- coding: utf-8 -*-
"""
Created on Fri Mar  5 13:33:24 2021

@author: anton
"""
import numpy as np
import math
import matplotlib.pyplot as plt

    
def bode(f,a,b,N):
    if (N-1)%4 != 0:
        print("N is not a multiple of the form 4p+1")
        return None
    s=0
    h=(b-a)/N
    for i in range(0,N-3,4):
        integ=(7*f(a+i*h)+32*f(a+(i+1)*h)+12*f(a+(i+2)*h)+32*f(a+(i+3)*h)+7*f(a+(i+4)*h))
        integ*=(2*h)/45
        s+=integ
    return s



def theta(b,rmax,rmin):
    return np.subtract(bode(theta1,b+1e-6,rmax,N), bode(theta2,rmin+1e-6,rmax,N))



N=100001 # 
E=2 # 
V=-1
rmax=10

# first integrand term for b<r<r_max
theta1 = lambda r: 2*b*(r**(-2)*1/(np.sqrt(1-b**2/r**2)))
# second integrand term for r_min<r<r_max
theta2 = lambda r: 2*b*(r**(-2)*1/(np.sqrt(1-b**2/r**2-V/E)))



analytical = lambda b: 2*(np.arccos(b/rmax)-np.arccos(b/rmax*(E/(E-V))**0.5))



bvec=np.linspace(0.1,rmax,20)
rmin=bvec*np.sqrt(E/(E-V))
its=len(rmin[rmin<rmax])
bvec=bvec[0:its]
rmin=rmin[0:its]
sol=np.zeros(len(bvec))
error=np.zeros(len(bvec))
thetaPrime=np.zeros(len(bvec)-1)
crossSec=np.zeros(len(bvec)-1)
h=bvec[1]-bvec[0]

for i in range(0,its):
    b=bvec[i]
    sol[i]=(theta(b,rmax,rmin[i]))
    error[i]=(abs(sol[i]-analytical(b)))
    
    
for j in range(1,its):
    thetaPrime[j-1]=(sol[j]-sol[j-1])/((h))
    crossSec[j-1]=(bvec[j])/(math.sin(sol[j])*abs(thetaPrime[j-1]))
    

plt.figure(1)
plt.plot(bvec,analytical(bvec),label="analytical",marker='o')
plt.plot(bvec,sol,label="numerical",marker='*')
plt.xlabel("b")
plt.ylabel('Angle in radians')
plt.legend(loc='best')
plt.show()


plt.figure(2)
plt.plot(bvec,error,label="error",marker='*')
plt.xlabel("b")
plt.ylabel('Angle in radians')
plt.legend(loc='best')
plt.show()

plt.figure(3)
# plt.plot(bvec[1:its],thetaPrime,label="thetaPrime",marker='*')
plt.plot(bvec[1:its],crossSec,label="crossSection",marker='*')
plt.xlabel("b")
plt.ylabel('')
plt.legend(loc='best')
plt.show()

plt.figure(4)
plt.plot(bvec[1:its],thetaPrime,label="thetaPrime",marker='*')
# plt.plot(bvec[1:its],crossSec,label="crossSection",marker='*')
plt.xlabel("b")
plt.ylabel('')
plt.legend(loc='best')
plt.show()
