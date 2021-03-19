# -*- coding: utf-8 -*-
"""
Created on Wed Mar 17 15:34:45 2021

@author: anton
"""
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  5 13:33:24 2021

@author: anton
"""
import numpy as np
import math
import matplotlib.pyplot as plt

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


# first integrand term for b<r<r_max
theta1 = lambda r: 2*b*(r**(-2)*1/(np.sqrt(1-b**2/r**2)))

# numberical solution
def theta(b,E,rmax,rmin):
    return np.subtract(bode(theta1,b+1e-6,rmax,N), bode(theta2,rmin+1e-6,rmax,N))

# part 2
# LJ dependet on E
V=1
a=3.333
b=1
rmax=3*a
N=100001  
# Energy is varying
Evec = V*np.linspace(0.1,100,20)

# Leonard Jones potential
pot = lambda r: 4*V*((a/r)**12-(a/r)**6)

rminfunc = lambda r,E: 1-b**2/r**2-pot(r)/E
rmin = np.zeros(len(Evec))
sol = np.zeros(len(Evec))

 # Find the positive root of the function
for i in range(0,len(Evec)):
    rminfuncE = lambda r: 1-b**2/r**2-pot(r)/Evec[i]
    rmin[i] = bissection(0.1, rmax, rminfuncE, 100)
    E=Evec[i]
    theta2 = lambda r: 2*b*(r**(-2)*1/(np.sqrt(1-b**2/r**2-pot(r)/E)))
    sol[i]= theta(b,E, rmax, rmin[i])
    
  
    
# part 2 cross section
bvec=np.linspace(0.1,6.87368,20)
h=bvec[1]-bvec[0]
E=2   
 
theta2 = lambda r: 2*b*(r**(-2)*1/(np.sqrt(1-b**2/r**2-pot(r)/E)))
def theta(b,rmax,rmin):
    return np.subtract(bode(theta1,b+1e-6,rmax,N), bode(theta2,rmin+1e-6,rmax,N))


rminLJ=np.zeros(len(bvec))
solLJ=np.zeros(len(bvec))
for i in range(0,len(bvec)):
    b=bvec[i]
    rminfuncLJ = lambda r: 1-b**2/r**2-pot(r)/E
    rminLJ[i] = bissection(0.1, rmax, rminfuncLJ, 100)
    solLJ[i]=theta(b,rmax,rminLJ[i])
   
thetaPrimeLJ=np.zeros(len(bvec)-1)
crossSecLJ=np.zeros(len(bvec)-1)
for j in range(1,len(bvec)):
    thetaPrimeLJ[j-1]=((solLJ[j]-solLJ[j-1])/((h)))
    crossSecLJ[j-1]=(bvec[j])/(math.sin(solLJ[j])*abs(thetaPrimeLJ[j-1]))
    

    
plt.figure(1)
plt.plot(Evec,sol,label="Numerical solution",marker='*')
plt.xlabel("Energy E")
plt.ylabel('Deflection angle [radians]')
plt.legend(loc='best')
plt.show()

plt.figure(2)
plt.plot(bvec,solLJ,label="Numerical solution",marker='*')
plt.xlabel("b")
plt.ylabel('Deflection angle [radians]')
plt.legend(loc='best')
plt.show()

plt.figure(3)
plt.plot(bvec[1:100],crossSecLJ,label="crossSectionLJ",marker='*')
plt.xlabel("b")
plt.ylabel('')
plt.legend(loc='best')
plt.show()

plt.figure(47)
plt.plot(bvec[1:100],thetaPrimeLJ,label="thetaPrimeLJ",marker='*')
plt.xlabel("b")
plt.ylabel('')
plt.legend(loc='best')
plt.show()

