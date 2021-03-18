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

# part 1

# def real(b):
#     retrun 2*(np.arcsin(b/rmax) -np.pi/2-np.arccos(b/rmax*np.sqrt(E/(E-V))))


def theta(b,rmax,rmin):
    return np.subtract(bode(theta1,b+1e-6,rmax,N), bode(theta2,rmin+1e-6,rmax,N))

N=100001 # 
E=2 # 
V=1
rmax=10

# first integrand term for b<r<r_max
theta1 = lambda r: 2*b*(r**(-2)*1/(np.sqrt(1-b**2/r**2)))
# second integrand term for r_min<r<r_max
theta2 = lambda r: 2*b*(r**(-2)*1/(np.sqrt(1-b**2/r**2-V/E)))

# def theta(b,rmax,rmin):
#     return np.subtract(bode(theta1,b+1e-6,rmax,N), bode(theta2,rmin+1e-6,rmax,N))

analytical = lambda b: 2*(np.arccos(b/rmax)-np.arccos(b/rmax*(E/(E-V))**0.5))

# def real(b):
#     return 2*(np.arccos(b/rmax)-np.arccos(b/rmax*(E/(E-V))**0.5))
# real=np.vectorize(real)

bvec=np.linspace(0.1,rmax,20)
rmin=bvec*np.sqrt(E/(E-V))
its=len(rmin[rmin<rmax])
bvec=bvec[0:its]
rmin=rmin[0:its]
sol=[]
error=[]
thetaPrime=[]
crossSec=[]
h=bvec[1]-bvec[0]

for i in range(0,its):
    b=bvec[i]
    sol.append(theta(b,rmax,rmin[i]))
    error.append(abs(sol[i]-analytical(b)))
    
# part 2 cross section
    
for j in range(1,its):
    thetaPrime.append((sol[j]-sol[j-1])/((h)))
    crossSec.append(bvec[j]/(math.sin(sol[j]))*thetaPrime[j-1])
    

# plt.figure(1)
# plt.plot(bvec,analytical(bvec),label="analytical",marker='o')
# plt.plot(bvec,sol,label="numerical",marker='*')
# plt.xlabel("b")
# plt.ylabel('Angle in radians')
# plt.legend(loc='best')
# plt.show()


# plt.figure(2)
# plt.plot(bvec,error,label="error",marker='*')
# plt.xlabel("b")
# plt.ylabel('Angle in radians')
# plt.legend(loc='best')
# plt.show()

# plt.figure(3)
# # plt.plot(bvec[1:its],thetaPrime,label="thetaPrime",marker='*')
# plt.plot(bvec[1:its],crossSec,label="crossSection",marker='*')
# plt.xlabel("b")
# plt.ylabel('')
# plt.legend(loc='best')
# plt.show()


# part 2
# LJ dependet on E
V=1
a=3.333
b=1
rmax=3*a

# Energy is varying
Evec = V*np.linspace(0.1,100,20)

# Leonard Jones potential
pot = lambda r: 4*V*((a/r)**12-(a/r)**6)

rminfunc = lambda r,E: 1-b**2/r**2-pot(r)/E
rmin = np.zeros(len(Evec))

# Find the positive root of the function
for i in range(0,len(Evec)):
    rminfuncE = lambda r: 1-b**2/r**2-pot(r)/Evec[i]
    rmin[i] = bissection(0.1, rmax, rminfuncE, 100)
    
    
# numberical solution
def theta(b,E,rmax,rmin):
    return np.subtract(bode(theta1,b+1e-6,rmax,N), bode(theta2,rmin+1e-6,rmax,N))

sol = np.zeros(len(Evec))


for i in range(0,len(Evec)):
    E=Evec[i]
    theta2 = lambda r: 2*b*(r**(-2)*1/(np.sqrt(1-b**2/r**2-pot(r)/E)))
    sol[i]= theta(b,E, rmax, rmin[i])
  
########################
# part 2 cross section #
########################
L = 100
bvec=np.linspace(0.1,6.8,L)
E=2    
theta2 = lambda r: 2*b*(r**(-2)*1/(np.sqrt(1-b**2/r**2-pot(r)/E)))
def theta(b,rmax,rmin):
    return np.subtract(bode(theta1,b+1e-6,rmax,N), bode(theta2,rmin+1e-6,rmax,N))

# rminfuncB = lambda r,b: 1-b**2/r**2-pot(r)/E
rminLJ=np.zeros(len(bvec))
solLJ=np.zeros(len(bvec))
for i in range(0,len(bvec)):
    b=bvec[i]
    rminfuncLJ = lambda r: 1-b**2/r**2-pot(r)/E
    rminLJ[i] = bissection(0.1, rmax, rminfuncLJ, L)
    solLJ[i]=theta(b,rmax,rminLJ[i])
   
thetaPrimeLJ=np.zeros(len(bvec)-1)
crossSecLJ=np.zeros(len(bvec)-1)
for j in range(1,len(bvec)):
    thetaPrimeLJ[j-1]=((solLJ[j]-solLJ[j-1])/((h)))
    crossSecLJ[j-1]=(bvec[j]/(math.sin(solLJ[j]))*thetaPrimeLJ[j-1])
    

    
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
plt.plot(bvec[1:L],thetaPrimeLJ,label="thetaPrimeLJ",marker='*')
plt.plot(bvec[1:L],crossSecLJ,label="crossSectionLJ",marker='*')
plt.xlabel("b")
plt.ylabel('')
plt.legend(loc='best')
plt.show()