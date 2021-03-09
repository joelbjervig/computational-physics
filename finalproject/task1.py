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



from library import bode


# def real(b):
#     retrun 2*(np.arcsin(b/rmax) -np.pi/2-np.arccos(b/rmax*np.sqrt(E/(E-V))))


def theta(b,rmax,rmin):
    return np.subtract(bode(theta1,b+1e-5,rmax,N), bode(theta2,rmin+1e-5,rmax,N))


N=100001 # 
E=2 # 
V=1
rmax=10

# first integrand term for b<r<r_max
theta1 = lambda r: 2*b*(r**(-2)*1/(np.sqrt(1-b**2/r**2)))
# second integrand term for r_min<r<r_max
theta2 = lambda r: 2*b*(r**(-2)*1/(np.sqrt(1-b**2/r**2-V/E)))

analytical = lambda b: 2*(np.arcsin(b/rmax) -np.pi/2-np.arccos(b/rmax*np.sqrt(E/(E-V))))

def real(b):
    return 2*(np.arccos(b/rmax)-np.arccos(b/rmax*(E/(E-V))**0.5))
real=np.vectorize(real)

bvec=np.linspace(0.1,rmax,20)
rmin=bvec*np.sqrt(E/(E-V))
its=len(rmin[rmin<rmax])
bvec=bvec[0:its]
rmin=rmin[0:its]

sol=np.zeros(len(rmin))
error=np.zeros(len(rmin))

for i in range(0,its):
    b = bvec[i]
    sol[i] = theta(b,rmax,rmin[i])

error=(abs(sol-real(bvec)))




plt.plot(bvec,real(bvec),label="Analytical solution",marker='o')
plt.plot(bvec,sol,label="Bode method",marker='*')
plt.xlabel("Impact paramenter b")
plt.ylabel('Scattering angle [radians]')
plt.legend(loc='best')
plt.show()


plt.figure()
plt.plot(bvec,error,label="Error",marker='*')
plt.xlabel("Impact parameter b")
plt.ylabel('Error [radians]')
plt.legend(loc='best')
plt.show()
