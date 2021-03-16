import numpy as np
from numpy.random import rand
import matplotlib.pyplot as plt

##############
#Functions
##############

def initial(N):
    lattice = np.random.choice([-1,1],(N,N))
    return lattice


def metropolis(latt,temp):
    for i in range(N):
        for j in range(N):
             
                neigh = latt[(i+1)%N,j] + latt[i,(j+1)%N] + latt[(i-1)%N,j] + latt[i,(j-1)%N]
                dE = 2*latt[i,j]*neigh
                if dE < 0:
                    latt[i, j] *= -1
                elif rand() < np.exp(-dE/(kb*temp)):
                    latt[i, j] *= -1
    return latt


def energy(latt):
    #vectorial computation of the S_alpha*S_beta for the closest neighbours
    Sup = np.roll(latt, -1, axis=0)
    Sdown = np.roll(latt, 1, axis=0)
    Sleft = np.roll(latt, 1, axis=1)
    Sright = np.roll(latt, -1, axis=1)
    
    sum_neighbours= Sup+Sdown+Sleft+Sright
    
    H=-np.sum(sum_neighbours*latt)
    return H/4


def magnetization(config):
    mag = np.sum(config)
    return mag

################
#Main Code
###############

kb=1
J=1    
N= 32       
equilibrium = 300 
sweeps = 300

T= np.linspace(1, 5, 30); 
E,M,Cb,Khi = np.zeros(T.shape), np.zeros(T.shape), np.zeros(T.shape), np.zeros(T.shape)
#U=np.zeros((3,len(T)))

norm1=sweeps*N**2
norm2=sweeps**2*N**2

for n in range(len(T)):
    e1 = m1 = e2 = m2 = m4 =0
    
    config = initial(N)
    
    for i in range(equilibrium): #loops to equilibrate the system
        metropolis(config, T[n]) #metropolis algorithm

    for i in range(sweeps): 
        metropolis(config, T[n])
        
        En = energy(config) 
        Mag = magnetization(config) 

        e1+=En
        m1+=Mag
        m2+=Mag**2 
        e2+=En**2
        m4+=Mag**4

    E[n] = e1/norm1
    M[n] = m1/norm1
    Cb[n] = (e2/norm1 - e1**2/norm2)/(kb*T[n]**2)
    Khi[n] = (m2/norm1 - m1**2/norm2)/(kb*T[n])
    Ul=1-m4/(norm1**4*3*m2**2/norm2)
    if N==8:
        U[0,n]=Ul
    elif N==16:
        U[1,n]=Ul
    elif N==32:
        U[2,n]=Ul
        
        
##############
#Plots
##############


plt.figure(1)
plt.plot(T,abs(M),"co")
plt.xlabel("Temperature (T)") 
plt.ylabel("Order parameter")
plt.title("N="+str(N)+", sweeps = equilibriate sweeps ="+str(sweeps))
plt.savefig("Magnetization"+str(N)+".jpg")

plt.figure(2)
plt.plot(T,Khi,"bo")
plt.xlabel("Temperature (T)") 
plt.ylabel("Susceptibility (Khi)")
plt.title("N="+str(N)+", sweeps = equilibriate sweeps ="+str(sweeps))
plt.savefig("Susceptibility"+str(N)+".jpg")

plt.figure(3)
plt.plot(T,Cb,"ro")
plt.xlabel("Temperature (T)") 
plt.ylabel("Specific heat (Cb)")
plt.title("N="+str(N)+", sweeps = equilibriate sweeps ="+str(sweeps))
plt.savefig("Specheat"+str(N)+".jpg")


plt.figure(4)
plt.plot(T,U[0,:],"y-",label="N=8")
plt.plot(T,U[1,:],"g-",label="N=16")
plt.plot(T,U[2,:],"b-",label="N=32")
plt.legend()
plt.xlabel("Temperature (T)") 
plt.ylabel("Fourth order cumulant (Ul)")


