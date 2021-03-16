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
    N=latt.shape[0]
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

def ising(N,T):
    kb=1
    J=1      
    equilibrium = 1000 
    sweeps = 1000

    E,M,Cb,Khi,U = np.zeros(T.shape), np.zeros(T.shape), np.zeros(T.shape), np.zeros(T.shape), np.zeros(T.shape)
    

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
        U[n]=1-m4/(3*m2**2)
        
    return E,M,Cb,Khi,U

T= np.linspace(1, 5, 30)
N=np.array([8,16,32])

E = np.zeros((len(N),len(T)))
M,Cb,Khi,U = np.zeros(E.shape), np.zeros(E.shape), np.zeros(E.shape), np.zeros(E.shape)

for k in range(len(N)):
    E[k,:],M[k,:],Cb[k,:],Khi[k,:],U[k,:]=ising(N[k],T)     
 
##############
#Plots
##############


plt.figure(1)
for i in range(len(N)):
    plt.plot(T,abs(M[i,:]),label="N="+str(N[i]))
plt.legend()
plt.xlabel("Temperature (T)") 
plt.ylabel("Order parameter")
plt.savefig("Magnetization.jpg")
plt.title("sweeps = equilibriate sweeps ="+str(sweeps))

plt.figure(2)
for i in range(len(N)):
    plt.plot(T,Khi[i,:],label="N="+str(N[i]))
plt.legend()
plt.xlabel("Temperature (T)") 
plt.ylabel("Susceptibility (Khi)")
plt.title("sweeps = equilibriate sweeps ="+str(sweeps))
plt.savefig("Susceptibility.jpg")

plt.figure(3)
for i in range(len(N)):
    plt.plot(T,Cb[i,:],label="N="+str(N[i]))
plt.legend()
plt.xlabel("Temperature (T)") 
plt.ylabel("Specific heat (Cb)")
plt.title("sweeps = equilibriate sweeps ="+str(sweeps))
plt.savefig("Specheat.jpg")

plt.figure(4)
for i in range(len(N)):
    plt.plot(T,U[i,:],label="N="+str(N[i]))
plt.legend()
plt.xlabel("Temperature (T)") 
plt.ylabel("Fourth order cumulant")
plt.title("sweeps = equilibriate sweeps ="+str(sweeps))
plt.savefig("Cumulant.jpg")


