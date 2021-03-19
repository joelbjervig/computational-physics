import numpy as np
from numpy.random import rand
import matplotlib.pyplot as plt

def initial(N):
    lattice = np.random.choice([-1,1],(N,N))
    return lattice


def heatbath(latt,temp):
    N=latt.shape[0]
    for i in range(N):
        for j in range(N):
             
                neigh = latt[(i+1)%N,j] + latt[i,(j+1)%N] + latt[(i-1)%N,j] + latt[i,(j-1)%N]
                arg=2*J*neigh/(kb*temp)
                p=np.exp(arg)/(1+np.exp(arg))

                if rand()<p:
                    latt[i,j]=1
                else:
                    latt[i,j]=-1
    return latt

def metropolis(latt,temp):
    N=latt.shape[0]
    for i in range(N):
        for j in range(N):
             
                neigh = latt[(i+1)%N,j] + latt[i,(j+1)%N] + latt[(i-1)%N,j] + latt[i,(j-1)%N]
                dE = 2*J*latt[i,j]*neigh+2*latt[i,j]*B
                if dE < 0:
                    latt[i, j] *= -1
                elif rand() < np.exp(-dE/(kb*temp)):
                    latt[i, j] *= -1
    return latt

def energy(latt,B):
    #vectorial computation of the S_alpha*S_beta for the closest neighbours
    Sup = np.roll(latt, -1, axis=0)
    Sdown = np.roll(latt, 1, axis=0)
    Sleft = np.roll(latt, 1, axis=1)
    Sright = np.roll(latt, -1, axis=1)
    
    sum_neighbours= Sup+Sdown+Sleft+Sright
    
    H=-J*np.sum(sum_neighbours*latt)-B*np.sum(latt)
    return H


def magnetization(config):
    mag = np.sum(config)
    return mag
  
  
  def ising(N,T,B,algo):

    E,M,Cb,Khi,U = np.zeros(T.shape), np.zeros(T.shape), np.zeros(T.shape), np.zeros(T.shape), np.zeros(T.shape)
    

    norm1=sweeps*N**2
    norm2=sweeps**2*N**2

    for n in range(len(T)):
        e1 = m1 = e2 = m2 = m4 =0
    
        config = initial(N)
    
        for i in range(equilibrium): #loops to equilibrate the system
            config=algo(config, T[n]) 
            
        for i in range(sweeps): 
            config=algo(config, T[n])
        
            En = energy(config,B) 
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
  
  
  kb=1
J=1
B=0
equilibrium = 300
sweeps = 300

T= np.linspace(1, 3.5, 50)
N=16


Em = np.zeros(T.shape)
Mm,Cbm,Khim,Um = np.zeros(T.shape), np.zeros(T.shape), np.zeros(T.shape), np.zeros(T.shape)


Eh = np.zeros(T.shape)
Mh,Cbh,Khih,Uh = np.zeros(T.shape), np.zeros(T.shape), np.zeros(T.shape), np.zeros(T.shape)


Em,Mm,Cbm,Khim,Um=ising(N,T,B,metropolis)
Eh,Mh,Cbh,Khih,Uh=ising(N,T,B,heatbath)

plt.figure(1)
plt.plot(T,abs(Mm),'ro',label="metropolis")
plt.plot(T,abs(Mh),'co',label="heatbath")
plt.legend()
plt.xlabel("Temperature (T)") 
plt.ylabel("Order parameter")
plt.title("sweeps = equilibriate sweeps ="+str(sweeps)+", N="+str(N))
plt.savefig("MagnetP4b.jpg")

plt.figure(2)
plt.plot(T,Khim,'ro',label="metropolis")
plt.plot(T,Khih,'co',label="heatbath")
plt.legend()
plt.xlabel("Temperature (T)") 
plt.ylabel("Susceptibility (Khi)")
plt.title("sweeps = equilibriate sweeps ="+str(sweeps)+", N="+str(N))
plt.savefig("SuceptP4b.jpg")

plt.figure(3)
plt.plot(T,Cbm,'ro',label="metropolis")
plt.plot(T,Cbh,'co',label="heatbath")
plt.legend()
plt.xlabel("Temperature (T)") 
plt.ylabel("Specific heat (Cb)")
plt.title("sweeps = equilibriate sweeps ="+str(sweeps)+", N="+str(N))
plt.savefig("SpeheatP4b.jpg")
