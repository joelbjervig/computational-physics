#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  5 16:09:37 2021

@author: joelbjervig
"""

import numpy as np
import matplotlib.pyplot as plt
from math import *

def initial(N):
    lattice = np.random.choice([-1,1],(N,N))
    return lattice

def H(latt,J):
    #vectorial computation of the S_alpha*S_beta for the closest neighbours
    Sup = np.roll(latt, -1, axis=0)
    Sdown = np.roll(latt, 1, axis=0)
    Sleft = np.roll(latt, 1, axis=1)
    Sright = np.roll(latt, -1, axis=1)
    
    sum_neighbours= Sup+Sdown+Sleft+Sright
    
    H=-J*np.sum(sum_neighbours*latt)
    return H/4

def Energy(latt,k,l):
    N=latt.shape[0]
    s=latt[k,l]
    neighb = latt[(k+1)%N,l] + latt[k,(l+1)%N] + latt[(k-1)%N,l] + latt[k,(l-1)%N]
    return 2*neighb*s

def metropolis(latt,T):
    Nx,Ny=latt.shape
    for i in range(Nx):
        for j in range(Ny):
            dE=Energy(latt,i,j)
            
            if dE<=0:
                latt[i,j]*=-1
                
            else:
                r=np.random.uniform(0,1)
                if r<exp(-dE/(kb*T)):
                    latt[i,j]*=-1
    return latt

def magnetization(latt):
    return np.sum(latt)


N=8
J=1
kb=1
T=np.linspace(1.6,3.2,70)

relax=100

orderparam,Khi,Cb= np.zeros(T.shape),np.zeros(T.shape),np.zeros(T.shape)

for n in range(len(T)):
    lattice= initial(N)
    
    for i in range(relax):
        metropolis(lattice,T[n])
        
    for i in range(ite):
        absM=E=M=E2=M2=0
        
        lattice = metropolis(lattice,T[n])
        ener = H(lattice,J)
        mag = magnetization(lattice)
        
        absM+=abs(mag)
        E+=ener
        M+=mag
        E2+=ener**2
        M2+=mag**2
        
    orderparam[n]=absM/(N*N) 
    Khi[n]=(M2/ite - (M/ite)**2)/T[n]
    Cb[n]=(E2/ite - (E/ite)**2)/(kb*T[n]**2)
