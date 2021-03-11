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
    Sup = np.roll(latt, -1, axis=0)*latt
    Sdown = np.roll(latt, 1, axis=0)*latt
    Sleft = np.roll(latt, 1, axis=1)*latt
    Sright = np.roll(latt, -1, axis=1)*latt
    
    sum_neighbours= Sup+Sdown+Sleft+Sright
    
    H=-J*np.sum(sum_neighbours)/4
    return H

def Energy(latt,k,l):
    N=latt.shape[0]
    s=latt[k,l]
    neighb = latt[(k+1)%N,l] + latt[k,(l+1)%N] + latt[(k-1)%N,l] + latt[k,(l-1)%N]
    return 2*neighb*s

def metropolis(latt,T):
    Nx,Ny=latt.shape
    for i in range(Nx):
        for j in range(Ny):
            k = np.random.randint(0,Nx)
            l = np.random.randint(0,Ny)
            
            Ei=Energy(latt,k,l)
            latt[k,l]*=-1
            Ef=Energy(latt,k,l)
            
            dE=Ef-Ei
            
            if dE<0:
                continue
            else:
                r=np.random.uniform(0,1)
                if r<exp(-dE/(kb*T)):
                    
                    latt[k,l]*=-1
    return latt

def magnetization(latt):
    return np.mean(latt)

N=8
J=1
kb=1
T=np.linspace(1.5,3,50)

ite=100

orderparam,Khi,Cb= np.zeros(T.shape),np.zeros(T.shape),np.zeros(T.shape)

for n in range(len(T)):
    lattice= initial(N)
    
    for i in range(ite):
        metropolis(lattice,T[n])
    for i in range(ite):
        absM=E=M=E2=M2=0
        
        metropolis(lattice,T[n])
        ener = H(lattice,J)
        mag = magnetization(lattice)
        
        absM+=abs(mag)
        E+=ener
        M+=mag
        E2+=ener**2
        M2+=mag**2
        
    orderparam[n]=absM/norm #independant of the size of the lattice
    Khi[n]=M2/ite - (M/ite)**2
    Cb[n]=E2/ite - (E/ite)**2
    
    
