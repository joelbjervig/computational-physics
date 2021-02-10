#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb  8 14:16:41 2021

@author: joelbjervig Thomas Herard
"""
import numpy as np
from math import *
import matplotlib.pyplot as plt

import sys
sys.path.insert(1, '/Users/joelbjervig/documents/universitet/kurser/pagaende/computational physics/computational-physics/homeExc')

from library import numerov
from library import search



# guess trial eigenvalue as something small
epsilon1 = 10**(-8)
epsilon2 = 10**(-6)
#k_init = 0.001 
x0 = 0.01
dx = 0.1
dk = 0.1
y0 = 0
y1 = y0+epsilon2

x = np.zeros(100)

f = lambda x:sin(pi*x)
k = lambda x:0
S = lambda x:0


while():
    