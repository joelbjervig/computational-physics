#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 12 09:37:12 2021

@author: joelbjervig
"""

import numpy as np
from math import *
import matplotlib.pyplot as plt

#from library import

# const def
a = 4;
r_max = 30

n = 10
N = 4*n

# function def
phi_GT  = lambda r: -1/sqrt(2*pi)*(np.exp(-a*r))
phi_LT  = lambda r: 1/sqrt(2*pi)*(np.exp(a*r)-np.exp(-a*r))
S       = lambda r: -1/2*r*np.exp(-r)

for 