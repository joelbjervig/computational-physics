#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  3 10:01:20 2021

@author: Joel Bjervig, Thomas Herrad, Anton Palm
"""

from math import *
from library import bode

# ODE
dydt = lambda p,p
dpdt = lambda y, -4*pi*y


# exacto solution of problem
p_exact = lambda t, c1*cos(2*pi*t)- 2*pi*c2*sin(2*pi*t)
y_exact = lambda t, c2*cos(2*pi*t) + c1/(2*pi)*sin(2*pi*t)

