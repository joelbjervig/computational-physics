#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  4 15:50:03 2021

@author: Joel Bjervig, Anton Palm, Thomas Herard
"""
phiExact = lambda r: 1-0.5*(r+2)*exp(-r)
source = lambda r: -0.5*r*exp(-r)
chargeDensity = lambda r: (1/(8*pi))*exp(-r)
f = lambda r: -4*pi*r*(1/(8*pi))*exp(-r)