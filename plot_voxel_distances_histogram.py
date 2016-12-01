#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  1 13:02:35 2016

@author: lukas
"""
import numpy as np
import matplotlib.pyplot as pl
from mpl_toolkits.mplot3d import axes3d

all_sizes = np.genfromtxt('all_sizes.txt')
sizes = pl.hist(all_sizes)
pl.plot(sizes[0])