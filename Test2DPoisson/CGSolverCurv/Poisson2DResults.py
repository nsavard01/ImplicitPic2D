#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 11 13:42:42 2022

@author: nicolas
"""

import numpy as np
import matplotlib
matplotlib.use('Qt5Agg')
import matplotlib.pyplot as plt
import math
import scipy
import scipy.optimize as opt
from scipy.stats import cosine
import scipy.sparse as sp
import time
import sys


# Define constants
eps_0 = scipy.constants.epsilon_0
c = scipy.constants.c
m_e = scipy.constants.m_e
m_p = scipy.constants.m_p
mu_0 = scipy.constants.mu_0
k_boltz = scipy.constants.k
e = scipy.constants.e

Length = 0.05
Width = 0.05

def sinFunc(numNodes, length, del_x):
    l = np.arange(1, numNodes+1) #like fortran style
    grid = length * ((l-1)/(numNodes-1) - (1/(numNodes-1) - del_x/length) * np.sin(2 * np.pi * (l-1) / (numNodes-1)) /np.sin(2 * np.pi / (numNodes-1)))
    return grid


# numNodes_x = 1025
# numNodes_y = 129
# GSRes = np.fromfile('finalSol.dat')
# gridX = np.linspace(0, Length, numNodes_x)
# gridY = np.linspace(0, Width, numNodes_y)
# grid2DX, grid2DY = np.meshgrid(gridX, gridY, indexing = 'ij')
# plt.figure()
# plt.pcolormesh(grid2DX, grid2DY, np.reshape(GSRes, (numNodes_x, numNodes_y), order = 'F'))
# plt.colorbar()
#

numNodes = np.fromfile('NumNodes.dat', dtype=int)
gridX = np.fromfile('gridX.dat')
gridY = np.fromfile('gridY.dat')
GSRes = np.fromfile('finalSol.dat')
# gridX = np.linspace(0, Length, numNodes_x)
# gridY = np.linspace(0, Width, numNodes_y)
grid2DX, grid2DY = np.meshgrid(gridX, gridY, indexing = 'ij')
plt.figure()
plt.pcolormesh(grid2DX, grid2DY, np.reshape(GSRes, (numNodes[0], numNodes[1]), order = 'F'))
plt.colorbar()



