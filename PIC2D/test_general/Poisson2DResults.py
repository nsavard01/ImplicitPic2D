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
def format_coord(x, y):
    xarr = X[0,:]
    yarr = Y[:,0]
    if ((x > xarr.min()) & (x <= xarr.max()) &
        (y > yarr.min()) & (y <= yarr.max())):
        col = np.searchsorted(xarr, x)-1
        row = np.searchsorted(yarr, y)-1
        z = Z[row, col]
        return f'x={x:1.4f}, y={y:1.4f}, z={z:1.4f}   [{row},{col}]'
    else:
        return f'x={x:1.4f}, y={y:1.4f}'


numNodes = np.fromfile('NumNodes.dat', dtype=int)
gridX = np.fromfile('gridX.dat')
gridY = np.fromfile('gridY.dat')
numDiag = np.fromfile('NumDiag.dat', dtype = int)[0]
# source_term = np.fromfile('source_term.dat')
# E_x = np.fromfile('E_x.dat')
# E_y = np.fromfile('E_y.dat')

grid2DX, grid2DY = np.meshgrid(gridX, gridY, indexing = 'ij')
Length = gridX[-1] - gridX[0]
Width = gridY[-1] - gridY[0]



stage = 1

numNodes_x = int((numNodes[0] + (2**(stage-1)-1)) / 2**(stage-1))
numNodes_y = int((numNodes[1] + (2**(stage-1)-1)) / 2**(stage-1))
gridX = gridX[0:numNodes[0]:2**(stage-1)]
gridY = gridY[0:numNodes[1]:2**(stage-1)]
grid2DX, grid2DY = np.meshgrid(gridX, gridY, indexing = 'ij')

plt.figure()
for i in range(numDiag):

    sol_1 = np.fromfile('finalSol' + str(i+1) + '.dat')
    plt.pcolormesh(grid2DX, grid2DY, np.reshape(sol_1, (numNodes_x, numNodes_y), order = 'F'), shading = 'nearest')
    cb = plt.colorbar()
    plt.pause(0.08)
    cb.remove()
    plt.cla()


# plt.figure()
# plt.pcolormesh(grid2DX, grid2DY, np.reshape(E_x, (numNodes_x, numNodes_y), order = 'F'), shading = 'nearest')
# plt.colorbar()
#
# plt.figure()
# plt.pcolormesh(grid2DX, grid2DY, np.reshape(E_y, (numNodes_x, numNodes_y), order = 'F'), shading = 'nearest')
# plt.colorbar()

# test = np.reshape(sol_1, (numNodes[0], numNodes[1]), order = 'F')[0:numNodes[0]:2, 0:numNodes[1]:2] - np.reshape(GSRes, (numNodes_x, numNodes_y), order = 'F')
# plt.figure()
# plt.pcolormesh(grid2DX, grid2DY, test, shading = 'nearest')
# plt.colorbar()

# test = np.fromfile('test.dat')
# grid2DX, grid2DY = np.meshgrid(gridX[0:numNodes[0]:2], gridY[0:numNodes[1]:2], indexing = 'ij')
# other = np.reshape(GSRes, (numNodes[0], numNodes[1]), order = 'F')
#
# plt.figure()
# plt.pcolormesh(grid2DX, grid2DY, np.absolute(np.reshape(test, (int((numNodes[0]+1)/2), int((numNodes[1]+1)/2)), order = 'F') - other[0:numNodes[0]:2, 0:numNodes[1]:2]), norm = 'log')
# plt.colorbar()



