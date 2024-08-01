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

# numNodes = 145
# GSRes = np.fromfile('Res_1.dat')
# gridX = np.linspace(0, Length, numNodes)
# gridY = np.linspace(0, Length, numNodes)
# grid2DX, grid2DY = np.meshgrid(gridX, gridY, indexing = 'ij')
# plt.figure()
# plt.pcolormesh(grid2DX, grid2DY, np.reshape(GSRes, (numNodes, numNodes), order = 'F'))
# plt.colorbar()
#
# numNodes = 73
# GSRes = np.fromfile('Res_1_2.dat')
# gridX = np.linspace(0, Length, numNodes)
# gridY = np.linspace(0, Length, numNodes)
# grid2DX, grid2DY = np.meshgrid(gridX, gridY, indexing = 'ij')
# plt.figure()
# plt.pcolormesh(grid2DX, grid2DY, np.reshape(GSRes, (numNodes, numNodes), order = 'F'))
# plt.colorbar()
#
# numNodes = 73
# GSRes = np.fromfile('Res_2.dat')
# gridX = np.linspace(0, Length, numNodes)
# gridY = np.linspace(0, Length, numNodes)
# grid2DX, grid2DY = np.meshgrid(gridX, gridY, indexing = 'ij')
# plt.figure()
# plt.pcolormesh(grid2DX, grid2DY, np.reshape(GSRes, (numNodes, numNodes), order = 'F'))
# plt.colorbar()
#
# numNodes = 37
# GSRes = np.fromfile('Res_2_3.dat')
# gridX = np.linspace(0, Length, numNodes)
# gridY = np.linspace(0, Length, numNodes)
# grid2DX, grid2DY = np.meshgrid(gridX, gridY, indexing = 'ij')
# plt.figure()
# plt.pcolormesh(grid2DX, grid2DY, np.reshape(GSRes, (numNodes, numNodes), order = 'F'))
# plt.colorbar()
#
# numNodes = 37
# GSRes = np.fromfile('Sol_3.dat')
# gridX = np.linspace(0, Length, numNodes)
# gridY = np.linspace(0, Length, numNodes)
# grid2DX, grid2DY = np.meshgrid(gridX, gridY, indexing = 'ij')
# plt.figure()
# plt.pcolormesh(grid2DX, grid2DY, np.reshape(GSRes, (numNodes, numNodes), order = 'F'))
# plt.colorbar()
#
# numNodes = 73
# GSRes = np.fromfile('Sol_2_3.dat')
# gridX = np.linspace(0, Length, numNodes)
# gridY = np.linspace(0, Length, numNodes)
# grid2DX, grid2DY = np.meshgrid(gridX, gridY, indexing = 'ij')
# plt.figure()
# plt.pcolormesh(grid2DX, grid2DY, np.reshape(GSRes, (numNodes, numNodes), order = 'F'))
# plt.colorbar()
#
# numNodes = 73
# GSRes = np.fromfile('Res_2.dat')
# gridX = np.linspace(0, Length, numNodes)
# gridY = np.linspace(0, Length, numNodes)
# grid2DX, grid2DY = np.meshgrid(gridX, gridY, indexing = 'ij')
# plt.figure()
# plt.pcolormesh(grid2DX, grid2DY, np.reshape(GSRes, (numNodes, numNodes), order = 'F'))
# plt.colorbar()

numNodes_x = 1025
numNodes_y = 129
GSRes = np.fromfile('finalSol.dat')
gridX = np.linspace(0, Length, numNodes_x)
gridY = np.linspace(0, Width, numNodes_y)
grid2DX, grid2DY = np.meshgrid(gridX, gridY, indexing = 'ij')
plt.figure()
plt.pcolormesh(grid2DX, grid2DY, np.reshape(GSRes, (numNodes_x, numNodes_y), order = 'F'))
plt.colorbar()



