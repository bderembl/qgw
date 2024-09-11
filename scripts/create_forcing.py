#!/usr/bin/env python


import matplotlib.pyplot as plt
import numpy as np
from scipy.io import netcdf_file
import qgutils as qg

plt.ion()

dirOUT = './'

# physical constants
N = 512
nl = 3
L0 = 4000e3 # size domain
tau0 = 2e-5 # wind magnitude 
dh = np.array([300, 700, 4000]) # layer thickness
f0 = 1e-4
beta = 2e-11
h_shelf = 1000 # topo anomaly in m
w_shelf = 100e3  # width of the shelf in m
# zapiola parameters
h_zap = 700
w_zap = 400e3
x_zap = 1000e3
y_zap = 1000e3


si_x = N+1;
si_y = N+1;

# create grid
x = np.linspace(0,L0,si_x)
y = np.linspace(0,L0,si_y)

xx,yy = np.meshgrid(x,y)

# create forcing
q_forcing = -tau0/dh[0]*2*np.pi/L0*np.sin(2*np.pi*yy/L0);

# modify forcing
q_forcing[int(N/2):,:] = 0 # only subtropical gyre


# topography

# sides
def shelf(x,d):
  return (1-np.exp(-x**2/(2*d**2)))


topo = h_shelf*(1-shelf(xx,w_shelf)*shelf(L0-xx,w_shelf)*shelf(yy,w_shelf)*shelf(L0-yy,w_shelf))

topo += h_zap*np.exp(-((xx - x_zap)**2 + (yy - y_zap)**2)/(2*w_zap**2))



# Save output
qg.write_nc("forcing.nc",{'topo':topo},timeDim=True)
