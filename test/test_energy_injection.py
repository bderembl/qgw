import numpy as np
import matplotlib.pyplot as plt
from scipy.io import netcdf_file
import sys

plt.ion()

dir0 = "outdir_"
if len(sys.argv) > 1:
  irun = sys.argv[1]
  dir0 = dir0 + str(irun).zfill(4) + '/'

f = netcdf_file(dir0 + 'vars.nc','r')

# read model parameters
exec(open(dir0 + "params.in").read())

Delta = Lx/Nx

psi = f.variables['psi'][:].copy()
q = f.variables['q'][:].copy()
time = f.variables['time'][:].copy()

energy = -np.sum(psi*q, axis = (-2,-1))*Delta**2/Lx**2
enstrophy = np.sum(q**2, axis = (-2,-1))*Delta**2/Lx**2


plt.figure()
plt.plot(time, energy, label = "observed")
plt.plot(time, sigma_f**2*time, label = "expected") 
plt.legend()
plt.title('energy input')
plt.ylabel('Energy')
plt.xlabel('time')

plt.figure()
plt.plot(time, enstrophy, label = "observed")
plt.plot(time, time*(2*np.pi)**2*k_f**2*sigma_f**2)
plt.title('enstrophy input')
plt.ylabel('Enstrophy')
plt.xlabel('time')

f.close()
