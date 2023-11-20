import numpy as np
import matplotlib.pyplot as plt
from diagnostics_pkg import operators
import xarray as xr

ds = xr.open_dataset('outdir_0001/vars.nc')

nconv = 0
Nt = len(ds['time']) - nconv
N = len(ds['x'])

L = ds['x'][-1].values
psi = ds['psi'][nconv:,:,:].values
q = ds['q'][nconv:,:,:].values
time = ds['time'][nconv:].values


kf = 10
epsilon = 1

# delta = L/(N - 1) # dx
# ke = operators.energy(psi, q, delta)
# plt.plot(time, ke)

delta = ds['x'][-1] - ds['x'][-2]
delta = delta.values
energy = -np.sum(psi*q, axis = (-2,-1))*delta**2/L**2
plt.title('energy input')
plt.ylabel('Energy')
plt.xlabel('time')
plt.plot(time, energy, label = "observed")
plt.plot(time, epsilon*time, label = "expected") # so the energy injection rate is 5*10e-6
plt.legend()
plt.show()

enstrophy = np.sum(q**2, axis = (-2,-1))*delta**2/L**2
plt.title('enstrophy input')
plt.ylabel('Enstrophy')
plt.xlabel('time')
plt.plot(time, enstrophy, label = "observed")
plt.plot(time, time*(2*np.pi)**2*kf**2*epsilon)

# There is a slight offset, but frankly I don't care enough to dig here. It might be
# just due to the randomness of the forcing, and as the random number generator I use
# outputs the same series every time I can't test it quickly.