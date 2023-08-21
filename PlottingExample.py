
import matplotlib.pyplot as plt
import DirectInterpRZPotential


from galpy.potential import MWPotential2014
import numpy as np
import astropy.units as u
import matplotlib.pyplot as plt
import matplotlib

import JeansPDESolver

h, rho = JeansPDESolver.PDESolve(cut_off_increasing=True)

#Density Plot
plt.imshow(rho[:,:], interpolation='nearest', extent =[-20, 20, -20, 20], norm=matplotlib.colors.LogNorm(), cmap='Set1',  origin='lower')
plt.xlabel('R (kpc)')
plt.ylabel('z (kpc)')
plt.title("SIDM Density")
plt.colorbar()
plt.show()


R_grid = np.linspace(-20,20,len(rho))*u.kpc
z_grid = np.linspace(-20,20,len(rho))*u.kpc

#Creating Interpolation Object
ip = DirectInterpRZPotential.DirectInterpRZPotential(directPotentialGrid=-h*(150*u.km/u.s)**2, directRGrid=R_grid.value, directzGrid=z_grid.value, interpPot=True, logR=False, fallbackPotential=MWPotential2014[2], interpDens=True, directDensGrid=rho)

test_R_grid = np.linspace(0,19.9,400)*u.kpc
test_z_grid = np.linspace(0,19.9,400)*u.kpc

#Density Interpolation
interpolated_dens = np.zeros((len(test_R_grid), len(test_z_grid)))

for i,R in enumerate(test_R_grid):
    for j,z in enumerate(test_z_grid):
        interpolated_dens[j][i] = ip.dens_eval(R.value, z.value)

plt.imshow(interpolated_dens, interpolation='nearest', extent =[test_R_grid.min().value, test_R_grid[:].max().value, test_z_grid.min().value, test_z_grid.max().value], cmap='hsv',  origin='lower', norm=matplotlib.colors.LogNorm())
plt.colorbar()
plt.title("Interpolated Density")
plt.show()

# R_Force plot
new_evaluated_R_force = ip.Rforce(test_R_grid,test_z_grid)
plt.imshow(new_evaluated_R_force, interpolation='nearest', extent =[test_R_grid.min().value, test_R_grid[:].max().value, test_z_grid.min().value, test_z_grid.max().value], cmap='RdBu',  origin='lower', )
plt.colorbar()
plt.show()


# v_circ plot (versus actual vcirc from MWPotential2014)
from galpy.potential import vcirc

plt.plot(test_R_grid[1:], vcirc(MWPotential2014, test_R_grid[1:], quantity=True, ro=8,vo=220))
        # MWPotential2014[2].vcirc(test_R_grid[1:], quantity=True, ro=10,vo=100))
plt.plot(test_R_grid, ip.vcirc(test_R_grid))

plt.show()

