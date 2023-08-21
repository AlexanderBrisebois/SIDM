import matplotlib.pyplot as plt
import astropy.units as u
import numpy as np
from galpy.potential import MWPotential2014
from galpy.potential import NFWPotential

import matplotlib
from astropy.constants import G


def PDESolve(R_min=-20,R_max=20,z_min=-20,z_max=20,dark_matter_pot=MWPotential2014[2], milky_way=MWPotential2014[0:2],
             grid_size=70,vel_disp=150*u.km/u.s,rho_0=0.05*u.M_sun/u.pc**3,r_0 = 2.7 * u.kpc, cut_off_increasing=False, error=1e-6):
    # Returns the dimensionless gravitational potential h, and the density

    # Note: The iteration over the grid points continues until the provided error threshold is reached, or until the
    # maximum step-by-step change across the grid begins to increase.  If the variable cut_off_increasing is set to
    # False, then the iteration will only stop when the error threshold is reached.

    a0 = (4 * np.pi * G * rho_0 * r_0**2 / vel_disp**2).to('')
    c=11.5


    h = np.zeros((grid_size,grid_size))




    R_grid = np.linspace(R_min,R_max,grid_size)*u.kpc / r_0
    z_grid = np.linspace(z_min,z_max,grid_size)*u.kpc / r_0
    for i,R in enumerate(R_grid):
        for j,z in enumerate(z_grid):
            h[j][i] = np.log(dark_matter_pot.dens(R*r_0,z*r_0, quantity=True, ro=8,vo=220) / rho_0)


    delta = R_grid[1] - R_grid[0]

    mw_density = np.zeros((grid_size,grid_size)) *  u.M_sun/u.pc**3

    for pot in milky_way:
        for i in range(grid_size):
            mw_density[i] += pot.dens(np.abs(R_grid*r_0), np.abs(z_grid[i]*r_0), quantity=True, ro=8, vo=220)

    prev_diff = 100

    if cut_off_increasing==True:
        for i in range(100000):
            h_last = h.copy()

            h[1:-1, 1:-1] = 0.0002 * ((h[1:-1, 2:] + h[1:-1, :-2] - 2 * h[1:-1, 1:-1]) / delta ** 2 + 1 / (R_grid[1:-1]) * (
                        h[1:-1, 2:] - h[1:-1, :-2]) / (2 * delta) + (
                                                  h[2:, 1:-1] + h[:-2, 1:-1] - 2 * h[1:-1, 1:-1]) / delta ** 2 + a0 * (
                                                  mw_density[1:-1, 1:-1] / rho_0 + np.e ** (h[1:-1, 1:-1]))) + h[1:-1, 1:-1]

            diff = np.max(abs(h - h_last))
            print(diff, i)
            if diff < error or prev_diff < diff:
                print(np.max(abs(h - h_last)))
                print(f'Solution converged in {i} steps.')
                break
            prev_diff = diff

    else:
        for i in range(100000):
            h_last = h.copy()

            h[1:-1, 1:-1] = 0.0002 * (
                        (h[1:-1, 2:] + h[1:-1, :-2] - 2 * h[1:-1, 1:-1]) / delta ** 2 + 1 / (R_grid[1:-1]) * (
                        h[1:-1, 2:] - h[1:-1, :-2]) / (2 * delta) + (
                                h[2:, 1:-1] + h[:-2, 1:-1] - 2 * h[1:-1, 1:-1]) / delta ** 2 + a0 * (
                                mw_density[1:-1, 1:-1] / rho_0 + np.e ** (h[1:-1, 1:-1]))) + h[1:-1, 1:-1]

            diff = np.max(abs(h - h_last))
            print(diff, i)
            if diff < error:
                print(np.max(abs(h - h_last)))
                print(f'Solution converged in {i} steps.')
                break



    rho_new = rho_0 * np.e**h
    return h, rho_new

# h, rho = PDESolve()

# plt.imshow(rho[:,:], interpolation='nearest', extent =[-20, 20, -20, 20], norm=matplotlib.colors.LogNorm(), cmap='Set1',  origin='lower')
# plt.xlabel('R (kpc)')
# plt.ylabel('z (kpc)')
# plt.title("SIDM Density")
# plt.colorbar()
# plt.show()
