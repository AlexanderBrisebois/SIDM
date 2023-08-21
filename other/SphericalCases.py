from galpy.potential.mwpotentials import McMillan17
import numpy as np
from astropy.constants import M_sun, G
from galpy.potential import MWPotential2014
from galpy import potential
import astropy.units as u
from galpy.potential import evaluatePotentials



def non_sph_density(rho_0, potential, r, z, vel_disp):
    from galpy.potential import evaluatePotentials

    phi_0 = 0
    phi_rz = 0
    if potential==McMillan17[::2]:
        phi_0 += evaluatePotentials(potential, 0.000000001*u.kpc, 0.00000001*u.kpc, quantity=True)
        phi_rz += evaluatePotentials(potential, r, z, quantity=True)
    else:
        phi_0 += evaluatePotentials(potential, 0.000000001*u.kpc, 0.00000001*u.kpc, quantity=True, ro=8,vo=220)
        phi_rz += evaluatePotentials(potential, r, z, quantity=True, ro=8,vo=220)

    rho = rho_0 * np.e**((phi_0 - phi_rz) / (vel_disp)**2)

    #index = min(range(len(rho)), key=lambda i: abs(rho[i] - rho_0/2))
    #print("core radius" , r[index], "at index:", index, f'(veloctiy dispersion = {vel_disp})')
    #print(rho[index])
    return rho.value

milky_way = MWPotential2014[0:2]


def density(rho_0, potential, r, vel_disp):
    #print(potential(r,0, quantity=True), potential(0,0, quantity=True))
    #print(potential(r,0,quantity=True)[0] - potential(0,0, quantity=True))
    rho = rho_0 * np.e**((potential(0,0,quantity=True) - potential(r,0,quantity=True)) / (vel_disp)**2)

    index = min(range(len(rho)), key=lambda i: abs(rho[i] - rho_0/2))
    print("core radius" , r[index], "at index:", index, f'(veloctiy dispersion = {vel_disp})')
    #print(rho[index])
    return rho

def rforce(r):
    return G * mass(r) / r**2


def mass(r):
    mass_bulge = milky_way[0].mass(r, ro=8,vo=220, use_physical=False)
    mass_disk = milky_way[1].mass(r, ro=8, vo=220, use_physical=False)
    total_mass = mass_bulge + mass_disk
    return total_mass