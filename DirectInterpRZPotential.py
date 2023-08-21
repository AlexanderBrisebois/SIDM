import copy
import ctypes
import ctypes.util
from functools import wraps

import numpy
from numpy.ctypeslib import ndpointer
from scipy import interpolate

# from ..util import _load_extension_libs, multi
# from ..util.conversion import physical_conversion
# from .Potential import Potential
import numpy as np


import numpy
from galpy.potential import interpRZPotential
import astropy.units as u

class DirectInterpRZPotential(interpRZPotential):
    def __init__(self, RZPot=None,rgrid=(numpy.log(0.01),numpy.log(20.),101),
                 zgrid=(0.,1.,101),logR=True,
                 interpPot=False,interpRforce=False,interpzforce=False,
                 interpDens=False,
                 interpvcirc=False,
                 interpdvcircdr=False,
                 interpepifreq=False,interpverticalfreq=False,
                 ro=None,vo=None,
                 use_c=False,enable_c=False,zsym=True,
                 numcores=None,
                 densityGrid=None,
                 directPotentialGrid=None,
                 directRGrid=None,
                 directzGrid=None,
                 fallbackPotential=None,
                 directDensGrid=None):
        # if not roSet: self._roSet= False
        # if not voSet: self._voSet= False
        self._origPot= RZPot
        self._rgrid= numpy.linspace(*rgrid)
        self._logR= logR
        if self._logR:
            self._rgrid= numpy.exp(self._rgrid)
            self._logrgrid= numpy.log(self._rgrid)
        self._zgrid= numpy.linspace(*zgrid)
        self._interpPot= interpPot
        self._interpRforce= interpRforce
        self._interpzforce= interpzforce
        self._interpDens= interpDens
        self._interpvcirc= interpvcirc
        self._interpdvcircdr= interpdvcircdr
        self._interpepifreq= interpepifreq
        self._interpverticalfreq= interpverticalfreq
        #self._enable_c= enable_c*ext_loaded
        # self.hasC= self._enable_c
        self._zsym= zsym
        self._densityGrid = densityGrid
        self._directPotentialGrid = directPotentialGrid
        self._directRGrid = directRGrid
        self.directzGrid = directzGrid
        self._fallbackPotential = fallbackPotential
        self.directDensGrid = directDensGrid

        #potential interpolation
        self._potInterp= interpolate.RectBivariateSpline(self._directRGrid,
                                                                 self.directzGrid,
                                                                 self._directPotentialGrid,
                                                                 kx=3,ky=3,s=0.)
        self.used_fallback = False
                     
        if interpDens:
            self.densInterp = interpolate.RectBivariateSpline(self._directRGrid,
                                                                     self.directzGrid,
                                                                     self.directDensGrid,
                                                                     kx=3,ky=3,s=0.)


    def _evaluate(self,R,z,phi=0.,t=0.):
        from galpy.potential import evaluatePotentials
        if self._interpPot:
            out= numpy.empty(R.shape)
            indx= (R >= self._directRGrid[0])*(R <= self._directRGrid[-1])\
                *(z <= self.directzGrid[-1])*(z >= self.directzGrid[0])
            if numpy.sum(indx) > 0:
                if self._logR:
                    out[indx]= self._potInterp.ev(z[indx],numpy.log(R[indx]))
                else:
                    out[indx]= self._potInterp.ev(z[indx],R[indx])

            if numpy.sum(True^indx) > 0:
                if self.used_fallback==False:
                    print( "Fallback potential used rather than interpolation, because requested coordinates are outside the provided grid.")
                out[True^indx]= evaluatePotentials(self._fallbackPotential, R[True^indx]*u.kpc,z[True^indx]*u.kpc, quantity=True, ro=8,vo=220)
                used_fallback = True
            return out
        else:
            print("Fallback potential used, rather than interpolation.")
            return evaluatePotentials(self._fallbackPotential,R,z)

    def Rforce(self, r_grid, z_grid):
        grid = np.zeros((len(z_grid),len(r_grid)))*u.km/u.s/u.Myr

        r_grid = r_grid.to(u.kpc)
        z_grid = z_grid.to(u.kpc)

        delta = r_grid[1] - r_grid[0]
        for i,r in enumerate(r_grid):
            for j,z in enumerate(z_grid):
                grid[j][i] = (((self._evaluate((r+delta).value,z.value) - self._evaluate((r).value,z.value)) / (delta.to(u.km))) * u.km**2/u.s**2).to(u.km/u.s/u.Myr)
        return -grid

    def vcirc(self, R_grid):
        R_force = self.Rforce(R_grid,np.array([0])*u.kpc)
        v_c = R_grid[:] * -R_force[0,:]
        v_c = np.sqrt(v_c).to(u.km/u.s)
        return v_c

    def dens_eval(self, R,z,phi=0.,t=0.):
        if self._interpDens:
            out= numpy.empty(R.shape)
            indx= True#(R >= self._directRGrid[0])*(R <= self._directRGrid[-1])\
                #*(z <= self.directzGrid[-1])*(z >= self.directzGrid[0])
            if numpy.sum(indx) > 0:
                if self._logR:
                    out[indx]= self.densInterp.ev(numpy.log(z[indx]),R[indx])
                else:
                    out[indx]= self.densInterp.ev(z[indx],R[indx])
            return out
        else:
            pass


