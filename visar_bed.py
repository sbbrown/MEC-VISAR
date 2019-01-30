#!/usr/bin/env python
'''
Name:          visar_bed.py
License:       MIT
Version:       2.0
Last modified: 30 Jan. 2019 (SBB)
Authors:       Akel Hashim          (ahashim@slac.stanford.edu)
               Bob Nagler           (bnagler@slac.stanford.edu)
               Shaughnessy Brown    (sbbrown@slac.stanford.edu)
Description:   Creates visar bed object to calcluate and store intrinsic 
               properties for later calculations
'''

'''Imports -----------------------------------------------------------------'''
from numpy import *

'''Class Definitions -------------------------------------------------------''' 
class VISAR_Bed():
    ''' Class defining visar bed. Calculates VPF, etalons, etc.
        All SI units except in status where they are printed with units
    '''
    
    def __init__(self,etal_h,landa=532e-9,n=1.46071,delta=0.0318,theta=11.31):
        '''Initializes visar bed.

           etalon_h : etalon thickness
           landa    : wavelength. Standard for mec visar is 532nm
           n        : index of refraction at landa. Standard is 1.46071 (UVFS)
           delta    : chromatic dispersion. Standard is 0.0318 (UVFS)
           theta    : angle (in degrees) of beams on visar. Standard is 11.31
        '''
        
        self._h     = etal_h            # self._h in meters
        self._landa = landa  
        self._n     = n
        self._delta = delta
        self._theta = theta/180.0*pi    # theta in radians
        self._c     = 299792458         # speed of light in m/s

    def __call__(self):
        self.status()

    def delay(self):
        '''Finds temporal delay between two arms of interferometer'''
        
        t0 = 2*self._h*(self._n-1/self._n)/self._c
        
        # Correction factor for non-normal incidence
        angle_correction = 1.0/(cos(arcsin(sin(self._theta/2.0)/self._n))) 
        
        # Returns t0 in seconds
        return t0*angle_correction 

    def vpf0(self):
        '''Finds VFP. Uncorrected for windows or fast lens'''
        
        tau = self.delay()
        
        return self._landa/(2*tau*(1+self._delta))
