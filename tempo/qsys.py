#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun 25 14:07:34 2023

The Quantum System (Qsys) class for representing a quantum system of one or more spins.

@author: hakalas
"""

#
# Import Python Packages
#
from qutip import identity, jmat, tensor, basis, qobj
import os
import numpy as np

from tempo.exceptions import *
from collections.abc import Iterable

class Qsys:
    """
    A class for representing a system of one or more spins and their spin operators.
    
    The object is constructed by passing a tuple of the dimensions of each spin's Hilbert space. For example, for a coupled system of a spin-1 and spin-1/2, use ``dimensions = (3,2)``.
    
    Also includes methods for obtaining the spin operators Sx, Sy, and Sz, as well as a list of the basis vectors for each spin.
    
    Parameters
    ----------
    dimensions : tuple
        A tuple of the dimensions of each spin's Hilbert space.
        The length of the tuple is the number of spins in the system.
        
    Attributes
    ----------
    dimensions: tuple of int
        A tuple of the dimensions of each spin's Hilbert space.
        The length of the tuple is the number of spins in the system.
    numparticles : int
        Number of spins in the system.
    stot : list of float
        A list of the spin quantum numbers of each spin. 
    Sx : list of `qutip.Qobj`
        A list of the Sx spin operators of each spin.
    Sy : list of `qutip.Qobj`
        A list of the Sy spin operators of each spin.
    Sz : list of `qutip.Qobj`
        A list of the Sz spin operators of each spin.
    """
    
    def __init__(self, dimensions: tuple):
        """
        Qsys constructor.

        Parameters
        ----------
        dimensions: tuple of int
            A tuple of the dimensions of each spin's Hilbert space.
            The length of the tuple is the number of spins in the system.
        """
        self._dimensions = dimensions

        # Input check
        if not isinstance(self.dimensions, tuple) or len(dimensions) == 0:
            raise TEMPO_ImproperInputException("Dimensions is not a tuple, or is empty")

        for c_idx, obj in enumerate(self.dimensions):
            if not isinstance(obj, int):
                raise TEMPO_ImproperInputException(f"Dimensions includes a non-int object at index {c_idx}")
            if obj <= 0:
                raise TEMPO_ImproperInputException("Invalid Dimension")
            

        self._numparticles = len(self._dimensions)
        
        self._stot = [(self._dimensions[i]-1)/2 for i in np.arange(self._numparticles)]
        
        self._Sx = self.oprs('x')
        self._Sy = self.oprs('y')
        self._Sz = self.oprs('z')
        
    def oprs(self, axis: str) -> Iterable[qobj.Qobj]:
        """
        Spin operators of each spin along `axis`.
        
        Parameters
        ----------
        axis : str {'x', 'y', 'z'}
            Axis along which to create spin operators; `'x'`, `'y'`, or `'z'`.
            
        Returns
        ----------
        oprs : list of `qutip.Qobj`
            List of (`qutip.Qobj`) spin operators along `axis` in the same order as spins in the `dimensions` tuple. 
        """
        if axis != "x" and axis != "y" and axis != "z":
            raise TEMPO_ImproperInputException("Axis for operator creation is not x,y,z")

        oprs = []
        idlist = [identity(self._dimensions[i]) for i in np.arange(self._numparticles)]
        
        for i in np.arange(self._numparticles):
            idlist[i] = jmat(self._stot[i], axis)
            oprs.append(tensor(idlist))
            idlist[i] = identity(self._dimensions[i])
        
        #oprs.append(oprs[0]+oprs[1])
        
        return oprs
    
    def basisstates(self, particlenum: int) -> Iterable[qobj.Qobj]:
        """
        Basis state vectors of spin at index `particlenum` in the `dimensions` tuple. 
        
        Parameters
        ----------
        particlenum : int
            Index of spin for which to generate basis vectors.
        
        Returns
        ----------
        basisstates : list of `qutip.Qobj`
            List of (`qutip.Qobj`) basis state vectors. 
            Length of list is the dimension of the spin's Hilbert space.
        """
        if not castable(particlenum, int):
            raise TEMPO_ImproperInputException("particlenum is not an integer")
        if particlenum >= len(self._dimensions) or particlenum < -1 * len(self._dimensions):
            raise TEMPO_ImproperInputException(f"Spin index {particlenum} was given a value outside the range -1 * len(dimensions) --> len(dimensions) - 1.")

        particlenum = int(particlenum)

        basisstates = []
        
        idlist = [identity(self._dimensions[i]) for i in np.arange(self._numparticles)]
        
        for i in np.arange(self._dimensions[particlenum]):
            idlist[particlenum] = basis(self._dimensions[particlenum], i)*basis(self._dimensions[particlenum], i).dag()
            basisstates.append(tensor(idlist))
            idlist[particlenum] = identity(self._dimensions[particlenum])
        
        return basisstates
    
    @property
    def dimensions(self):
        return self._dimensions
    
    @dimensions.setter
    def dimensions(self, dimensions):
        self._dimensions = dimensions
        self._numparticles = len(self._dimensions)
        self._stot = [(self._dimensions[i]-1)/2 for i in np.arange(self._numparticles)]
        self._Sx = self.oprs('x')
        self._Sy = self.oprs('y')
        self._Sz = self.oprs('z')

    
    @dimensions.deleter 
    def dimensions(self):
        del self._dimensions
    
    @property
    def numparticles(self):
        return self._numparticles
    
    @numparticles.setter
    def numparticles(self, num):
        self._numparticles = num
     
    @numparticles.deleter 
    def numparticles(self):
        del self._numparticles
    
    @property
    def stot(self):
        return self._stot
    
    @stot.setter
    def stot(self, ls):
        self._stot = ls
    
    @stot.deleter 
    def stot(self):
        del self._stot
    
    @property
    def Sx(self):
        return self._Sx
    
    @Sx.setter
    def Sx(self, Sx):   
        if type(Sx) == qobj.Qobj:
            self._Sx = Sx
       
    @Sx.deleter 
    def Sx(self):
        del self._Sx
    
    @property
    def Sy(self):
        return self._Sy
    
    @Sy.setter
    def Sy(self, Sy): 
        if type(Sy) == qobj.Qobj:
             self._Sy = Sy
       
    @Sy.deleter 
    def Sy(self):
        del self._Sy
    
    @property
    def Sz(self):
        return self._Sz
    
    @Sz.setter       
    def Sz(self, Sz): 
        if type(Sz) == qobj.Qobj:
             self._Sz = Sz

    @Sz.deleter   
    def Sz(self):
        del self._Sz
