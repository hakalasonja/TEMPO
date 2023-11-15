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
    
    def __init__(self, dimensions):
        """
        Qsys constructor.
        """
        self.setdimensions(dimensions)
        self._numparticles = len(self._dimensions)
        
        self._stot = [(self._dimensions[i]-1)/2 for i in np.arange(self._numparticles)]
        
        self._Sx = self.oprs('x')
        self._Sy = self.oprs('y')
        self._Sz = self.oprs('z')
        
    def oprs(self, axis):
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
        oprs = []
        idlist = [identity(self._dimensions[i]) for i in np.arange(self._numparticles)]
        
        for i in np.arange(self._numparticles):
            idlist[i] = jmat(self._stot[i], axis)
            oprs.append(tensor(idlist))
            idlist[i] = identity(self._dimensions[i])
        
        #oprs.append(oprs[0]+oprs[1])
        
        return oprs
    
    def basisstates(self, particlenum):
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
        basisstates = []
        
        idlist = [identity(self._dimensions[i]) for i in np.arange(self._numparticles)]
        
        for i in np.arange(self._dimensions[particlenum]):
            idlist[particlenum] = basis(self._dimensions[particlenum], i)*basis(self._dimensions[particlenum], i).dag()
            basisstates.append(tensor(idlist))
            idlist[particlenum] = identity(self._dimensions[particlenum])
        
        return basisstates
    
    def getdimensions(self):
        return self._dimensions
    
    def getnumparticles(self):
        return self._numparticles
    
    def getstot(self):
        return self._stot
    
    def getSx(self):
        return self._Sx
    
    def getSy(self):
        return self._Sy
    
    def getSz(self):
        return self._Sz
    
    def setdimensions(self, dimensions):
        
        # dimensions must be a tuple
        
        if type(dimensions) == tuple:
            self._dimensions = dimensions
            self._numparticles = len(self._dimensions)
            self._stot = [(self._dimensions[i]-1)/2 for i in np.arange(self._numparticles)]
            self._Sx = self.oprs('x')
            self._Sy = self.oprs('y')
            self._Sz = self.oprs('z')
        else:
            raise TypeError("Dimensions must be a tuple")
            
    def setnumparticles(self, num):
        
        if type(num) == int:
            self._numparticles = num
        else:
            raise TypeError("Number of particles must be an integer")
            
    def setstot(self, ls):
        
        if type(ls) == list:
            self._stot = ls
        else: 
            raise TypeError("Maximum spins must be stored in a list")
            
    def setSx(self, Sx):
        
        if type(Sx) == qobj.Qobj:
            self._Sx = Sx
            
    def setSy(self, Sy):
         
        if type(Sy) == qobj.Qobj:
             self._Sy = Sy
             
    def setSz(self, Sz):
         
        if type(Sz) == qobj.Qobj:
             self._Sz = Sz
            
    def deldimensions(self):
        del self._dimensions
        
    def delnumparticles(self):
        del self._numparticles
        
    def delstot(self):
        del self._stot
        
    def delSx(self):
        del self._Sx
        
    def delSy(self):
        del self._Sy
        
    def delSz(self):
        del self._Sz
    
    dimensions = property(getdimensions, setdimensions, deldimensions)
    Sx = property(getSx, setSx, delSx)
    Sy = property(getSy, setSy, delSy)
    Sz = property(getSz, setSz, delSz)
    numparticles = property(getnumparticles, setnumparticles, delnumparticles)

        
