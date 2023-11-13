#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  5 22:01:39 2023

The Hamiltonian class for representing Hamiltonian operators.

@author: hakalas
"""

import qutip
import os
import numpy as np

# quality of life: if someone wants just a single matrix, should they still define a function like lambda x: x????
# if user doesn't input any function, just do lambda x: x by default
# can pass Hmats as qobj if only one matrix
# if Hmats is dict but only one entry, unwrap with list(Hmats.keys())[0] // DONE

# add type checks // done?
# add setters and getters for each var // DONE

# do we want wrapper functions for diagonalizing & eigenvalues?

# consider "add" method for adding together two Hamiltonians 

# add H as instance variable w/ getters, setters (getter is getH)

class Hamiltonian():
    """
    A class for representing time-independent Hamiltonian operators.
    
    The Hamiltonian object is composed of a scalar coefficient and an operator. The constructor takes a dictionary of scalar parameters and a dictionary of operators, along with a callback function. The Python function tells the Hamiltonian how to combine the scalars and operators into a single term. For example, if the Hamiltonian term is :math:`a/bC\cdot D`, then 
    
    .. code-block:: python
        
        def func(Hmats, Hpars):
            return Hpars['a']/Hpars['b']*Hmats['C']*Hmats['D']
            
        Hpars = {'a': a, 'b': b}
        Hmats = {'C': C, 'D': D}
        
    Parameters
    ----------
    Hmats : dict of `qutip.Qobj`, `qutip.Qobj`, or `np.ndarray`
        Dictionary of all operators that make up the final Hamiltonian.
        If the Hamiltonian has no scalar coefficient and is only composed of one operator, the operator may be passed directly as a `qutip.Qobj` or a `np.ndarray`. In this case `Hpars` and `func` are not needed as inputs.
    Hpars: dict, optional
        Dictionary of scalars that make up the coefficient of the final Hamiltonian.
    func : function, optional
        Callback function that takes the operators and scalars to combine.
        
    Attributes
    ----------
    Hmats : dict of `qutip.Qobj`, `qutip.Qobj`, or `np.ndarray`
        Dictionary of all operators that make up the final Hamiltonian.
    Hpars: dict of float
        Dictionary of scalars that make up the coefficient of the final Hamiltonian.
    func : function
        Callback function that takes the operators and scalars to combine.  
    H : `qutip.Qobj`
        Time-independent Hamiltonian with a scalar coefficient part and an operator part.
    """
    
    def __init__(self, Hmats = {}, Hpars = {}, func = None):
        """
        Hamiltonian constructor.
        """
        self.setHmats(Hmats)
        self.setHpars(Hpars)
        self.setHfunc(func)
    
    def getHmats(self):
        return self._Hmats
    
    def getHpars(self):
        return self._Hpars
    
    def getHfunc(self):
        return self._Hfunc
    
    def getH(self, mats = None, pars = None):
        """
        Time-independent Hamiltonian term as a `qutip.Qobj`.
        
        Parameters
        ----------
        mats : dict of `qutip.Qobj`, optional
            Dictionary of operators to pass to the callback function. If None, then the dictionary passed in the constructor is used.
        pars : dict of float, optional
            Dictionary of scalars to pass to the callback function. If None, then the dictionary passed in the constructor is used.
            
        Returns
        ----------
        ham : `qutip.Qobj`
            Time-independent Hamiltonian with a scalar coefficient part and an operator part.
        """
        # Check if user has assigned a Hfunc;
        # if not, if Hmats only has a single matrix, the default function returns that matrix.
        # If Hfunc is unspecified and Hmats has more than one term, raise error.
        # add option to pass a single matrix for mats here too. in that case, it shouldn't raise an error when no func was defined.
        
        mats = self._Hmats if mats is None else mats
        pars = self._Hpars if pars is None else pars
        
        if self._Hfunc == None:
            raise TypeError("Hamiltonian function has not been specified")
        else:
            ham = self._Hfunc(mats, pars)
            return ham
        
    def setHmats(self, Hmats):
        
        # Check that Hmats is a dictionary.
        # If Hmats is not a dictionary but a Qobj, wrap Hmats in a dictionary.
        # can also be array
        
        if type(Hmats) == dict:
            self._Hmats = Hmats
        elif type(Hmats) == qutip.qobj.Qobj:
            self._Hmats = Hmats
            self._Hfunc = lambda mat, pars: mat
        elif type(Hmats) == ndarray:
            self._Hmats = Hmats
            self._Hfunc = lambda mat, pars: mat
        else: 
            raise TypeError("Matrices must be in a dictionary")
        
    def setHpars(self, Hpars):
        
        if type(Hpars) == dict:
            self._Hpars = Hpars
        elif Hpars == None:
            self._Hpars = {}
        else:
            raise TypeError("Parameters must be in a dictionary")

    def setHfunc(self, Hfunc):
        
        if callable(Hfunc):
            self._Hfunc = Hfunc
        else:
            raise TypeError("Must be a function")
            
    def setH(self, H):
        if type(H) == qutip.qobj.Qobj:
            self._H = H
        else:
            raise TypeError("The full Hamiltonian must be a QuTiP Qobj object.")
        
    def delHmats(self):
        del self._Hmats
        
    def delHpars(self):
        del self._Hpars
        
    def delHfunc(self):
        del self._Hfunc
        
    def delH(self):
        del self._H

    Hmats = property(getHmats, setHmats, delHmats)
    Hpars = property(getHpars, setHpars, delHpars)
    Hfunc = property(getHfunc, setHfunc, delHfunc)
    H = property(getH, setH, delH)