#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 10 18:20:45 2023

The Pulsetype class to create time-dependent pulse models. 

@author: hakalas
"""

import qutip
import os
import qsys
import ham
import numpy as np

# add type checks
# add setters and getters for each var 
# do checks on datatypes of each var in setters
# --> pulsepars only accepts LIST right now, change that later to also accept np array (and others?)


# pulsetype will probably change to not accept numerical values in pulsepars, but only a string of keys (var names) for
# pulsefunc; pulse object will input numerical values and evaluate the coefficient (?)
# can leave evalcoeff in pulsetype, but it will have to get numerical values as an input

class Pulsetype():
    """
    A class for defining types of pulses.
    
    The Pulsetype class is used to create models for different types of pulses or time-dependent Hamiltonians. The Pulse class is then used to create individual pulse instances. 
    
    The Pulsetype object is composed of a time-dependent coefficient and an operator. The constructor takes the operator, a list of scalar parameter names, and a callback function. The function, which calculates the time-dependent coefficient of the operator at a time t, must have a signature ``f(t: float, args: dict) -> float``, for example
    
    .. code-block:: python
        
        def func(t, args):
            return args['amp']*np.cos(2*np.pi*args['freq']*t)
            
    In this example the list of parameter names passed would be ``keys = ['amp', 'freq']``. The dictionary `args` is not part of the pulsetype class.
    
    Parameters 
    ----------
    pulsemat : :obj:`ham.Hamiltonian` or `qutip.Qobj`
        Hamiltonian operator for pulse.
    pulseparkeys: list of str
        List of names of scalar parameters used in `pulsefunc`. Names should be strings.
    pulsefunc : function
        Callback function that calculates the time-dependent coefficient using a time `t` and a dictionary of `args`.
        
    Attributes
    ----------
    pulsemat : `qutip.Qobj`
        Hamiltonian operator for pulse.
    pulseparkeys: list of str
        List of names of scalar parameters used in ``pulsefunc``. Names should be strings.
    pulsefunc : function
        Callback function that calculates the time-dependent coefficient using a time `t` and a dictionary of `args`.
    """
    
    def __init__(self, pulsemat = None, pulseparkeys = None, pulsefunc = None):
        """
        Pulsetype constructor.
        """
        self.setpulsemat(pulsemat) # Hamiltonian or Qobj input, but self._pulsemat is a Qobj
        self.setpulseparkeys(pulseparkeys) # list of strings (key names for dictionary of parameters)
        self.setpulsefunc(pulsefunc) # function
            
    def evalcoeff(self, t, pars):
        """
        Pulse coefficient at time `t`. 
        
        Parameters
        ----------
        t : float
            The time for which to evaluate the coefficient.
        pars : dict of float
            Parameters to use in `pulsefunc` to evaluate coefficient.
            
        Returns
        ----------
        coeff : float
            Pulse coefficient at time `t`.
        """
        # pars must be a dictionary!
        # ! check that input type = output type (array, float, etc)
        # t can be tuple, numpy array, list, float... but ppl need to make sure that their own function returns the right type 
        coeff = self._pulsefunc(t, pars)
        return coeff
        
    def getpulsemat(self):
        return self._pulsemat
    
    def getpulseparkeys(self):
        return self._pulseparkeys
    
    def getpulsefunc(self):
        return self._pulsefunc
    
    def setpulsemat(self, pulsemat):
        
        # pulsemat can be a Hamiltonian object or a quantum object
        # add also array option
        
        if type(pulsemat) == ham.Hamiltonian:
            self._pulsemat = pulsemat.getH()
        elif type(pulsemat) == qutip.qobj.Qobj:
            self._pulsemat = pulsemat
        elif pulsemat == None:
            self._pulsemat = pulsemat
        else: 
            print(type(pulsemat))
            raise TypeError("Matrix must be a Hamiltonian object or a Quantum object")
        
    def setpulseparkeys(self, pulseparkeys):
        
        # Check that pulseparkeys is a list
        
        if type(pulseparkeys) == list:
            self._pulseparkeys = pulseparkeys
        elif pulseparkeys == None:
            self._pulseparkeys = []
        else:
            raise TypeError("Keys must be in a list of strings")
        
    def setpulsefunc(self, pulsefunc):
        
        # check that Hfunc is a function
        
        if callable(pulsefunc):
            self._pulsefunc = pulsefunc
        else:
            raise TypeError("Must be a function")
        
    def delpulsemat(self):
        del self._pulsemat
        
    def delpulseparkeys(self):
        del self._pulseparkeys
        
    def delpulsefunc(self):
        del self._pulsefunc
    
    pulsemat = property(getpulsemat, setpulsemat, delpulsemat)
    pulseparkeys = property(getpulseparkeys, setpulseparkeys, delpulseparkeys)
    pulsefunc = property(getpulsefunc, setpulsefunc, delpulsefunc)

    
