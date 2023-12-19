#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 10 18:20:45 2023

The Pulserecipe class to create time-dependent pulse models. 

@author: hakalas

"""
 
from qutip import qobj
import os
from tempo.hamiltonian import Hamiltonian
import numpy as np

class Pulse_recipe():
    """
    A class for defining types of pulses.
    
    The Pulserecipe class is used to create models, or recipes, for different types of pulses or time-dependent Hamiltonians. The Pulse class is then used to create individual pulse instances. 
    
    The Pulserecipe object is composed of a time-dependent coefficient and an operator. The constructor takes the operator, a list of scalar parameter names, and a callback function. The function, which calculates the time-dependent coefficient of the operator at a time t, must have a signature ``f(t: float, args: dict) -> float``, for example
    
    .. code-block:: python
        
        def func(t, args):
            return args['amp']*np.cos(2*np.pi*args['freq']*t)
            
    In this example the list of parameter names passed would be ``keys = ['amp', 'freq']``. The dictionary `args` is not part of the Pulserecipe class; it is incorporated in the `Puls`' class, where one can evaluate the numerical value of the callback function based on a particular dictionary of values. The callback function should not include checks for whether the pulse is on or off; this can be achieved by defining the particular pulse's starttime and endtime in the `Pulse` class.
    
    Parameters 
    ----------
    op : :obj:`ham.Hamiltonian` or `qutip.Qobj`
        Hamiltonian operator for pulse.
    param_keys: list of str
        List of names of scalar parameters used in `coeff_func`. Names should be strings.
    coeff_func : function
        Callback function that calculates the time-dependent coefficient using a time `t` and a dictionary of `args`.
        
    Attributes
    ----------
    op : `qutip.Qobj`
        Hamiltonian operator for pulse.
    param_keys: list of str
        List of names of scalar parameters used in `coeff_func`. Names should be strings.
    coeff_func : function
        Callback function that calculates the time-dependent coefficient using a time `t` and a dictionary of `args`.
    """
    
    def __init__(self, op = None, param_keys = None, coeff_func = None):
        """
        Pulserecipe constructor.
        """
        self.op = op # Hamiltonian or Qobj input, but self._op is a Qobj
        self.param_keys = param_keys # list of strings (key names for dictionary of parameters)
        self.coeff_func = coeff_func # function
            
    def eval_coeff(self, t, params):
        """
        Pulse coefficient at time `t`. 
        
        Parameters
        ----------
        t : float
            The time for which to evaluate the coefficient.
        pars : dict of float
            Parameters to use in `pulsefunc` to evaluate coefficient.
            
        Returns
        -------
        coeff : float
            Pulse coefficient at time `t`.
        """
        # pars must be a dictionary!
        # ! check that input type = output type (array, float, etc)
        # t can be tuple, numpy array, list, float... but ppl need to make sure that their own function returns the right type 
        coeff = self._coeff_func(t, params)
        return coeff
     
    @property
    def op(self):
        return self._op
    
    @op.setter
    def op(self, op):
        
        # op can be a Hamiltonian object or a quantum object
        
        if type(op) == Hamiltonian:
            self._op = op.H
        elif type(op) == qobj.Qobj:
            self._op = op
        elif op == None:
            self._op = op
        else: 
            print(type(op))
            raise TypeError("Matrix must be a Hamiltonian object or a Quantum object")
    
    @op.deleter
    def op(self):
        del self._op
      
    @property
    def param_keys(self):
        return self._param_keys
    
    @param_keys.setter
    def param_keys(self, param_keys):
        
        if type(param_keys) == list:
            self._param_keys = param_keys
        elif param_keys == None:
            self._param_keys = []
        else:
            raise TypeError("Keys must be in a list of strings")
    
    @param_keys.deleter     
    def param_keys(self):
        del self._param_keys
    
    @property
    def coeff_func(self):
        return self._coeff_func
    
    @coeff_func.setter   
    def coeff_func(self, coeff_func):
        
        if callable(coeff_func):
            self._coeff_func = coeff_func
        else:
            raise TypeError("Must be a function")
        
    @coeff_func.deleter    
    def coeff_func(self):
        del self._coeff_func
    