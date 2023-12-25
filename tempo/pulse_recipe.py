#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 10 18:20:45 2023

The Pulse_recipe class to create time-dependent pulse models. 

@author: hakalas

"""
 
from qutip import qobj
import os
from tempo.hamiltonian import Hamiltonian
import numpy as np

class Pulse_recipe():
    """
    A class for defining types of pulses.
    
    The Pulse_recipe class is used to create models, or recipes, for different types of pulses or time-dependent Hamiltonians. The `pulse.Pulse` class is then used to create individual pulse instances. 
    
    The Pulse_recipe object is composed of a time-dependent coefficient and an operator. The constructor takes the operator, a list of parameter names, and a callback function. The function, which calculates the time-dependent coefficient of the operator at a time t, must have a signature ``f(t: float, args: dict) -> float``, for example
    
    .. code-block:: python
        
        def func(t, args):
            return args['amp']*np.cos(2*np.pi*args['freq']*t)
            
    In this example the list of parameter names passed would be ``keys = ['amp', 'freq']``. The dictionary `args` is not part of the Pulse_recipe class; it is incorporated in the `pulse.Pulse` class, where one can evaluate the numerical value of the callback function based on a particular dictionary of values. The callback function should not include checks for whether the pulse is on or off; this can be achieved by defining the particular pulse's start time and end time in the `pulse.Pulse` class.
    
    Parameters 
    ----------
    ham : :obj:`hamiltonian.Hamiltonian`
        Hamiltonian instance that describes the operator part of the pulse.
    param_keys: list of str
        List of names of parameters used in `coeff_func`. Names should be strings.
    coeff_func : function
        Callback function that calculates the time-dependent coefficient using a time `t` and a dictionary of `args`.
        
    Attributes
    ----------
    ham : :obj:`hamiltonian.Hamiltonian`
        Hamiltonian instance that describes the operator part of the pulse.
    param_keys: list of str
        List of names of parameters used in `coeff_func`. Names should be strings.
    coeff_func : function
        Callback function that calculates the time-dependent coefficient using a time `t` and a dictionary of `args`.
    """
    
    def __init__(self, ham = None, param_keys = None, coeff_func = None):
        """
        Pulse_recipe constructor.
        """
        self.ham = ham # Hamiltonian or Qobj
        self.param_keys = param_keys # list of strings (key names for dictionary of parameters)
        self.coeff_func = coeff_func # function
            
    def eval_coeff(self, t, params):
        """
        Pulse coefficient at time `t`. 
        
        Parameters
        ----------
        t : float
            The time for which to evaluate the coefficient.
        params : dict of float
            Parameters to use in `coeff_func` to evaluate coefficient.
            
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
    def ham(self):
        return self._ham
    
    @ham.setter
    def ham(self, ham):
        
        if isinstance(ham, Hamiltonian):
            self._ham = ham
        elif ham == None:
            self._ham = ham
        else: 
            print(type(ham))
            raise TypeError("ham must be a Hamiltonian instance")
    
    @ham.deleter
    def ham(self):
        del self._ham
      
    @property
    def param_keys(self):
        return self._param_keys
    
    @param_keys.setter
    def param_keys(self, param_keys):
        
        if isinstance(param_keys, list):
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
    
