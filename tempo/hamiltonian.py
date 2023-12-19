#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  5 22:01:39 2023

The Hamiltonian class for representing Hamiltonian operators.

@author: hakalas

"""

from qutip import qobj
import os
import numpy as np

# quality of life: if someone wants just a single matrix, should they still define a function like lambda x: x????
# if user doesn't input any function, just do lambda x: x by default
# can pass Hmats as qobj if only one matrix

# do we want wrapper functions for diagonalizing & eigenvalues?

# consider "add" method for adding together two Hamiltonians 

class Hamiltonian():
    """
    A class for representing time-independent Hamiltonian operators.
    
    The Hamiltonian object is composed of a numerical coefficient and an operator. The constructor takes a dictionary of parameters and a dictionary of operators, along with a callback function. The Python function tells the Hamiltonian how to combine the parameters and operators into a single term. For example, if the Hamiltonian term is :math:`a/bC\cdot D`, then 
    
    .. code-block:: python
        
        def func(ops, params):
            return params['a']/params['b']*ops['C']*ops['D']
            
        ops = {'C': C, 'D': D}
        params = {'a': a, 'b': b}
        
    Parameters
    ----------
    ops : dict of `qutip.Qobj` or `qutip.Qobj`
        Dictionary of all operators that make up the final Hamiltonian.
        If the Hamiltonian has no specified numerical coefficient and is only composed of one operator, the operator may be passed directly as a `qutip.Qobj`. In this case `params` and `func` are not needed as inputs.
    params: dict, optional
        Dictionary of numerical parameters that make up the coefficient of the final Hamiltonian.
    func : function, optional
        Callback function that takes the operators and parameters to combine.
        
    Attributes
    ----------
    ops : dict of `qutip.Qobj` or `qutip.Qobj`
        Dictionary of all operators that make up the final Hamiltonian.
    params: dict of float
        Dictionary of scalars that make up the coefficient of the final Hamiltonian.
    func : function
        Callback function that takes the operators and scalars to combine.  
    H : `qutip.Qobj`
        Time-independent Hamiltonian with a scalar coefficient part and an operator part evaluated via func.
    """
    
    def __init__(self, ops = {}, params = {}, func = None):
        """
        Hamiltonian constructor.
        """
        self.ops = ops
        self.params = params
        self.func = func
    
    @property
    def ops(self):
        return self._ops
    
    @ops.setter
    def ops(self, ops):
        
        # Check that ops is a dictionary or Qobj
        # if Qobj, assume the singular operator is returned in func
        
        if type(ops) == dict:
            self._ops = ops
        elif type(ops) == qobj.Qobj:
            self._ops = ops
            self._func = lambda op, pars: op
        else: 
            raise TypeError("Matrices must be in a dictionary")
    
    @ops.deleter
    def ops(self):
        del self._ops
    
    @property
    def params(self):
        return self._params
    
    
    @params.setter
    def params(self, params):
        
        if type(params) == dict:
            self._params = params
        elif params == None:
            self._params = {}
        else:
            raise TypeError("Parameters must be in a dictionary")
    
    @params.deleter
    def params(self):
        del self._params
      
    @property
    def func(self):
        return self._func
    
    @func.setter
    def func(self, func):
        
        if callable(func):
            self._func = func
        else:
            raise TypeError("Must be a function")
    
    @func.deleter
    def func(self):
        del self._func
    
    @property
    def H(self, ops = None, params = None):
        """
        Time-independent Hamiltonian term as a `qutip.Qobj`.
        
        Parameters
        ----------
        ops : dict of `qutip.Qobj`, optional
            Dictionary of operators to pass to the callback function. If None, then the dictionary passed in the constructor is used.
        params : dict of float, optional
            Dictionary of scalars to pass to the callback function. If None, then the dictionary passed in the constructor is used.
            
        Returns
        -------
        ham : `qutip.Qobj`
            Time-independent Hamiltonian with a scalar coefficient part and an operator part.
        """
        # Check if user has assigned a func;
        # if not, if ops only has a single matrix, the default function returns that matrix.
        # If func is unspecified and ops has more than one term, raise error.
        # add option to pass a single matrix for mats here too. in that case, it shouldn't raise an error when no func was defined.
        
        ops = self._ops if ops is None else ops
        params = self._params if params is None else params
        
        if self._func == None:
            raise TypeError("Hamiltonian function has not been specified")
        else:
            ham = self._func(ops, params)
            return ham
    
    @H.setter
    def H(self, H):
        if type(H) == qobj.Qobj:
            self._H = H
        else:
            raise TypeError("The full Hamiltonian must be a QuTiP Qobj object.")    
    
    @H.deleter
    def H(self):
        del self._H
