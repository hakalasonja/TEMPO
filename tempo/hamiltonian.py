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
        
        def func(ops, op_params):
            return op_params['a']/op_params['b']*ops['C']*ops['D']
            
        ops = {'C': C, 'D': D}
        op_params = {'a': a, 'b': b}
        
    Parameters
    ----------
    ops : dict of :obj:`qutip.Qobj` or :obj:`qutip.Qobj`
        Dictionary of all operators that make up the final Hamiltonian.
        If the Hamiltonian has no specified numerical coefficient and is only composed of one operator, the operator may be passed directly as a :obj:`qutip.Qobj`. In this case `params` and `func` are not needed as inputs.
    op_params: dict, optional
        Dictionary of numerical parameters that make up the coefficient of the final Hamiltonian.
    func : function, optional
        Callback function that takes the operators and parameters to combine.
        
    Attributes
    ----------
    ops : dict of :obj:`qutip.Qobj` or :obj:`qutip.Qobj`
        Dictionary of all operators that make up the final Hamiltonian.
    op_params: dict of float
        Dictionary of scalars that make up the coefficient of the final Hamiltonian.
    func : function
        Callback function that takes the operators and scalars to combine.  
    H : :obj:`qutip.Qobj`
        Time-independent Hamiltonian with a scalar coefficient part and an operator part evaluated via func.
    """
    
    def __init__(self, ops = {}, op_params = {}, func = None):
        """
        Hamiltonian constructor.
        """
        self.ops = ops
        
        if not isinstance(ops, qobj.Qobj):
            self.op_params = op_params
            self.func = func
    
    @property
    def ops(self):
        return self._ops
    
    @ops.setter
    def ops(self, ops):
        
        # Check that ops is a dictionary or Qobj
        # if Qobj, assume the singular operator is returned in func
        
        if isinstance(ops, dict):
            self._ops = ops
        elif isinstance(ops, qobj.Qobj):
            self._ops = ops
            self._op_params = {}
            self._func = lambda op, params: op
        else: 
            raise TypeError("Matrices must be in a dictionary")
    
    @ops.deleter
    def ops(self):
        del self._ops
    
    @property
    def op_params(self):
        return self._params
    
    
    @op_params.setter
    def op_params(self, op_params):
        
        if isinstance(op_params, dict):
            self._op_params = op_params
        elif op_params == None:
            self._op_params = {}
        else:
            raise TypeError("Parameters must be in a dictionary")
    
    @op_params.deleter
    def op_params(self):
        del self._op_params
      
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
    def H(self):
        """
        Time-independent Hamiltonian term as a :obj:`qutip.Qobj`.
        
        Parameters
        ----------
        ops : dict of :obj:`qutip.Qobj`, optional
            Dictionary of operators to pass to the callback function. If None, then the dictionary `ops` passed in the constructor is used.
        op_params : dict of float, optional
            Dictionary of numerical parameters to pass to the callback function. If None, then the dictionary `op_params` passed in the constructor is used.
            
        Returns
        -------
        H : :obj:`qutip.Qobj`
            Time-independent Hamiltonian with a numerical coefficient part and an operator part.
        """
        # Check if user has assigned a func;
        # if not, if ops only has a single matrix, the default function returns that matrix.
        # If func is unspecified and ops has more than one term, raise error.
        # add option to pass a single matrix for ops here too. in that case, it shouldn't raise an error when no func was defined.
        
        if self._func == None:
            raise ValueError("Hamiltonian function has not been specified")
        else:
            H = self._func(self._ops, self._op_params)
            return H
    
    @H.setter
    def H(self, H):
        if isinstance(H, qobj.Qobj):
            self._H = H
        else:
            raise TypeError("The full Hamiltonian must be a QuTiP Qobj object.")    
    
    @H.deleter
    def H(self):
        del self._H
