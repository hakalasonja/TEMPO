#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 17 19:43:07 2023

Create time-dependent pulses.

@author: hakalas

"""

import os
import numpy as np
from tempo.pulse_recipe import Pulse_recipe
from tempo.hamiltonian import Hamiltonian

class Pulse():
    """
    A class for creating pulses or time-dependent Hamiltonian objects.
    
    One of the attributes of :obj:`Pulse` is `recipe`, which is an instance of :obj:`pulse_recipe.Pulse_recipe`, but there is a clear distinction between the two. Pulse recipe can be summarized as a blueprint (a recipe), while individual pulses are created according to the blueprint. Any number of pulses can be created with the same pulse recipe, but they can differ in the numerical values for the input parameters of the time-dependent coefficient. 
    
    A dictionary of parameters is passed in the constructor and used to evaluate the time-dependent coefficient of the pulse. The dictionary should contain all of the entries listed in the pulse recipe's `recipe.param_keys` list of labels.  
    
    Parameters
    ----------
    recipe : :obj:`pulse_recipe.Pulse_recipe`
        Pulse_recipe instance to provide a model for the kind of pulse.
    start_time : float
        Start time of pulse.
    duration : float
        Duration of pulse. 
    params : dict of float
        Dictionary of parameters to be passed in the coefficient function of `recipe`. All entries of `recipe`'s `param_keys` must be included in the dictionary. 
        
    Attributes
    ----------
    recipe : :obj:`pulse_recipe.Pulse_recipe`
        :obj:`pulse_recipe.Pulse_recipe` instance to provide a model for the kind of pulse.
    start_time : float
        Start time of pulse.
    duration : float
        Duration of pulse. 
    end_time : float
        End time of pulse.
    coeff_params : dict of float
        Dictionary of parameters to be passed in the time-dependent coefficient function `pulse_recipe.coeff_func`.
    ham : :obj:`hamiltonian.Hamiltonian`
        The operator part of the time-dependent term.
    """
    
    def __init__(self, recipe, start_time = 0, duration = 0, coeff_params = None):
        """
        Pulse constructor.
        """
        # check if self.pulsefunc == None then set it to lambda t, args: 1 
        
        self.recipe = recipe
        self._start_time = start_time
        self.duration = duration
        self._coeff_params = {}
        
        if coeff_params == None:
            if len(self._recipe.param_keys) != 0:
                self.update_params(self._recipe.param_keys, [0]*len(self._recipe.param_keys))
        else:
            self._coeff_params = coeff_params
           
        self.ham = recipe.ham

        
    def update_params(self, keys, values):
        """
        Add key-value pairs to `coeff_params`. `keys` and `values` should be the same length. Elements are paired by index. 
        
        Parameters
        ----------
        keys : list of str
            List of parameter labels.
        values : list of float
            List of parameter values.
        """
        for i in range(len(keys)):
            self._coeff_params.update({keys[i]: values[i]})
        
    def eval_coeff(self, t, args = {}):
        """
        Evaluate coefficient of pulse at time `t`. If `args` is not given, then `coeff_params` is used as parameter input for coefficient function. If the pulse is off at time t, the coefficient is 0.
        
        Parameters
        ----------
        t : float
            The time for which to evaluate the coefficient.
        args : dict of float, optional
            Parameters to use in `recipe.coeff_func` to evaluate coefficient. If not given, `coeff_params` is used.
            
        Returns
        -------
        pulsecoeff : float
            Pulse coefficient at time `t`.
        """
        # params must be a dictionary
        # tlist can be a list or a scalar (only scalar right now)

        if len(args) == 0:
            coeff_params = self._coeff_params
        else:
            coeff_params = args
        
        pulsecoeff = 0
        
        if t >= self._start_time and t <= self._end_time:
            f = self._recipe.coeff_func
            pulsecoeff = f(t, coeff_params)
        
        return pulsecoeff
    
    def serial_eval_coeff(self, t, args = {}):
        """
        Evaluate coefficient of pulse at time `t`. Differs from `eval_coeff` in that the start- and endtimes of the pulse are not taken into account; the coefficient will be evaluated whether the pulse is on or off. This method is called by the method `evolver.serial_evolve` in the :obj:`evolver.Evolver` class. 
        
        Parameters
        ----------
        t : float
            The time for which to evaluate the coefficient.
        args : dict of float, default = {}
            Always an empty dictionary. `coeff_params` is used to evaluate the coefficient via `recipe.coeff_func`
            
        Returns
        -------
        coeff : float
            Pulse coefficient at time `t`.
        """ 
        f = self._recipe.coeff_func
        pulsecoeff = f(t, self._coeff_params)
        return pulsecoeff

    @property
    def recipe(self):
        return self._recipe
    
    @recipe.setter
    def recipe(self, recipe):
        if isinstance(recipe, Pulse_recipe):
            self._recipe = recipe
        else:
            raise TypeError('pulserecipe must be a member of the Pulserecipe class')

    @recipe.deleter
    def recipe(self):
        del self._recipe
    
    @property
    def start_time(self):
        return self._start_time
    
    @start_time.setter
    def start_time(self, start_time):
        self._start_time = start_time
        
    @start_time.deleter
    def start_time(self):
        del self._start_time
    
    @property
    def duration(self):
        return self._duration
    
    @duration.setter
    def duration(self, duration):
        self._duration = duration
        self._end_time = self._start_time+duration
        
    @duration.deleter
    def duration(self):
        del self._duration
    
    @property
    def end_time(self):
        return self._end_time
    
    @end_time.setter
    def end_time(self, end_time):
        
        if end_time < self._start_time:
            raise ValueError("Start time must be before end time")
        else:
            self._end_time = end_time
            self._duration = end_time-self._start_time
     
    @end_time.deleter
    def end_time(self):
        del self._end_time
    
    @property
    def coeff_params(self):
        return self._coeff_params

    @coeff_params.setter
    def coeff_params(self, coeff_params):

        if isinstance(coeff_params, dict):
            
            for key in self._recipe.param_keys:
                if key not in coeff_params.keys():
                    raise KeyError("{key} value not specified".format(key = key))
            
            self.update_params(list(coeff_params.keys()), list(coeff_params.values()))
            
        else: 
            raise TypeError("Pulse coefficient parameters must be in a dictionary")
    
    @coeff_params.deleter
    def coeff_params(self):
        del self._coeff_params
        
    @property
    def ham(self):
        return self._ham
    
    @ham.setter
    def ham(self, ham):
        
        if isinstance(ham, Hamiltonian):        
            self._ham = ham
        else:
            raise TypeError("Pulse Hamiltonian must be a Hamiltonian instance")
            
    @ham.deleter
    def ham(self):
        del self._ham
