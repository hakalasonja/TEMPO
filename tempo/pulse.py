#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Create time-dependent pulses.

@author: hakalas

"""

import os
import numpy as np
import qutip
from tempo.pulse_recipe import Pulse_recipe
from tempo.hamiltonian import Hamiltonian
from tempo.exceptions import *

from collections.abc import Iterable

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
    coeff_params : dict of float
        Dictionary of parameters to be passed in the time-dependent coefficient function `pulse_recipe.coeff_func`.
    """
    
    def __init__(self, recipe: Pulse_recipe, start_time: float = 0, duration: float = 0, coeff_params: dict = None):
        """
        Pulse constructor.
        """
        # check if self.pulsefunc == None then set it to lambda t, args: 1 
        
        if recipe is None:
            raise TEMPO_NullValueException("Recipe for pulse is null")
        if not isinstance(recipe, Pulse_recipe):
            raise TEMPO_ImproperInputException("Recipe is not a pulse_recipe.Pulse_recipe")
        if not castable(start_time, float):
            try: 
                start_time = float(start_time)
            except Exception:
                raise TEMPO_ImproperInputException("Start time is not a float")
        if not castable(start_time, float):
            try: 
                duration = float(duration)
            except Exception:
                raise TEMPO_ImproperInputException("Duration time is not a float")     
        if coeff_params is not None and type(coeff_params) != type({}):
            raise TEMPO_ImproperInputException("Coefficient parameters are not a dictionary")        

        self.recipe = recipe
        self._start_time = start_time
        self.duration = duration
        self._end_time = start_time + duration # --> Ensure that end time is set before possibly being called.
        self._coeff_params = {}
        
        if coeff_params == None:
            if len(self._recipe.param_keys) != 0:
                self.update_params(self._recipe.param_keys, [0]*len(self._recipe.param_keys))
        else:
            self._coeff_params = coeff_params
           
        self.ham = recipe.ham

        
    def update_params(self, keys: Iterable[str], values: Iterable[float]):
        """
        Add key-value pairs to `coeff_params`. `keys` and `values` should be the same length. Elements are paired by index. 
        
        Parameters
        ----------
        keys : list of str
            List of parameter labels.
        values : list of float
            List of parameter values.
        """
        if not isinstance(keys, Iterable):
            raise TEMPO_ImproperInputException("Keys are not an iterable")
        if not isinstance(values, Iterable):
            raise TEMPO_ImproperInputException("Values are not an iterable")
        if type(keys) == type({}):
            if sum([1 if not isinstance(x, str) else 0 for x in keys.values()]) != 0:
                raise TEMPO_ImproperInputException("Keys include a non string")
        else:
            if sum([1 if not isinstance(x, str) else 0 for x in keys]) != 0:
                raise TEMPO_ImproperInputException("Keys include a non string")
        if type(values) == type({}):
            if sum([1 if not castable(x, float) else 0 for x in values.values()]) != 0:
                raise TEMPO_ImproperInputException("Values include a non float")  
        else:
            if sum([1 if not isinstance(x, float) else 0 for x in values]) != 0:
                raise TEMPO_ImproperInputException("Values include a non float")  


        for i in range(len(keys)):
            self._coeff_params.update({keys[i]: values[i]})
        
    def eval_coeff(self, t: float, args: dict = {}) -> float:
        """
        Evaluate coefficient of pulse at time `t`. If `args` is not given, then `coeff_params` is used as parameter input for coefficient function. If the pulse is off at time t, the coefficient is 0.
        Warning that using args will mean coeff_params is entirely ignored this execution.

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
        if not castable(t, float):
            try:
                t = float(t)
            except Exception:
                raise TEMPO_ImproperInputException("Time cannot be cast to float")
        if type(args) != type({}):
            raise TEMPO_ImproperInputException("Args were not provided in dictionary form")
            
        if len(args) == 0:
            coeff_params = self._coeff_params
        else:
            coeff_params = args
        
        pulsecoeff = 0
        
        if t >= self._start_time and t <= self._end_time:
            f = self._recipe.coeff_func
            pulsecoeff = f(t, coeff_params)
        
        return pulsecoeff
    
    def serial_eval_coeff(self, t: float, args: dict = {}) -> float:
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
        if len(args) > 0:
            print_warning("Warning: Args is non-empty but these will NOT be used for serial evaluation.")

        f = self._recipe.coeff_func
        pulsecoeff = f(t, self._coeff_params)
        return pulsecoeff

    @property
    def recipe(self):
        """
        Hide on documentation page
        :meta private:
        """
        return self._recipe
    
    @recipe.setter
    def recipe(self, recipe):
        if isinstance(recipe, Pulse_recipe):
            self._recipe = recipe
        else:
            raise TEMPO_ImproperInputException('pulserecipe must be a member of the Pulserecipe class')

    @recipe.deleter
    def recipe(self):
        del self._recipe
    
    @property
    def start_time(self):
        """
        Hide on documentation page
        :meta private:
        """
        return self._start_time
    
    @start_time.setter
    def start_time(self, start_time):
        if not castable(start_time, float):
            try:
                start_time = float(start_time)
            except Exception:
                raise TEMPO_ImproperInputException("Time could not be cast to float")

        self._start_time = start_time
        
    @start_time.deleter
    def start_time(self):
        del self._start_time
    
    @property
    def duration(self):
        """
        Hide on documentation page
        :meta private:
        """
        return self._duration
    
    @duration.setter
    def duration(self, duration):
        self._duration = duration

        if not castable(duration, float):
            try:
                duration = float(duration)
            except Exception:
                raise TEMPO_ImproperInputException("Duration could not be cast to float")

        self._end_time = self._start_time+duration
        
    @duration.deleter
    def duration(self):
        del self._duration
    
    @property
    def end_time(self):
        """
        Hide on documentation page
        :meta private:
        """
        return self._end_time
    
    @end_time.setter
    def end_time(self, end_time):
        
        if end_time < self._start_time:
            raise TEMPO_ImproperInputException("Start time must be before end time")
        elif end_time < self._start_time + self._duration:
            # Note, the above if statement ensures that end time is already definitely greater than start time. 
            print_warning("Warning: End time is set before the duration would dictate. May cause undefined behavior.")
        else:
            self._end_time = end_time
            self._duration = end_time-self._start_time
     
    @end_time.deleter
    def end_time(self):
        del self._end_time
    
    @property
    def coeff_params(self):
        """
        Hide on documentation page
        :meta private:
        """
        return self._coeff_params

    @coeff_params.setter
    def coeff_params(self, coeff_params):

        if isinstance(coeff_params, dict):
            
            for key in self._recipe.param_keys:
                if key not in coeff_params.keys():
                    raise KeyError("{key} value not specified".format(key = key))
            
            self.update_params(list(coeff_params.keys()), list(coeff_params.values()))
            
        else: 
            raise TEMPO_ImproperInputException("Pulse coefficient parameters must be in a dictionary")
    
    @coeff_params.deleter
    def coeff_params(self):
        del self._coeff_params
        
    @property
    def ham(self):
        """
        Hide on documentation page
        :meta private:
        """
        return self._ham
    
    @ham.setter
    def ham(self, ham):
        
        if isinstance(ham, Hamiltonian):        
            self._ham = ham
        else:
            raise TEMPO_ImproperInputException("Pulse Hamiltonian must be a Hamiltonian instance")
            
    @ham.deleter
    def ham(self):
        del self._ham
