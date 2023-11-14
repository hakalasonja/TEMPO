#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 17 19:43:07 2023

Create time-dependent pulses.

@author: hakalas
"""

import qutip
import os
from tempo.qsys import Qsys
from tempo.ham import Hamiltonian
import numpy as np
from tempo.pulsetype import Pulsetype
 
# add update method: update dictionary using input pulseparams and the string keys from pulsetp
# update only in init
# in pulsepars setter, check that all keys of input dict are included in string keys
# initialize pulsepars dictionary with default values (0)

class Pulse(pulsetype.Pulsetype):
    """
    A class for creating pulses or time-dependent Hamiltonian objects.
    
    Pulse extends Pulsetype, but there is a clear distinction between the two. Pulsetype can be summarized as a blueprint, while individual Pulses are created according to the blueprint. Any number of pulses can be created with the same pulsetype, each with different numerical values for the input parameters of the pulse coefficient. Because of this, a :obj:`pulsetype` object must be passed to the Pulse constructor. 
    
    A dictionary of parameters is also passed in the constructor and used to evaluate the time-dependent coefficient of the pulse. The dictionary should contain all of the entries listed in the pulsetype's `pulseparkeys` list of labels.  
    
    Parameters
    ----------
    pulsetp : :obj:`pulsetype.Pulsetype`
        Pulsetype object to provide a model for the kind of pulse.
    starttime : float
        Start time of pulse.
    duration : float
        Duration of pulse. 
    pulsepars : dict of float
        Dictionary of parameters to be passed in the coefficient function of `pulsetp`. All entries of `pulsetp`'s `pulseparkeys` must be included in the dictionary. 
        
    Attributes
    ----------
    pulsemat : `qutip.Qobj`
        Hamiltonian operator for pulse.
    pulseparkeys: list of str
        List of names of scalar parameters used in `pulsefunc`. Names are strings.
    pulsefunc : function
        Callback function that calculates the time-dependent coefficient using a time `t` and a dictionary of `args`. 
    starttime : float
        Start time of pulse.
    duration : float
        Duration of pulse. 
    endtime : float
        End time of pulse.
    pulsepars : dict of float
        Dictionary of parameters to be passed in the coefficient function `pulsefunc`.
    """
    
    def __init__(self, pulsetp, starttime = 0, duration = 0, pulsepars = None):
        """
        Pulse constructor.
        """
        # check if self.pulsefunc == None then set it to lambda t, args: 1 
        super().__init__(pulsetp.getpulsemat(), pulsetp.getpulseparkeys(), pulsetp.getpulsefunc())
        
        self._starttime = starttime
        self.setduration(duration)
        self._pulsepars = {}
        self.setendtime(starttime+duration)
        
        if pulsepars == None:
            if len(self._pulseparkeys) != 0:
                self.updatepars(self._pulseparkeys, [0]*len(self._pulseparkeys))
        else:
            self.setpulsepars(pulsepars)
        
    def updatepars(self, keys, values):
        """
        Add key-value pairs to `pulsepars`. `keys` and `values` should be the same length. Elements are paired by index. 
        
        Parameters
        ----------
        keys : list of str
            List of parameter labels.
        values : list of float
            List of parameter values.
        """
        for i in range(len(keys)):
            self._pulsepars.update({keys[i]: values[i]})
        
    def evalcoeff(self, t, args = {}):
        """
        Evaluate coefficient of pulse at time `t`. Overrides :obj:`Pulsetype`'s method. If `args` is None then `pulsepars` is used as parameter input for coefficient function. 
        
        Parameters
        ----------
        t : float
            The time for which to evaluate the coefficient.
        args : dict of float
            Parameters to use in `pulsefunc` to evaluate coefficient. If None, `pulsepars` is used.
            
        Returns
        ----------
        pulsecoeff : float
            Pulse coefficient at time `t`.
        """
        # pulsepars must be a dictionary
        # tlist can be a list or a scalar (only scalar right now)

        if len(args) == 0:
            pars = self._pulsepars
        else:
            pars = args
        
        pulsecoeff = 0
        
        if t >= self._starttime and t <= self._endtime:
            pulsecoeff = self._pulsefunc(t, pars)
        
        return pulsecoeff
    
    def serial_evalcoeff(self, t, args = {}):
        return self._pulsefunc(t, self._pulsepars)
        
    
#     def checkpulsefunc(self, t = 0):
        
#         # read about type hints
        
#         if type(t) == type(self.evalcoeff(t)):
#             return 1
#         else:
#             raise TypeError("")

    def getstarttime(self):
        return self._starttime
    
    def getduration(self):
        return self._duration
    
    def getendtime(self):
        return self._endtime
    
    def getpulsepars(self):
        return self._pulsepars
    
    def setstarttime(self, starttime):
        self._starttime = starttime
        
    def setduration(self, duration):
        self._duration = duration
        self._endtime = self._starttime+duration
        
    def setendtime(self, endtime):
        
        if endtime < self._starttime:
            raise ValueError("Start time must be before end time")
        else:
            self._endtime = endtime
            self._duration = endtime-self._starttime
        
    def setpulsepars(self, pulsepars):
        
        # finish exception 
        
        if type(pulsepars) == dict:
            
            for key in self._pulseparkeys:
                if key not in pulsepars.keys():
                    raise KeyError("{key} value not specified".format(key = key))
            
            self.updatepars(list(pulsepars.keys()), list(pulsepars.values()))
            
        else: 
            raise TypeError("Pulse parameters must be in a dictionary")
        
    def delstarttime(self):
        del self._starttime
        
    def delduration(self):
        del self._duration
        
    def delendtime(self):
        del self._endtime
        
    def delpulsepars(self):
        del self._pulsepars
        
    starttime = property(getstarttime, setstarttime, delstarttime)
    duration = property(getduration, setduration, delduration)
    endtime = property(getendtime, setendtime, delendtime)
    pulsepars = property(getpulsepars, setpulsepars, delpulsepars)
    
