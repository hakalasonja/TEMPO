#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 24 20:02:47 2023

Sequence (list) of pulses.

@author: hakalas

"""

from qutip import qobj
import os
from tempo.hamiltonian import Hamiltonian
from tempo.pulse import Pulse
import numpy as np

# look into making input types more general; list --> array_like (consider iterable)

class Pulse_sequence():
    """
    Class for creating sequences or lists of pulses.
    
    A Pulsesequence object stores all pulses applied to a system during a single simulation. If there is a static (time-independent) Hamiltonian that applies throughout the entire simulation, this should be passed to the constructor separately. 
    
    The pulses may be added all at once in a list in the constructor, or they may be added one by one after the initialization of the pulse sequence.
    
    Parameters
    ----------
    pulses : list of :obj:`pulse.Pulse` 
        List of :obj:`pulse.Pulse` objects that make up the pulse sequence. Pulses do not have to be ordered in any way.
    Hstat : :obj:`hamiltonian.Hamiltonian` or `qutip.Qobj`, optional
        Time-independent Hamiltonian that applies at all times in the simulation.
    
    Attributes
    ----------
    pulses : list of :obj:`pulse.Pulse`
        List of :obj:`pulse.Pulse` objects that make up the pulse sequence. 
    Hstat : `qutip.Qobj`
        Time-independent Hamiltonian that applies at all times in the simulation.
    """
    
    def __init__(self, pulses = None, Hstat = None):
        
        self.pulses = pulses
        self.Hstat = Hstat
        
    def add_pulse(self, pls):
        """
        Add a pulse or list of pulses to the sequence.
        
        Parameters
        ----------
        pls : :obj:`pulse.Pulse` or list of :obj:`pulse.Pulse`
            Pulse(s) to be added. 
        """
        if type(pls) == Pulse:
            self._pulses.append(pls)
        elif type(pls) == list:
            for p in pls: 
                if type(p) == Pulse:
                    self._pulses.append(p)
                else:
                    raise TypeError("All pulses in list must be Pulse objects")
        else: 
            raise TypeError("Pulses must be Pulse objects. Multiple pulses must be in a list") 
               
    def remove_pulse(self, pls):
        """
        Remove a pulse or list of pulses from the sequence. If there are multiple copies of the same pulse in the sequence, then only the first instance is removed.
        
        Parameters
        ----------
        pls : :obj:`pulse.Pulse` or list of :obj:`pulse.Pulse`
            Pulse(s) to be removed. 
        """
        if type(pls) == Pulse:
            self._pulses.remove(pls)
        elif type(pls) == list:
            for p in pls: 
                if type(p) == Pulse:
                    self._pulses.remove(p)
                else:
                    raise TypeError("All pulses in list must be Pulse objects")
        else: 
            raise TypeError("Pulses must be Pulse objects. Multiple pulses must be in a list") 
        
    def timing_info(self):
        """
        Get a dictionary with the start time, duration, and end time of the entire pulse sequence. Start time is the start time of the earliest pulse in the sequence, and end time is accordingly the end time of the last-ending pulse in the sequence. Note that the simulation itself may span a longer time duration than this if the times at which the state of the system is evaluated extend outside of the pulse sequence time range; for example, there may be a duration of time at the beginning or end of the simulation where only the static Hamiltonian applies. 
        """
        start_time = self._pulses[0].start_time
        duration = self._pulses[0].duration
        end_time = self._pulses[0].end_time
        
        for i in np.arange(len(self._pulses)):
            if self._pulses[i].start_time < start_time:
                start_time = self._pulses[i].start_time
                duration = end_time-start_time
            if self._pulses[i].end_time > end_time: 
                end_time = self._pulses[i].end_time
                duration = end_time-start_time
        
        # return a dictionary
        return {'Starttime': start_time, 'Duration': duration, 'Endtime': end_time}
     
    @property
    def pulses(self):
        return self._pulses
    
    @pulses.setter
    def pulses(self, pulses):
        if type(pulses) == list:
            self._pulses = pulses
        elif type(pulses) == Pulse:
            self._pulses = [pulses]
        elif pulses == None:
            self._pulses = []
        else: 
            raise TypeError("Pulses must be in a list")
    
    @pulses.deleter
    def pulses(self):
        del self._pulses
    
    @property
    def Hstat(self):
        return self._Hstat
    
    @Hstat.setter        
    def Hstat(self, Hstat):
        
        if Hstat == None:
            self._Hstat = Hstat
        elif type(Hstat) == Hamiltonian:
            self._Hstat = Hstat.H
        elif type(Hstat) == qobj.Qobj:
            self._Hstat = Hstat
        else: 
            print(type(Hstat))
            raise TypeError("Operator must be a Hamiltonian object or a Quantum object")
    
    @Hstat.deleter
    def Hstat(self):
        del self._Hstat
        
        
