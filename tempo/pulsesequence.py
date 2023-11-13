#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 24 20:02:47 2023

Sequence (list) of pulses.

@author: hakalas
"""

import qutip
import os
import ham
import pulse
import numpy as np

# TODO: option to add multiple pulses at once in addpulse method
# look into making input types more general; list --> array_like (consider iterable)

class Pulsesequence():
    """
    Class for creating sequences or lists of pulses.
    
    A Pulsesequence object stores all pulses applied to a system during a single simulation. If there is a static (time-independent) Hamiltonian that applies throughout the entire simulation, this can be passed to the constructor separately. 
    
    The pulses may be added all at once in a list in the constructor, or they may be added one by one after the initialization of the pulse sequence.
    
    Parameters
    ----------
    pulsels : list of :obj:`pulse.Pulse` 
        List of :obj:`pulse.Pulse` objects that make up the pulse sequence. Pulses do not have to be ordered in any way.
    statham : :obj:`ham.Hamiltonian` or `qutip.Qobj`, optional
        Time-independent Hamiltonian that applies at all times in the simulation.
    
    Attributes
    ----------
    pulsels : list of :obj:`pulse.Pulse`
        List of :obj:`pulse.Pulse` objects that make up the pulse sequence. 
    statham : `qutip.Qobj`
        Time-independent Hamiltonian that applies at all times in the simulation.
    """
    
    def __init__(self, pulsels = None, statham = None):
        
        self.setpulsels(pulsels)
        self.setstatham(statham)
        
    def addpulse(self, pls):
        """
        Add a pulse or list of pulses to the sequence.
        
        Parameters
        ----------
        pulse : :obj:`pulse.Pulse` or list of :obj:`pulse.Pulse`
            Pulse(s) to be added. 
        """
        if type(pls) == pulse.Pulse:
            self._pulsels.append(pls)
        elif type(pls) == list:
            if len(pls) == 1:
                if type(pls[0]) == pulse.Pulse:
                    self._pulsels.append(pls[0])
                else:
                    raise TypeError("Pulse must be a Pulse object")
            else: 
                for p in pls: 
                    if type(p) == pulse.Pulse:
                        self._pulsels.append(p)
                    else:
                        raise TypeError("All pulses in list must be Pulse objects")
        else: 
            raise TypeError("Pulses must be Pulse objects. Multiple pulses must be in a list") 
               
    def printtiminginfo(self):
        """
        Print the start time, duration, and end time of the entire pulse sequence. Start time is the start time of the earliest pulse in the sequence, and end time is accordingly the end time of the last-ending pulse in the sequence. Note that the simulation itself may be longer than this with empty or static pulses outside of the pulse sequence time range.
        """
        starttime = self._pulsels[0].getstarttime()
        duration = self._pulsels[0].getduration()
        endtime = self._pulsels[0].getendtime()
        
        for i in np.arange(len(self._pulsels)):
            if self._pulsels[i].getstarttime() < starttime:
                starttime = self._pulsels[i].getstarttime()
                duration = endtime-starttime
            if self._pulsels[i].getendtime() > endtime: 
                endtime = self._pulsels[i].getendtime()
                duration = endtime-starttime
        
        print("Starttime: ", starttime, ", Duration: ", duration, ", Endtime: ", endtime)
        
    def getpulsels(self):
        return self._pulsels
    
    def getstatham(self):
        return self._statham
    
    def setpulsels(self, pulsels):
        if type(pulsels) == list:
            self._pulsels = pulsels
        elif type(pulsels) == pulse.Pulse:
            self._pulsels = [pulsels]
        elif pulsels == None:
            self._pulsels = []
        else: 
            raise TypeError("Pulses must be in a list")
            
    def setstatham(self, statham):
        
        if statham == None:
            self._statham = statham
        elif type(statham) == ham.Hamiltonian:
            self._statham = statham.getH()
        elif type(statham) == qutip.qobj.Qobj:
            self._statham = statham
        else: 
            print(type(statham))
            raise TypeError("Operator must be a Hamiltonian object or a Quantum object")
            
    def delpulsels(self):
        del self._pulsels
        
    def delstatham(self):
        del self._statham
              
    pulsels = property(getpulsels, setpulsels, delpulsels)  
    statham = property(getstatham, setstatham, delstatham)
        
