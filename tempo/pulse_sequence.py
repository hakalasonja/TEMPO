#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Sequence (list) of pulses.

@author: hakalas

"""

from qutip import qobj
from tempo.hamiltonian import Hamiltonian
from tempo.pulse import Pulse
import numpy as np

from typing import Union
from collections.abc import Iterable
from tempo.exceptions import *


class Pulse_sequence():
    """
    Class for creating sequences or lists of pulses.
    
    A Pulse sequence object stores all pulses applied to a system during a single simulation. If there is a static (time-independent) Hamiltonian that applies throughout the entire simulation, this should be passed to the constructor separately. 
    
    The pulses may be added all at once in a list in the constructor, or they may be added one by one after the initialization of the pulse sequence.
    
    Parameters
    ----------
    pulses : list of :obj:`pulse.Pulse` 
        List of :obj:`pulse.Pulse` instances that make up the pulse sequence. Pulses do not have to be ordered in any way.
    Hstat : :obj:`hamiltonian.Hamiltonian` or :obj:`qutip.Qobj`, optional
        Time-independent Hamiltonian that applies at all times in the simulation.
    
    Attributes
    ----------
    pulses : list of :obj:`pulse.Pulse`
        List of :obj:`pulse.Pulse` instances that make up the pulse sequence. 
    Hstat : :obj:`qutip.Qobj`
        Time-independent Hamiltonian that applies at all times in the simulation.
    """
    
    def __init__(self, pulses: Iterable[Pulse] = None, Hstat: Union[Hamiltonian, qobj.Qobj] = None):
        """
        Constructor for storing pulse sequences (lists of pulses).

        Parameters
        ----------
        pulses : list of :obj:`pulse.Pulse` 
            List of :obj:`pulse.Pulse` instances that make up the pulse sequence. Pulses do not have to be ordered in any way.
        Hstat : :obj:`hamiltonian.Hamiltonian` or :obj:`qutip.Qobj`, optional
            Time-independent Hamiltonian that applies at all times in the simulation.
        """
        if pulses is not None and not isinstance(pulses, Iterable):
            raise TEMPO_ImproperInputException("Pulses is not an iterable")
        if pulses is not None:
            if type(pulses) == type({}):
                if sum([1 if not isinstance(x, Pulse) else 0 for x in pulses.values()]) != 0:
                    raise TEMPO_ImproperInputException("Pulses includes a non Pulse object")
            else:
                if sum([1 if not isinstance(x, Pulse) else 0 for x in pulses]) != 0:
                    raise TEMPO_ImproperInputException("Pulses includes a non Pulse object")
        if Hstat is not None:
            if not isinstance(Hstat, qobj.Qobj) and not isinstance(Hstat, Hamiltonian):
                raise TEMPO_ImproperInputException("The provided time-independent Hamiltonian is not a Qobj or Hamiltonian")

        self.pulses = pulses
        self.Hstat = Hstat
        
    def add_pulse(self, pls: Union[Pulse, Iterable[Pulse]]):
        """
        Add a pulse or list of pulses to the sequence.
        
        Parameters
        ----------
        pls : :obj:`pulse.Pulse` or list of :obj:`pulse.Pulse`
            Pulse(s) to be added. 
        """
        if isinstance(pls, Pulse):
            self._pulses.append(pls)
        elif isinstance(pls, list):
            for p in pls: 
                if isinstance(p, Pulse):
                    self._pulses.append(p)
                else:
                    raise TEMPO_ImproperInputException("All pulses in list must be Pulse instances")
        else: 
            raise TEMPO_ImproperInputException("Pulses must be Pulse instances. Multiple pulses must be in a list") 
               
    def remove_pulse(self, pls: Union[Pulse, Iterable[Pulse]]):
        """
        Remove a pulse or list of pulses from the sequence. If there are multiple copies of the same pulse in the sequence, then only the first instance is removed.
        
        Parameters
        ----------
        pls : :obj:`pulse.Pulse` or list of :obj:`pulse.Pulse`
            Pulse(s) to be removed. 
        """
        if isinstance(pls, Pulse):
            self._pulses.remove(pls)
        elif isinstance(pls, list):
            for p in pls: 
                if isinstance(p, Pulse):
                    self._pulses.remove(p)
                else:
                    raise TypeError("All pulses in list must be Pulse instances")
        else: 
            raise TypeError("Pulses must be Pulse instances. Multiple pulses must be in a list") 
        
    def timing_info(self) -> dict:
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
        """
        Hide on documentation page
        
        :meta private:

        """
        return self._pulses
    
    @pulses.setter
    def pulses(self, pulses):
        if isinstance(pulses, list):
            self._pulses = pulses
        elif isinstance(pulses, Pulse):
            self._pulses = [pulses]
        elif pulses == None:
            self._pulses = []
        else: 
            raise TEMPO_ImproperInputException("Pulses must be in a list")
            # Checking for individual pulses will not be repeated. It's done when the pulses are added...
    
    @pulses.deleter
    def pulses(self):
        del self._pulses
    
    @property
    def Hstat(self):
        """
        Hide on documentation page
        
        :meta private:

        """
        return self._Hstat
    
    @Hstat.setter        
    def Hstat(self, Hstat):
        if Hstat == None:
            self._Hstat = Hstat
        elif isinstance(Hstat, Hamiltonian):
            self._Hstat = Hstat.H
        elif isinstance(Hstat, qobj.Qobj):
            self._Hstat = Hstat
        else: 
            raise TEMPO_ImproperInputException("Operator must be a Hamiltonian instance or a QuTiP Qobj instance")
    
    @Hstat.deleter
    def Hstat(self):
        del self._Hstat
        
        
