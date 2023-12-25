#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 24 20:10:22 2023

Class to implement time-evolution of a state under applied pulses.

@author: hakalas

"""

from qutip import mesolve, solver, qobj, expect
import os
from tempo.hamiltonian import Hamiltonian
from tempo.pulse_sequence import Pulse_sequence
import numpy as np
from math import isclose

class Evolver():
    """
    Time-evolution of a state under a pulse sequence.
    
    The Evolver object stores information about a particular pulse sequence simulation. Its `evolve` method evolves the given initial state `state_init` when the pulses in `pulse_seq` are applied. `evolve` relies on QuTiP's `qutip.mesolve` method and the inputs `c_ops`, `e_ops`, and `opts`, if defined, are directly passed to `mesolve`. For these, see `mesolve` documentation https://qutip.org/docs/latest/apidoc/functions.html#module-qutip.mesolve. 
    
    Parameters 
    ----------
    state_init : `qutip.Qobj`
        Initial state vector or density matrix.
    pulse_seq : :obj:`pulse_sequence.Pulse_sequence`
        Sequence of pulses to be applied to initial state.
    tlist : list of float
        List of times at which to evaluate the system's state or the expectation value of the operator(s) in `e_ops`. 
    Hstat : :obj:`hamiltonian.Hamiltonian` or `qutip.Qobj`, optional
        Time-independent Hamiltonian that applies at all times in the simulation.
    c_ops : list of `qutip.Qobj`
        List of collapse operators.
    e_ops : list of `qutip.Qobj`
        List of operators for which to evaluate expectation values.
    opts : `qutip.Options`
        Options for `qutip.mesolve`.
    
    Attributes
    ----------
    state_init : `qutip.Qobj`
        Initial state vector or density matrix.
    pulse_seq : :obj:`pulse_sequence.Pulse_sequence`
        Sequence of pulses to be applied to initial state.
    tlist : list of float
        List of times at which to evaluate the system's state or the expectation value of the operator(s) in `e_ops`. 
    Hstat : `qutip.Qobj`
        Time-independent Hamiltonian that applies at all times in the simulation.
    c_ops : list of `qutip.Qobj`
        List of collapse operators.
    e_ops : list of `qutip.Qobj`
        List of operators for which to evaluate expectation values.
    opts : `qutip.Options`
        Options for `mesolve`.
    """
    
    def __init__(self, state_init = None, tlist = None, pulse_seq = None, c_ops=None, e_ops=None, opts = None):
        
        self._state_init = state_init
        self._tlist = tlist
        self.pulse_seq = pulse_seq
        self.Hstat = self._pulse_seq.Hstat
        self._c_ops = c_ops if c_ops is not None else []
        self._e_ops = e_ops if e_ops is not None else []
        self._opts = opts
        
    def generate_H(self, pulse_seq):
        """
        Put Hamiltonian terms (static Hamiltonian and pulses) in a list together. The static Hamiltonian is first. Time-dependent (pulse) terms are [operator, callback function] pairs. The resulting list is in the correct format to be passed directly into mesolve as `H`.
        
        Returns
        -------
        H : list of `qutip.Qobj` or list of [`qutip.Qobj`, function]
            List of `qutip.Qobj` objects and [`qutip.Qobj`, callback function] pairs.
            List of Hamiltonian terms and their coefficients at each point in time.
        """
        H = []
        is_Hstat = False
        
        
        if pulse_seq.Hstat != None:
            H.append(pulse_seq.Hstat)
            is_Hstat = True
        for pulse in pulse_seq.pulses:
            H.append([pulse.ham.H, pulse.eval_coeff])
            
        if len(pulse_seq.pulses) == 0 and is_Hstat:
            H = H[0]
        
        return H
    
    def serial_generate_H(self, pulselist, safe = False):
        H = []

        if safe:
            for pulse in pulselist:
                H.append([pulse.ham.H, pulse.eval_coeff])
        else:
            for pulse in pulselist:
                H.append([pulse.ham.H, pulse.serial_eval_coeff])

        return H
    
    def evolve(self, state_init = None, tlist = None, pulse_seq = None, c_ops = None, e_ops = None, opts = None, method = 'serial', t_rtol = 1e-8):
        """
        Evolve `state_init` using given pulse sequence. Collapse operators `c_ops` may be passed. Return a `qutip.Result` object which stores output for each timestamp given in `tlist`. If `e_ops` is left as None, `qutip.Result` contains the state vector of the system at each timestamp in its `states` attribute. Otherwise, the expectation values of the operator(s) listed in `e_ops` are stored in `qutip.Result.expect` in a 2-dimensional list. 
        
        The parameter `evolve_method` may take one of three values: 'regular', 'serial', or 'serial_safe'. These describe different ways of solving the evolution of the state, but ultimately all use `qutip.mesolve` at their base. The two serial options significantly reduce runtime for evolving systems with many pulses, especially if the total number of pulses is much larger than the number of overlapping pulses at any given instant. 'serial' is the fastest of the three, but 'serial_safe' is slightly more true to the plain `qutip.mesolve` (denoted by 'regular') in terms of state accuracy. One can additionally control the relative error tolerance in integration by changing the parameter `t_rtol`, which applies to 'serial' and 'serial_safe'. 
        
        Parameters
        ----------
        state_init : `qutip.Qobj`
            Initial state vector or density matrix.
        tlist : list of float
            List of times at which to evaluate the system's state or the expectation value of the operator(s) in `e_ops`. 
        pulse_seq : `pulse_sequence.Pulse_sequence` 
            Pulse sequence object containing the pulses for this simulation.
        c_ops : list of `qutip.Qobj`
            List of collapse operators.
        e_ops : list of `qutip.Qobj`
            List of operators for which to evaluate expectation values.
        opts : `qutip.Options`
            Options for `mesolve`.
        method : str {'regular', 'serial', 'serial_safe'}, default = 'serial'
            Choice of which solving method to use.
        t_rtol : float, default = 1e-8
            Relative time-step tolerance. Determines which relative time differences can be resolved in `serial` and `serial_safe` evolve methods.
            
        Returns
        -------
        result : `qutip.Result`
            System state or expectation value at `tlist` timestamps.
        """
        if state_init == None:
            if self._state_init == None:
                raise AttributeError("Initial state has not been defined")
            else:
                state_init = self._state_init
        if tlist is None:
            if self._tlist is None:
                raise AttributeError("tlist has not been defined")
            else:
                tlist = self._tlist
        if pulse_seq == None:
            if self._pulse_seq == None:
                raise AttributeError('Pulse sequence has not been defined')
            else:
                pulse_seq = self._pulse_seq
        if c_ops == None:
            c_ops = self._c_ops
        if e_ops == None:
            e_ops = self._e_ops
        if opts == None:
            opts = self._opts
                
        if method == 'regular':
            H = self.generate_H(pulse_seq)
            result = mesolve(H, state_init, tlist, c_ops, e_ops, args = None, options = opts)
        elif method == 'serial':
            result = self.serial_evolve(pulse_seq, state_init, tlist, c_ops, e_ops, opts, t_rtol)
        elif method == 'serial_safe':
            result = self.serial_safe_evolve(pulse_seq, state_init, tlist, c_ops, e_ops, opts, t_rtol)
        else: 
            raise ValueError('Solving method must be one of the following: \'regular\', \'serial\', or \'serial_safe\' ')
            
        return result
            
    def serial_evolve(self, pulse_seq, state_init, tlist, c_ops, e_ops, opts, t_rtol):

#         timediff = np.diff(tlist)
#         if np.min(timediff) < t_tol:
#             raise ValueError("tlist values too close together; try decreasing t_rel_tol or changing tlist")

        for i in np.arange(len(tlist)-1):
            if isclose(tlist[i], tlist[i+1], rel_tol=t_rtol):
                raise ValueError("tlist values too close together; try decreasing t_rtol or changing tlist")
        
        ps_init = pulse_seq
    
        ps_starts = sorted(ps_init.pulses, key = lambda p: p.start_time)
        H_starts = self.serial_generate_H(ps_starts)
        i = 0

        
        ps_ends = sorted(ps_init.pulses, key = lambda p: p.end_time)
        H_ends = self.serial_generate_H(ps_ends)
        j = 0

        to_user_states = [state_init]
        #to_user_expect = [[expect(op, state_init)] for op in e_ops]
        state_init = state_init

        # sort pulse start and end times in order to create intervals; within each interval, the pulses that fall in it are ON for the entire interval's duration
        segment_breaks_unsorted = set([p.start_time for p in ps_starts] + [p.end_time for p in ps_ends] + [tlist[0]] + [tlist[-1]])
        segment_breaks = sorted(list(segment_breaks_unsorted))
        
        l = 0
        while l < len(segment_breaks):
            if segment_breaks[l] < tlist[0]:
                segment_breaks.remove(segment_breaks[l])
                l += 1
            else: 
                break

        segment_H = []
        
        is_Hstat = True
        if pulse_seq.Hstat == None:
            is_Hstat = False
        else:
            segment_H = [pulse_seq.Hstat]

        num_pulses = len(ps_init.pulses)

        k = 0


        # loop through segments
        for n in np.arange(len(segment_breaks)-1):
            segment_start = segment_breaks[n]
            if segment_start < tlist[0]:
                continue
            if segment_start >= tlist[-1]:
                break
            segment_end = segment_breaks[n+1]

            segment_tlist = [segment_start]

            while i < num_pulses and ps_starts[i].start_time < segment_end:
                segment_H.append(H_starts[i])
                i += 1

            while j < num_pulses and ps_ends[j].end_time <= segment_start:
                segment_H.remove(H_ends[j])
                j += 1

            k += 1
                
            while k < len(tlist) and not isclose(tlist[k], segment_end, rel_tol = t_rtol) and tlist[k] < segment_end: # as long as tlist value is not too close to segment end...
                segment_tlist.append(tlist[k])
                k += 1
                
            
            segment_tlist.append(segment_end)
            
            if len(segment_H) == 0 or isclose(segment_start, segment_end, rel_tol= t_rtol):
                timestamps = len(segment_tlist)
                states = [state_init]*timestamps
                #exp = np.repeat(expect(e_ops, state_init), timestamps, axis = 1)
                
            else:
                if len(segment_H) == 1 and is_Hstat:
                    res = mesolve(segment_H[0], state_init, segment_tlist, c_ops = c_ops, options = opts)
                else:
                    res = mesolve(segment_H, state_init, segment_tlist, c_ops = c_ops, options = opts)
                states = res.states
                #exp = expect(e_ops, states)
                state_init = states[-1] # for the next iteration, this is the initial state

            if k < len(tlist) and not isclose(tlist[k], segment_end, rel_tol = t_rtol): # if we're close but above, ignore; if we're close but below, ignore; and if we're not close, delete the last state. Not close implies that tlist[k] is necessarily above segment_end, so we delete states[-1] which is the value at segment_end, something the user is not interested in. 
                # cannot be both close and above + close and below in the same segment
                states = states[:-1]
                #exp = np.delete(exp, -1, 1)
                k -= 1
            elif k < len(tlist)-1 and isclose(tlist[k+1], segment_end, rel_tol = t_rtol): 
                states.append(states[-1])
                #exp = np.hstack((exp, exp[:, -1][:, np.newaxis]))
                k += 1
            
            states = states[1:] # delete initial state every time; if it is in tlist, it will be reported by the previous iteration
            #exp = np.delete(exp, 0, 1)

            to_user_states += states
            
            is_e_ops = False
            if isinstance(e_ops, qobj.Qobj):
                to_user_expect = [expect(e_ops, to_user_states)]
                is_e_ops = True
            elif callable(e_ops):
                to_user_expect = [] 
                for t_ind in np.arange(len(tlist)): 
                    to_user_expect.append(e_ops(t_ind, to_user_states[t_ind])) 
                    
                is_e_ops = True
            elif e_ops == None or len(e_ops) == 0:
                to_user_expect = []
            else:
                to_user_expect = []
                for op_ind in np.arange(len(e_ops)):
                    op = e_ops[op_ind]
                    if not isinstance(op, qobj.Qobj) and callable(op):
                        to_user_expect.append(np.zeros(len(tlist), dtype = complex))
                        for t_ind in np.arange(len(tlist)):
                            to_user_expect[op_ind][t_ind] = op(t_ind, to_user_states[t_ind])
                    else:
                        to_user_expect.append(expect(op, to_user_states))
                is_e_ops = True
                
             
            if is_e_ops:
                to_user_states = []
                    
            
        to_user = solver.Result()
        to_user.solver = res.solver
        to_user.times = tlist
        to_user.states = to_user_states
        to_user.expect = to_user_expect
        to_user.num_expect = res.num_expect
        to_user.num_collapse = res.num_collapse
        to_user.ntraj = res.ntraj
        to_user.col_times = res.col_times
        to_user.col_which = res.col_which
        
        return to_user

    def serial_safe_evolve(self, pulse_seq, state_init, tlist, c_ops = None, e_ops = None, opts = None, t_rtol = 1e-8):
        
        for i in np.arange(len(tlist)-1):
            if isclose(tlist[i], tlist[i+1], rel_tol=t_rtol):
                raise ValueError("tlist values too close together; try decreasing t_rtol or changing tlist")
        
        ps_init = pulse_seq
    
        ps_starts = sorted(ps_init.pulses, key = lambda p: p.start_time)

        ps_ends = sorted(ps_init.pulses, key = lambda p: p.end_time)
        
        segment_breaks_unsorted = set([p.start_time for p in ps_starts] + [p.end_time for p in ps_ends] + [tlist[0]] + [tlist[-1]])
        segment_breaks = sorted(list(segment_breaks_unsorted))
        
        l = 0
        while l < len(segment_breaks):
            if segment_breaks[l] < tlist[0]:
                segment_breaks.remove(segment_breaks[l])
                l += 1
            else: 
                break
        
        n = 0
        k = 0 
        while n < len(segment_breaks)-1:
            increment = True
            
            if isclose(segment_breaks[n], segment_breaks[n+1], rel_tol = t_rtol):
                if tlist[0] == segment_breaks[n+1] or tlist[-1] == segment_breaks[n+1]:
                    segment_breaks.remove(segment_breaks[n])
                elif n+1 < len(segment_breaks)-1:
                    segment_breaks.remove(segment_breaks[n+1])
                    increment = False
                
            while k < len(tlist) and tlist[k] < segment_breaks[n+1]:
                if isclose(tlist[k], segment_breaks[n], rel_tol = t_rtol) and tlist[k] != segment_breaks[n]:
                    segment_breaks.remove(segment_breaks[n])
                    increment = False
                    break
                else:
                    k += 1
                    
            if k >= len(tlist):
                break
            if k > 0:
                k -= 1
            if increment:
                n += 1
                
        H_starts = self.serial_generate_H(ps_starts, safe = True)
        i = 0
        H_ends = self.serial_generate_H(ps_ends, safe = True)
        j = 0
        
        to_user_states = [state_init]
        #to_user_expect = [[expect(op, state_init)] for op in e_ops]
        state_init = state_init
        
        segment_H = []
        
        is_Hstat = True
        if pulse_seq.Hstat == None:
            is_Hstat = False
        else:
            segment_H = [pulse_seq.Hstat]

        num_pulses = len(ps_init.pulses)

        k = 0

        # loop through segments
        for n in np.arange(len(segment_breaks)-1):
            segment_start = segment_breaks[n]
            if segment_start < tlist[0]:
                continue
            elif segment_start >= tlist[-1]:
                break
            segment_end = segment_breaks[n+1]
            segment_tlist = [segment_start]

            while i < num_pulses and ps_starts[i].start_time < segment_end:
                segment_H.append(H_starts[i])
                i += 1

            while j < num_pulses and ps_ends[j].end_time <= segment_start:
                segment_H.remove(H_ends[j])
                j += 1

            k += 1
                
            while k < len(tlist) and tlist[k] < segment_end:
                segment_tlist.append(tlist[k])
                k += 1
                
              
            
            segment_tlist.append(segment_end)
            
            if len(segment_H) == 0:
                timestamps = len(segment_tlist)
                states = [state_init]*timestamps
                #exp = [[expect(op, state_init)]*timestamps for op in e_ops]
            
            else:
                if len(segment_H) == 1 and is_Hstat:
                    res = mesolve(segment_H[0], state_init, segment_tlist, c_ops = c_ops, options = opts)
                    
                else:
                    res = mesolve(segment_H, state_init, segment_tlist, c_ops = c_ops, options = opts)
                states = res.states
                #exp = [[expect(op, state) for op in e_ops] for state in states]
                state_init = states[-1] # for the next iteration, this is the initial state

            if k < len(tlist) and tlist[k] != segment_end: 
                states = states[:-1]
                #exp = [row[:-1] for row in exp]
                k -= 1
            
            states = states[1:] # delete initial state every time; if it is in tlist, it will be reported by the previous iteration
            #exp = [row[1:] for row in exp]

            to_user_states += states
            
            is_e_ops = False
            if isinstance(e_ops, qobj.Qobj):
                to_user_expect = [expect(e_ops, to_user_states)]
                is_e_ops = True
            elif callable(e_ops):
                to_user_expect = [] 
                for t_ind in np.arange(len(tlist)): 
                    to_user_expect.append(e_ops(t_ind, to_user_states[t_ind])) 
                    
                is_e_ops = True
            elif e_ops == None or len(e_ops) == 0:
                to_user_expect = []
            else:
                to_user_expect = []
                for op_ind in np.arange(len(e_ops)):
                    op = e_ops[op_ind]
                    if not isinstance(op, qobj.Qobj) and callable(op):
                        to_user_expect.append(np.zeros(len(tlist), dtype = complex))
                        for t_ind in np.arange(len(tlist)):
                            to_user_expect[op_ind][t_ind] = op(t_ind, to_user_states[t_ind])
                    else:
                        to_user_expect.append(expect(op, to_user_states))
                is_e_ops = True
                
             
            if is_e_ops:
                to_user_states = []
            
        to_user = solver.Result()
        to_user.solver = res.solver
        to_user.times = tlist
        to_user.states = to_user_states
        to_user.expect = to_user_expect
        to_user.num_expect = res.num_expect
        to_user.num_collapse = res.num_collapse
        to_user.ntraj = res.ntraj
        to_user.col_times = res.col_times
        to_user.col_which = res.col_which
        
        return to_user
    
    @property
    def state_init(self):
        return self._state_init
    
    @state_init.setter
    def state_init(self, state_init):
        self._state_init = state_init
        
    @state_init.deleter
    def state_init(self):
        del self._state_init
        
    @property
    def tlist(self):
        return self._tlist
    
    @tlist.setter
    def tlist(self, tlist):
        self._tlist = tlist 
        
    @tlist.deleter
    def tlist(self):
        del self._tlist
    
    @property
    def pulse_seq(self):
        return self._pulse_seq
    
    @pulse_seq.setter
    def pulse_seq(self, pulse_seq):
        if isinstance(pulse_seq, Pulse_sequence):
            self._pulse_seq = pulse_seq
        else:
            raise TypeError('Pulse sequence must be a Pulse_sequence object')
        
    @pulse_seq.deleter
    def pulse_seq(self):
        del self._pulse_seq
    
    @property
    def Hstat(self):
        return self._Hstat
    
    @Hstat.setter
    def Hstat(self, Hstat):
        if isinstance(Hstat, Hamiltonian):
            self._Hstat = Hstat.H
        elif isinstance(Hstat, qobj.Qobj):
            self._Hstat = Hstat
        else: 
            print(type(Hstat))
            raise TypeError("Static Hamiltonian operator must be a Hamiltonian object or a Quantum object")
    
    @Hstat.deleter
    def Hstat(self):
        del self._Hstat
    
    @property
    def c_ops(self):
        return self._c_ops
    
    @c_ops.setter
    def c_ops(self, c_ops):
        self._c_ops = c_ops
        
    @c_ops.deleter
    def c_ops(self):
        del self._c_ops
    
    @property
    def e_ops(self):
        return self._e_ops
    
    @e_ops.setter
    def e_ops(self, e_ops):
        self._e_ops = e_ops
    
    @e_ops.deleter
    def e_ops(self):
        del self._e_ops
   
    @property
    def opts(self):
        return self._opts

    @opts.setter
    def opts(self, opts):
        self._opts = opts
    
    @opts.deleter
    def opts(self):
        del self._opts
            
