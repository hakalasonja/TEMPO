#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 24 20:10:22 2023

Class to implement time-evolution of a state under applied pulses.

@author: hakalas
"""

from qutip import mesolve, solver, qobj, expect
import os
from tempo.ham import Hamiltonian
from tempo.pulse import Pulse
import numpy as np
from math import isclose

# incorporate static ham input // DONE
# checks for evolver static ham arg, if not then checks for it in pulseseq, if not then no stat ham
# getters setters // DONE
# put opts in notebook, not here as defaults, since different systems can need wildly different things // DONE

#TODO:
# what to put for qutip parameter documentation like c_ops, result
# make generatehamls consistent in the sense that hams should be in input options as well, not just pulling instance variables
# Hamls as instance variable or no? It's not used in serial case --> remove Hamls option bceause the user could just use qutip mesolve on their own if they want it // DONE
# serial_evolve as separate method from evolve? OR evolve with serial=True option? In the latter case, what to do about the fact that regular evolve can take Hamls input, but serial cannot use Hamls? And what if user passes in both to regular (which one takes priority)? If pulseseq takes priority, make sure to add pulseseq input to generateHamls as it's called in evolve() --> hamls not a problem anymore (see above), but make sure to incorporate user input pulsesequence in serial case // DONE
# update docstrings  
# what to do about static hamiltonian input in evolve() and serial? currently not an input to either, so the only way to get it to evolver is through pulsesequence or Hamls input or instance var pulseseq, and the only way to get it to serial is through pulsesequence input --> make tlist, in_state, pulseseq optional in constructor, //DONE// then require at least one input of tlist and in_state in evolve() call //DONE//. evolve() inputs override instance variables. //DONE// put inputs in in_state, tlist, ps order. //DONE// remove statham input from constructor, //DONE// pull statham from pulsesequence //DONE
# how to add expectation values to Result object? Unclear what to do in the case where e_ops are callback functions, or what to do if user inputs none (don't want it to save [[][][][]] in Result.expect) --> why is this a problem? look into mesolve evolve function to see how it looks different in opr vs callback func case
# pick a lane: interval or segment?
# what do about e_ops: if user inputs e_ops, no states are generated for result object. We need at least the final state of each interval to mesolve the next interval. Should we just always do mesolve with no expectation operators and calculate the exp. values ourselves after obtaining the state? or should we do mesolve with the oprs and then an extra mesolve to obtain the state at the end of the interval? 
# hwo to deal with tiny intervals like [1, 1.00000000002] that come from rounding error when creating tlist or pulse timings? --> add tolerance t_tol, default 1e-8; numpy should have a function to compare two values in a fractional difference. 
# simplify expect instance variable if e_ops = None; don't want it to save [[],[],[]] but whatever mesolve would normally do

class Evolver():
    """
    Time-evolution of a state under a pulse sequence.
    
    The Evolver object stores information about a particular pulse sequence simulation. Its `evolve()` method evolves the given initial state `state_init` when the pulses in `pulsesequence` are applied. `evolve()` relies on QuTiP's `mesolve` method and the inputs `c_ops`, `e_ops`, and `opts`, if defined, are directly passed to `mesolve`. For these, see `mesolve` documentation https://qutip.org/docs/latest/apidoc/functions.html#module-qutip.mesolve. 
    
    If `e_ops` is left as None, `evolve()` returns the state vector of the system at each timestamp defined in `tlist`. Otherwise, the expectation value of the operator(s) listed in `e_ops` is returned.
    
    If a static (time-independent) Hamiltonian applies at all times in `tlist`, the Hamiltonian may be passed to the constructor as `statham`. 
    
    Parameters 
    ----------
    state_init : `qutip.Qobj`
        Initial state vector or density matrix.
    pulsesequence : :obj:`pulsesequence.Pulsesequence`
        Sequence of pulses to be applied to initial state.
    tlist : list of float
        List of times at which to evaluate the system's state or the expectation value of the operator(s) in `e_ops`. 
    statham : :obj:`ham.Hamiltonian` or `qutip.Qobj`, optional
        Time-independent Hamiltonian that applies at all times in the simulation.
    c_ops : list of `qutip.Qobj`
        List of collapse operators.
    e_ops : list of `qutip.Qobj`
        List of operators for which to evaluate expectation values.
    opts : `qutip.Options`
        Options for `mesolve`.
    
    Attributes
    ----------
    state_init : `qutip.Qobj`
        Initial state vector or density matrix.
    pulsesequence : :obj:`pulsesequence.Pulsesequence`
        Sequence of pulses to be applied to initial state.
    tlist : list of float
        List of times at which to evaluate the system's state or the expectation value of the operator(s) in `e_ops`. 
    statham : `qutip.Qobj`
        Time-independent Hamiltonian that applies at all times in the simulation.
    c_ops : list of `qutip.Qobj`
        List of collapse operators.
    e_ops : list of `qutip.Qobj`
        List of operators for which to evaluate expectation values.
    opts : `qutip.Options`
        Options for `mesolve`.
    """
    
    def __init__(self, state_init = None, tlist = None, pulsesequence = None, c_ops=None, e_ops=None, opts = None):
        
        self._stateinit = state_init
        self._tlist = tlist
        self._pulsesequence = pulsesequence
        self._statham = self._pulsesequence.getstatham()
        self._c_ops = c_ops if c_ops is not None else []
        self._e_ops = e_ops if e_ops is not None else []
        self._opts = opts
        
    def generateHamls(self, pulsesequence = None):
        """
        Put Hamiltonian terms (static Hamiltonian and pulses) in a list together. The static Hamiltonian is first. Time-dependent (pulse) terms are [operator, callback function] pairs. `pulsesequence` and `statham` attributes are used to obtain terms.
        
        Returns
        ----------
        Hls : list of `qutip.Qobj` or list of [`qutip.Qobj`, function]
            List of `qutip.Qobj` objects and [`qutip.Qobj`, callback function] pairs.
            List of Hamiltonian terms and their coefficients at each point in time.
        """
        Hls = []
        
        if pulsesequence == None:
            if self._statham != None:
                Hls.append(self._statham)
            elif self._pulsesequence.getstatham() != None:
                Hls.append(self._pulsesequence.getstatham())

            for pulse in self._pulsesequence.getpulsels():
                Hls.append([pulse.getpulsemat(), pulse.evalcoeff])
        else:
            if pulsesequence.getstatham() != None:
                Hls.append(pulsesequence.getstatham())
            for pulse in pulsesequence.getpulsels():
                Hls.append([pulse.getpulsemat(), pulse.evalcoeff])
                     
        return Hls
    
    def serial_generateHamls(self, pulsels):
        Hls = []

        for pulse in pulsels:
            Hls.append([pulse.getpulsemat(), pulse.serial_evalcoeff])

        return Hls
   
    def serial_generateHamls_alt(self, pulsels):
        Hls = []

        for pulse in pulsels:
            Hls.append([pulse.getpulsemat(), pulse.evalcoeff])

        return Hls
    
    def evolve(self, state_init = None, tlist = None, pulsesequence = None, c_ops = None, e_ops = None, opts = None, serial = True, t_rel_tol = 1e-8):
        """
        Evolve `state_init` using given Hamiltonians `Hamls` and collapse operators `c_ops`. Return output at each time given in `tlist`. Output may be state vector or expectation value depending on if `e_ops` was specified.
        
        Parameters
        ----------
        state_init : `qutip.Qobj`
            Initial state vector or density matrix.
        Hamls : list of `qutip.Qobj` or list of [`qutip.Qobj`, function] 
            List of Hamiltonian terms, both time-independent and time-dependent. Time-dependent terms are [operator, callback function] pairs. 
        tlist : list of float
            List of times at which to evaluate the system's state or the expectation value of the operator(s) in `e_ops`. 
        c_ops : list of `qutip.Qobj`
            List of collapse operators.
        e_ops : list of `qutip.Qobj`
            List of operators for which to evaluate expectation values.
        opts : `qutip.Options`
            Options for `mesolve`.
        
        Returns
        ----------
        result : `qutip.Result`
            System state or expectation value at 'tlist' timestamps.
        """
        
        if state_init == None:
            if self._stateinit == None:
                raise AttributeError("Initial state has not been defined")
            else:
                state_init = self._stateinit
        if tlist is None:
            if self._tlist is None:
                raise AttributeError("tlist has not been defined")
            else:
                tlist = self._tlist
        if pulsesequence == None:
            pulsesequence = self._pulsesequence
        if c_ops == None:
            c_ops = self._c_ops
        if e_ops == None:
            e_ops = self._e_ops
        if opts == None:
            opts = self._opts
                
        if serial == False:
            if pulsesequence != None:
                Hamls = self.generateHamls(pulsesequence)
            else:
                Hamls = self.generateHamls()
            result = mesolve(Hamls, state_init, tlist, c_ops, e_ops, args = None, options = opts)
        else: 
            result = self.serial_evolve(pulsesequence, state_init, tlist, c_ops, e_ops, opts, t_rel_tol)
        
        return result
            
    def serial_evolve(self, pulseseq, state_init, tlist, c_ops, e_ops, opts, t_rel_tol):
        
#         timediff = np.diff(tlist)
#         if np.min(timediff) < t_tol:
#             raise ValueError("tlist values too close together; try decreasing t_rel_tol or changing tlist")

        for i in np.arange(len(tlist)-1):
            if isclose(tlist[i], tlist[i+1], rel_tol=t_rel_tol):
                raise ValueError("tlist values too close together; try decreasing t_rel_tol or changing tlist")
        
        ps_init = pulseseq
    
        ps_starts = sorted(ps_init.getpulsels(), key = lambda p: p.getstarttime())
        Ham_starts = self.serial_generateHamls(ps_starts)
        i = 0

        
        ps_ends = sorted(ps_init.getpulsels(), key = lambda p: p.getendtime())
        Ham_ends = self.serial_generateHamls(ps_ends)
        j = 0

        to_user_states = [state_init]
        to_user_expect = [[expect(op, state_init)] for op in e_ops]
        state_init = state_init

        # only include tlist[0] if earliest pulse starts after tlist[0] 
        # sort pulse start and end times in order to create intervals; within each interval, the pulses that fall in it are ON for the entire interval's duration
        times = sorted(list(set([p.getstarttime() for p in ps_starts] + [p.getendtime() for p in ps_ends] + [tlist[0]] + [tlist[-1]])))
        # remove times before tlist[0]
        # split this up into multiple lines

        Hamlist_final = []
        is_statham = True
        if pulseseq.getstatham() == None:
            is_statham = False
        else:
            Hamlist_final = [pulseseq.getstatham()]

        len_hams = len(ps_init.getpulsels())

        k = 0

        badcount = 0

        # loop through segments
        for n in np.arange(len(times)-1):
            segment_start = times[n]
            if segment_start < tlist[0]:
                continue
            if segment_start >= tlist[-1]:
                break
            segment_end = times[n+1]

            segment_tlist = [segment_start]

            while i < len_hams and ps_starts[i].getstarttime() < segment_end:
                Hamlist_final.append(Ham_starts[i])
                i += 1

            while j < len_hams and ps_ends[j].getendtime() <= segment_start:
                Hamlist_final.remove(Ham_ends[j])
                j += 1

#             while k < len(tlist) and isclose(tlist[k], segment_start, rel_tol = t_rel_tol):
            k += 1
                
            while k < len(tlist) and not isclose(tlist[k], segment_end, rel_tol = t_rel_tol) and tlist[k] < segment_end: # as long as tlist value is not too close to segment end...
                segment_tlist.append(tlist[k])
                k += 1
                
            
            segment_tlist.append(segment_end)
#             print(segment_tlist)
            
            if len(Hamlist_final) == 0 or (value := isclose(segment_start, segment_end, rel_tol= t_rel_tol)):
                # remove this error for now
#                 if value:
#                     badcount += 1
#                     if badcount == 3:
#                         raise ValueError("Multiple pulse timings too close together; try decreasing t_rel_tol or changing starttimes or durations")
#                 else:
#                     badcount = 0
                timestamps = len(segment_tlist)
                states = [state_init]*timestamps
                exp = [[expect(op, state_init)]*timestamps for op in e_ops]
                
            elif len(Hamlist_final) == 1:
#                 badcount = 0
                H = Hamlist_final[0] if is_statham else Hamlist_final # can put this in elif condition
                res = mesolve(H, state_init, segment_tlist, c_ops, e_ops, args = None, options = opts)
                states = res.states
                exp = res.expect
                state_init = states[-1] # for the next iteration, this is the initial state
            
            else:
#                 badcount = 0
                res = mesolve(Hamlist_final, state_init, segment_tlist, c_ops, e_ops, args = None, options = opts)
                states = res.states # avoid duplicate code with above by making a boolean tomesolve 
                exp = res.expect
                state_init = states[-1] # for the next iteration, this is the initial state

            if k < len(tlist) and not isclose(tlist[k], segment_end, rel_tol = t_rel_tol): # if we're close but above, ignore; if we're close but below, ignore; and if we're not close, delete the last state. Not close implies that tlist[k] is necessarily above segment_end, so we delete states[-1] which is the value at segment_end, something the user is not interested in. 
                # cannot be both close and above + close and below in the same segment
                states = states[:-1]
                exp = [row[:-1] for row in exp]
                k -= 1
            elif k < len(tlist)-1 and isclose(tlist[k+1], segment_end, rel_tol = t_rel_tol): 
                states.append(states[-1])
                exp = [row + row[-1] for row in exp]
                k += 1
            
            states = states[1:] # delete initial state every time; if it is in tlist, it will be reported by the previous iteration
            exp = [row[1:] for row in exp]

            to_user_states += states
            to_user_expect = [to_user_expect[i]+exp[i] for i in np.arange(len(exp))]
            
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

    def serial_evolve_alt(self, pulseseq, state_init, tlist, c_ops = None, e_ops = None, opts = None, t_rel_tol = 1e-8):
        
        for i in np.arange(len(tlist)-1):
            if isclose(tlist[i], tlist[i+1], rel_tol=t_rel_tol):
                raise ValueError("tlist values too close together; try decreasing t_rel_tol or changing tlist")
        
        ps_init = pulseseq
    
        ps_starts = sorted(ps_init.getpulsels(), key = lambda p: p.getstarttime())

        ps_ends = sorted(ps_init.getpulsels(), key = lambda p: p.getendtime())
        
        times = sorted(list(set([p.getstarttime() for p in ps_starts] + [p.getendtime() for p in ps_ends] + [tlist[0]] + [tlist[-1]])))
        #print('times: ', times)
        
        n = 0
        k = 0 
        while n < len(times)-1:
            increment = True
            
            if isclose(times[n], times[n+1], rel_tol = t_rel_tol):
                if tlist[0] == times[n+1] or tlist[-1] == times[n+1]:
                    #print('remove 1: ', times[n])
                    times.remove(times[n])
                elif n+1 < len(times)-1:
                    #print('remove 2: ', times[n+1])
                    times.remove(times[n+1])
                    increment = False
                
            while k < len(tlist) and tlist[k] < times[n+1]:
                if isclose(tlist[k], times[n], rel_tol = t_rel_tol) and tlist[k] != times[n]:
                    #print('remove 3: ', times[n])
                    times.remove(times[n])
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
            
        #print('final: ', times)
                
        Ham_starts = self.serial_generateHamls_alt(ps_starts)
        i = 0
        Ham_ends = self.serial_generateHamls_alt(ps_ends)
        j = 0
        
        to_user_states = [state_init]
        to_user_expect = [[expect(op, state_init)] for op in e_ops]
        state_init = state_init
        
        Hamlist_final = []
        is_statham = True
        if pulseseq.getstatham() == None:
            is_statham = False
        else:
            Hamlist_final = [pulseseq.getstatham()]

        len_hams = len(ps_init.getpulsels())

        k = 0

        # loop through segments
        for n in np.arange(len(times)-1):
            segment_start = times[n]
            if segment_start < tlist[0]:
                continue
            elif segment_start >= tlist[-1]:
                break
            segment_end = times[n+1]
            segment_tlist = [segment_start]

            while i < len_hams and ps_starts[i].getstarttime() < segment_end:
                Hamlist_final.append(Ham_starts[i])
                i += 1

            while j < len_hams and ps_ends[j].getendtime() <= segment_start:
                Hamlist_final.remove(Ham_ends[j])
                j += 1

            k += 1
                
            while k < len(tlist) and tlist[k] < segment_end:
                segment_tlist.append(tlist[k])
                k += 1
                
            segment_tlist.append(segment_end)
            
            if len(Hamlist_final) == 0:
                timestamps = len(segment_tlist)
                states = [state_init]*timestamps
                exp = [[expect(op, state_init)]*timestamps for op in e_ops]
                
            elif len(Hamlist_final) == 1:
                H = Hamlist_final[0] if is_statham else Hamlist_final 
                res = mesolve(H, state_init, segment_tlist, c_ops, e_ops, args = None, options = opts)
                states = res.states
                exp = res.expect
                state_init = states[-1] # for the next iteration, this is the initial state
            
            else:
                res = mesolve(Hamlist_final, state_init, segment_tlist, c_ops, e_ops, args = None, options = opts)
                states = res.states
                exp = res.expect
                state_init = states[-1] # for the next iteration, this is the initial state

            if k < len(tlist) and tlist[k] != segment_end: 
                states = states[:-1]
                exp = [row[:-1] for row in exp]
                k -= 1
            
            states = states[1:] # delete initial state every time; if it is in tlist, it will be reported by the previous iteration
            exp = [row[1:] for row in exp]

            to_user_states += states
            to_user_expect = [to_user_expect[i]+exp[i] for i in np.arange(len(exp))]
            
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
    
    def getstateinit(self):
        return self._stateinit
    
    def getpulsesequence(self):
        return self._pulsesequence
    
    def gettlist(self):
        return self._tlist
    
    def getstatham(self):
        return self._statham
    
    def getc_ops(self):
        return self._c_ops
    
    def gete_ops(self):
        return self._e_ops
    
    def getopts(self):
        return self._opts
    
    def setstateinit(self, state_init):
        self._stateinit = state_init    
        
    def setpulsesequence(self, pulsesequence):
        self._pulsesequence = pulsesequence
        
    def settlist(self, tlist):
        self._tlist = tlist   
        
    def setstatham(self, statham):
        
        if type(statham) == ham.Hamiltonian:
            self._statham = statham.getH()
        elif type(statham) == qobj.Qobj:
            self._statham = statham
        else: 
            print(type(statham))
            raise TypeError("Static Hamiltonian operator must be a Hamiltonian object or a Quantum object")
    
    def setc_ops(self, c_ops):
        self._c_ops = c_ops
    
    def sete_ops(self, e_ops):
        self._e_ops = e_ops
        
    def setopts(self, opts):
        self._opts = opts
    
#     def setHamls(self, Hamls):
#         if type(Hamls) in [list, np.ndarray]:
#             self._Hamls = Hamls
#         else:
#             raise TypeError("Hamiltonians must be in a list")
    
    def delstateinit(self):
        del self._stateinit
    
    def delpulsesequence(self):
        del self._pulsesequence

    def deltlist(self):
        del self._tlist
        
    def delstatham(self):
        del self._statham
        
    def delc_ops(self):
        del self._c_ops
        
    def dele_ops(self):
        del self._e_ops
        
    def delopts(self):
        del self._opts
            
#     def delHamls(self):
#         del self._Hamls
    
    stateinit = property(getstateinit, setstateinit, delstateinit)
    pulsesequence = property(getpulsesequence, setpulsesequence, delpulsesequence)
    tlist = property(gettlist, settlist, deltlist)
    statham = property(getstatham, setstatham, delstatham)
    c_ops = property(getc_ops, setc_ops, delc_ops)
    e_ops = property(gete_ops, sete_ops, dele_ops)
    opts = property(getopts, setopts, delopts)
#     Hamls = property(getHamls, setHamls, delHamls)
    
        
