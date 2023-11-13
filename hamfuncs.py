#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 10 19:52:34 2023

@author: hakalas
"""

from qutip import *
import os
import qsys
import ham
import numpy as np

def ZFSfunc(Hmats, Hpars):
    return Hpars['coeff']*Hpars['ZFSconst']*Hmats['Sz']*Hmats['Sz']

def ZeeNVfunc(Hmats, Hpars):
    return Hpars['coeff']*Hpars['gammaNV']*tensor(dotproduct(Hpars['Bfield'], Hmats['S1']), identity(Hpars['P2dims']))

def ZeeNucfunc(Hmats, Hpars):
    return Hpars['coeff']*Hpars['gammaNuc']*tensor(identity(Hpars['P1dims']), dotproduct(Hpars['Bfield'], Hmats['S2']))

def HFfunc(Hmats, Hpars):
    return Hmats['A'][2,2] * Hmats['Sz'] * Hmats['Iz'] + Hmats['A'][0,0] * (Hmats['Sx']*Hmats['Ix'] + Hmats['Sy']*Hmats['Iy'])

def ZFS(qsystem):
    
    # assumes that ZFS is on particle 1 in the qsystem tuple "dimensions"
    
    ZFSpars = {'coeff': 2*np.pi, 'ZFSconst': 2.87e3}
    ZFSmats = {'Sz': qsystem.getSz()[0]}
    
    return ham.Hamiltonian(ZFSmats, ZFSpars, ZFSfunc)

def Zeeman(qsystem, Bfield, nuc = False, NV = 15):
    
    # assumes that in qsystem, the particle dimensions are listed in the order (electron, nucleus)
    
    Zeepars = {'coeff': -2*np.pi, 'Bfield': Bfield}
    Zeemats = {}
    
    if nuc:
        Zeepars['gammaNuc'] = -431.6e-6 if NV == 15 else 307.7e-6 
        Zeepars['P1dims'] = qsystem.getdimensions()[0]
        Zeemats['S2'] = jmat(qsystem.getstot()[1])
    else:
        Zeepars['gammaNV'] = -2.8025
        Zeepars['P2dims'] = qsystem.getdimensions()[1]
        Zeemats['S1'] = jmat(qsystem.getstot()[0])
        
    Zeefunc = ZeeNucfunc if nuc else ZeeNVfunc
        
    return ham.Hamiltonian(Zeemats, Zeepars, Zeefunc)

def Hyperfine(qsystem, NV = 15):
    
    # assumes that in qsystem, the particle dimensions are listed in the order (electron, nucleus)
    
    A_N15 = np.array([[3.65,0,0],[0,3.65,0],[0,0,3.03]]) # Nitrogen-15 HF tensor
    A_N14 = np.array([[-2.62,0,0],[0,-2.62,0],[0,0,-2.162]]) # Nitrogen-14 HF tensor
    
    HFmats = {'A': A_N15 if NV == 15 else A_N14, 'Sz': qsystem.getSz()[0], 'Iz': qsystem.getSz()[1], 
              'Sx': qsystem.getSx()[0], 'Ix': qsystem.getSx()[1], 'Sy': qsystem.getSy()[0], 'Iy': qsystem.getSy()[1]}
    
    return ham.Hamiltonian(HFmats, None, HFfunc)
    
def dotproduct(vecV, vecU):
    #
    # Dot product between vector V = (Vx,Vy,Vz) and vector U = (Ux,Uy,Uz)
    #
    return sum([Vcomp*Ucomp for Vcomp, Ucomp in zip(vecV, vecU)])      
    
