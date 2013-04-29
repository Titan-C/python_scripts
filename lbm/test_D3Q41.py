# -*- coding: utf-8 -*-
"""
Created on Wed Apr 24 10:04:00 2013

@author: Oscar Najera

Tests for D3Q41 lattice
"""

import lb_D3Q41 as lb
import numpy as np
from test_D2Q9 import ev_U, U_sol, anSolution, reynolds

#intial conditions
size=[10,3,3]
tau =0.65
r=1
F=np.array([1e-6,0,0])

def test_consPBC():
    """Test if on periodic boundary conditions density and fluid velocity are conserved"""
    U=[0.05,0.05,0.05]
#    import pdb; pdb.set_trace()
    lat = lb.lattice(size,U,tau)
    for step in range(300):
       lat.updateMacroVariables()
       lat.collision()
       lat.streamming()

    print np.abs(lat.U-U), np.abs(lat.rho-r)
    U=np.array(U)
    assert ( np.abs(lat.U-U) < 3e-15 ).all()
    assert ( np.abs(lat.rho-r) < 3e-15 ).all()

def test_consArbiWalls():
    """ Test if with top & bottom walls density is conserved"""
    U=[0.08,0.0,0.]
    lat = lb.lattice(size,U,tau)
    wall=np.zeros(size)
    wall[[0,1,-2,-1]]=1
    lat.setWalls(wall)
#    import pdb; pdb.set_trace()
    for step in range(300):
       lat.updateMacroVariables()
       lat.collision()
       lat.streamming()
       lat.onWallBBNS_BC(wall)
    print lat.rho-1.
    assert ( np.abs(lat.rho - r) < 4e-14 ).all()