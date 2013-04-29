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
size=[3,3,3]
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

    print lat.U, lat.rho
    U=np.array(U)
    assert ( np.abs(lat.U-U) < 3e-15 ).all()
    assert ( np.abs(lat.rho-r) < 3e-15 ).all()
