# -*- coding: utf-8 -*-
"""
Created on Wed Apr 24 10:04:00 2013

@author: Oscar Najera

Tests for D3Q19 lattice
"""

import lb_D3Q19
import numpy as np
from test_D2Q9 import ev_U, U_sol, anSolution, reynolds

#intial conditions
size=[10,10,10]
tau =0.503
r=1
F=np.array([1e-6,0,0])

def test_consPBC():
    """Test if on periodic boundary conditions density and fluid velocity are conserved"""
    U=[0.1,0.1,0.1]
    lat = lb_D3Q19.lattice(size,U,tau)
    for step in range(300):
       lat.updateMacroVariables()
       lat.collision()
       lat.streamming()

    print lat.U, lat.rho
    U=np.array(U)
    assert ( np.abs(lat.U-U) < 3e-15 ).all()
    assert ( np.abs(lat.rho-r) < 3e-15 ).all()

def test_consPBCforced():
    """Test if on periodic boundary conditions density and fluid velocity are conserved, when sistem is forced"""
    lat = lb_D3Q19.lattice(size,[0.,0.,0.],tau)
    for step in range(300):
       lat.updateMacroVariables(F)
       lat.collision()
       lat.force(F)
       lat.streamming()
    print lat.rho-1.
    assert ( np.abs(lat.rho - r) < 4e-14 ).all()

def test_consArbiWalls():
    """ Test if with top & bottom walls density is conserved"""
    U=[0.1,0.0,0.]
    lat = lb_D3Q19.lattice(size,U,tau)
    wall=np.zeros(size)
    wall[[0,-1]]=1
    lat.setWalls(wall)
    for step in range(300):
       lat.updateMacroVariables()
       lat.collision()
       lat.streamming()
       lat.onWallBBNS_BC(wall)
    print lat.rho-1.
    assert ( np.abs(lat.rho - r) < 4e-14 ).all()


def test_poiseuille():
    """Test poiseuille Flow"""
    size=[21,1,1]
    lat = lb_D3Q19.lattice(size,[0.,0.,0.],0.65)
    wall=np.zeros(size)
    wall[[0,-1]]=1
    lat.setWalls(wall)
    for i in range(18980):
       lat.collision()
       lat.force(F)
       lat.streamming()
       lat.onWallBBNS_BC(wall)
       lat.updateMacroVariables(F)
            
    an= anSolution(F[0],size[0],lat.nu)
    print lat.U[:,0,0,0]
    print an
    print 'diff num - analytic: ', np.abs(lat.U[:,0,0,0]-an)
    
    assert ( np.abs(lat.U[:,0,0,0]-an) < 2.25e-6 )[1:-1].all()
