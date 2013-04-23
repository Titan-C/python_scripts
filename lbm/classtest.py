# -*- coding: utf-8 -*-
"""
Test over 2D LBM
"""

import lb_D2Q9 
import numpy as np

#intial conditions
size=[21,3]
tau =0.505
r=1
F=[1e-6,0]

def test_consPBC():
    """Test if on periodic boundary conditions density and fluid velocity are conserved"""
    
    lat = lb_D2Q9.lattice(size,[0.1,0.],tau)
    for step in range(300):
       lat.updateMacroVariables()
       lat.collision()
       lat.streamming()

    magnitudes = lat.ux, lat.uy, lat.rho
    U=[0.1, 0., 1]
    for i in range(len(magnitudes)):
        assert ( np.abs(magnitudes[i] - U[i]) < 1e-15 ).all()

def test_consPBCforced():
    """Test if on periodic boundary conditions density and fluid velocity are conserved, when sistem is forced"""
    lat = lb_D2Q9.lattice(size,[0.1,0.1],tau)
    for step in range(300):
       lat.updateMacroVariables(F)
       lat.collision()
       lat.force(F)
       lat.streamming()

    assert ( np.abs(lat.rho - r) < 1e-12 ).all()

def test_consArbiWalls():
    """ Test if on top & bottom walls density is conserved"""
    lat = lb_D2Q9.lattice(size,[0.,0.1],tau)
    wall=np.zeros(size)
    wall[:,[0,-1]]=1
    lat.setWalls(wall)
    for step in range(2000):
       lat.updateMacroVariables()
       lat.collision()
       lat.streamming()
       lat.onWallBBNS_BC(wall)
    assert ( np.abs(lat.rho - r) < 1e-12 ).all()

def ev_U(F,L,nu): return 0.5*F * (L/2.)**2/nu
def U_sol(x,U,L): return U*(1-(2.*x/L)**2)

def anSolution(F,L,nu):
    U0=ev_U(F,L-2,nu)
    y=np.arange(L)-L/2
    return U_sol(y,U0,L-2)

def reynolds(F,L,D,nu):
    u=ev_U(F,L,nu)
    re=u*2/3.*D/nu
    return u,re
    
def test_poiseuille():
    lat = lb_D2Q9.lattice(size,[0.,0.],0.503)
    uxold=1
    step = 1
    wall=np.zeros(size)
    wall[[0,-1]]=1
    lat.setWalls(wall)
    while ( np.abs(uxold - lat.ux) > 1e-14 )[wall==0].all():
       lat.collision()
       lat.force(F)
       lat.streamming()
       lat.onWallBBNS_BC(wall)
       uxold=lat.ux
       lat.updateMacroVariables(F)
       step +=1
            
    an= anSolution(F[0],size[0],lat.nu)
    print lat.ux[:,1]
    print an
    print 'diff num - analytic: ', np.abs(lat.ux[:,1]-an)
    
    assert ( np.abs(lat.ux[:,1]-an) < 1e-5 )[1:-1].all()
