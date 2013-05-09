# -*- coding: utf-8 -*-
"""
Created on Wed Apr 24 10:04:00 2013

@author: Oscar Najera

Tests for D3Q41 lattice
"""

import lb_D3Q41 as lb
import numpy as np
from test_D2Q9 import ev_U, U_sol, reynolds
from pylab import show,plot


#intial conditions
size=[10,10,10]
tau =0.503
r=1
F=np.array([0,1e-7,0])
#def test_weights():
#    metrics=[[2.,1.,1.],[1.,2.,1.],[1.,1.,2.],[.5,1,1],[1,.5,1],[1,1,.5],[1,1,1]]
#    for metric in metrics:
#        lat = lb.lattice(size,[0,0,0],tau,metric)
#        print 'for metric: ', metric
#        assert ( np.abs(lb.sum(lat.w)-1) < 5e-15)
#        assert ( (np.abs(lb.sum(lat.E,axis=0)) < 5e-15).all() )
#
#
#def test_Macro():
#    """Test if Macroscopic variables can be generated from distributions functions"""
#    U=[0.05,0.05,0.05]
#    metrics=[[2.,1.,1.],[1.,2.,1.],[1.,1.,2.],[.5,1,1],[1,.5,1],[1,1,.5],[1,1,1]]
#    for metric in metrics:
#        lat = lb.lattice(size,U,tau,metric)
#        lat.updateMacroVariables()
#        print 'for metric: ', metric
#        print 'rho_err= ',np.max(np.abs(lat.rho-r))
#        print 'U_err0= ',np.abs(lat.U-U)
#        U=np.array(U)
#        assert ( np.abs(lat.rho-r) < 5e-15 ).all()
#        assert ( np.abs(lat.U-U) < 5e-15 ).all()
#
def test_consPBC():
    """Test if on periodic boundary conditions density and fluid velocity are conserved"""
    U=[0.05,0.05,0.05]
    size=[20,20,20]
    metrics=[[1,1,1],[.5,1,1],[1,.5,1],[1,1,.5],[2.,1.,1.],[1.,2.,1.],[1.,1.,2.]]
    for metric in metrics:
#        import pdb; pdb.set_trace()
        lat = lb.lattice(size,U,tau,metric)
        for step in range(700):
           lat.updateMacroVariables()
           lat.collision()
           lat.streamming()
        print 'for metric: ', metric
        print 'rho_err= ',np.max(np.abs(lat.rho-r))
        print 'U_err0= ',np.max(np.abs(lat.U-U))
        U=np.array(U)
        assert ( np.abs(lat.rho-r) < 5e-13 ).all()
        assert ( np.abs(lat.U-U) < 6e-13 ).all()
        plot((lat.U-U)[:,5,5,0],'x-');show()
#
#def test_consPBCforced():
#    """Test if on periodic boundary conditions density is conserved, when sistem is forced"""
#    lat = lb.lattice(size,[0.,0.,0.],tau)
#    for step in range(300):
#       lat.updateMacroVariables(F)
#       lat.collision()
#       lat.force(F)
#       lat.streamming()
#    print 'rho_err= ',np.max(np.abs(lat.rho-r))
#    print 'U=',lat.U[1,2,1,:]
#    assert ( np.abs(lat.rho - r) < 5e-10 ).all()
##
def test_consArbiWalls():
    """ Test if with top & bottom walls density is conserved"""
    U=[0.05,0.0,0.]
    size=[21,7,7]
    wall=np.zeros(size)
    wall[[0,1,2,-3,-2,-1]]=1
    metrics=[[2.,1.,1.],[1.,2.,1.],[1.,1.,2.],[.5,1,1],[1,.5,1],[1,1,.5],[1,1,1]]
    for metric in metrics:
        lat = lb.lattice(size,U,tau,metric)
        lat.setWalls(wall)
        for step in range(500):
           lat.updateMacroVariables()
           if step%20==0: plot(lat.U[:,5,5,0],'x-')#;import pdb;pdb.set_trace()
           lat.collision()
           lat.streamming()
           lat.onWallBBNS_BC(wall)
        show()
        print 'for metric: ', metric
        print 'rho_err= ',np.max(np.abs(lat.rho-r))
        assert ( np.abs(lat.rho - r) < 4e-15 ).all()

def anSolution(F,L,nu):
    U0=ev_U(F,L-4,nu)
    y=np.arange(L)-L/2
    return y,U_sol(y,U0,L-4)

#def test_poiseuille():
#    """Test poiseuille Flow"""
#    size=[3,3,63]
#    lat = lb.lattice(size,[0.,0.,0.],0.65)
#    wall=np.zeros(size)
#    wall[:,:,[0,1,-2,-1]]=1
#    lat.setWalls(wall)
#    
#    for i in range(20000):
#       lat.collision()
#       lat.force(F)
#       lat.streamming()
#       lat.onWallBBNS_BC(wall)
#       lat.updateMacroVariables(F)
#       if i%50==0: plot(lat.U[0,0,:,1],'x-'); show()#;import pdb;pdb.set_trace()
#
#    y,an= anSolution(F[0],size[1],lat.nu)
#    print lat.U[:,0,0,0]
#    print an
#    print 'diff num - analytic: ', np.abs(lat.U[0,0,:,0]-0.5*an)
#    print 'err_max= ',np.max(np.abs(lat.U[0,0,:,0]-0.5*an)[2:-2])
#    plot(y,lat.U[0,0,:,0],'x-')
#    plot(y,0.5*an)
#    show()
#    assert ( np.abs(lat.U[:,0,0,0]-an) < 7e-6 )[2:-2].all()