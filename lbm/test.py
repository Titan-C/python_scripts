# -*- coding: utf-8 -*-
"""
Test over 2D LBM
"""
import lbm_kernel as lb
import numpy as np
from pylab import show,plot
#intial conditions
size=[33,5]
tau =1
r=1
F=[1e-6,0]

def init(LatticeSize,U,r):
    """Unitiate dicrete distribution functions for a given LatticeSize
       and velocity"""
    ux  = U[0]*np.ones(LatticeSize)
    uy  = U[1]*np.ones(LatticeSize)
    rho = r*np.ones(LatticeSize)
    return lb.eqdistributions(ux,uy,rho)

def conservationPBC():
    """Test if on periodic boundary conditions density
       and fluid velocity are conserved"""
    for i in range(2):
        U = [0.1,0.1]
        fd = init(size,U,r)
        for step in range(200):
            ux,uy,rho = lb.eqmacroVariables(fd)
            lb.collision(ux,uy,rho,fd,tau)
            lb.streamming(fd)

        magnitudes = lb.eqmacroVariables(fd)
        U.append(r)
        for i in range(len(magnitudes)):
            assert ( np.abs(magnitudes[i] - U[i]) < 1e-13 ).all()
        print 'Conservation of density and velocites on periodic BC: [OK]'

def conservationPBCforced():
    """Test if on periodic boundary conditions density
       and fluid velocity are conserved, when sistem is forced"""
    for i in range(2):
        U = [1,1*i]
        fd = init(size,U,r)
        for step in range(200):
            ux,uy,rho = lb.eqmacroVariables(fd,F)
            lb.collision(ux,uy,rho,fd,tau)
            lb.force(ux,uy,rho,fd,tau,F)
            lb.streamming(fd)

        ux,uy,rho = lb.eqmacroVariables(fd)
        assert ( np.abs(rho - r) < 1e-13 ).all()
        print 'Conservation of density and velocites on periodic BC and forced: [OK]'

def conservationArbiWalls():
    """ Test if on top & bottom walls density is conserved"""
    U = [0.0,0.1]
    fd = init(size,U,r)
    wall=np.zeros(size)
    wall[:,[0,-1]]=1
    lb.setWalls(fd,wall)
    for step in range(300):
        ux,uy,rho = lb.eqmacroVariables(fd)
        lb.collision(ux,uy,rho,fd,tau)
        lb.streamming(fd)
        lb.onWallBBNS_BC(fd,wall)
    ux,uy,rho = lb.eqmacroVariables(fd)
    assert ( np.abs(rho - r) < 1e-12 ).all()
    print 'Conservation of density on periodic Flow restrained arbitrary walls: [OK]'

def ev_nu(tau): return (tau-0.5)/3.
def ev_U(F,L,nu): return 0.5*F * (L/2.)**2/nu
def U_sol(x,U,L): return U*(1-(2.*x/L)**2)
def anSolution(F,L,tau):
    nu = ev_nu(tau)
    U0=ev_U(F,L-2,nu)
    y=np.arange(L)-L/2
    return U_sol(y,U0,L-2)

def conservationArbiWallsforced():
    """ Test if on top & bottom walls density is conserved"""
    U = [0.,0.]
    fd = init(size,U,r)
    wall=np.zeros(size)
    wall[[0,-1]]=1
    lb.setWalls(fd,wall)

    for step in range(3200):
        ux,uy,rho = lb.eqmacroVariables(fd,F)
        lb.collision(ux,uy,rho,fd,tau)
        lb.force(ux,uy,rho,fd,tau,F)
        lb.streamming(fd)
        lb.onWallBBNS_BC(fd,wall)
    ux,uy,rho = lb.eqmacroVariables(fd,F)
    assert ( np.abs(rho - r) < 1e-12 ).all()
    print 'Conservation of density on forced periodic Flow and Upper and lower walls: [OK]'
    an= anSolution(F[0],size[0],tau)
    assert ( np.abs(ux[:,1]-an) < 5e-6)[1:-1].all()
    print 'Parabolic velocity profile on Poiseuille Flow on Upper and lower walls: [OK]'

def main():
    conservationPBC()
    conservationPBCforced()
    conservationArbiWalls()
    conservationArbiWallsforced()
        
if __name__ == "__main__":
    main()