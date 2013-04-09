# -*- coding: utf-8 -*-
"""
Lattice Boltzmann Method test exercise
Fluid Flow inside periodic slit with external force field

"""

import numpy as np

def eqdistributions(ux,uy,rho):
    """Evaluate Equilibrium distribution functions"""
    ux2,uy2,u2 = velocities(ux,uy)

    f0 = 4.0/9.0    *rho*(1                           -1.5*u2)
    f1 = 1.0/9.0    *rho*(1 +3*ux     +4.5*ux2        -1.5*u2)
    f2 = 1.0/9.0    *rho*(1 +3*uy     +4.5*uy2        -1.5*u2)
    f3 = 1.0/9.0    *rho*(1 -3*ux     +4.5*ux2        -1.5*u2)
    f4 = 1.0/9.0    *rho*(1 -3*uy     +4.5*uy2        -1.5*u2)
    f5 = 1.0/36.0 *rho*(1 +3*(ux+uy)  +4.5*(ux+uy)**2 -1.5*u2)
    f6 = 1.0/36.0 *rho*(1 +3*(-ux+uy) +4.5*(-ux+uy)**2-1.5*u2)
    f7 = 1.0/36.0 *rho*(1 +3*(-ux-uy) +4.5*(-ux-uy)**2-1.5*u2)
    f8 = 1.0/36.0 *rho*(1 +3*(ux-uy)  +4.5*(ux-uy)**2 -1.5*u2)

    return np.array([f0,f1,f2,f3,f4,f5,f6,f7,f8])

def streamming(fd):
    """Stream distributions funtions along al discrete velocity directions
       keep periodic boundary conditions on all edges"""
    fd[1] = np.roll(fd[1], 1,axis=1)  #x >
    fd[2] = np.roll(fd[2],-1,axis=0)  #y ^
    fd[3] = np.roll(fd[3],-1,axis=1)  #x <
    fd[4] = np.roll(fd[4], 1,axis=0)  #y v
    fd[5] = np.roll( np.roll(fd[5],-1,axis=0) , 1,axis=1) #y ^,x >
    fd[6] = np.roll( np.roll(fd[6],-1,axis=0) ,-1,axis=1) #y ^,x <
    fd[7] = np.roll( np.roll(fd[7], 1,axis=0) ,-1,axis=1) #y v,x <
    fd[8] = np.roll( np.roll(fd[8], 1,axis=0) , 1,axis=1) #y v,x >

def topbottomWalls(fd):
    """Bounceback no-slip boundary conditions for top and bottom walls"""
    #up down, directions 2-4
    fd[4][0],fd[2][-1] = fd[2][-1], fd[4][0]
#    #diagonal 5-7
    fd[7][0],fd[5][-1] = np.roll( fd[5][-1],-1 ), np.roll(fd[7][0],1)
#    #diagonal 6-8
    fd[8][0], fd[6][-1] = np.roll( fd[6][-1],1 ),  np.roll(fd[8][0],-1)

def collision(ux,uy,rho,fd,tau):
    feq = eqdistributions(ux,uy,rho)
    for a in range(len(feq)):
        fd[a] -=  1.0/tau*(fd[a] - feq[a])


def init(LatticeSize,U,r):
    """Unitiate dicrete distribution functions for a given LatticeSize
       and velocity"""
    ux  = U[0]*np.ones(LatticeSize)
    uy  = U[1]*np.ones(LatticeSize)
    rho = r*np.ones(LatticeSize)
    return eqdistributions(ux,uy,rho)

#Common operations
def macroVariables(fd):
    """Returns the macroscopic variables density and velocity.(rho, ux,uy)"""
    rho = sum(fd)
    ux  = (fd[1]-fd[3] + (fd[5]-fd[6]) + (fd[8]-fd[7]))/rho
    uy  = (fd[2]-fd[4] + (fd[5]+fd[6]) - (fd[8]+fd[7]))/rho
    return ux,uy,rho

def velocities(ux,uy):
    ux2 =ux**2
    uy2 =uy**2
    u2  =ux2+uy2
    return ux2,uy2,u2
    
def preasuredrop(tau,endVel,rho,fd):
    """Preasure drop as implemented by Succi"""
    nu=(tau-0.5)/3.0
    F =8*nu*endVel*rho/(6.0*len(rho)**2)
    fd[[1,5,8]]+=F
    fd[[3,6,7]]-=F