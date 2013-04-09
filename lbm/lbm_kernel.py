# -*- coding: utf-8 -*-
"""
Lattice Boltzmann Method test exercise
Fluid Flow inside periodic slit with external force field

"""

import numpy as np
import argparse
from pylab import quiver,show,plot
from scipy.optimize import curve_fit

#Input parsing
parser = argparse.ArgumentParser(description='LBM flow in slit')
parser.add_argument('-s','--size', metavar='n', type=int,
                    nargs=2, default=[32,128], help='Domain size')
parser.add_argument('-U','--velocity',metavar='U',type=float,
                    nargs=1, default=0.1, help='Magnitude of desired final longitudinal macroscopic velocity')
parser.add_argument('-t','--tau',metavar='t',type=float,
                    nargs=1, default=0.8, help='Relaxation time')
parser.add_argument('-T','--test')

#initiate dicrete distribution functions
def init(LatticeSize):

    ux  = 3*np.ones(LatticeSize)#+0.3*(np.random.random_sample(LatticeSize)-0.5)
    uy  = np.zeros(LatticeSize)#+0.2*(np.random.random_sample(LatticeSize)-0.5)
    rho = np.ones(LatticeSize)#*(1+0.1*np.sin(np.arange(LatticeSize[1])*2*np.pi/LatticeSize[1]))
    return eqdistributions(rho,ux,uy)

#Evaluate Equilibrium distribution functions
def eqdistributions(rho,ux,uy):
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

#Streaming
def streamming(fd):
    fd[1] = np.roll(fd[1], 1,axis=1)  #x >
    fd[2] = np.roll(fd[2],-1,axis=0)  #y ^
    fd[3] = np.roll(fd[3],-1,axis=1)  #x <
    fd[4] = np.roll(fd[4], 1,axis=0)  #y v
    fd[5] = np.roll( np.roll(fd[5],-1,axis=0) , 1,axis=1) #y ^,x >
    fd[6] = np.roll( np.roll(fd[6],-1,axis=0) ,-1,axis=1) #y ^,x <
    fd[7] = np.roll( np.roll(fd[7], 1,axis=0) ,-1,axis=1) #y v,x <
    fd[8] = np.roll( np.roll(fd[8], 1,axis=0) , 1,axis=1) #y v,x >
    boundary(fd)

def collision(rho,ux,uy,fd,tau):
    feq = eqdistributions(rho,ux,uy)
    for a in range(len(feq)):
        fd[a] -=  1.0/tau*(fd[a] - feq[a])

def boundary(fd):
    #up down, directions 2-4
    fd[4][0],fd[2][-1] = fd[2][-1], fd[4][0]
#    #diagonal 5-7
    fd[7][0],fd[5][-1] = np.roll( fd[5][-1],-1 ), np.roll(fd[7][0],1)
#    #diagonal 6-8
    fd[8][0], fd[6][-1] = np.roll( fd[6][-1],1 ),  np.roll(fd[8][0],-1)

#iteration step
def iteration(endVel,tau,fd):    
    streamming(fd)
    rho, ux,uy = macroVariables(fd)
    collision(rho,ux,uy,fd,tau)
    force(tau,endVel,rho,fd)
#    quiver(ux*rho,uy*rho)
#    show()

#Common operations
def macroVariables(fd):
    rho = sum(fd)
    ux  = (fd[1]-fd[3] + (fd[5]-fd[6]) + (fd[8]-fd[7]))/rho
    uy  = (fd[2]-fd[4] + (fd[5]+fd[6]) - (fd[8]+fd[7]))/rho
    return rho,ux,uy

def velocities(ux,uy):
    ux2 =ux**2
    uy2 =uy**2
    u2  =ux2+uy2
    return ux2,uy2,u2
    
def force(tau,endVel,rho,fd):
    nu=(tau-0.5)/3.0
    F =8*nu*endVel*rho/(6.0*len(rho)**2)
    fd[[1,5,8]]+=F
    fd[[3,6,7]]-=F

def parabola(x,a,b,c):
    return a*x**2+b*x+c


def main(args):
   
    LatticeSize=args.size
    fd = init(LatticeSize)
    for x in range(300):
        iteration(args.velocity,args.tau,fd)
#        if x%20 == 0:
#            rho, ux,uy = macroVariables(fd)
#            plot(ux[:,2])
            
    
    rho, ux,uy = macroVariables(fd)
#    quiver(ux*rho,uy*rho)
    popt, pcov = curve_fit(parabola, np.arange(len(ux[:,10])), ux[:,10])
    print popt
    y=np.arange(len(rho))
    U=parabola(y,popt[0],popt[1],popt[2])
    print y,U
    plot(y,U)
    plot(ux[:,2])
#    print np.abs(rho-1)#<=2e-13
    show()

        
if __name__ == "__main__":
    args=parser.parse_args()
    main(args)