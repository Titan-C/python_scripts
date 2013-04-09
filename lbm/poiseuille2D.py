# -*- coding: utf-8 -*-
"""
Poiseuille Flow 2D
"""

import argparse
import lbm_kernel as lb
import numpy as np
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

def boundaries(fd):
    #up down, directions 2-4
    fd[4][0],fd[2][-1] = fd[2][-1], fd[4][0]
#    #diagonal 5-7
    fd[7][0],fd[5][-1] = np.roll( fd[5][-1],-1 ), np.roll(fd[7][0],1)
#    #diagonal 6-8
    fd[8][0], fd[6][-1] = np.roll( fd[6][-1],1 ),  np.roll(fd[8][0],-1)
#iteration step
def iteration(endVel,tau,fd):    
    lb.streamming(fd)
    boundaries(fd)
    rho, ux,uy = lb.macroVariables(fd)
    lb.collision(rho,ux,uy,fd,tau)
    lb.force(tau,endVel,rho,fd)
#    quiver(ux*rho,uy*rho)
#    show()    

def parabola(x,a,b,c):
    return a*x**2+b*x+c


def main(args):
   
    LatticeSize=args.size
    fd = lb.init(LatticeSize)
    for x in range(300):
        lb.iteration(args.velocity,args.tau,fd)
#        if x%20 == 0:
#            rho, ux,uy = macroVariables(fd)
#            plot(ux[:,2])
            
    
    rho, ux,uy = lb.macroVariables(fd)
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