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
                    nargs=2, default=[10,20], help='Domain size')
parser.add_argument('-U','--velocity',metavar='U',type=float,
                    nargs=2, default=[0.1,0.1], help='Magnitude of desired final longitudinal macroscopic velocity')
parser.add_argument('-t','--tau',metavar='t',type=float,
                    nargs=1, default=0.8, help='Relaxation time')
parser.add_argument('-T','--test')


def parabola(x,a,b,c):
    return a*x**2+b*x+c


def main(args):
    fd = lb.init(args.size,args.velocity,1)
    for x in range(50):
        ux,uy,rho = lb.macroVariables(fd)
        quiver(ux*rho,uy*rho)
        show()
        lb.collision(ux,uy,rho,fd,args.tau)
        lb.streamming(fd)
        lb.topbottomWalls(fd)            
    
    ux,uy,rho = lb.macroVariables(fd)
    
    popt, pcov = curve_fit(parabola, np.arange(len(ux[:,10])), ux[:,10])
    print popt
    y=np.arange(len(rho))
    U=parabola(y,popt[0],popt[1],popt[2])
    print y,U
#    plot(y,U)
#    plot(ux[:,2])
##    print np.abs(rho-1)#<=2e-13
    show()

        
if __name__ == "__main__":
    args=parser.parse_args()
    main(args)