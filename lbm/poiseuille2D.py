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
                    nargs=2, default=[33,32], help='Domain size')
parser.add_argument('-F','--force',metavar='U',type=float,
                    nargs=2, default=[1e-6,0.], help='External force field')
parser.add_argument('-U','--velocity',metavar='U',type=float,
                    nargs=2, default=[0.,0.], help='Magnitude of desired final longitudinal macroscopic velocity')
parser.add_argument('-t','--tau',metavar='t',type=float,
                    nargs=1, default=1.2, help='Relaxation time')
parser.add_argument('-T','--test')


def U_sol(x,U,L):
    return U*(1-(x/L)**2)

def ev_nu(tau): return (tau-0.5)/3.
def ev_U(F,L,nu): return 0.5*F * (L/2.)**2/nu

def anSolution(F,L,tau):
    L2= L/2
    nu = ev_nu(tau)
    U0=ev_U(F,L,nu)
    y=np.arange(1.*L)-L2
    return U0*(1-(y/L2)**2)

def main(args):
    size = args.size
    U = args.velocity
    ux  = U[0]*np.ones(size)
    uy  = U[1]*np.ones(size)
    rho = np.ones(size)
    fd = lb.eqdistributions(ux,uy,rho)
    uxold=1
    step = 1
    while ( np.abs(uxold - ux) > 1e-14 ).all():
        lb.collision(ux,uy,rho,fd,args.tau)
        lb.force(ux,uy,rho,fd,args.tau,args.force)
        lb.streamming(fd)
        lb.topbottomWallsBBNS_BC(fd)
        uxold=ux
        ux,uy,rho = lb.eqmacroVariables(fd,args.force)
        step +=1
    print step, 'iter'
            
    
    ux,uy,rho = lb.eqmacroVariables(fd)
    domain = np.arange(len(ux[:,10]))- size[0]/2
    popt, pcov = curve_fit(U_sol,domain , ux[:,10])
    print ux[:,10],popt
    U=U_sol(domain,popt[0],popt[1])
    print domain,U
    print 'diff num - fit: ', np.abs(ux[:,10]-U)
    an= anSolution(args.force[0],size[0],args.tau)
    print an
    print 'diff num - analytic: ', np.abs(ux[:,10]-an)
    plot(domain,ux[:,10])
    plot(domain,U)
    show()
    
    return ux

        
if __name__ == "__main__":
    args=parser.parse_args()
    main(args)