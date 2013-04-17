# -*- coding: utf-8 -*-
"""
Poiseuille Flow 2D
"""

import argparse
import lbm_kernel as lb
import numpy as np
from pylab import show,plot

#Input parsing
parser = argparse.ArgumentParser(description='LBM flow in slit')
parser.add_argument('-s','--size', metavar='n', type=int,
                    nargs=2, default=[31,10], help='Domain size')
parser.add_argument('-F','--force',metavar='U',type=float,
                    nargs=2, default=[1e-6,0.], help='External force field')
parser.add_argument('-U','--velocity',metavar='U',type=float,
                    nargs=2, default=[0.,0.], help='Magnitude of desired final longitudinal macroscopic velocity')
parser.add_argument('-t','--tau',metavar='t',type=float,
                    nargs=1, default=[1.1], help='Relaxation time')
parser.add_argument('-T','--test')



def ev_nu(tau): return (tau-0.5)/3.
def ev_U(F,L,nu): return 0.5*F * (L/2.)**2/nu
def U_sol(x,U,L): return U*(1-(2.*x/L)**2)

def anSolution(F,L,tau):
    nu = ev_nu(tau)
    U0=ev_U(F,L-2,nu)
    y=np.arange(L)-L/2
    return U_sol(y,U0,L-2)

def reynolds(F,L,D,tau):
    u=ev_U(F,L,ev_nu(tau))
    re=u*2/3.*D/ev_nu(tau)
    return u,re

def main(args):
    size = args.size
    U = args.velocity
    ux  = U[0]*np.ones(size)
    uy  = U[1]*np.ones(size)
    rho = np.ones(size)
    fd = lb.eqdistributions(ux,uy,rho)
    uxold=1
    step = 1
    wall=np.zeros(size)
    wall[[0,-1]]=1
    lb.setWalls(fd,wall)
    ux,uy,rho = lb.eqmacroVariables(fd,[0.,0.])
    while ( np.abs(uxold - ux) > 1e-10 )[wall==0].all():
        lb.collision(ux,uy,rho,fd,args.tau[0])
        lb.force(ux,uy,rho,fd,args.tau[0],args.force)
        lb.streamming(fd)
        lb.onWallBBNS_BC(fd,wall)
        uxold=ux
        ux,uy,rho = lb.eqmacroVariables(fd,args.force)
        step +=1
    print step, 'iter'
            
    ux,uy,rho = lb.eqmacroVariables(fd,args.force)
    domain = np.arange(len(ux[:,1]))- size[0]/2
    print ux[:,1]
    an= anSolution(args.force[0],size[0],args.tau[0])
    print an
    print 'diff num - analytic: ', np.abs(ux[:,1]-an)
    plot(domain,ux[:,1],'x-')
    plot(domain,an)
    show()
    
    return ux

        
if __name__ == "__main__":
    args=parser.parse_args()
    main(args)
