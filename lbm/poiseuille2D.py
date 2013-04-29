# -*- coding: utf-8 -*-
"""
Poiseuille Flow 2D
"""

import argparse
import lb_D2Q9 as lb
import numpy as np
from pylab import show,plot
from classtest import ev_U, U_sol, anSolution, reynolds

#Input parsing
parser = argparse.ArgumentParser(description='LBM flow in slit')
parser.add_argument('-s','--size', metavar='n', type=int,
                    nargs=2, default=[31,1], help='Domain size')
parser.add_argument('-F','--force',metavar='U',type=float,
                    nargs=2, default=[1e-6,0.], help='External force field')
parser.add_argument('-U','--velocity',metavar='U',type=float,
                    nargs=2, default=[0.,0.], help='Magnitude of desired final longitudinal macroscopic velocity')
parser.add_argument('-t','--tau',metavar='t',type=float,
                    nargs=1, default=[0.65], help='Relaxation time')
parser.add_argument('-T','--test')

def main(args):
    size = args.size
    U = args.velocity
    lat = lb.lattice(args.size,args.velocity,args.tau[0])
    F=np.array(args.force)
    wall=np.zeros(args.size)
    wall[[0,-1]]=1
    lat.setWalls(wall)
    lat.updateMacroVariables(F)

    uxold=1
    step = 1
    while ( np.abs(uxold - lat.ux) > 1e-12 )[wall==0].all():
        lat.collision()
        lat.force(F)
        lat.streamming()
        lat.onWallBBNS_BC(wall)
        uxold=lat.ux
        lat.updateMacroVariables(F)
        step +=1
    print step, 'iter'
    domain = np.arange(len(lat.ux))- size[0]/2
    print lat.ux
    an= anSolution(args.force[0],size[0],lat.nu)
    print an
    print 'diff num - analytic: ', np.abs(lat.ux.reshape(size[0])-an)
    plot(domain,lat.ux,'x-')
    plot(domain,an)
    show()

    return lat.ux

        
if __name__ == "__main__":
    args=parser.parse_args()
    main(args)
