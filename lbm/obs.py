# -*- coding: utf-8 -*-
"""
Poiseuille Flow 2D
"""

import argparse
import lbm_kernel as lb
import numpy as np
from pylab import show,plot,imshow

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
parser.add_argument('-r','--radius',metavar='r',type=float,
		    nargs=1, default=[20], help='Radius of obstacle')
parser.add_argument('-T','--test')

def createobstacle(d):
    X=np.array([[x for x in range(d+1)] for y in range(d+1)])
    Y=np.array([[y for x in range(d+1)] for y in range(d+1)])
    x =X-d/2.
    y =Y-d/2.
    R = np.sqrt(x**2+y**2) <= d/2.
    wall=np.zeros([3*(d+1),8*d])
    wall[d+1:-d-1,d+3:2*d+4]=R
    wall[[0,-1]]=1 
    return wall

def main(args):
    wall = createobstacle(40)
    U = args.velocity
    ux  = U[0]*np.ones(wall.shape)
    uy  = U[1]*np.ones(wall.shape)
    rho = np.ones(wall.shape)
    fd = lb.eqdistributions(ux,uy,rho)
    lb.setWalls(fd,wall)
    for x in range(300):
        
        ux,uy,rho = lb.eqmacroVariables(fd,[0.,0.])
        imshow(np.sqrt(ux**2+uy**2))        
        show()	
        lb.collision(ux,uy,rho,fd,args.tau[0])
        lb.force(ux,uy,rho,fd,args.tau[0],args.force)
        lb.streamming(fd)
        lb.onWallBBNS_BC(fd,wall)

        
if __name__ == "__main__":
    args=parser.parse_args()
    main(args)
