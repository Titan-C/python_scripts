# -*- coding: utf-8 -*-
"""
Poiseuille Flow 2D
"""

import argparse
import lbm_kernel as lb
import numpy as np
import matplotlib.pyplot as plt
from poiseuille2D import anSolution

#Input parsing
parser = argparse.ArgumentParser(description='LBM flow in slit')
parser.add_argument('-d','--diameter', metavar='D', type=int,
                    nargs=2, default=[50], help='Diameter of disk obstacle')
parser.add_argument('-F','--force',metavar='F',type=float,
                    nargs=2, default=[1e-6,0.], help='External force field')
parser.add_argument('-t','--tau',metavar='t',type=float,
                    nargs=1, default=[0.535], help='Relaxation time')
parser.add_argument('-i','--ifile',help='input file with macroscopic config')

def createobstacle(d,L):
    X=np.array([[x for x in range(d+1)] for y in range(d+1)])
    Y=np.array([[y for x in range(d+1)] for y in range(d+1)])
    x =X-d/2.
    y =Y-d/2.
    R = np.sqrt(x**2+y**2) <= d/2.
    wall=np.zeros([L,9*d])
    wall[25:26+d,d+3:2*d+4]=R
    wall[[0,-1]]=1 
    return wall
    

def main(args):
    fig = plt.figure()
    axU = fig.add_subplot(211)
    axR = fig.add_subplot(212)

    wall = createobstacle(args.diameter[0],90)
    ux  = np.zeros(wall.shape)
    an=anSolution(args.force[0],wall.shape[0],args.tau[0])
    ux[:,0]=an
    uy  = np.zeros(wall.shape)
    rho = np.ones(wall.shape)
    fd = lb.eqdistributions(ux,uy,rho)
    lb.setWalls(fd,wall)
    for i in range(1000):
        ux,uy,rho = lb.eqmacroVariables(fd,args.force)
        ux[:,0]=an
        if i%30 == 0:
            axU.cla()
            axU.imshow(np.sqrt(ux**2+uy**2))
            fname = '%08d.png'%i
            dyUy,dxUy=np.gradient(uy)
            dyUx,dxUx=np.gradient(ux)
            vort=np.abs(dxUy-dyUx)
            axR.cla()
            axR.imshow(vort)
            fig.savefig('tot'+fname)
        
        lb.collision(ux,uy,rho,fd,args.tau[0])
        lb.force(ux,uy,rho,fd,args.tau[0],args.force)
        lb.streamming(fd)
        fd[:,:,-1]=np.copy(fd[:,:,-2]) #gradient to cero on outlet
        lb.onWallBBNS_BC(fd,wall)
    ux,uy,rho = lb.eqmacroVariables(fd,args.force)
    np.savez('out.npz',ux,uy,rho,wall,an)

if __name__ == "__main__":
    args=parser.parse_args()
    main(args)
