# -*- coding: utf-8 -*-
"""
Poiseuille Flow 2D
"""

import argparse
import lb_D2Q9
from classtest import anSolution
import numpy as np
import matplotlib.pyplot as plt

#Input parsing
parser = argparse.ArgumentParser(description='LBM flow in slit')
parser.add_argument('-d','--diameter', metavar='D', type=int,
                    nargs=2, default=[25], help='Diameter of disk obstacle')
parser.add_argument('-F','--force',metavar='F',type=float,
                    nargs=2, default=[1e-6,0.], help='External force field')
parser.add_argument('-t','--tau',metavar='t',type=float,
                    nargs=1, default=[0.503], help='Relaxation time')
parser.add_argument('-i','--ifile',help='input file with macroscopic config')

def createobstacle(d,L):
    X=np.array([[x for x in range(d+1)] for y in range(d+1)])
    Y=np.array([[y for x in range(d+1)] for y in range(d+1)])
    x =X-d/2.
    y =Y-d/2.
    R = np.sqrt(x**2+y**2) <= d/2.
    wall=np.zeros([L,10*d])
    wall[42:43+d,d+3:2*d+4]=R
    wall[[0,-1]]=1 
    return wall
    

def main(args):
    fig = plt.figure()
    axU = fig.add_subplot(211)
    axR = fig.add_subplot(212)

    wall = createobstacle(args.diameter[0],111)
    lat = lb_D2Q9.lattice(wall.shape,[0.,0.],args.tau[0])
    lat.ux[:,0]=anSolution(5e-8,111,lat.nu) #set intel macroscopic eq vel
    lat.fd = lat.eqdistributions()#New equilibrium
    lat.setWalls(wall)
    finlet = np.copy( lat.fd[:,:,:1]) #get inlet macroscopic vel
    lat.updateMacroVariables()

    for i in range(400000):
        if i%50 == 0:
            axU.cla()
            axU.imshow(np.sqrt(lat.ux**2+lat.uy**2),vmin=0,vmax=0.2)
            fname = '%08d.png'%i
            dyUy,dxUy=np.gradient(lat.uy)
            dyUx,dxUx=np.gradient(lat.ux)
            vort=np.abs(dxUy-dyUx)
            axR.cla()
            axR.imshow(vort,vmin=0,vmax=0.1)
            fig.savefig('tot'+fname)
        
        lat.collision()
        lat.streamming()
        lat.onWallBBNS_BC(wall)
        lat.fd[:,:,-1]=np.copy(lat.fd[:,:,-2]) #gradient to cero on outlet
        lat.fd[:,:,:1]=np.copy(finlet)
        lat.updateMacroVariables()
    np.savez('out.npz',lat.ux,lat.uy,lat.rho,wall)

if __name__ == "__main__":
    args=parser.parse_args()
    main(args)
