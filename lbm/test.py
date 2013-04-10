# -*- coding: utf-8 -*-
"""
Test over 2D LBM
"""
import lbm_kernel as lb
import numpy as np
from pylab import quiver,show,plot
#intial conditions
size=[64,64]
tau =1
r=1
F=[0.01,0]

def conservationPBC():
    """Test if on periodic boundary conditions density
       and fluid velocity are conserved"""
    for i in range(2):
        U = [2,1*i]
        fd = lb.init(size,U,r)
        for step in range(200):
            ux,uy,rho = lb.eqmacroVariables(fd)
            lb.collision(ux,uy,rho,fd,tau)
            lb.streamming(fd)

        magnitudes = lb.eqmacroVariables(fd)
        U.append(r)
        for i in range(len(magnitudes)):
            assert ( np.abs(magnitudes[i] - U[i]) < 1e-13 ).all()
        print 'Conservation of density and velocites on periodic BC: [OK]'
        
def conservationPBCforced():
    """Test if on periodic boundary conditions density
       and fluid velocity are conserved, when sistem is forced"""
    for i in range(2):
        U = [1,1*i]
        fd = lb.init(size,U,r)
        for step in range(200):
            ux,uy,rho = lb.eqmacroVariables(fd,F)
            lb.collision(ux,uy,rho,fd,tau)
            lb.force(ux,uy,rho,fd,tau,F)
            lb.streamming(fd)

        ux,uy,rho = lb.eqmacroVariables(fd)
        assert ( np.abs(rho - r) < 1e-13 ).all()
        print 'Conservation of density and velocites on periodic BC and forced: [OK]'

def conservationTBWalls():
    """ Test if on top & bottom walls density is conserved"""
    for i in range(1):
        U = [0.1,0.1*i]
        fd = lb.init(size,U,r)
        for step in range(2400):
            ux,uy,rho = lb.eqmacroVariables(fd)
            lb.collision(ux,uy,rho,fd,tau)
            lb.streamming(fd)
            lb.topbottomWallsBBNS_BC(fd)
        ux,uy,rho = lb.eqmacroVariables(fd)
        assert ( np.abs(rho - r) < 1e-12 ).all()
        plot(ux[:,10])
        show()
        
        print 'Conservation of density on periodic Flow and Upper and lower walls: [OK]'
def conservationTBWallsforced():
    """ Test if on top & bottom walls density is conserved"""
    for i in range(1):
        U = [0.1,0.1*i]
        fd = lb.init(size,U,r)
        for step in range(2800):
            ux,uy,rho = lb.eqmacroVariables(fd,F)
            lb.collision(ux,uy,rho,fd,tau)
            lb.force(ux,uy,rho,fd,tau,F)
            lb.streamming(fd)
            lb.topbottomWallsBBNS_BC(fd)
        ux,uy,rho = lb.eqmacroVariables(fd)
        assert ( np.abs(rho - r) < 1e-12 ).all()
        plot(ux[:,10])
        show()
        
        print 'Conservation of density on forced periodic Flow and Upper and lower walls: [OK]'


def main():
#    conservationPBC()
#    conservationPBCforced()
    conservationTBWalls()
    conservationTBWallsforced()
        
if __name__ == "__main__":
    main()