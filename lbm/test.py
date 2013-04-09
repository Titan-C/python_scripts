# -*- coding: utf-8 -*-
"""
Test over 2D LBM
"""
import lbm_kernel as lb
import numpy as np

def conservationPBC():
    """Test if on periodic boundary conditions density
       and fluid velocity are conserved"""
    size=[10,20]
    tau =1
    r=1
    U = [0.1,0.2]
    fd = lb.init(size,U,r)
    for step in range(200):
        lb.streamming(fd)
        ux,uy,rho = lb.macroVariables(fd)
        lb.collision(ux,uy,rho,fd,tau)
    
    magnitudes = lb.macroVariables(fd)
    U.append(r)
    for i in range(len(magnitudes)):
        assert ( np.abs(magnitudes[i] - U[i]) < 1e-13 ).all()
    print 'Conservation of density and velocites on periodic BC: [OK]'



def main():
    conservationPBC()

        
if __name__ == "__main__":
    main()