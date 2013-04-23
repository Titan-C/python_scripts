# -*- coding: utf-8 -*-
"""
Created on Tue Apr 23 13:45:43 2013

@author: Oscar Najera
Lattice Boltzmann 2D en D2Q9 class
"""
import numpy as np

class lattice:
    def __init__(self,LatticeSize,tau,r=1):
        """Unitiate dicrete distribution functions for a given LatticeSize
           and velocity"""
        self.ux  = np.zeros(LatticeSize)
        self.uy  = np.zeros(LatticeSize)
        self.rho = r*np.ones(LatticeSize)
        self.fd  = self.eqdistributions()
        self.tau = tau
        self.nu  = (tau-0.5)/3.

    def eqdistributions(self):
        """Evaluate Equilibrium distribution functions"""
        rho,ux,uy = self.rho, self.ux, self.uy
        ux2,uy2,u2 = self.velocities()
    
        f0 = 4.0/9.0    *rho*(1                           -1.5*u2)
        f1 = 1.0/9.0    *rho*(1 +3*ux     +4.5*ux2        -1.5*u2)
        f2 = 1.0/9.0    *rho*(1 +3*uy     +4.5*uy2        -1.5*u2)
        f3 = 1.0/9.0    *rho*(1 -3*ux     +4.5*ux2        -1.5*u2)
        f4 = 1.0/9.0    *rho*(1 -3*uy     +4.5*uy2        -1.5*u2)
        f5 = 1.0/36.0 *rho*(1 +3*(ux+uy)  +4.5*(ux+uy)**2 -1.5*u2)
        f6 = 1.0/36.0 *rho*(1 +3*(-ux+uy) +4.5*(-ux+uy)**2-1.5*u2)
        f7 = 1.0/36.0 *rho*(1 +3*(-ux-uy) +4.5*(-ux-uy)**2-1.5*u2)
        f8 = 1.0/36.0 *rho*(1 +3*(ux-uy)  +4.5*(ux-uy)**2 -1.5*u2)
    
        return np.array([f0,f1,f2,f3,f4,f5,f6,f7,f8])
    
    def velocities(self):
        ux2 = self.ux**2
        uy2 = self.uy**2
        u2  =ux2+uy2
        return ux2,uy2,u2

    #Common operations
    def updateMacroVariables(self,F=[0.,0.]):
        """Returns the macroscopic variables density and velocity.(rho, ux,uy)
           Equilibrium velocity includes forcing"""
        fd= self.fd
        
        self.rho = fd[0]+fd[1]+fd[2]+fd[3]+fd[4]+fd[5]+fd[6]+fd[7]+fd[8]
        self.rho[self.rho<1e-13]=1
        self.ux  = (fd[1]-fd[3] + (fd[5]-fd[6]) + (fd[8]-fd[7]))/self.rho
        self.uy  = (fd[2]-fd[4] + (fd[5]+fd[6]) - (fd[8]+fd[7]))/self.rho
        #Additional forcing
        self.ux  += 0.5*F[0]/self.rho
        self.uy  += 0.5*F[1]/self.rho

    def force(self,F):
        """Adds microscopic forcing term"""
        cf = (1-0.5/self.tau)*3.0
    
        F0 = 4.0/9.0  *(   -self.ux*F[0]  -    self.uy *F[1] )
        F1 = 1.0/9.0  *( (1-self.ux)*F[0] -    self.uy *F[1] + 3*self.ux*F[0] )
        F2 = 1.0/9.0  *(   -self.ux *F[0] + (1-self.uy)*F[1] + 3*self.uy*F[1] )
        F3 = 1.0/9.0  *((-1-self.ux)*F[0] -    self.uy *F[1] + 3*self.ux*F[0] )
        F4 = 1.0/9.0  *(   -self.ux *F[0] +(-1-self.uy)*F[1] + 3*self.uy*F[1] )
    
        F5 = 1.0/36.0 *( (1-self.ux)*F[0] + (1-self.uy)*F[1] + 3*( self.ux+self.uy)*( F[0]+F[1]) )
        F6 = 1.0/36.0 *((-1-self.ux)*F[0] + (1-self.uy)*F[1] + 3*(-self.ux+self.uy)*(-F[0]+F[1]) )
        F7 = 1.0/36.0 *((-1-self.ux)*F[0] +(-1-self.uy)*F[1] + 3*(-self.ux-self.uy)*(-F[0]-F[1]) )
        F8 = 1.0/36.0 *( (1-self.ux)*F[0] +(-1-self.uy)*F[1] + 3*( self.ux-self.uy)*( F[0]-F[1]) )
    
        self.fd+= cf*np.array([F0,F1,F2,F3,F4,F5,F6,F7,F8])

#LBM steps
    def streamming(self):
        """Stream distributions funtions along al discrete velocity directions
           keep periodic boundary conditions on all edges"""
        self.fd[1] = np.roll(self.fd[1], 1,axis=1)  #x >
        self.fd[2] = np.roll(self.fd[2],-1,axis=0)  #y ^
        self.fd[3] = np.roll(self.fd[3],-1,axis=1)  #x <
        self.fd[4] = np.roll(self.fd[4], 1,axis=0)  #y v
        self.fd[5] = np.roll( np.roll(self.fd[5],-1,axis=0) , 1,axis=1) #y ^,x >
        self.fd[6] = np.roll( np.roll(self.fd[6],-1,axis=0) ,-1,axis=1) #y ^,x <
        self.fd[7] = np.roll( np.roll(self.fd[7], 1,axis=0) ,-1,axis=1) #y v,x <
        self.fd[8] = np.roll( np.roll(self.fd[8], 1,axis=0) , 1,axis=1) #y v,x >

    def collision(self):
        """Modify distribution funtions according to collision term"""
        feq = self.eqdistributions()
        self.fd -=  1.0/self.tau*(self.fd - feq)

##boundaries handling
    def setWalls(self,walls):
        """Anulates fd where a wall(solid) is declared"""
        self.fd[:,walls==1]=0

    def onWallBBNS_BC(self,walls):
        """Bounceback no-slip boundary conditions for arbitrary walls"""
        Ny,Nx=walls.shape
        solid=np.array(np.where(walls==1))
        #generate direction vectors, mantain periodic boundary conditions
        goright=solid[1]+1
        goright[goright==Nx]=0
        goleft=solid[1]-1
        goleft[goleft==-1]+=Nx
        goup=solid[0]-1
        goup[goup==-1]+=Ny
        godown=solid[0]+1
        godown[godown==Ny]=0
        #colision to right,left, direction 1-3
        self.fd[1,solid[0],goright] = self.fd[3,solid[0],solid[1]]
        self.fd[3,solid[0],goleft] = self.fd[1,solid[0],solid[1]]
        #colision top,botto direction 2-4
        self.fd[2,goup,solid[1]] = self.fd[4,solid[0],solid[1]]
        self.fd[4,godown,solid[1]] = self.fd[2,solid[0],solid[1]]
        #diagonal 5-7
        self.fd[5,goup,goright] = self.fd[7,solid[0],solid[1]]
        self.fd[7,godown,goleft] = self.fd[5,solid[0],solid[1]]
        #diagonal 6-8
        self.fd[6,goup,goleft] = self.fd[8,solid[0],solid[1]]
        self.fd[8,godown,goright] = self.fd[6,solid[0],solid[1]]
        
        #clean solid
        self.fd[:,solid[0],solid[1]]=0
