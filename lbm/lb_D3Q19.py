# -*- coding: utf-8 -*-
"""
Created on Tue Apr 23 13:45:43 2013

@author: Oscar Najera
Lattice Boltzmann 3D en D3Q19 class
"""
import numpy as np

class lattice:
    def __init__(self,LatticeSize,U,tau,r=1):
        """Unitiate dicrete distribution functions for a given LatticeSize
           and velocity"""
        self.ux = U[0]*np.ones(LatticeSize)
        self.uy = U[1]*np.ones(LatticeSize)
        self.uz = U[1]*np.ones(LatticeSize)
        self.rho = r*np.ones(LatticeSize)
        self.fd  = self.eqdistributions()
        self.tau = tau
        self.nu  = (tau-0.5)/3.

    def eqdistributions(self):
        """Evaluate Equilibrium distribution functions"""
        rho,ux,uy,uz = self.rho, self.ux, self.uy, self.uz
        ux2,uy2,uz2,u2 = self.velocities()
        #no shift
        f0 = 1.0/3.0    *rho*(1                           -1.5*u2)
        #mayor axes
        f1 = 1.0/18.0    *rho*(1 +3*ux     +4.5*ux2        -1.5*u2)
        f2 = 1.0/18.0    *rho*(1 -3*ux     +4.5*ux2        -1.5*u2)
        f3 = 1.0/18.0    *rho*(1 +3*uy     +4.5*uy2        -1.5*u2)
        f4 = 1.0/18.0    *rho*(1 -3*uy     +4.5*uy2        -1.5*u2)
        f5 = 1.0/18.0    *rho*(1 +3*uz     +4.5*uz2        -1.5*u2)
        f6 = 1.0/18.0    *rho*(1 -3*uz     +4.5*uz2        -1.5*u2)
        #Diagonal axes on planes
        f7 = 1.0/36.0 *rho*(1 +3*(ux+uy)  +4.5*(ux+uy)**2  -1.5*u2)
        f8 = 1.0/36.0 *rho*(1 +3*(-ux-uy) +4.5*(-ux-uy)**2 -1.5*u2)
        f9 = 1.0/36.0 *rho*(1 +3*(ux-uy)  +4.5*(ux-uy)**2  -1.5*u2)
        f10= 1.0/36.0 *rho*(1 +3*(-ux+uy) +4.5*(-ux+uy)**2 -1.5*u2)
        
        f11= 1.0/36.0 *rho*(1 +3*(ux+uz)  +4.5*(ux+uz)**2  -1.5*u2)
        f12= 1.0/36.0 *rho*(1 +3*(-ux-uz) +4.5*(-ux-uz)**2 -1.5*u2)
        f13= 1.0/36.0 *rho*(1 +3*(ux-uz)  +4.5*(ux-uz)**2  -1.5*u2)
        f14= 1.0/36.0 *rho*(1 +3*(-ux+uz) +4.5*(-ux+uz)**2 -1.5*u2)
        
        f15= 1.0/36.0 *rho*(1 +3*(uy+uz)  +4.5*(uy+uz)**2  -1.5*u2)
        f16= 1.0/36.0 *rho*(1 +3*(-uy-uz) +4.5*(-uy-uz)**2 -1.5*u2)
        f17= 1.0/36.0 *rho*(1 +3*(uy-uz)  +4.5*(uy-uz)**2  -1.5*u2)
        f18= 1.0/36.0 *rho*(1 +3*(-uy+uz) +4.5*(-uy+uz)**2 -1.5*u2)
    
        return np.array([f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13,f14,f15,f16,f17,f18])
    
    def velocities(self):
        ux2 = self.ux**2
        uy2 = self.uy**2
        uz2 = self.uz**2
        u2  =ux2+uy2+uz2
        return ux2,uy2,uz2,u2

    #Common operations
    def updateMacroVariables(self,F=[0.,0.,0.]):
        """Returns the macroscopic variables density and velocity.(rho, ux,uy)
           Equilibrium velocity includes forcing"""
        f= self.fd
        
        self.rho = np.sum(f,axis=0)
        self.rho[self.rho<1e-13]=1
        self.ux = ( (f[1]-f[2]) + (f[7]-f[8])   + (f[9] -f[10]) + (f[11]-f[12]) + (f[13]-f[14]) )/self.rho
        self.uy = ( (f[3]-f[4]) + (f[7]-f[8])   - (f[9] -f[10]) + (f[15]-f[16]) + (f[17]-f[18]) )/self.rho
        self.uz = ( (f[5]-f[6]) + (f[11]-f[12]) - (f[13]-f[14]) + (f[15]-f[16]) - (f[17]-f[18]) )/self.rho
        #Additional forcing
        self.ux  += 0.5*F[0]/self.rho
        self.uy  += 0.5*F[1]/self.rho
        self.uz  += 0.5*F[1]/self.rho

#LBM steps
    def streamming(self):
        """Stream distributions funtions along al discrete velocity directions
           keep periodic boundary conditions on all edges"""
        
        self.fd[1] = np.roll(self.fd[1], 1,axis=2)  #x >
        self.fd[2] = np.roll(self.fd[2],-1,axis=2)  #x <
        self.fd[3] = np.roll(self.fd[3], 1,axis=1)  #y v
        self.fd[4] = np.roll(self.fd[4],-1,axis=1)  #y ^
        self.fd[5] = np.roll(self.fd[5], 1,axis=0)  #z /^
        self.fd[6] = np.roll(self.fd[6],-1,axis=0)  #z v/
        
        self.fd[7] = np.roll( np.roll(self.fd[7], 1,axis=2) , 1,axis=1) #x >,y v
        self.fd[8] = np.roll( np.roll(self.fd[8],-1,axis=2) ,-1,axis=1) #x <,y ^
        self.fd[9] = np.roll( np.roll(self.fd[9], 1,axis=2) ,-1,axis=1) #x >,y ^
        self.fd[10]= np.roll( np.roll(self.fd[10],-1,axis=2), 1,axis=1) #x <,y v
        
        self.fd[11] = np.roll( np.roll(self.fd[11], 1,axis=2) , 1,axis=0) #x >,z /^
        self.fd[12] = np.roll( np.roll(self.fd[12],-1,axis=2) ,-1,axis=0) #x <,z v/
        self.fd[13] = np.roll( np.roll(self.fd[13], 1,axis=2) ,-1,axis=0) #x >,z v/
        self.fd[14] = np.roll( np.roll(self.fd[14],-1,axis=2) , 1,axis=0) #x <,z /^
        
        self.fd[15] = np.roll( np.roll(self.fd[15], 1,axis=1) , 1,axis=0) #y v,z /^
        self.fd[16] = np.roll( np.roll(self.fd[16],-1,axis=1) ,-1,axis=0) #y ^,z v/
        self.fd[17] = np.roll( np.roll(self.fd[17], 1,axis=1) ,-1,axis=0) #y v,z v/
        self.fd[18] = np.roll( np.roll(self.fd[17],-1,axis=1) , 1,axis=0) #y ^,z /^

    def collision(self):
        """Modify distribution funtions according to collision term"""
        feq = self.eqdistributions()
        self.fd -=  1.0/self.tau*(self.fd - feq)

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
