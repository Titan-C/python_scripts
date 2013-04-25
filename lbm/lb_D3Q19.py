# -*- coding: utf-8 -*-
"""
Created on Tue Apr 23 13:45:43 2013

@author: Oscar Najera
Lattice Boltzmann 3D en D3Q19 class
"""
from numpy import array, sum, ones, roll, where

class lattice:
    def __init__(self,LatticeSize,U,tau,r=1):
        """Unitiate dicrete distribution functions for a given LatticeSize
           and velocity"""
        self.ux = U[0]*ones(LatticeSize)
        self.uy = U[1]*ones(LatticeSize)
        self.uz = U[2]*ones(LatticeSize)
        self.rho = r*ones(LatticeSize)
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
    
        return array([f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13,f14,f15,f16,f17,f18])
    
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
        
        self.rho = sum(f,axis=0)
        self.rho[self.rho<1e-13]=1
        self.ux = ( (f[1]-f[2]) + (f[7]-f[8])   + (f[9] -f[10]) + (f[11]-f[12]) + (f[13]-f[14]) )/self.rho
        self.uy = ( (f[3]-f[4]) + (f[7]-f[8])   - (f[9] -f[10]) + (f[15]-f[16]) + (f[17]-f[18]) )/self.rho
        self.uz = ( (f[5]-f[6]) + (f[11]-f[12]) - (f[13]-f[14]) + (f[15]-f[16]) - (f[17]-f[18]) )/self.rho
        #Additional forcing
        self.ux  += 0.5*F[0]/self.rho
        self.uy  += 0.5*F[1]/self.rho
        self.uz  += 0.5*F[2]/self.rho

#LBM steps
    def streamming(self):
        """Stream distributions funtions along al discrete velocity directions
           keep periodic boundary conditions on all edges"""
        
        self.fd[1] = roll(self.fd[1], 1,axis=2)  #x >
        self.fd[2] = roll(self.fd[2],-1,axis=2)  #x <
        self.fd[3] = roll(self.fd[3], 1,axis=1)  #y v
        self.fd[4] = roll(self.fd[4],-1,axis=1)  #y ^
        self.fd[5] = roll(self.fd[5], 1,axis=0)  #z /^
        self.fd[6] = roll(self.fd[6],-1,axis=0)  #z v/
        
        self.fd[7] = roll( roll(self.fd[7], 1,axis=2) , 1,axis=1) #x >,y v
        self.fd[8] = roll( roll(self.fd[8],-1,axis=2) ,-1,axis=1) #x <,y ^
        self.fd[9] = roll( roll(self.fd[9], 1,axis=2) ,-1,axis=1) #x >,y ^
        self.fd[10]= roll( roll(self.fd[10],-1,axis=2), 1,axis=1) #x <,y v
        
        self.fd[11] = roll( roll(self.fd[11], 1,axis=2) , 1,axis=0) #x >,z /^
        self.fd[12] = roll( roll(self.fd[12],-1,axis=2) ,-1,axis=0) #x <,z v/
        self.fd[13] = roll( roll(self.fd[13], 1,axis=2) ,-1,axis=0) #x >,z v/
        self.fd[14] = roll( roll(self.fd[14],-1,axis=2) , 1,axis=0) #x <,z /^
        
        self.fd[15] = roll( roll(self.fd[15], 1,axis=1) , 1,axis=0) #y v,z /^
        self.fd[16] = roll( roll(self.fd[16],-1,axis=1) ,-1,axis=0) #y ^,z v/
        self.fd[17] = roll( roll(self.fd[17], 1,axis=1) ,-1,axis=0) #y v,z v/
        self.fd[18] = roll( roll(self.fd[18],-1,axis=1) , 1,axis=0) #y ^,z /^

    def collision(self):
        """Modify distribution funtions according to collision term"""
        feq = self.eqdistributions()
        self.fd -=  1.0/self.tau*(self.fd - feq)

    def force(self,F):
        """Adds microscopic forcing term"""
        ux,uy,uz = self.ux, self.uy, self.uz
        cf = (1-0.5/self.tau)*3.0

        F0 = 1./3.  *(   -ux *F[0] -    uy *F[1]  -    uz *F[2] )

        F1 = 1./18. *( (1-ux)*F[0] -    uy *F[1]  -    uz *F[2] + 3*ux*F[0])
        F2 = 1./18. *((-1-ux)*F[0] -    uy *F[1]  -    uz *F[2] + 3*ux*F[0])
        F3 = 1./18. *(   -ux *F[0] - (1-uy)*F[1]  -    uz *F[2] + 3*uy*F[0])
        F4 = 1./18. *(   -ux *F[0] -(-1-uy)*F[1]  -    uz *F[2] + 3*uy*F[0])
        F5 = 1./18. *(   -ux *F[0] -    uy *F[1]  - (1-uz)*F[2] + 3*uz*F[0])
        F6 = 1./18. *(   -ux *F[0] -    uy *F[1]  -(-1-uz)*F[2] + 3*uz*F[0])

        F7 = 1./36. *( (1-ux)*F[0] - (1-uy)*F[1]  -    uz *F[2] + 3*(ux+uy)*(F[0]+F[1]) )
        F8 = 1./36. *((-1-ux)*F[0] -(-1-uy)*F[1]  -    uz *F[2] + 3*(ux+uy)*(F[0]+F[1]) )
        F9 = 1./36. *( (1-ux)*F[0] -(-1-uy)*F[1]  -    uz *F[2] + 3*(ux-uy)*(F[0]-F[1]) )
        F10= 1./36. *((-1-ux)*F[0] - (1-uy)*F[1]  -    uz *F[2] + 3*(ux-uy)*(F[0]-F[1]) )

        F11= 1./36. *( (1-ux)*F[0] -    uy *F[1]  - (1-uz)*F[2] + 3*(ux+uz)*(F[0]+F[2]) )
        F12= 1./36. *((-1-ux)*F[0] -    uy *F[1]  -(-1-uz)*F[2] + 3*(ux+uz)*(F[0]+F[2]) )
        F13= 1./36. *( (1-ux)*F[0] -    uy *F[1]  -(-1-uz)*F[2] + 3*(ux-uz)*(F[0]-F[2]) )
        F14= 1./36. *((-1-ux)*F[0] -    uy *F[1]  - (1-uz)*F[2] + 3*(ux-uz)*(F[0]-F[2]) )

        F15= 1./36. *(   -ux *F[0] - (1-uy)*F[1]  - (1-uz)*F[2] + 3*(uy+uz)*(F[1]+F[2]) )
        F16= 1./36. *(   -ux *F[0] -(-1-uy)*F[1]  -(-1-uz)*F[2] + 3*(uy+uz)*(F[1]+F[2]) )
        F17= 1./36. *(   -ux *F[0] - (1-uy)*F[1]  -(-1-uz)*F[2] + 3*(uy-uz)*(F[1]-F[2]) )
        F18= 1./36. *(   -ux *F[0] -(-1-uy)*F[1]  - (1-uz)*F[2] + 3*(uy-uz)*(F[1]-F[2]) )

        self.fd+= cf*array([F0,F1,F2,F3,F4,F5,F6,F7,F8,F9,F10,F11,F12,F13,F14,F15,F16,F17,F18])

##boundaries handling
    def setWalls(self,walls):
        """Anulates fd where a wall(solid) is declared"""
        self.fd[:,walls==1]=0

    def onWallBBNS_BC(self,walls):
        """Bounceback no-slip boundary conditions for arbitrary walls"""
        Nz,Ny,Nx=walls.shape
        so=array(where(walls==1)) #so stands for Solid Object
        #generate direction vectors, mantain periodic boundary conditions
        xp=so[2]+1
        xp[xp==Nx]=0
        xm=so[2]-1
        xm[xm==-1]+=Nx
        
        yp=so[1]+1
        yp[yp==Ny]=0
        ym=so[1]-1
        ym[ym==-1]+=Ny
        
        zp=so[0]+1
        zp[zp==Nz]=0
        zm=so[0]-1
        zm[zm==-1]+=Nz

        #Collisions: Mayor Axis
        self.fd[1,so[0],so[1],xp   ] = self.fd[2,so[0],so[1],so[2]]
        self.fd[2,so[0],so[1],xm   ] = self.fd[1,so[0],so[1],so[2]]
        self.fd[3,so[0],yp   ,so[2]] = self.fd[4,so[0],so[1],so[2]]
        self.fd[4,so[0],ym   ,so[2]] = self.fd[3,so[0],so[1],so[2]]
        self.fd[5,zp   ,so[1],so[2]] = self.fd[6,so[0],so[1],so[2]]
        self.fd[6,zm   ,so[1],so[2]] = self.fd[5,so[0],so[1],so[2]]
        #Collisions: Diagonal Axis
        self.fd[7, so[0],yp   ,xp   ] = self.fd[8,so[0],so[1],so[2]]
        self.fd[8, so[0],ym   ,xm   ] = self.fd[7,so[0],so[1],so[2]]
        self.fd[9, so[0],ym   ,xp   ] = self.fd[10,so[0],so[1],so[2]]
        self.fd[10,so[0],yp   ,xm   ] = self.fd[9,so[0],so[1],so[2]]
        
        self.fd[11,zp   ,so[1],xp   ] = self.fd[12,so[0],so[1],so[2]]
        self.fd[12,zm   ,so[1],xm   ] = self.fd[11,so[0],so[1],so[2]]
        self.fd[13,zm   ,so[1],xp   ] = self.fd[14,so[0],so[1],so[2]]
        self.fd[14,zp   ,so[1],xm   ] = self.fd[13,so[0],so[1],so[2]]
        
        self.fd[15,zp   ,yp   ,so[2]] = self.fd[16,so[0],so[1],so[2]]
        self.fd[16,zm   ,ym   ,so[2]] = self.fd[15,so[0],so[1],so[2]]
        self.fd[17,zm   ,yp   ,so[2]] = self.fd[18,so[0],so[1],so[2]]
        self.fd[18,zp   ,ym   ,so[2]] = self.fd[17,so[0],so[1],so[2]]

        #clean solid
        self.fd[:,so[0],so[1],so[2]]=0
