# -*- coding: utf-8 -*-
"""
Created on Tue Apr 23 13:45:43 2013

@author: Oscar Najera
Lattice Boltzmann 3D en D3Q19 class
"""
from numpy import array, sum, ones,zeros, roll, where,outer,tensordot,dot

class lattice:
    def __init__(self,LatticeSize,U,tau,r=1):
        """Unitiate dicrete distribution functions for a given LatticeSize
           and velocity"""
        self.Nz,self.Ny,self.Nx=LatticeSize
        self.U = array(U*self.Nz*self.Ny*self.Nx).reshape(self.Nz,self.Ny,self.Nx,3)
        self.E = array([[0,0,0],[1,0,0],[-1,0,0],[0,1,0],[0,-1,0],[0,0,1],[0,0,-1],
                        [1,1,0],[-1,-1,0],[1,-1,0],[-1,1,0],
                        [1,0,1],[-1,0,-1],[1,0,-1],[-1,0,1],
                        [0,1,1],[0,-1,-1],[0,1,-1],[0,-1,1]])
        self.w=array([1./3., 1./18.,1./18.,1./18.,1./18.,1./18.,1./18.,
                 1./36.,1./36.,1./36.,1./36.,1./36.,1./36.,
                 1./36.,1./36.,1./36.,1./36.,1./36.,1./36.])
        self.rho = r*ones(LatticeSize)
        
        self.fd  = self.eqdistributions()
        self.tau = tau
        self.nu  = (tau-0.5)/3.

    def eqdistributions(self):
        """Evaluate Equilibrium distribution functions"""

        f=self.w.reshape(19,self.Nz,self.Ny,self.Nx)*self.rho
        eu=tensordot(self.E,self.U,axes=(1,3))
        u2=sum(self.U**2,axis=3)
        f*=(1+3*eu+4.5*eu**2-1.5*u2)

        return f

    #Common operations
    def updateMacroVariables(self,F=zeros(3)):
        """Returns the macroscopic variables density and velocity.(rho, ux,uy)
           Equilibrium velocity includes forcing"""
        
        self.rho = sum(self.fd,axis=0)
        self.rho[self.rho<1e-13]=1
        self.U = tensordot(self.fd,self.E,axes=(0,0))/self.rho.reshape(self.Nz,self.Ny,self.Nx,1)
        #Additional forcing
        self.U += 0.5*F/self.rho.reshape(self.Nz,self.Ny,self.Nx,1)

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
        cf = (1-0.5/self.tau)*3.0*self.w.reshape(19,1,1,1)
        self.fd+= cf*( dot(self.E.reshape(19,1,1,1,3)-self.U, F) + 3 * tensordot(self.E,self.U,axes=(1,3)) * dot(self.E,F).reshape(19,1,1,1))

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
