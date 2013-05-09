# -*- coding: utf-8 -*-
"""
Created on Tue Apr 23 13:45:43 2013

@author: Oscar Najera
Lattice Boltzmann 3D en D3Q19 class
"""
from numpy import array, sum, ones,zeros, roll, where,tensordot,dot,sqrt,eye,trace,diag
from numpy.linalg import inv

class lattice:
    def __init__(self,LatticeSize,U,tau,d=[1.,1.,1.],r=1):
        """Unitiate dicrete distribution functions for a given LatticeSize
           and velocity"""
        self.Nz,self.Ny,self.Nx=LatticeSize
        self.U = array([U[0],U[1],U[2]]*self.Nz*self.Ny*self.Nx).reshape(self.Nz,self.Ny,self.Nx,3)
        self.E = array([[0,0,0],[1,0,0],[-1,0,0],[0,1,0],[0,-1,0],[0,0,1],[0,0,-1],
                        [1,1,0],[-1,-1,0],[1,-1,0],[-1,1,0],
                        [1,0,1],[-1,0,-1],[1,0,-1],[-1,0,1],
                        [0,1,1],[0,-1,-1],[0,1,-1],[0,-1,1],
                        [1,1,1],[-1,-1,-1],[1,1,-1],[-1,-1,1],
                        [1,-1,1],[-1,1,-1],[1,-1,-1],[-1,1,1],
                        [3,0,0],[-3,0,0],[0,3,0],[0,-3,0],[0,0,3],[0,0,-3],
                        [3,3,3],[-3,-3,-3],[3,3,-3],[-3,-3,3],
                        [3,-3,3],[-3,3,-3],[3,-3,-3],[-3,3,3] ],dtype=float)
        metric=[inv(diag(d))]
        self.g = array([metric]*self.Nz*self.Ny*self.Nx).reshape(self.Nz,self.Ny,self.Nx,3,3)
        w0=[2.*(5045-1507*sqrt(10))/2025.]
        w1_6=6*[37./5/sqrt(10)-91./40.]
        w7_18=12*[(55-17*sqrt(10))/50.]
        w19_26=8*[(233*sqrt(10)-730)/1600.]
        w27_32=6*[(295-92*sqrt(10))/16200.]
        w33_40=8*[(130-41*sqrt(10))/129600.]
        self.w=array(w0+w1_6+w7_18+w19_26+w27_32+w33_40)
        self.rho = r*ones(LatticeSize)
        
        self.wp=self.w.reshape(41,1,1,1)*self.rho
        self.ge=tensordot(self.E,self.g,axes=(1,4))
        self.ege=sum(self.E.reshape(41,1,1,1,3)*self.ge,axis=4)
        self.e2=sum(self.E**2,axis=1).reshape(41,1,1,1)
        self.gii=trace(self.g,axis1=3,axis2=4)
        self.cs2=1-sqrt(2./5.)
        self.fd  = self.eqdistributions()
        
        self.tau = tau
        self.nu  = (tau-0.5)*self.cs2

    def eqdistributions(self):
        """Evaluate Equilibrium distribution functions"""
        cs2=self.cs2
        eu=tensordot(self.E,self.U,axes=(1,3))
        u2=sum(self.U**2,axis=3)
        uge=sum(self.U*self.ge,axis=4)
#        import pdb;pdb.set_trace()

        f =self.wp*(5./2.*cs2 + 2*eu + 0.5*self.ege - 0.5*self.e2
        + 0.5*eu**2/cs2 - 0.5*self.gii*cs2 - 0.5*u2 + eu**3/cs2**2/6.
        - 0.5*eu*u2/cs2 + 0.5*eu*(self.ege-self.e2)/cs2
        - 0.5*eu*(self.gii-3) - uge)/cs2

        return f

    #Common operations
    def updateMacroVariables(self,F=zeros(3)):
        """Returns the macroscopic variables density and velocity.(rho, ux,uy)
           Equilibrium velocity includes forcing"""
        
        self.rho = sum(self.fd,axis=0)
        self.rho[self.rho<1e-13]=1
        self.U = tensordot(self.fd,self.E,axes=(0,0))/self.rho.reshape(self.Nz,self.Ny,self.Nx,1)
        #Additional forcing
#        self.U += 0.5*F/self.rho.reshape(self.Nz,self.Ny,self.Nx,1)

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

        self.fd[19] = roll( roll( roll(self.fd[19], 1,axis=2) , 1,axis=1) , 1,axis=0) #x >,y v,z /^
        self.fd[20] = roll( roll( roll(self.fd[20],-1,axis=2) ,-1,axis=1) ,-1,axis=0) #x <,y ^,z v/
        self.fd[21] = roll( roll( roll(self.fd[21], 1,axis=2) , 1,axis=1) ,-1,axis=0) #x >,y v,z v/
        self.fd[22] = roll( roll( roll(self.fd[22],-1,axis=2) ,-1,axis=1) , 1,axis=0) #x <,y ^,z /^

        self.fd[23] = roll( roll( roll(self.fd[23], 1,axis=2) ,-1,axis=1) , 1,axis=0) #x >,y ^,z /^
        self.fd[24] = roll( roll( roll(self.fd[24],-1,axis=2) , 1,axis=1) ,-1,axis=0) #x <,y v,z v/
        self.fd[25] = roll( roll( roll(self.fd[25], 1,axis=2) ,-1,axis=1) ,-1,axis=0) #x >,y ^,z v/
        self.fd[26] = roll( roll( roll(self.fd[26],-1,axis=2) , 1,axis=1) , 1,axis=0) #x <,y v,z /^

        self.fd[27] = roll(self.fd[27], 3,axis=2)  #3x >
        self.fd[28] = roll(self.fd[28],-3,axis=2)  #3x <
        self.fd[29] = roll(self.fd[29], 3,axis=1)  #3y v
        self.fd[30] = roll(self.fd[30],-3,axis=1)  #3y ^
        self.fd[31] = roll(self.fd[31], 3,axis=0)  #3z /^
        self.fd[32] = roll(self.fd[32],-3,axis=0)  #3z v/

        self.fd[33] = roll( roll( roll(self.fd[33], 3,axis=2) , 3,axis=1) , 3,axis=0) #3x >,3y v,3z /^
        self.fd[34] = roll( roll( roll(self.fd[34],-3,axis=2) ,-3,axis=1) ,-3,axis=0) #3x <,3y ^,3z v/
        self.fd[35] = roll( roll( roll(self.fd[35], 3,axis=2) , 3,axis=1) ,-3,axis=0) #3x >,3y v,3z v/
        self.fd[36] = roll( roll( roll(self.fd[36],-3,axis=2) ,-3,axis=1) , 3,axis=0) #3x <,3y ^,3z /^

        self.fd[37] = roll( roll( roll(self.fd[37], 3,axis=2) ,-3,axis=1) , 3,axis=0) #3x >,3y ^,3z /^
        self.fd[38] = roll( roll( roll(self.fd[38],-3,axis=2) , 3,axis=1) ,-3,axis=0) #3x <,3y v,3z v/
        self.fd[39] = roll( roll( roll(self.fd[39], 3,axis=2) ,-3,axis=1) ,-3,axis=0) #3x >,3y ^,3z v/
        self.fd[40] = roll( roll( roll(self.fd[40],-3,axis=2) , 3,axis=1) , 3,axis=0) #3x <,3y v,3z /^



    def collision(self):
        """Modify distribution funtions according to collision term"""
        feq = self.eqdistributions()
        self.fd -=  1.0/self.tau*(self.fd - feq)

    def force(self,Fext):
        """Adds microscopic forcing term"""

        F=Fext.reshape(1,1,1,1,3)
        cs2=self.cs2
        eF=sum(self.E.reshape(41,1,1,1,3)*F,axis=4)
        eu=tensordot(self.E,self.U,axes=(1,3))
        uF=sum(self.U*F,axis=4)
        geF=sum(self.ge*F,axis=4)
        u2=sum(self.U**2,axis=3)

        self.fd+= self.wp*( 3.5*eF + eu*eF/cs2 - uF
            + 0.5*self.ege*eF/cs2 - geF - 0.5*self.gii*eF
            - 0.5*self.e2*eF/cs2
            + 0.5*eu*eu*eF/cs2**2 - eu*eF/cs2 - 0.5*u2*eF/cs2)/cs2

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
        x3p=so[2]+3
        x3p[x3p>=Nx]-=Nx
        xm=so[2]-1
        xm[xm==-1]+=Nx
        x3m=so[2]-3
        x3m[x3m<0]+=Nx
        
        yp=so[1]+1
        yp[yp==Ny]=0
        y3p=so[1]+3
        y3p[y3p>=Ny]-=Ny
        ym=so[1]-1
        ym[ym==-1]+=Ny
        y3m=so[1]-3
        y3m[y3m<0]+=Ny
        
        zp=so[0]+1
        zp[zp==Nz]=0
        z3p=so[0]+3
        z3p[z3p>=Nz]-=Nz
        zm=so[0]-1
        zm[zm==-1]+=Nz
        z3m=so[0]-3
        z3m[z3m<0]+=Nz

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
        #Collision: Cube vertices
        self.fd[19,zp   ,yp   ,xp   ] = self.fd[20,so[0],so[1],so[2]]
        self.fd[20,zm   ,ym   ,xm   ] = self.fd[19,so[0],so[1],so[2]]
        self.fd[21,zm   ,yp   ,xp   ] = self.fd[22,so[0],so[1],so[2]]
        self.fd[22,zp   ,ym   ,xm   ] = self.fd[21,so[0],so[1],so[2]]

        self.fd[23,zp   ,ym   ,xp   ] = self.fd[24,so[0],so[1],so[2]]
        self.fd[24,zm   ,yp   ,xm   ] = self.fd[23,so[0],so[1],so[2]]
        self.fd[25,zm   ,ym   ,xp   ] = self.fd[26,so[0],so[1],so[2]]
        self.fd[26,zp   ,yp   ,xm   ] = self.fd[25,so[0],so[1],so[2]]
        #Collision: Skip 3 cube faces
        self.fd[27,so[0],so[1],x3p   ] = self.fd[28,so[0],so[1],so[2]]
        self.fd[28,so[0],so[1],x3m   ] = self.fd[27,so[0],so[1],so[2]]
        self.fd[29,so[0],y3p   ,so[2]] = self.fd[30,so[0],so[1],so[2]]
        self.fd[30,so[0],y3m   ,so[2]] = self.fd[29,so[0],so[1],so[2]]
        self.fd[31,z3p   ,so[1],so[2]] = self.fd[32,so[0],so[1],so[2]]
        self.fd[32,z3m   ,so[1],so[2]] = self.fd[31,so[0],so[1],so[2]]
        #Collision: Skip 3 cube vertices
        self.fd[33,z3p   ,y3p   ,x3p   ] = self.fd[34,so[0],so[1],so[2]]
        self.fd[34,z3m   ,y3m   ,x3m   ] = self.fd[33,so[0],so[1],so[2]]
        self.fd[35,z3m   ,y3p   ,x3p   ] = self.fd[36,so[0],so[1],so[2]]
        self.fd[36,z3p   ,y3m   ,x3m   ] = self.fd[35,so[0],so[1],so[2]]

        self.fd[37,z3p   ,y3m   ,x3p   ] = self.fd[38,so[0],so[1],so[2]]
        self.fd[38,z3m   ,y3p   ,x3m   ] = self.fd[37,so[0],so[1],so[2]]
        self.fd[39,z3m   ,y3m   ,x3p   ] = self.fd[40,so[0],so[1],so[2]]
        self.fd[40,z3p   ,y3p   ,x3m   ] = self.fd[39,so[0],so[1],so[2]]

        #clean solid
        self.fd[:,so[0],so[1],so[2]]=0
