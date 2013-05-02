#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>

#define N 41
#define Lx 128
#define Ly 1
#define Lz 1
#define Dim 3
#define dx 0.001
#define radio 1000.0
#define fa 0.1

using namespace std;

const double c_s2 = 1.0 - sqrt(2.0/5.0); 
const double invc_s2 = 1.0/(1.0 - sqrt(2.0/5.0)); 
const double invcs = 1.0/sqrt(c_s2);
const double dt = dx;

int tmax,tmaxold=0; int Tini=0;

inline int delta(int a, int b){
  if(a==b)
    return 1;
  else
    return 0;
}
double metric(int a, int b, int i, int j, int k){
  double r=0; double g[Dim][Dim];
  int ik, jk, kk;

  for(ik=0;ik<Dim;ik++)
    for(jk=0;jk<Dim;jk++)
      g[ik][jk]=0;
 
  g[0][0]=2.0;
  g[1][1]=1.0;//*( i + radio )*( i + radio )*dx*dx;
  g[2][2]=1.0;
 
  return (g[a][b]);
}
double Gamma(int a, int b, int c, int i, int j, int k){
  double r=0; double ggmm[Dim][Dim][Dim];
  int ik, jk, kk;

  for(ik=0;ik<Dim;ik++)
    for(jk=0;jk<Dim;jk++)
      for(kk=0;kk<Dim;kk++)
	ggmm[ik][jk][kk]=0;

  ggmm[1][0][1] = 1.0/(dx*( i + radio ));
  ggmm[1][1][0] = 1.0/(dx*( i + radio ));
  ggmm[0][1][1] = -1.0*(dx*i + dx*radio);

  return (ggmm[a][b][c]);
}
double trazag(int i, int j, int k){
  double r=0; 

  r = metric(0,0,i,j,k) + metric(1,1,i,j,k) + metric(2,2,i,j,k);
 
  return (r);
}
double scalarg(int i, int j, int k){
  double r=0; 

  r = metric(0,0,i,j,k)*( metric(1,1,i,j,k)*metric(2,2,i,j,k) - metric(1,2,i,j,k)*metric(2,1,i,j,k)  ) - metric(0,1,i,j,k)*( metric(1,0,i,j,k)*metric(2,2,i,j,k) - metric(1,2,i,j,k)*metric(2,0,i,j,k)  ) + metric(0,2,i,j,k)*( metric(1,0,i,j,k)*metric(2,1,i,j,k) - metric(1,1,i,j,k)*metric(2,0,i,j,k)  );
 
  return (r);
}
double invmetric(int a, int b, int i, int j, int k){
  double r=0; double g[Dim][Dim];
  int ik, jk, kk;
  
  g[0][0] = metric(2,2,i,j,k)*metric(1,1,i,j,k) - metric(2,1,i,j,k)*metric(1,2,i,j,k);
  g[0][1] = -metric(2,2,i,j,k)*metric(0,1,i,j,k) + metric(2,1,i,j,k)*metric(0,2,i,j,k);
  g[0][2] = metric(1,2,i,j,k)*metric(0,1,i,j,k) - metric(1,1,i,j,k)*metric(0,2,i,j,k);

  g[1][0] = -metric(2,2,i,j,k)*metric(1,0,i,j,k) + metric(2,0,i,j,k)*metric(1,2,i,j,k);
  g[1][1] = metric(2,2,i,j,k)*metric(0,0,i,j,k) - metric(2,0,i,j,k)*metric(0,2,i,j,k);
  g[1][2] = -metric(1,2,i,j,k)*metric(0,0,i,j,k) + metric(1,0,i,j,k)*metric(0,2,i,j,k);

  g[2][0] = metric(2,1,i,j,k)*metric(1,0,i,j,k) - metric(2,0,i,j,k)*metric(1,1,i,j,k);
  g[2][1] = -metric(2,1,i,j,k)*metric(0,0,i,j,k) + metric(2,0,i,j,k)*metric(0,1,i,j,k);
  g[2][2] = metric(1,1,i,j,k)*metric(0,0,i,j,k) - metric(1,0,i,j,k)*metric(0,1,i,j,k); 
  
  return (g[a][b]/scalarg(i,j,k));

}
double trazainvg(int i, int j, int k){
  double r=0; 

  r = invmetric(0,0,i,j,k) + invmetric(1,1,i,j,k) + invmetric(2,2,i,j,k);

  return (r);
}
double densidad(int i, int j, int k){
  double g;
  
  g=1.0;//*( i + radio )*( i + radio )*dx*dx;
  
  return (g);
}
double Voxi(int i, int j, int k){
    
    return (0.005);
}
double Voyi(int i, int j, int k){
    
    return (0);
}
double Vozi(int i, int j, int k){
    
    return (0);
}
class Cplasma{
  double fi[Lx][Ly][Lz][N], fip[Lx][Ly][Lz][N]; double taoi,waoi,rho[Lx][Ly][Lz],Vm[Lx][Ly][Lz][Dim]; int v[N][Dim]; double w[N]; bool Astro[Lx][Ly][Lz];
public:
  void CalculoMacro(int t);
  double CalculoFeq(int n,int i, int j, int k,int t);
  double CalculoFeqBound(int n,int i, int j, int k,int t);
  double Getrho(int i, int j, int k){return ( rho[i][j][k] );};
  double GetVx(int i, int j, int k){return (Vm[i][j][k][0]);};
  double GetVy(int i, int j, int k){return (Vm[i][j][k][1]);};
  double GetVz(int i, int j, int k){return (Vm[i][j][k][2]);};
  double GetP(int i, int j, int k){return (rho[i][j][k]*c_s2);};
  double GetT(int i, int j, int k){return c_s2;};
  void Inicio(void);
  void Evolucion(int t);
  double Fx(int n, int i, int j, int k);
  double Fy(int n, int i, int j, int k);
  double Fz(int n, int i, int j, int k);
  double Tfuerza(int n, int i, int j, int k);
  void test_V();
  double test_w();
  bool test_macro();
};
bool Cplasma::test_macro(){
  CalculoMacro(0);
  double vin[3];
  vin[0]=Voxi(0,0,0);
  vin[1]=Voyi(0,0,0);
  vin[2]=Vozi(0,0,0);
  double diff;
  bool r=true;
  for(int i=0;i<Lx;i++){
    for(int j=0;j<Ly;j++){
      for(int k=0;k<Lz;k++){
	for(int m=0;m<3;m++){
	  diff=abs(Vm[i][j][k][m]-vin[m]);
	if (diff>5e-14){
	  cout<<diff;
	  r=false;
	  break;
	}}}}}
  return r;
}

  
void Cplasma::test_V(){
  int vtes[3];
  vtes[0]=vtes[1]=vtes[2]=0;
  for(int n=0;n<41;n++){
    for(int m=0;m<3;m++){
      vtes[m]+=v[n][m];
    }}
  cout<<"added velocities micro: ";
  for(int m=0;m<3;m++)
    cout<<vtes[m]<<" ";
  cout<<endl;
}
double Cplasma::test_w(){
  double sumw=0;
  for(int n=0;n<41;n++)
    sumw+=w[n];
  return sumw;
}
  
void Cplasma::Inicio(void){
  int n,p,h,i,j,k,a,b,c,d; double t,g, gh, ghh, ghhh;

  w[0] = 2.0*(5045.0 - 1507.0*sqrt(10))/2025.0;
  w[1] = (55.0 - 17.0*sqrt(10))/50.0;
  w[2] = (55.0 - 17.0*sqrt(10))/50.0;
  w[3] = (55.0 - 17.0*sqrt(10))/50.0;
  w[4] = (55.0 - 17.0*sqrt(10))/50.0;
  w[5] = (55.0 - 17.0*sqrt(10))/50.0;
  w[6] = (55.0 - 17.0*sqrt(10))/50.0;
  w[7] = (55.0 - 17.0*sqrt(10))/50.0;
  w[8] = (55.0 - 17.0*sqrt(10))/50.0;
  w[9] = (55.0 - 17.0*sqrt(10))/50.0;
  w[10] = (55.0 - 17.0*sqrt(10))/50.0;
  w[11] = (55.0 - 17.0*sqrt(10))/50.0;
  w[12] = (55.0 - 17.0*sqrt(10))/50.0;

  w[13] = 37.0/(5.0*sqrt(10)) - 91.0/40.0;
  w[14] = 37.0/(5.0*sqrt(10)) - 91.0/40.0;
  w[15] = 37.0/(5.0*sqrt(10)) - 91.0/40.0;
  w[16] = 37.0/(5.0*sqrt(10)) - 91.0/40.0;
  w[17] = 37.0/(5.0*sqrt(10)) - 91.0/40.0;
  w[18] = 37.0/(5.0*sqrt(10)) - 91.0/40.0;

  w[19] = (233.0*sqrt(10) - 730.0)/1600.0;
  w[20] = (233.0*sqrt(10) - 730.0)/1600.0;
  w[21] = (233.0*sqrt(10) - 730.0)/1600.0;
  w[22] = (233.0*sqrt(10) - 730.0)/1600.0;
  w[23] = (233.0*sqrt(10) - 730.0)/1600.0;
  w[24] = (233.0*sqrt(10) - 730.0)/1600.0;
  w[25] = (233.0*sqrt(10) - 730.0)/1600.0;
  w[26] = (233.0*sqrt(10) - 730.0)/1600.0;

  w[27] = (295.0 - 92.0*sqrt(10))/16200.0;
  w[28] = (295.0 - 92.0*sqrt(10))/16200.0;
  w[29] = (295.0 - 92.0*sqrt(10))/16200.0;
  w[30] = (295.0 - 92.0*sqrt(10))/16200.0;
  w[31] = (295.0 - 92.0*sqrt(10))/16200.0;
  w[32] = (295.0 - 92.0*sqrt(10))/16200.0;

  w[33] = (130.0 - 41.0*sqrt(10))/129600.0;
  w[34] = (130.0 - 41.0*sqrt(10))/129600.0;
  w[35] = (130.0 - 41.0*sqrt(10))/129600.0;
  w[36] = (130.0 - 41.0*sqrt(10))/129600.0;
  w[37] = (130.0 - 41.0*sqrt(10))/129600.0;
  w[38] = (130.0 - 41.0*sqrt(10))/129600.0;
  w[39] = (130.0 - 41.0*sqrt(10))/129600.0;
  w[40] = (130.0 - 41.0*sqrt(10))/129600.0;

  
  v[0][0] = 0;
  v[0][1] = 0;
  v[0][2] = 0;

  v[1][0] = 1;
  v[1][1] = 1;
  v[1][2] = 0;

  v[2][0] = -1;
  v[2][1] = 1;
  v[2][2] = 0;

  v[3][0] = -1;
  v[3][1] = -1;
  v[3][2] = 0;

  v[4][0] = 1;
  v[4][1] = -1;
  v[4][2] = 0;

  v[5][0] = 1;
  v[5][1] = 0;
  v[5][2] = 1;

  v[6][0] = -1;
  v[6][1] = 0;
  v[6][2] = 1;

  v[7][0] = -1;
  v[7][1] = 0;
  v[7][2] = -1;

  v[8][0] = 1;
  v[8][1] = 0;
  v[8][2] = -1;

  v[9][0] = 0;
  v[9][1] = 1;
  v[9][2] = 1;

  v[10][0] = 0;
  v[10][1] = -1;
  v[10][2] = 1;

  v[11][0] = 0;
  v[11][1] = -1;
  v[11][2] = -1;

  v[12][0] = 0;
  v[12][1] = 1;
  v[12][2] = -1;

  v[13][0] = -1;
  v[13][1] = 0;
  v[13][2] = 0;

  v[14][0] = 1;
  v[14][1] = 0;
  v[14][2] = 0;

  v[15][0] = 0;
  v[15][1] = 1;
  v[15][2] = 0;

  v[16][0] = 0;
  v[16][1] = -1;
  v[16][2] = 0;

  v[17][0] = 0;
  v[17][1] = 0;
  v[17][2] = -1;

  v[18][0] = 0;
  v[18][1] = 0;
  v[18][2] = 1;

  v[19][0] = 1;
  v[19][1] = 1;
  v[19][2] = 1;

  v[20][0] = -1;
  v[20][1] = 1;
  v[20][2] = 1;

  v[21][0] = 1;
  v[21][1] = -1;
  v[21][2] = 1;

  v[22][0] = -1;
  v[22][1] = -1;
  v[22][2] = 1;

  v[23][0] = 1;
  v[23][1] = 1;
  v[23][2] = -1;

  v[24][0] = -1;
  v[24][1] = 1;
  v[24][2] = -1;

  v[25][0] = 1;
  v[25][1] = -1;
  v[25][2] = -1;

  v[26][0] = -1;
  v[26][1] = -1;
  v[26][2] = -1;

  v[27][0] = -3;
  v[27][1] = 0;
  v[27][2] = 0;

  v[28][0] = 3;
  v[28][1] = 0;
  v[28][2] = 0;

  v[29][0] = 0;
  v[29][1] = 3;
  v[29][2] = 0;

  v[30][0] = 0;
  v[30][1] = -3;
  v[30][2] = 0;

  v[31][0] = 0;
  v[31][1] = 0;
  v[31][2] = -3;

  v[32][0] = 0;
  v[32][1] = 0;
  v[32][2] = 3;  
  
  v[33][0] = 3;
  v[33][1] = 3;
  v[33][2] = 3;

  v[34][0] = -3;
  v[34][1] = 3;
  v[34][2] = 3;

  v[35][0] = 3;
  v[35][1] = -3;
  v[35][2] = 3;

  v[36][0] = -3;
  v[36][1] = -3;
  v[36][2] = 3;

  v[37][0] = 3;
  v[37][1] = 3;
  v[37][2] = -3;

  v[38][0] = -3;
  v[38][1] = 3;
  v[38][2] = -3;

  v[39][0] = 3;
  v[39][1] = -3;
  v[39][2] = -3;

  v[40][0] = -3;
  v[40][1] = -3;
  v[40][2] = -3;
  
  for(i=0;i<Lx;i++)
    for(j=0;j<Ly;j++)
      for(k=0;k<Lz;k++){
	
	if( i<=0 || i>=Lx-1 ){
	  Astro[i][j][k]=1;
	}else{
	  Astro[i][j][k]=0;
	}
	
      }
  
    taoi=1.0;
    waoi=1.0/taoi;

    for(i=0;i<Lx;i++)
      for(j=0;j<Ly;j++)
	for(k=0;k<Lz;k++){
  
	  for(n=0;n<N;n++){

	    g=v[n][0]*Voxi(i,j,k)+v[n][1]*Voyi(i,j,k)+v[n][2]*Vozi(i,j,k);
	    gh = v[n][0]*v[n][0]*invmetric(0,0,i,j,k) + v[n][1]*v[n][0]*invmetric(1,0,i,j,k) + v[n][0]*v[n][1]*invmetric(0,1,i,j,k) + v[n][1]*v[n][1]*invmetric(1,1,i,j,k) + v[n][1]*v[n][2]*invmetric(1,2,i,j,k) + v[n][2]*v[n][1]*invmetric(2,1,i,j,k) + v[n][2]*v[n][2]*invmetric(2,2,i,j,k) + v[n][2]*v[n][0]*invmetric(2,0,i,j,k) + v[n][0]*v[n][2]*invmetric(0,2,i,j,k);
	    ghh = v[n][0]*v[n][0] + v[n][1]*v[n][1] + v[n][2]*v[n][2];
	    ghhh = v[n][0]*Voxi(i,j,k)*invmetric(0,0,i,j,k) + v[n][1]*Voxi(i,j,k)*invmetric(1,0,i,j,k) + v[n][0]*Voyi(i,j,k)*invmetric(0,1,i,j,k) + v[n][1]*Voyi(i,j,k)*invmetric(1,1,i,j,k) + v[n][1]*Vozi(i,j,k)*invmetric(1,2,i,j,k) + v[n][2]*Voyi(i,j,k)*invmetric(2,1,i,j,k) + v[n][2]*Vozi(i,j,k)*invmetric(2,2,i,j,k) + v[n][2]*Voxi(i,j,k)*invmetric(2,0,i,j,k) + v[n][0]*Vozi(i,j,k)*invmetric(0,2,i,j,k);

	    fi[i][j][k][n]=w[n]*densidad(i,j,k)*( 1 + g*invc_s2 + 0.5*gh*invc_s2 - 0.5*ghh*invc_s2 - 0.5*trazainvg(i,j,k) + 1.5 + 0.5*g*g*(invc_s2*invc_s2) - 0.5*(Voxi(i,j,k)*Voxi(i,j,k)+Voyi(i,j,k)*Voyi(i,j,k)+Vozi(i,j,k)*Vozi(i,j,k))*invc_s2  + g*g*g*invc_s2*invc_s2*invc_s2/6.0  - 0.5*g*(Voxi(i,j,k)*Voxi(i,j,k)+Voyi(i,j,k)*Voyi(i,j,k)+Vozi(i,j,k)*Vozi(i,j,k))*(invc_s2*invc_s2) + 0.5*g*(gh - ghh)*invc_s2*invc_s2 - 0.5*g*( trazainvg(i,j,k) - 3 )*invc_s2 - ghhh*invc_s2 + g*invc_s2 );
	    
	  }
	  
	}
    
}
void Cplasma::CalculoMacro(int t){
  int n, i,j,k, a, b;
  
  for(i=0;i<Lx;i++)
    for(j=0;j<Ly;j++)
      for(k=0;k<Lz;k++){
	
	Vm[i][j][k][0]=0;
	Vm[i][j][k][1]=0;
	Vm[i][j][k][2]=0;
	rho[i][j][k]=fi[i][j][k][0];
	
	for(n=1;n<N;n++){
	  rho[i][j][k]+=fi[i][j][k][n];
	  Vm[i][j][k][0]+=fi[i][j][k][n]*v[n][0];
	  Vm[i][j][k][1]+=fi[i][j][k][n]*v[n][1];
	  Vm[i][j][k][2]+=fi[i][j][k][n]*v[n][2];
	}
	
	Vm[i][j][k][0]/=rho[i][j][k];
	Vm[i][j][k][1]/=rho[i][j][k];
	Vm[i][j][k][2]/=rho[i][j][k];

 	if(Astro[i][j][k]==1){
 	  Vm[i][j][k][0]=0;
 	  Vm[i][j][k][1]=0;
 	  Vm[i][j][k][2]=0;
 	}
	
//  	Vm[i][j][k][0]=0;
//  	Vm[i][j][k][2]=0;

      }
//  
//  for(i=0;i<Lx;i++)
//    for(j=0;j<Ly;j++)
//      for(k=0;k<Lz;k++){
//
//	rho[Lx/2][j][k]=0.909175;
//	// 	rho[0][j][k]=0.993;
//	
// 	rho[1][j][k]=(4.0*rho[2][j][k]-rho[3][j][k])/3.0;
// 	rho[Lx-2][j][k]=(4.0*rho[Lx-3][j][k]-rho[Lx-4][j][k])/3.0;
//
//   	rho[0][j][k]=(4.0*rho[1][j][k]-rho[2][j][k])/3.0;
//	rho[Lx-1][j][k]=(4.0*rho[Lx-2][j][k]-rho[Lx-3][j][k])/3.0;
//
//      }

}
inline double Cplasma::CalculoFeq(int n, int i, int j ,int k,int t){
  double r=0,g, gh, ghh, ghhh, vx,vy,vz;
  
  vx=Vm[i][j][k][0];
  vy=Vm[i][j][k][1];
  vz=Vm[i][j][k][2];

  g=v[n][0]*vx+v[n][1]*vy+v[n][2]*vz;
  gh = v[n][0]*v[n][0]*invmetric(0,0,i,j,k) + v[n][1]*v[n][0]*invmetric(1,0,i,j,k) + v[n][0]*v[n][1]*invmetric(0,1,i,j,k) + v[n][1]*v[n][1]*invmetric(1,1,i,j,k) + v[n][1]*v[n][2]*invmetric(1,2,i,j,k) + v[n][2]*v[n][1]*invmetric(2,1,i,j,k) + v[n][2]*v[n][2]*invmetric(2,2,i,j,k) + v[n][2]*v[n][0]*invmetric(2,0,i,j,k) + v[n][0]*v[n][2]*invmetric(0,2,i,j,k);
  ghh = v[n][0]*v[n][0] + v[n][1]*v[n][1] + v[n][2]*v[n][2];
  ghhh = v[n][0]*vx*invmetric(0,0,i,j,k) + v[n][1]*vx*invmetric(1,0,i,j,k) + v[n][0]*vy*invmetric(0,1,i,j,k) + v[n][1]*vy*invmetric(1,1,i,j,k) + v[n][1]*vz*invmetric(1,2,i,j,k) + v[n][2]*vy*invmetric(2,1,i,j,k) + v[n][2]*vz*invmetric(2,2,i,j,k) + v[n][2]*vx*invmetric(2,0,i,j,k) + v[n][0]*vz*invmetric(0,2,i,j,k);
  
  r=w[n]*rho[i][j][k]*( 1 + 0.5*gh*invc_s2 - 0.5*trazainvg(i,j,k) + g*invc_s2 + 0.5*g*g*(invc_s2*invc_s2) - 0.5*(vx*vx+vy*vy+vz*vz)*invc_s2 + 1.5 - 0.5*ghh*invc_s2 + g*g*g*(invc_s2*invc_s2*invc_s2)/6.0  - 0.5*g*( vx*vx+vy*vy+vz*vz )*(invc_s2*invc_s2) + 0.5*g*(gh - ghh)*(invc_s2*invc_s2) - 0.5*g*( trazainvg(i,j,k) - 3 )*invc_s2 - ghhh*invc_s2 + g*invc_s2 );
  
  return r;
  
}
void Cplasma::Evolucion(int t){
  int i,j,k,ia,ja,ka,ik,jk,kk, n,nop; double g, aux1[N];

  CalculoMacro(t);
  
  for(i=0;i<Lx;i++)
    for(j=0;j<Ly;j++)
      for(k=0;k<Lz;k++){
 	for(n=0;n<N;n++){
	  fip[i][j][k][n]=fi[i][j][k][n]-waoi*(fi[i][j][k][n]-CalculoFeq(n,i,j,k,t));// + dt*Tfuerza(n, i, j, k);
 	}
      }	
  
  for(i=0;i<Lx;i++)
    for(j=0;j<Ly;j++)
      for(k=0;k<Lz;k++)
 	for(n=0;n<N;n++){
	  if(n==1) nop=3; if(n==2) nop=4; if(n==3) nop=1; if(n==4) nop=2; if(n==5) nop=7; if(n==7) nop=5; if(n==6) nop=8; if(n==8) nop=6; if(n==9) nop=11; if(n==11) nop=9; if(n==10) nop=12; if(n==12) nop=10;if(n==13) nop=14; if(n==14) nop=13; if(n==15) nop=16; if(n==16) nop=15; if(n==17) nop=18; if(n==18) nop=17; if(n==0) nop=0; if(n==19) nop=26; if(n==26) nop=19; if(n==20) nop=25; if(n==25) nop=20; if(n==24) nop=21; if(n==21) nop=24; if(n==22) nop=23; if(n==23) nop=22; if(n==27) nop=28; if(n==28) nop=27; if(n==29) nop=30; if(n==30) nop=29;  if(n==31) nop=32; if(n==32) nop=31; if(n==33) nop=40; if(n==40) nop=33; if(n==34) nop=39; if(n==39) nop=34; if(n==35) nop=38; if(n==38) nop=35;  if(n==36) nop=37; if(n==37) nop=36;
	  
 	  ia=(i+v[n][0]+Lx)%Lx; 
 	  ja=(j+v[n][1]+Ly)%Ly;
 	  ka=(k+v[n][2]+Lz)%Lz;

    	  ia=i+v[n][0]; 

//  	  if( ia > Lx-1 || ia < 0 )
//  	    ia = i;
	
//	  if(ia>=0 && ia<=Lx-1)
//	    fi[i][j][k][nop]=fip[ia][ja][ka][nop];
//	  else
//	    fi[i][j][k][nop]=CalculoFeq(nop,i,j,k,t);
        if(ia>=0 && ia<=Lx-1)
          fi[i][j][k][nop]=fip[ia][ja][ka][nop];
        else
          fi[i][j][k][nop]=fip[i][j][k][n];


	}

  
//   for(i=0;i<Lx;i++)
//     for(j=0;j<Ly;j++)
//       for(k=0;k<Lz;k++){
	
//   	for(n=0;n<N;n++)
//   	  aux1[n] = fi[i][j][k][n];
	
//   	for(n=0;n<N;n++){
//   	  if(n==1) nop=3; if(n==2) nop=4; if(n==3) nop=1; if(n==4) nop=2; if(n==5) nop=7; if(n==7) nop=5; if(n==6) nop=8; if(n==8) nop=6; if(n==9) nop=11; if(n==11) nop=9; if(n==10) nop=12; if(n==12) nop=10;if(n==13) nop=14; if(n==14) nop=13; if(n==15) nop=16; if(n==16) nop=15; if(n==17) nop=18; if(n==18) nop=17; if(n==0) nop=0; if(n==19) nop=26; if(n==26) nop=19; if(n==20) nop=25; if(n==25) nop=20; if(n==24) nop=21; if(n==21) nop=24; if(n==22) nop=23; if(n==23) nop=22; if(n==27) nop=28; if(n==28) nop=27; if(n==29) nop=30; if(n==30) nop=29;  if(n==31) nop=32; if(n==32) nop=31; if(n==33) nop=40; if(n==40) nop=33; if(n==34) nop=39; if(n==39) nop=34; if(n==35) nop=38; if(n==38) nop=35;  if(n==36) nop=37; if(n==37) nop=36;
	  
//    	  ia=(i+v[n][0]+Lx)%Lx; 
//  	  ja=(j+v[n][1]+Ly)%Ly;
//   	  ka=(k+v[n][2]+Lz)%Lz;
	  
//   	  if( Astro[i][j][k]==1 && Astro[ia][ja][ka]==0  )
//   	    fi[i][j][k][n]=aux1[nop];
  	  
//   	}
	
//       }
    
}
inline double Cplasma::Fx(int n, int i, int j, int k){
  double r;
  int ik, jk;
  
  r=0;
  for(ik=0;ik<Dim;ik++)
    for(jk=0;jk<Dim;jk++)
      r += -Gamma(0,ik,jk,i,j,k)*v[n][ik]*v[n][jk];
  
  return r;
}
inline double Cplasma::Fy(int n, int i, int j, int k){
  double r;
  int ik, jk;
  
  r=0;
  for(ik=0;ik<Dim;ik++)
    for(jk=0;jk<Dim;jk++)
      r += -Gamma(1,ik,jk,i,j,k)*v[n][ik]*v[n][jk];
 
  return r + fa/rho[i][j][k];

}
inline double Cplasma::Fz(int n, int i, int j, int k){
  double r;
  int ik, jk;
  
  r=0;
  for(ik=0;ik<Dim;ik++)
    for(jk=0;jk<Dim;jk++)
      r += -Gamma(2,ik,jk,i,j,k)*v[n][ik]*v[n][jk];
  
  return r;
}
double Cplasma::Tfuerza(int n, int i, int j, int k){
  double r=0,g,gh,ghh, ghhh, ghhhh, ghhhhh, vx,vy,vz, thirdorder=0, faux[Dim];
  int a,b,c,d;
  double aH[Dim][Dim][Dim];
  double aHe[Dim][Dim][Dim];
  
  vx=Vm[i][j][k][0];
  vy=Vm[i][j][k][1];
  vz=Vm[i][j][k][2];


  g=v[n][0]*vx+v[n][1]*vy+v[n][2]*vz;
  gh=v[n][0]*Fx(n,i,j,k)+v[n][1]*Fy(n,i,j,k)+v[n][2]*Fz(n,i,j,k);
  ghh=vx*Fx(n,i,j,k)+vy*Fy(n,i,j,k)+vz*Fz(n,i,j,k);  
  ghhh = v[n][0]*v[n][0]*invmetric(0,0,i,j,k) + v[n][1]*v[n][0]*invmetric(1,0,i,j,k) + v[n][0]*v[n][1]*invmetric(0,1,i,j,k) + v[n][1]*v[n][1]*invmetric(1,1,i,j,k) + v[n][1]*v[n][2]*invmetric(1,2,i,j,k) + v[n][2]*v[n][1]*invmetric(2,1,i,j,k) + v[n][2]*v[n][2]*invmetric(2,2,i,j,k) + v[n][2]*v[n][0]*invmetric(2,0,i,j,k) + v[n][0]*v[n][2]*invmetric(0,2,i,j,k);
  ghhhh = v[n][0]*v[n][0] + v[n][1]*v[n][1] + v[n][2]*v[n][2];
  ghhhhh = v[n][0]*Fx(n,i,j,k)*invmetric(0,0,i,j,k) + v[n][1]*Fx(n,i,j,k)*invmetric(1,0,i,j,k) + v[n][0]*Fy(n,i,j,k)*invmetric(0,1,i,j,k) + v[n][1]*Fy(n,i,j,k)*invmetric(1,1,i,j,k) + v[n][1]*Fz(n,i,j,k)*invmetric(1,2,i,j,k) + v[n][2]*Fy(n,i,j,k)*invmetric(2,1,i,j,k) + v[n][2]*Fz(n,i,j,k)*invmetric(2,2,i,j,k) + v[n][2]*Fx(n,i,j,k)*invmetric(2,0,i,j,k) + v[n][0]*Fz(n,i,j,k)*invmetric(0,2,i,j,k);

  
  //Third order in Hermite for the force term
//   faux[0]=Fx(n,i,j,k)*invcs;
//   faux[1]=Fy(n,i,j,k)*invcs;
//   faux[2]=Fz(n,i,j,k)*invcs;

//   for(a=0;a<Dim;a++)
//     for(b=0;b<Dim;b++)
//       for(c=0;c<Dim;c++){
// 	aH[a][b][c] = ( invmetric(a,b,i,j,k)  -  delta(a,b) )*Vm[i][j][k][c]*invcs + ( invmetric(b,c,i,j,k)  -  delta(b,c) )*Vm[i][j][k][a]*invcs + ( invmetric(a,c,i,j,k)  -  delta(a,c) )*Vm[i][j][k][b]*invcs + Vm[i][j][k][a]*Vm[i][j][k][b]*Vm[i][j][k][c]*invcs*invcs*invcs; 
// 	aHe[a][b][c] = 0;
// 	for(d=0;d<Dim;d++)
// 	  aHe[a][b][c] += ( (invcs*invcs*invcs*invcs)*v[n][a]*v[n][b]*v[n][c]*v[n][d] - invcs*invcs*( v[n][a]*v[n][b]*delta(c,d) + v[n][b]*v[n][c]*delta(d,a) + v[n][c]*v[n][d]*delta(a,b) + v[n][d]*v[n][a]*delta(b,c) + v[n][d]*v[n][b]*delta(a,c) + v[n][a]*v[n][c]*delta(b,d)   ) + ( delta(a,b)*delta(c,d) + delta(a,c)*delta(b,d) +  delta(a,d)*delta(b,c)  ) )*faux[d]; 
//       }
  
//   thirdorder=0;

//   for(a=0;a<Dim;a++)
//     for(b=0;b<Dim;b++)
//       for(c=0;c<Dim;c++)
// 	thirdorder+=aH[a][b][c]*aHe[a][b][c];
  
//   thirdorder*=w[n]*rho[i][j][k]/6.0;
  //-----------------------------------------
  
  r=w[n]*rho[i][j][k]*(  gh/c_s2 + g*gh/(c_s2*c_s2) - ghh/c_s2 + 0.5*ghhh*gh/(c_s2*c_s2) - ghhhhh/c_s2 - 0.5*gh*trazainvg(i,j,k)/c_s2 - 0.5*ghhhh*gh/(c_s2*c_s2)  + 2.5*gh/c_s2  + 0.5*gh*g*g/(c_s2*c_s2*c_s2) - g*ghh/(c_s2*c_s2) - 0.5*(vx*vx+vy*vy+vz*vz)*gh/(c_s2*c_s2)  ) + thirdorder;
   
  return r;

}
int main(){
  int a,t,i,j,k,n;
  Cplasma plasma;
  char filename[4][20];
  system("rm *.dat");    
  plasma.Inicio();
  FILE *X[4];
  
  cout<<"Tests\n";
  plasma.test_V();
  cout<<"\n Sum W "<<plasma.test_w();
  cout<<"recreate macro";
  if (plasma.test_macro())
    cout<<"sucess";
  else
    cout<<"fail";
  do{
    
    //	cout << "Introduzca tiempo de parada:";
    cin >> tmax;
    
    for(t=tmaxold;t<tmax;t++){
      plasma.Evolucion(t);
      //      cout << t+1 << "\t" << endl;
     
      if(t%500==0){
	
	sprintf(filename[0], "rho_%d.dat",t + 1000000);
	sprintf(filename[1], "vx_%d.dat",t + 1000000);
	sprintf(filename[2], "vy_%d.dat",t + 1000000);
	sprintf(filename[3], "vz_%d.dat",t + 1000000);
	
	for(a=0;a<4;a++)
	  X[a]=fopen(filename[a],"w");
	
	k=0;
	
	for(i=0;i<Lx;i++){
	  for(j=0;j<Ly;j++){
	    fprintf(X[0],"%.7e \t",plasma.Getrho(i,j,k));
	    fprintf(X[1],"%.7e \t",plasma.GetVx(i,j,k)*plasma.Getrho(i,j,k));
	    fprintf(X[2],"%.7e \t",plasma.GetVy(i,j,k)*plasma.Getrho(i,j,k));
	    fprintf(X[3],"%.7e \t",plasma.GetVz(i,j,k)*plasma.Getrho(i,j,k));
	    
	  }
	  
	  for(a=0;a<4;a++)
	    fprintf(X[a],"\n");
	  
	}
	
	for(a=0;a<4;a++)
	  fclose(X[a]);
	
	cout << "set view 0,0,1; set xlabel '" << t+1000000 << "' ; set size ratio -1; set pm3d; unset surface; plot [][] 'rho_" << t+1000000 << ".dat' w l;" << endl; 
	
	
      }
    }
    
    tmaxold=tmax;
    
  }while(1);
  
  return 0;
}


    
    
