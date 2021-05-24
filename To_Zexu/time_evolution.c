#include "global_variables.h"
#include "time_evolution.h"
#include "pfr_cufft1.h"
#include "misc.h"
#include "inhom_v3.h"

double  cr,nr,mur,hr;
double g2;
cufftDoubleComplex ck,nk,muk;
float randnum[L][M][N];
void GetEintC(int timestep) //this function interpret the governing equation
{
  int g_i,g_j,g_k,gijk;
  int x_i,x_j,x_k,xijk;
  int l;
  int i;
  double n;
  char filename[500];
  //Evolving composition c(k)
  
//  cufftcr3k(&c_k[0][0][0],&c_r[0][0][0]);
/*
#pragma acc data present(c_r,hc_r)
#pragma acc parallel loop collapse(3) private(cr,x_i,x_j,x_k)
for_xijk
//cr=c_r[x_i][x_j][x_k];
//hc_r[x_i][x_j][x_k]=cr*cr*(3-2*cr);
if(c_r[x_i][x_j][x_k]>=0.5)
continue;
//hc_r[x_i][x_j][x_k]=1/3*c_r[x_i][x_j][x_k]-1.767;
efor_xijk
*/
/*
  compute_mur_c(mu_cr);
  cufftrc3k(&mu_cr[0][0][0],&mu_ck[0][0][0]);
#pragma omp parallel for private(g_i,g_j,g_k,gijk,l,temp1)
#pragma acc data present(c_k,g_mod2,mu_ck)
#pragma acc parallel loop collapse(3) private(muk,g2,ck,g_i,g_j,g_k)
//#pragma acc kernels
  fgijk
	ck=c_k[g_i][g_j][g_k];
  	g2=g_mod2[g_i][g_j][g_k];
        muk=mu_ck[g_i][g_j][g_k];
  	c_k[g_i][g_j][g_k].x=(ck.x-Mc*g2*muk.x*dt)/            \
						(1+kc*Mc*g2*g2*dt);
  	c_k[g_i][g_j][g_k].y=(ck.y-Mc*g2*muk.y*dt)/            \
						(1+kc*Mc*g2*g2*dt);
  efgijk
  //--------------------------------------------------
*/
  //Evolving eta eta(k)
  compute_mur_eta(mu_r);
  //  #pragma acc update host(mu_r)
//		sprintf(filename,"output_files/mur_t%d.vtk",timestep);
//		write_vtk(&mu_r[0][0][0][0],filename,"etar");

  for(i=0;i<MartV;i++){
  cufftrc3k(&mu_r[i][0][0][0],&mu_k[i][0][0][0]); 
#pragma omp parallel for private(g_i,g_j,g_k,l)
#pragma acc data present(eta_k,g_mod2,mu_k)
#pragma acc parallel loop collapse(3) private(muk,g2,nk,g_i,g_j,g_k)
fgijk
	nk=eta_k[i][g_i][g_j][g_k];
  	g2=g_mod2[g_i][g_j][g_k];
    muk=mu_k[i][g_i][g_j][g_k];
	eta_k[i][g_i][g_j][g_k].x=(nk.x-Meta*(muk.x)*dt)/(1+Meta*keta*g2*dt);
	eta_k[i][g_i][g_j][g_k].y=(nk.y-Meta*(muk.y)*dt)/(1+Meta*keta*g2*dt);
  efgijk
  cufftcr3k(&eta_k[i][0][0][0],&eta_r[i][0][0][0]);
}
#pragma acc parallel loop collapse(3) private(x_i,x_j,x_k)
  for_xijk
if(eta_r[4][x_i][x_j][x_k]>=0.5)
{
eta_r[0][x_i][x_j][x_k]=0;
eta_r[1][x_i][x_j][x_k]=0;
eta_r[2][x_i][x_j][x_k]=0;
eta_r[3][x_i][x_j][x_k]=0;
}
  efor_xijk
  
#ifdef NOISE_ON
  if(timestep<=NRAND)
  {
//if(L!=2 && M!=2)//for 3D case
//{
  //Adding noise 
#pragma acc data create(randnum) present(eta_r)
  for(i=0;i<MartV;i++)
 // if(i>3)
//  continue;
		{
curand_generate((float*)acc_deviceptr(&randnum[0][0][0]),L*M*N);
#pragma acc parallel loop collapse(3) private(x_i,x_j,x_k)
  for_xijk
	eta_r[i][x_i][x_j][x_k]+=Namp*(randnum[x_i][x_j][x_k]-0.5);
  efor_xijk
}
//}
/*
else if(L==2 && M!=2)//for 2D case
{
  //Adding noise
#pragma acc data create(randnum) present(eta_r)
		{
curand_generate((float*)acc_deviceptr(&randnum[0][0][0]),L*M*N);
#pragma acc parallel loop collapse(3) private(x_i,x_j,x_k)
  for_xijk
	if(x_i==0)
	{eta_r[x_i][x_j][x_k]+=Namp*(randnum[x_i][x_j][x_k]-0.5);
	eta_r[1][x_j][x_k]=eta_r[x_i][x_j][x_k];
	}
  efor_xijk
		}//end of data region
}
else if(L==2 && M==2)//for 1D case
{
  //Adding noise
  for_xijk
	if(x_i==0 && x_j==0)
	{eta_r[x_i][x_j][x_k]+=Namp*(unirand()-0.5);
	eta_r[0][1][x_k]=eta_r[x_i][x_j][x_k];
	eta_r[1][0][x_k]=eta_r[x_i][x_j][x_k];
	eta_r[1][1][x_k]=eta_r[x_i][x_j][x_k];
	}
  efor_xijk

}
  }
  else
  {
#pragma acc data create(randnum) present(eta_r)
		{
curand_generate((float*)acc_deviceptr(&randnum[0][0][0]),L*M*N);
#pragma acc parallel loop collapse(3) private(x_i,x_j,x_k)
  for_xijk
	eta_r[x_i][x_j][x_k]+=Namp2*(randnum[x_i][x_j][x_k]-0.5);
  efor_xijk
}
 */ 
  }
  
#endif  

//cut off eta field
#pragma acc data present(eta_r)
  for(i=0;i<4;i++)
  {
#pragma acc parallel loop collapse(3) private(x_i,x_j,x_k)
  fxijk
  	eta_r[i][x_i][x_j][x_k]=(eta_r[i][x_i][x_j][x_k]>1.0)?1.0:eta_r[i][x_i][x_j][x_k];
  	eta_r[i][x_i][x_j][x_k]=(eta_r[i][x_i][x_j][x_k]<0.0)?0.0:eta_r[i][x_i][x_j][x_k];
  efxijk
  cufftrc3k(&eta_r[i][0][0][0],&eta_k[i][0][0][0]);	 
}
/*
#pragma acc data present(c_r,eta_r)
  for(i=0;i<4;i++)
  {
#pragma acc parallel loop collapse(3) private(x_i,x_j,x_k)
for_xijk
//cr=c_r[x_i][x_j][x_k];
//hc_r[x_i][x_j][x_k]=cr*cr*(3-2*cr);
if(c_r[x_i][x_j][x_k]>=0.5)
eta_r[i][x_i][x_j][x_k]=0;
//hc_r[x_i][x_j][x_k]=1/3*c_r[x_i][x_j][x_k]-1.767;
efor_xijk
}
*/
#ifdef MISFIT_COUPLED_WITH_ETA
 //Updating h(eta)
// printf("misfit coupled with eta\n");
//#pragma omp parallel for private(x_i,x_j,x_k,xijk,n)
#pragma acc data present(heta_r,eta_r)
#pragma acc parallel loop collapse(3) private(x_i,x_j,x_k,n)
for_xijk
heta_r[x_i][x_j][x_k]=0;
  for(i=0;i<V;i++)
  {
	n=eta_r[i][x_i][x_j][x_k];
	heta_r[x_i][x_j][x_k]+=n*n*(3-2*n);
//	heta_r[x_i][x_j][x_k]+=eta_r[i][x_i][x_j][x_k];
}
efor_xijk
// printf("misfit coupled with eta end\n");
  cufftrc3k(&heta_r[0][0][0],&heta_k[0][0][0]);
 
#endif

//#ifdef MISFIT_COUPLED_WITH_CONC
//Updating h(c)
//cufftrc3k(&hc_r[0][0][0],&hc_k[0][0][0]);	
//#endif
/*
#ifdef MISFIT_COUPLED_WITH_ETA
//Updating h(c)
#pragma acc data present(eta_r,hc_r)
#pragma acc parallel loop collapse(3) private(nr,x_i,x_j,x_k)
for_xijk
nr=eta_r[x_i][x_j][x_k];
hc_r[x_i][x_j][x_k]=nr*nr*(3-2*nr);
//hc_r[x_i][x_j][x_k]=1/3*c_r[x_i][x_j][x_k]-1.767;
efor_xijk
cufftrc3k(&hc_r[0][0][0],&hc_k[0][0][0]);	
#endif
*/
}//end of GetEintC
/*
static inline void compute_mur_c(double mu[L][M][N])
{
int x_i,x_j,x_k,xijk;	
int g_i,g_j,g_k,gijk;
int j;
//inhom_solver();// to compute the local stress

#pragma omp parallel for private(x_i,x_j,x_k,xijk,j)
#pragma acc data present(eta_r,c_r,mu,Sig,E_oo,Sij0)
#pragma acc parallel loop collapse(3) private(nr,hr,cr,x_i,x_j,x_k,j)
for_xijk
	hr=0;
	for(j=0;j<V;j++)
	{
	nr=eta_r[j][x_i][x_j][x_k];
	hr+=nr*nr*(3-2*nr);
	}
	cr=c_r[x_i][x_j][x_k];
	//mu from chemical free energy
	  
		mu[x_i][x_j][x_k]=hr*A1* \
	( cr*(2+cr*(-6+cr*4)) + B1*6*cr*(cr-1) ) \
		+ (1-hr)*A2*6*cr*(1-cr);

		mu[x_i][x_j][x_k]=0;

#ifdef MISFIT_COUPLED_WITH_CONC
		voigt e_i;//elastic strain
		voigt sigma;// local stress
			real hrp=6*cr*(cr-1);
			//real hrp=1/3;
		V6_loop{sigma[mi]=Sig[xijk][mi];}
		V6_loop{
				e_i[mi]=0.0;
				V6p_loop{e_i[mi]+=Sij0[mi][mip]*sigma[mip];}
			}
		real mu_elas=0.0;
		V6_loop{
				mu_elas+=(-sigma[mi]*E_oo[mi]);
				V6p_loop{mu_elas+=(-0.5*Cij2[mi][mip]*e_i[mi]*e_i[mip]);}
		}
		mu_elas*=hrp;
		mu[x_i][x_j][x_k]+=mu_elas;
#endif

efor_xijk
}//end of compute_mu_c()
*/
static inline void compute_mur_eta(double mu[V][L][M][N])
{
int x_i,x_j,x_k,xijk;	
int g_i,g_j,g_k,gijk;
int i;
float temp;
float xNi;
float e1,e2,e3,e4,h,dh,hs,dg;
//float x_Ni,x_Pt;

   inhom_solver();// to compute the local stress
//Computing mu_eta
//for(i=0;i<V;i++)
//{
#pragma acc enter data create(mu,conc)
#pragma acc data present(eta_r,mu,Sig,Sij0,te_p,conc)
#pragma acc parallel loop collapse(3) private(xNi,nr,e2,e3,e4,h,dh,x_i,x_j,x_k,i,temp)
for_xijk
//	cr=c_r[x_i][x_j][x_k];
        temp=0;
xNi=conc[x_i][x_j][x_k]-0.022;
//xNi=0.5;
for(i=0;i<V;i++)
{
//if(i>3)
//continue;
	nr=eta_r[i][x_i][x_j][x_k];
	temp+=nr*nr;
	}

for(i=0;i<V;i++)
{
//x_Ni=0.503;
//if(i>3)
//continue;
//x_Pt=0.2;
	nr=eta_r[i][x_i][x_j][x_k];
	//	mu[i][x_i][x_j][x_k]=-0.1*nr;

/*
#ifdef ALLOW_APB
		mu[x_i][x_j][x_k]= \
			+6*nr*(1-fabs(nr))*(                       	   \
				A1*( cr*cr*(1+cr*(-2+cr)) + B1 )       \
				-(A2+A1*B1)*( cr*cr*(3-2*cr) )   \
				)                                  \
			+W*( nr*(2-6*fabs(nr)+4*nr*nr) ) ;
#else
*/
/*
		mu[x_i][x_j][x_k]= \
			+6*nr*(1-nr)*(                       	   \
				( cr*cr*(1+cr*(-2+cr)) + B1 )       \
				-(A2+A1*B1)*( cr*cr*(3-2*cr) )   \
				)                                  \
			+W*( nr*(2+nr*(-6+4*nr)) ) ;
*/
//#ifdef FREE_ENERGY_MARTENSITE
		mu[i][x_i][x_j][x_k]= \
               A1*nr*(0.0466*(T-(310-10000*(xNi-0.5)))-1.02*nr*nr+1.195*temp*temp);
//#endif


 //              A1*nr*(0.1-1.8213*nr*nr+2.1343*temp*temp);
//               A1*nr*(-0.1);

#ifdef MISFIT_COUPLED_WITH_ETA
		voigt e_i;//elastic strain
		voigt sigma;// local stress
			real hrp=6*nr*(1-nr);
			//real hrp=1/3;
		V6_loop{sigma[mi]=Sig[xijk][mi];}
		V6_loop{
				e_i[mi]=0.0;
				V6p_loop{e_i[mi]+=Sij0[mi][mip]*sigma[mip];}
			}
		real mu_elas=0.0;
		V6_loop{
				mu_elas+=(-sigma[mi]*te_p[i][mi]);
				V6p_loop{mu_elas+=(-0.5*Cij2[mi][mip]*e_i[mi]*e_i[mip]);}
		}
		mu_elas*=hrp;
		mu[i][x_i][x_j][x_k]+=mu_elas;
#endif
}
efor_xijk
//}
}//end of compute_mu_eta()

