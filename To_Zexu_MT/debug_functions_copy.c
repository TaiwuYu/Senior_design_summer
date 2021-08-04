//function declaration
#include "global_variables.h"
#include "debug_functions.h"

void cal_vfrac()
{
double vfrac,cavg,vfrac1,vfrac2;
int x_i,x_j,x_k,xijk;
#ifdef MISFIT_COUPLED_WITH_HETA
double *hetar=&heta_r[0][0][0];
#endif
double *cr=&c_r[0][0][0];
double *etar=&eta_r[0][0][0];
//Calculating volume fraction
 for(int i=0;i<ng*nv;i++)
 {
vfrac=0;cavg=0;vfrac1=0;vfrac2=0;
for_xijk
	vfrac+=fabs(eta_r[i][x_i][x_j][x_k]);
	vfrac1+=eta_r[i][x_i][x_j][x_k]*(eta_r[i][x_i][x_j][x_k]>0);
	vfrac2+=-eta_r[i][x_i][x_j][x_k]*(eta_r[i][x_i][x_j][x_k]<0);
	cavg+=cr[xijk];
efor_xijk
vfrac/=L*M*N;
vfrac1/=L*M*N;
vfrac2/=L*M*N;
cavg/=L*M*N;
printf("vf = v%d, %7.6lf (%7.3lf,%7.3lf) c_avg = %7.6lf\n",i,vfrac,vfrac1,vfrac2,cavg);
}
}//end of cal_vfrac()

double cal_Fchem()
{
int x_i,x_j,x_k,xijk;
int g_i,g_j,g_k,gijk;
double Fchem=0,Fbulk=0,Fgrad1=0,Fgrad2=0;
double f_o=0,f_d,g,h,c,n;
cufftDoubleComplex ck,nk;
//Bulk free energy
Fbulk=0;
for_xijk
	c=c_r[x_i][x_j][x_k];
 for(int i=0;i<ng*nv;i++)
 {
#ifdef ALLOW_APB
	n=fabs(eta_r[x_i][x_j][x_k]);
#else
	n=eta_r[i][x_i][x_j][x_k];
#endif
h+=n*n*(3-2*n);
g+=n*n*(1-n)*(1-n);
}
f_o=A1*(c*c*(1-c)*(1-c)+B1*(1-c*c*(3-2*c)));//free energy of ordered phase
f_d=A2*c*c*(3-2*c);//free energy of disordered phase
Fbulk+=h*f_o+(1-h)*f_d+W*g;
efor_xijk
Fbulk*=dx*dy*dz;

//Gradient free energy
//calculated in reciprocal space
Fgrad1=0;Fgrad2=0;
for_gijk
ck=c_k[g_i][g_j][g_k];	
 for(int i=0;i<ng*nv;i++)
 {
nk=eta_k[i][g_i][g_j][g_k];

Fgrad1+=g_mod2[g_i][g_j][g_k]*(ck.x*ck.x+ck.y*ck.y)*(1+(g_k!=0));
Fgrad2+=g_mod2[g_i][g_j][g_k]*(nk.x*nk.x+nk.y*nk.y)*(1+(g_k!=0));
}
efor_gijk
Fgrad1*=kc/2/(L*M*N);
Fgrad2*=keta/2/(L*M*N);

Fchem=Fbulk+Fgrad1+Fgrad2;
return Fchem;
}//end of cal_Fchem()
/*
double cal_Felastic()
{
int g_i,g_j,g_k,gijk;
double Eel=0;
cufftDoubleComplex ck,hk;
#ifdef MISFIT_COUPLED_WITH_HETA
Eel=0;
for_gijk
hk=heta_k[g_i][g_j][g_k];
Eel+=Bhh[g_i][g_j][g_k]*(hk.x*hk.x+hk.y*hk.y)*(1+(g_k!=0));
efor_gijk
Eel*=0.5/(L*M*N);
#endif
#ifdef MISFIT_COUPLED_WITH_CONC
Eel=0;
for_gijk
hk=hc_k[g_i][g_j][g_k];
Eel+=Bhh[g_i][g_j][g_k]*(hk.x*hk.x+hk.y*hk.y)*(1+(g_k!=0));
efor_gijk
Eel*=0.5/(L*M*N);
#endif
return Eel;
}
*/
double cal_sigma()
{
double sigma=0;
sigma=cal_Fchem()/2/(dx*dy*L*M);
return sigma;
}//end of sigma()
