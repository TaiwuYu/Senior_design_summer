#include "global_variables.h"
#include "debug_functions.h"
#include "pfr_cufft1.h"
#include "io.h"
#include "cal_kvec.h"
#include "initialization.h"
#include "time_evolution.h"
#include "curand_cuda.h"
#include "inhom_v3.h"
void arg_initialize(int argc, char *argv[]);

int main(int argc, char *argv[])
{
int timestep;
time_t t1,t2;
double tol=1e-10;
double Fchem,Felas;
double Fchem_prev,Felas_prev;
char dirname_output[]="output_files";
char filename_cin[]="cr.vtk";
char filename_etain[]="eta.vtk";
outdir=&dirname_output[0];
c_fin=&filename_cin[0];
eta_fin=&filename_etain[0];

printf("Starting the program ...\n");
arg_initialize(argc,argv);

//Initializing cufft
t1=time(NULL);
cufft_start();//initializing cufft
#ifdef NOISE_ON
curand_init();//initialization of curand
#endif
t2=time(NULL);
printf("Time taken to initialize cufft %lf sec \n",difftime(t2,t1));


//Creating a new directory to store data To store output data
if(!opendir(outdir))  mkdir(outdir,0777);


printf("total_step = %d\n",total_step);

printf("Initializing K^2 table ...\n");
grad();//Initializing k^2 table
printf("K^2 table initialized \n");


printf("Initializing the field ...\n");
init_field();//Initializing the field
printf("All the field variables are initialized \n");

printf("Starting the time evolution...\n");
#pragma acc data copy(c_r,eta_r,mu_r,hc_r,g_v,g_mod2,mu_cr,heta_r,conc)
#pragma acc data copy(c_k,eta_k,mu_k,hc_k,mu_ck,heta_k)
{
//Initializa inhomogeneous solver
inhom_init();
OutputMed(RESTART_TIMESTEP);//printing the initial configuration
Fchem_prev=cal_Fchem();
Felas_prev=cal_Felas();
printf("(t,Fchem,Felastic,Ftotal) %d %10lf %10lf %10lf \n",RESTART_TIMESTEP,Fchem_prev,Felas_prev,Fchem_prev+Felas_prev);
cal_vfrac();

t1=time(NULL);

for(timestep=RESTART_TIMESTEP+1;timestep<=total_step+RESTART_TIMESTEP;timestep++)
{	GetEintC(timestep);
	if( (timestep % save_config )==0 ){
     #pragma acc update host(c_r,eta_r,c_k,eta_k,hc_r,mu_r,Sig,heta_r,rSig,T_tr)
	OutputMed(timestep);
	 	Fchem=cal_Fchem();
		Felas=cal_Felas();
		printf("(t,Fchem,Felastic,Ftotal) %d %10lf %10lf %10lf \n",timestep,Fchem,Felas,Fchem+Felas);
	cal_vfrac();
	if(isnan(c_r[0][0][0]))
		{fprintf(stderr,"#Error:Solution diverged \n");break;}
#ifdef DEBUG_MODE_ON
	if( fabs(Fchem-Fchem_prev+Felas-Felas_prev)<tol*save_config )	
		{printf("The Free energy of the system converged within the tolerance of %lf \n",tol);break;}
#endif
	Fchem_prev=Fchem;
	Felas_prev=Felas;
	}//end of if condition
}//end ot time loop

t2=time(NULL);
}//end of data region

printf("\n");
	if( ( (timestep-1)%save_config )!=0 )
	{
		Fchem=cal_Fchem();
		Felas=cal_Felas();
	printf("(t,Fchem,Felastic,Ftotal) %d %10lf %10lf %10lf \n",timestep,Fchem,Felas,Fchem+Felas);
		cal_vfrac();
	OutputMed(timestep-1);
}
if(L==2 && M==2)
	printf("Surface energy = %lf\n",cal_sigma());

cufft_finish();
#ifdef NOISE_ON
curand_finish();//finish of curand()
#endif
inhom_finish();
printf("Total cpu time utilized for time evolution = %lf sec\n",difftime(t2,t1));


printf("#All is well \n");

return 0;
}//end of main()




