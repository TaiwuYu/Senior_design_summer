#include "global_variables.h"
#include "io.h"
#include "pfr_cufft1.h"
#include "initialization.h"
void initial_condition5(double r);
void usage(int argc, char *argv[]);
void init_field()
{
int x_i,x_j,x_k,xijk;
double p,c,n;
char filename[500];
#ifdef FROM_INPUT_FILE
#ifndef RESTART_FROM_END
//Reading for eta_r
	sprintf(filename,"%s",eta_fin);
	init_from_file(&eta_r[0][0][0],filename);
//Reading for c_r
	sprintf(filename,"%s",c_fin);
	init_from_file(&c_r[0][0][0],filename);
#else
//Reading for eta_r
	sprintf(filename,"%s/etar_rand3_t10000.vtk",outdir);
	init_from_file(&eta_r[0][0][0][0],filename);
//Reading for c_r
	sprintf(filename,"%s/cr_rand3_t10000.vtk",outdir);
	init_from_file(&c_r[0][0][0],filename);
#endif
#else
//	sprintf(filename,"%s/etar_3d_grow2_t1400.vtk",outdir);
//	init_from_file(&eta_r[0][0][0],filename);
//Reading for c_r
//	sprintf(filename,"%s/conc_dissolv2.vtk",outdir);
//	init_from_file(&c_r[0][0][0],filename);

//        sprintf(filename,"%s/etar_3d3_t1100.vtk",outdir);
//	init_from_file(&eta_r[0][0][0],filename);
//Reading for c_r
//	sprintf(filename,"%s/Cr_diss_new.vtk",outdir);
//	init_from_file(&c_r[0][0][0],filename);
	initial_condition4();
//	initial_condition5(8.0);
//        initial_particle(8,100);
//	initial_condition1();
#endif


#ifdef MISFIT_COUPLED_WITH_CONC
//Initializing h(c)
for_xijk
//c=c_r[x_i][x_j][x_k];
c=c_r[x_i][x_j][x_k];
hc_r[x_i][x_j][x_k]=1-c*c*(3-2*c);
//hc_r[x_i][x_j][x_k]=1/3*c_r[x_i][x_j][x_k]-1.767;
efor_xijk
#endif
/*
#ifdef MISFIT_COUPLED_WITH_ETA
for_xijk
//c=c_r[x_i][x_j][x_k];
n=eta_r[x_i][x_j][x_k];
hc_r[x_i][x_j][x_k]=n*n*(3-2*n);
//hc_r[x_i][x_j][x_k]=1/3*c_r[x_i][x_j][x_k]-1.767;
efor_xijk
#endif
*/
//Initializing stuffs in fourier space
cufftrc3(&c_r[0][0][0],&c_k[0][0][0]);
for(int i=0;i<nv*ng;i++)
{
cufftrc3(&eta_r[i][0][0][0],&eta_k[i][0][0][0]);
}
cufftrc3(&hc_r[0][0][0],&hc_k[0][0][0]);
}//end of init_field()

void arg_initialize(int argc, char *argv[])
{
int i;
for(i=1;i<argc;i++)
	{
	if(strcmp("-Nstep",argv[i])==0)
		{total_step=atoi(argv[++i]);printf("Total step = %d\n",total_step);}
	else if(strcmp("-outdir",argv[i])==0)
		{outdir=argv[++i];printf("outdir=%s\n",outdir);}
	else if(strcmp("-c_fin",argv[i])==0)
		{c_fin=argv[++i];printf("c_fin=%s\n",c_fin);}
	else if(strcmp("-n_fin",argv[i])==0)
		{eta_fin=argv[++i];printf("eta_fin=%s\n",eta_fin);}
	else if(strcmp("-save",argv[i])==0)
		{save_config=atoi(argv[++i]);printf("save_config=%d\n",save_config);}
	else if(strcmp("-dt",argv[i])==0)
		{dt=atof(argv[++i]);printf("dt=%lf\n",dt);}
	else if(strcmp("-h",argv[i])==0)
		{usage(argc,argv);exit(0);}
	else
	{fprintf(stderr,"#Error: the argument %s can't be recognised\n",argv[i]);usage(argc,argv);exit(0);}
	}
}//end of arg_initialize()

void usage(int argc, char *argv[])
{
printf("%s usage: \n",argv[0]);
printf("%s [option] [value] \n",argv[0]);
printf(" \n");
printf("DESCRIPTION: \n");
printf("-Nstep		Total number of timesteps \n");
printf("-outdir         Name of the output directory \n");
printf("-c_fin		Input file (conc) \n");
printf("-n_fin		Input file (eta)  \n");
printf("-save		Frequency of saving the configuration file\n");
printf("-dt		delta t for integration                   \n");
printf("-h		To display this help message\n");
printf("Example: \n");
printf("%s -Nstep %d -outdir %s -c_fin %s -n_fin %s -save %d -dt %lf\n",argv[0],total_step,outdir,c_fin,eta_fin,save_config,dt);
}//end of usage

