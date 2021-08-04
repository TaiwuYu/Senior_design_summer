#include "global_variables.h"
#include "io.h"
#include "pfr_cufft1.h"
#include "initialization.h"
#define CURR_LINE() {printf("%s at %d:\n",__FILE__,__LINE__);}
void initial_condition5(double r);
void usage(int argc, char *argv[]);
void init_field()
{
int x_i,x_j,x_k,xijk,i;
double p,c,n;
char filename[500];
#ifdef FROM_INPUT_FILE
/*#ifndef RESTART_FROM_END
//Reading for eta_r
	sprintf(filename,"%s",eta_fin);
	init_from_file(&eta_r[0][0][0],filename);
//Reading for c_r
	sprintf(filename,"%s",c_fin);
	init_from_file(&c_r[0][0][0],filename);
#else
*/
//Reading for eta_r
CURR_LINE(); 
	initial_condition1();
//	sprintf(filename,"%s/etar_rand3_t10000.vtk",outdir);
//	init_from_file(&eta_r[0][0][0],filename);
//Reading for c_r
CURR_LINE(); 
//	sprintf(filename,"input_files/etar_test_0_t10.vtk");
//	sprintf(filename,"input_files/eta_v1_hole40_09.vtk");
//	sprintf(filename,"input_files/etar_3d3_0_t1000.vtk");
//	init_from_file(&c_r[0][0][0],filename);
//	sprintf(filename,"input_files/cr_3d505_t12000_new.vtk");
//	sprintf(filename,"input_files/cr_3d503_t100000_new.vtk");
	sprintf(filename,"input_files/Ni_rand5_11.vtk");
	init_from_file(&conc[0][0][0],filename);
//	sprintf(filename,"input_files/Hf_new.vtk");
	sprintf(filename,"input_files/Hf_rand5_11.vtk");
	init_from_file(&conc2[0][0][0],filename);

//	sprintf(filename,"input_files/505_T55r_0_t10000_new.vtk");
//	sprintf(filename,"input_files/505_T45_pure_0_t40000_new.vtk");
//	sprintf(filename,"input_files/505_T110_0_t10000_new.vtk");
	sprintf(filename,"input_files/505_T10_multi_0_t2000_new.vtk");
//	sprintf(filename,"input_files/T70_nonequi_v1_01_new.vtk");
	init_from_file(&eta_r[0][0][0][0],filename);
//	sprintf(filename,"input_files/505_T45_pure_1_t40000_new.vtk");
//	sprintf(filename,"input_files/505_T55r_1_t10000_new.vtk");
//	sprintf(filename,"input_files/505_T110_1_t10000_new.vtk");
	sprintf(filename,"input_files/505_T10_multi_1_t2000_new.vtk");
//	sprintf(filename,"input_files/T70_nonequi_v2_01_new.vtk");
	init_from_file(&eta_r[1][0][0][0],filename);
//	sprintf(filename,"input_files/505_T110_2_t10000_new.vtk");
	sprintf(filename,"input_files/505_T10_multi_2_t2000_new.vtk");
//	sprintf(filename,"input_files/505_T45_pure_2_t40000_new.vtk");
//	sprintf(filename,"input_files/505_T55r_2_t10000_new.vtk");
//	sprintf(filename,"input_files/T70_nonequi_v3_01_new.vtk");
	init_from_file(&eta_r[2][0][0][0],filename);
	sprintf(filename,"input_files/505_T10_multi_3_t2000_new.vtk");
//	sprintf(filename,"input_files/T70_nonequi_v4_01_new.vtk");
//	sprintf(filename,"input_files/505_T110_3_t10000_new.vtk");
//	sprintf(filename,"input_files/505_T45_pure_3_t40000_new.vtk");
//	sprintf(filename,"input_files/505_T55r_3_t10000_new.vtk");
	init_from_file(&eta_r[3][0][0][0],filename);

//	sprintf(filename,"input_files/etar_3d505_0_t12000.vtk");
	sprintf(filename,"input_files/eta_v1_rand5_11.vtk");
	init_from_file(&eta_r[4][0][0][0],filename);
	sprintf(filename,"input_files/eta_v2_rand5_11.vtk");
	init_from_file(&eta_r[5][0][0][0],filename);
	sprintf(filename,"input_files/eta_v3_rand5_11.vtk");
	init_from_file(&eta_r[6][0][0][0],filename);
	sprintf(filename,"input_files/eta_v4_rand5_11.vtk");
	init_from_file(&eta_r[7][0][0][0],filename);
	sprintf(filename,"input_files/eta_v5_rand5_11.vtk");
	init_from_file(&eta_r[8][0][0][0],filename);
	sprintf(filename,"input_files/eta_v6_rand5_11.vtk");
	init_from_file(&eta_r[9][0][0][0],filename);

//	sprintf(filename,"input_files/pure_T310_0_t10000_new.vtk");
/*	sprintf(filename,"input_files/void_r40_T300_0_t10000_new.vtk");
	init_from_file(&eta_r[0][0][0][0],filename);
//	sprintf(filename,"input_files/pure_T310_1_t10000_new.vtk");
	sprintf(filename,"input_files/void_r40_T300_1_t10000_new.vtk");
	init_from_file(&eta_r[1][0][0][0],filename);
//	sprintf(filename,"input_files/pure_T310_2_t10000_new.vtk");
	sprintf(filename,"input_files/void_r40_T300_2_t10000_new.vtk");
	init_from_file(&eta_r[2][0][0][0],filename);
//	sprintf(filename,"input_files/pure_T310_3_t10000_new.vtk");
	sprintf(filename,"input_files/void_r40_T300_3_t10000_new.vtk");
	init_from_file(&eta_r[3][0][0][0],filename);
*/
CURR_LINE(); 

//#endif
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
//	initial_condition4();
//	initial_condition5(5.0);
//        initial_particle(15,40);
	initial_condition1();
#endif

//		sprintf(filename,"output_files/cr_init.vtk");
//		write_vtk2(&c_r[0][0][0],filename,"cr");

#ifdef MISFIT_COUPLED_WITH_CONC
//Initializing h(c)
for_xijk
//c=c_r[x_i][x_j][x_k];
c=c_r[x_i][x_j][x_k];
hc_r[x_i][x_j][x_k]=1-c*c*(3-2*c);
//hc_r[x_i][x_j][x_k]=1/3*c_r[x_i][x_j][x_k]-1.767;
efor_xijk
#endif

CURR_LINE(); 
for_xijk
c=c_r[x_i][x_j][x_k];
hc_r[x_i][x_j][x_k]=c*c*(3-2*c);
efor_xijk

#ifdef MISFIT_COUPLED_WITH_ETA
for_xijk
//c=c_r[x_i][x_j][x_k];
hetar[xijk]=0;
for(i=0;i<V;i++)
{
n=eta_r[i][x_i][x_j][x_k];
hetar[xijk]+=n*n*(3-2*n);
//hc_r[x_i][x_j][x_k]=1/3*c_r[x_i][x_j][x_k]-1.767;
}
efor_xijk
#endif

CURR_LINE(); 
//Initializing stuffs in fourier space
cufftrc3(&c_r[0][0][0],&c_k[0][0][0]);
for(i=0;i<V;i++)
cufftrc3(&eta_r[i][0][0][0],&eta_k[i][0][0][0]);

cufftrc3(&hc_r[0][0][0],&hc_k[0][0][0]);
cufftrc3(&heta_r[0][0][0],&heta_k[0][0][0]);
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

