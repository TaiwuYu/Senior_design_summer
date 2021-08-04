#include "global_variables.h"
#include "io.h"
#include "inhom_v3.h"
void write_vtk(double *rfield, char filename[], char fieldname[] );
void write_vtk2(double *rfield, char filename[], char fieldname[] );
double Sig_out[LMN];
//voigt total_strain[LMN];
double eta_out[V][L+2][M+2][N+2];

void OutputMed(int timestep)
{
char filename[500];
int i,x_i,x_j,x_k;
double sum;
#ifdef BIN_FILE
		//For producing binary files
		FILE *fp=NULL;
		sprintf(filename,"output_files/time_%d.bin",timestep);
		fp=fopen(filename,"wb");
		if(fp)
		{	fwrite(  &c_r[0][0][0],sizeof(double),LMN,fp);
			fwrite(&eta_r[0][0][0],sizeof(double),V*LMN,fp);
#ifdef MISFIT_COUPLED_WITH_HETA
			fwrite(&heta_r[0][0][0],sizeof(double),LMN,fp);
#endif
		fclose(fp);
		}
		else
		{fprintf(stderr,"#Error:Could open the file %s to write\n",filename);exit(0);}
#else 
for(i=0;i<V;i++)
for(x_i=0;x_i<L+2;x_i++) 
	for(x_j=0;x_j<M+2;x_j++) 
		for(x_k=0;x_k<N+2;x_k++) 
		{ eta_out[i][x_i][x_j][x_k]=0;	
		}

for(i=0;i<V;i++)
for(x_i=1;x_i<L+1;x_i++) 
	for(x_j=1;x_j<M+1;x_j++) 
		for(x_k=1;x_k<N+1;x_k++)
		{
		eta_out[i][x_i][x_j][x_k]=eta_r[i][x_i-1][x_j-1][x_k-1];
		}

		//For producing vtk files
//		sprintf(filename,"%s/cr_void2_t%d.vtk",outdir,timestep);
//		write_vtk2(&c_r[0][0][0],filename,"cr");
#pragma omp parallel for private(filename)
for(i=0;i<MartV;i++)
	{	sprintf(filename,"%s/test_%d_t%d.vtk",outdir,i,timestep);
		write_vtk(&eta_out[i][0][0][0],filename,"etar");
}
//		sprintf(filename,"%s/mur_t%d.vtk",outdir,timestep);
//		write_vtk(&mu_r[0][0][0][0],filename,"etar");
/*		sprintf(filename,"%s/mucr_t%d.vtk",outdir,timestep);
		write_vtk(&mu_cr[0][0][0],filename,"mucr");
		sprintf(filename,"%s/hcr_t%d.vtk",outdir,timestep);
		write_vtk(&hc_r[0][0][0],filename,"hcr");
*/
/*
#pragma acc update host(Sig)
                x_loop{
                Sig_out[pIDX]=Sig[pIDX][0];
                }
		sprintf(filename,"%s/Sig11_void1_t%d.vtk",outdir,timestep);
		write_vtk2(&Sig_out[0],filename,"Sig");
#pragma acc update host(Sig)
                x_loop{
                Sig_out[pIDX]=Sig[pIDX][3];
                }
		sprintf(filename,"%s/Sig23_void1_t%d.vtk",outdir,timestep);
		write_vtk2(&Sig_out[0],filename,"Sig");
#pragma acc update host(Sig)
                x_loop{
                Sig_out[pIDX]=Sig[pIDX][4];
                }
		sprintf(filename,"%s/Sig13_void1_t%d.vtk",outdir,timestep);
		write_vtk2(&Sig_out[0],filename,"Sig");
#pragma acc update host(Sig)
                x_loop{
                Sig_out[pIDX]=Sig[pIDX][5];
                }
		sprintf(filename,"%s/Sig12_void1_t%d.vtk",outdir,timestep);
		write_vtk2(&Sig_out[0],filename,"Sig");
*/
//#ifdef MISFIT_COUPLED_WITH_ETA
//		sprintf(filename,"%s/hetar_3d_t%d.vtk",outdir,timestep);
//		write_vtk(&heta_r[0][0][0],filename,"hetar");
//#endif

//	x_loop{ V6_loop{ total_strain[pIDX][mi]=rSig[pIDX][mi]+T_tr[pIDX][mi]; } }
//        x_loop{sum+=rSig[pIDX][0]+T_tr[pIDX][0];}
//	sum/=LMN;

//printf("total strain11 =%f \n",sum);

		printf("timestep=%d completed \n",timestep);
		fflush(stdout);
//		cal_vfrac();
#endif
}//end of OutputMed()

void write_vtk(double *rfield, char filename[], char fieldname[] )
{
FILE *fp=NULL;
int i,j,k;
fp=fopen(filename,"w");
if(!fp)
	{fprintf(stderr,"#Error: %s could not be opened for writing\n",filename);exit(0);}

//writing vtk header
fprintf(fp,"# vtk DataFile Version 2.0\n");
fprintf(fp,"VTK from C-program\n");
fprintf(fp,"ASCII\n");
fprintf(fp,"DATASET STRUCTURED_POINTS\n");
fprintf(fp,"DIMENSIONS %d %d %d\n",L+2,M+2,N+2);
fprintf(fp,"SPACING 1 1 1\n");
fprintf(fp,"ORIGIN 0 0 0\n");
fprintf(fp,"POINT_DATA %d\n",(L+2)*(M+2)*(N+2));
fprintf(fp,"SCALARS %s float 1\n",fieldname);
fprintf(fp,"LOOKUP_TABLE default\n");
for(k=0; k<N+2; k++)
	for(j=0; j<M+2; j++)
		for(i=0; i<L+2; i++)
			fprintf(fp,"%.6lf \n",rfield[i*(M+2)*(N+2)+j*(N+2)+k]);
fclose(fp);
}//end of write_vtk()

void write_vtk2(double *rfield, char filename[], char fieldname[] )
{
FILE *fp=NULL;
int i,j,k;
fp=fopen(filename,"w");
if(!fp)
	{fprintf(stderr,"#Error: %s could not be opened for writing\n",filename);exit(0);}

//writing vtk header
fprintf(fp,"# vtk DataFile Version 2.0\n");
fprintf(fp,"VTK from C-program\n");
fprintf(fp,"ASCII\n");
fprintf(fp,"DATASET STRUCTURED_POINTS\n");
fprintf(fp,"DIMENSIONS %d %d %d\n",L,M,N);
fprintf(fp,"SPACING 1 1 1\n");
fprintf(fp,"ORIGIN 0 0 0\n");
fprintf(fp,"POINT_DATA %d\n",(L)*(M)*(N));
fprintf(fp,"SCALARS %s float 1\n",fieldname);
fprintf(fp,"LOOKUP_TABLE default\n");
for(k=0; k<N; k++)
	for(j=0; j<M; j++)
		for(i=0; i<L; i++)
			fprintf(fp,"%.6lf \n",rfield[i*(M)*(N)+j*(N)+k]);
fclose(fp);
}//end of write_vtk()


void init_from_file(double *tpr, char *filename)
{
	int x_i,x_j,x_k,xijk,i;
	int l;
	FILE *fp=NULL;
	char ignore[40];
	int fL,fM,fN;
fp=fopen(filename,"r");
	if(fp)
	{
		l=strlen(filename);
		if(strcmp(&filename[l-4],".dat")==0)
		{	for_xijk
			fscanf(fp,"%lf",&tpr[xijk]);
			efor_xijk
		}
		else if(strcmp(&filename[l-4],".vtk")==0)
		{	for(i=0;i<4;i++)
				fgets(ignore,40,fp);//ignoring first four lines	
			fscanf(fp,"DIMENSIONS %d %d %d\n",&fL,&fM,&fN);
			for(i=0;i<5;i++)
				fgets(ignore,40,fp);//ignoring 5 lines
			if( (L!=fL) | (M!=fM) | (N!=fN) )
				{fprintf(stderr,"#Error:Input file size (%d,%d,%d) is different from the expected size (%d,%d,%d)\n",fL,fM,fN,L,M,N);exit(0);}
for(x_k=0; x_k<N; x_k++)
	for(x_j=0; x_j<M; x_j++)
		for(x_i=0; x_i<L; x_i++)
			fscanf(fp,"%lf",&tpr[x_i*M*N+x_j*N+x_k]);
		}
		else
		{fprintf(stderr,"#Error:Unknown input file format,  %s",filename);exit(0);}
		fclose(fp);
	}
	else
	{fprintf(stderr,"#Error: the file %s could not be opened for reading\n",filename);exit(0);
	}
}//end of init_from_file

