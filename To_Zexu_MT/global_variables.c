#include "global_variables.h"

//Field variables
double c_r[L][M][N];
double conc[L][M][N];
double conc2[L][M][N];
double eta_r[V][L][M][N];
double *etar=&eta_r[0][0][0][0];

//#ifdef MISFIT_COUPLED_WITH_ETA
double heta_r[L][M][N];
double *hetar=&heta_r[0][0][0];
//#endif
double mu_r[V][L][M][N];
double mu_cr[L][M][N];

double g_v[L][M][N/2+1][3];
double g_mod2[L][M][N/2+1];

cufftDoubleComplex eta_k[V][L][M][N/2+1];
cufftDoubleComplex *etak=&eta_k[0][0][0][0];

cufftDoubleComplex c_k[L][M][N/2+1];
cufftDoubleComplex heta_k[L][M][N/2+1];
cufftDoubleComplex *hetak=&heta_k[0][0][0];
cufftDoubleComplex mu_k[V][L][M][N/2+1];
cufftDoubleComplex mu_ck[L][M][N/2+1];

//to store h(c)=1-c*c*(3-2*c)
double hc_r[L][M][N];
cufftDoubleComplex hc_k[L][M][N/2+1];

int total_step=1000;
char *outdir=NULL;
char *c_fin=NULL;
char *eta_fin=NULL;
int save_config=100;
double dt=dt_value;

