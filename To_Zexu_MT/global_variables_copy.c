#include "global_variables.h"

//Field variables
double c_r[L][M][N];
double eta_r[L][M][N];
#ifdef MISFIT_COUPLED_WITH_HETA
double heta_r[L][M][N];
#endif
double mu_r[L][M][N];

double g_v[L][M][N/2+1][3];
double g_mod2[L][M][N/2+1];

cufftDoubleComplex c_k[L][M][N/2+1];
cufftDoubleComplex eta_k[L][M][N/2+1];
cufftDoubleComplex heta_k[L][M][N/2+1];
cufftDoubleComplex mu_k[L][M][N/2+1];

//to store h(c)=1-c*c*(3-2*c)
double hc_r[L][M][N];
cufftDoubleComplex hc_k[L][M][N/2+1];

int total_step=10;
char *outdir=NULL;
char *c_fin=NULL;
char *eta_fin=NULL;
int save_config=1;
double dt=dt_value;

