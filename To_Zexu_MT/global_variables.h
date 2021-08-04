#ifndef GLOBAL_VARIABLES_H
#define GLOBAL_VARIABLES_H
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <dirent.h>
#include <sys/stat.h>
#include <time.h>
#include <string.h>
#include <openacc.h>
#include <complex.h>
#include <cufft.h>
#include "constants.h"
#include "shortcut_wrappers.h"

//real field variables
extern double c_r[L][M][N];
extern double conc[L][M][N];//xNi
extern double conc2[L][M][N];//xHf
extern double eta_r[V][L][M][N];
extern double *etar;

extern double mu_r[V][L][M][N];
extern double mu_cr[L][M][N];

extern double g_v[L][M][N/2+1][3];
extern double g_mod2[L][M][N/2+1];

//complex field variables 
extern cufftDoubleComplex c_k[L][M][N/2+1];
extern cufftDoubleComplex eta_k[V][L][M][N/2+1];
extern cufftDoubleComplex *etak;

extern cufftDoubleComplex heta_k[L][M][N/2+1];
extern cufftDoubleComplex *hetak;

extern cufftDoubleComplex mu_k[V][L][M][N/2+1];
extern cufftDoubleComplex mu_ck[L][M][N/2+1];

//to store h(c)=1-c*c*(3-2*c)
extern double hc_r[L][M][N];
extern double heta_r[L][M][N];
extern double *hetar;

extern cufftDoubleComplex hc_k[L][M][N/2+1];

extern int total_step;
extern char *outdir;
extern char *c_fin;
extern char *eta_fin;
extern int save_config;
extern double dt;

#endif

