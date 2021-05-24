/*
 * =====================================================================================
 *
 *       Filename:  inhom.h
 *
 *    Description:  Phase-field inhomogeneous elasticity solver
 *
 *        Version:  1.0
 *        Created:  07/29/2016 13:20:03
 *       Revision:  none
 *       Compiler:  gcc
 *
 *
 * =====================================================================================
 */

#ifndef H_INHOM
#define H_INHOM

#include <math.h>
#include <cufft.h>
#include "constants.h"
/**********************
  types and structures
  *********************/
typedef double real;
typedef cufftDoubleComplex complx;
typedef real ten2nd[3][3];
typedef complx ten2ndk[3][3];
typedef real ten4th[3][3][3][3];
typedef real voigt[6];
typedef complx voigtk[6];
typedef real voigt66[6][6];


/**********************
Function to solve for stress
  *********************/
void inhom_init();
void inhom_finish();
void inhom_solver();
double cal_Felas();

/***************************
  loops, allocators, indices
 ***************************/
#define AllocMem(a, n, t) a = (t*)calloc((n),sizeof(t))

#define pIDX ((px*M + py)*N + pz)	/* index for local array in real space*/
#define kIDX ((kx*M + ky)*(N/2+1) + kz)	/* index for local array in k space*/
#define pxyz px][py][pz
#define kxyz kx][ky][kz

#define x_loop for(int px=0;px<L;px++)for(int py=0;py<M;py++)for(int pz=0;pz<N;pz++)
#define k_loop for(int kx=0;kx<L;kx++)for(int ky=0;ky<M;ky++)for(int kz=0;kz<N/2+1;kz++)

#define T2_loop for(int mi=0;mi<3;mi++)for(int mj=0;mj<3;mj++)
#define T2p_loop for(int mip=0;mip<3;mip++)for(int mjp=0;mjp<3;mjp++)
#define V6_loop for(int mi=0; mi<6; mi++)
#define V6p_loop for(int mip=0; mip<6; mip++)
#define C4_loop for(int mi=0;mi<3;mi++)for(int mj=0;mj<3;mj++)for(int mk=0;mk<3;mk++)for(int ml=0;ml<3;ml++)
#define C6_loop for(int mi=0;mi<6;mi++)for(int mj=0;mj<6;mj++)
#define field_loop for(int ip=0;ip<V;ip++)
#define field2_loop for(int ip=0;ip<V;ip++)for(int iq=0;iq<V;iq++)

/***************************
  units & constants
 ***************************/

#define RSQ2 0.70710678118654744

#define ERROR_REL_ELAST 1E-5
#define ERROR_ABS_ELAST 1E-6
#define TINY (1.0E-20)
#define PI 3.14159265359
#define TWO_PI 6.28318530718
#define EV_IN_J 1.60217733e-19
#define J_IN_EV (1.0/EV_IN_J)
#define MY_J_IN_EV 0.62415063631
#define BOLZ_IN_J__K 1.3806488e-23
#define BOLZ_IN_EV__K 8.6173324e-5
#define BOLZ BOLZ_IN_J__K
#define MY_BOLTZ    1.3806488

/*  Some math operations */
#define VSecNorm(v)	\
	(sqrt(Sqr(v[0]) + Sqr(v[1]) + Sqr(v[2])))
#define VScale(v, s)  \
	v[0] *= s,		\
	v[1] *= s,		\
	v[2] *= s
#define Sqr(x)	((x) * (x))

extern voigt E_oo;
extern voigt66 Sij0;
extern voigt66 Cij2;
extern voigt Sig[LMN];
extern voigt rSig[LMN];
extern voigt te_p[V]; 
extern voigt T_tr[LMN];
//extern voigt total_strain[LMN];//Output total strain


static inline void chg_basis_Kelvin_1(voigt V6, ten2nd T2)
{
        for(int i=0; i<3; i++) T2[i][i] = V6[i];
        T2[2][1] = T2[1][2] = V6[3]*RSQ2;
        T2[2][0] = T2[0][2] = V6[4]*RSQ2;
        T2[0][1] = T2[1][0] = V6[5]*RSQ2;
}
static inline void chg_basis_Kelvin_1k(voigtk V6, ten2ndk T2)
{
		//real
        for(int i=0; i<3; i++) T2[i][i].x = V6[i].x;
        T2[2][1].x = T2[1][2].x = V6[3].x*RSQ2;
        T2[2][0].x = T2[0][2].x = V6[4].x*RSQ2;
        T2[0][1].x = T2[1][0].x = V6[5].x*RSQ2;
		//imaginary
        for(int i=0; i<3; i++) T2[i][i].y = V6[i].y;
        T2[2][1].y = T2[1][2].y = V6[3].y*RSQ2;
        T2[2][0].y = T2[0][2].y = V6[4].y*RSQ2;
        T2[0][1].y = T2[1][0].y = V6[5].y*RSQ2;
}
static inline void chg_basis_Kelvin_2(ten2nd T2, voigt V6 )
{
        for(int i=0; i<3; i++) V6[i] = T2[i][i];
        V6[3] = T2[2][1]/RSQ2;
        V6[4] = T2[2][0]/RSQ2;
        V6[5] = T2[0][1]/RSQ2;
}
static inline void chg_basis_Kelvin_2k(ten2ndk T2, voigtk V6 )
{
        for(int i=0; i<3; i++) V6[i].x = T2[i][i].x;
        V6[3].x = T2[2][1].x/RSQ2;
        V6[4].x = T2[2][0].x/RSQ2;
        V6[5].x = T2[0][1].x/RSQ2;
        for(int i=0; i<3; i++) V6[i].y = T2[i][i].y;
        V6[3].y = T2[2][1].y/RSQ2;
        V6[4].y = T2[2][0].y/RSQ2;
        V6[5].y = T2[0][1].y/RSQ2;
}
static inline void chg_basis_Kelvin_3( voigt66 C2, ten4th T4)
{

        C4_loop{
            int p = (mi+1)*(mi==mj)+(1-(mi==mj))*(7-mi-mj);
            int q = (mk+1)*(mk==ml)+(1-(mk==ml))*(7-mk-ml);
            real t1 = ((real)(mi==mj)+(1.-(mi==mj))/RSQ2);
            real t2 = ((real)(mk==ml)+(1.-(mk==ml))/RSQ2);
            T4[mi][mj][mk][ml] = C2[p-1][q-1]/t1/t2;
        }
}
static inline void chg_basis_Kelvin_4(ten4th T4, voigt66 C2)
{
       C4_loop{
            int p = (mi+1)*(mi==mj)+(1-(mi==mj))*(7-mi-mj);
            int q = (mk+1)*(mk==ml)+(1-(mk==ml))*(7-mk-ml);
            real t1 = ((real)(mi==mj)+(1.-(mi==mj))/RSQ2);
            real t2 = ((real)(mk==ml)+(1.-(mk==ml))/RSQ2);
            C2[p-1][q-1] = t1*t2*T4[mi][mj][mk][ml];
        }
}
#endif


