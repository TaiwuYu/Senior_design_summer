//all the constants are derived from the literature
//J. Zhu, H.-H.Wu, Y. Wu et al, Acta Materialia (2021)
//The nomalization factor of chemical free energy Delta F=10e7
//
#ifndef CONSTANTS_H
 #define CONSTANTS_H
 // Constants for the HEA model
 
//#define BIN_FILE

#define FROM_INPUT_FILE
  //  	#define RESTART_FROM_END
		#define RESTART_TIMESTEP 0

#define SAVE_CONFIG 1000 

#define NOISE_ON

//#define ALLOW_APB

 
 //Constants in the free energy
 #define A1 (1.0)
 #define T (240.0) //temperature of MT
 //#define T0 (291.0) //temperature
 #define B1 (0.1)
 #define W	(0.0)
 #define V  (5)    //
 #define MartV (4) //4 variants for B2 to R transformation
 //Constants of the system size, discretization
 #define L (128)
   #define M (128)
      #define N (128)
      #define LMN (L*M*N)
#define HLMN (L*M*(N/2+1))
 #define dx (1.0)
 #define dy (1.0)
 #define dz (1.0)
 #define dt_value (0.01)

//Kinetic parameters and gradient coefficient
 #define Mc 	(1.0)
 #define Meta 	(5.0)
 #define kc    (0.0)
 #define keta  (0.094)

//Noise parameters
	#define Namp (0.08)
	#define Namp2 (0.0)
	#define NRAND (50)

//#define DEBUG_MODE_ON
//#define x_Ni (0.503)
//#define x_Pt (0.2)
#define s_app_inc (0.0)

#define MAXITER (100) // maximum iterations for the inhomogeneous solver
//#define MISFIT_COUPLED_WITH_CONC
#define MISFIT_COUPLED_WITH_ETA
#endif

  
 
 
 
