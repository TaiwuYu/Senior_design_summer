#ifndef TIME_EVOLUTION_H
#define TIME_EVOLUTION_H

void GetEintC(int timestep);//to evolve one time iteration
static inline void compute_mur_c(double mu[L][M][N]);//to compute diffusion potential of c
static inline void compute_mur_eta(double mu[ng*nv][L][M][N]);//to compute diffusion potential of eta
#endif
