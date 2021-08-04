#include "global_variables.h"
#include "misc.h"
//#include <iomanip>
#include <string.h>
#include <stdlib.h>

void initial_condition1();
void initial_condition2();
void initial_condition3();
void initial_condition4();
void initial_particle(int a, int b);
void init_homo();

void initial_condition1()
{// homogenous disordered phase with c=0.55 with a noise
double *cr=&c_r[0][0][0];
//double *etar[nv*ng]=&eta_r[ng*nv][0][0][0];
int x_i,x_j,x_k,xijk;
 init_homo();
        for(int i=0;i<V;i++)
	{
if(L!=2 && M!=2)
	{
	for_xijk	
	//		cr[xijk]=7.5+Namp*(unirand()-0.5);
			cr[xijk]=0.0;
		//	eta_r[i][x_i][x_j][x_k]=0+Namp*(unirand()-0.5);
			eta_r[i][x_i][x_j][x_k]=0;
	efor_xijk
	}
else if(L==2 && M!=2)//For 2D case
	{
	for_xijk	
			if(x_i==0)
		//	{cr[xijk]=6.0+Namp*unirand();
			{cr[xijk]=6.5+1.0*Namp*unirand();
			eta_r[i][x_i][x_j][x_k]=0+Namp*unirand();
		//	etar[xijk]=0.0;
			c_r[1][x_j][x_k]=cr[xijk];
			eta_r[i][1][x_j][x_k]=eta_r[i][x_i][x_j][x_k];
			}
	efor_xijk
	}
else if(L==2 && M==2)//For 1D case
	{
	for_xijk	
			cr[xijk]=6.5;
			if(x_i==0 && x_j==0)
			{
			eta_r[i][x_i][x_j][x_k]=0+Namp*(unirand()-0.5);
			eta_r[i][0][1][x_k]=eta_r[i][x_i][x_j][x_k];
			eta_r[i][1][0][x_k]=eta_r[i][x_i][x_j][x_k];
			eta_r[i][1][1][x_k]=eta_r[i][x_i][x_j][x_k];
			}
	efor_xijk
	}

	}

}//end of initial_condition1()

void initial_condition2()
{//Flat interface across z-direction between (c=0,eta=0) and (c=1,eta=1)
double *cr=&c_r[0][0][0];
//double *etar=&eta_r[0][0][0];
int x_i,x_j,x_k,xijk;

 init_homo();
        for(int i=0;i<V;i++)
	{
for_xijk	
		if(x_k<N/2)
		{
//		cr[xijk]=5.6;
//		etar[xijk]=0;
//		}
//		else 
//		{
		cr[xijk]=0;
		eta_r[0][x_i][x_j][x_k]=1;
		}

efor_xijk
}
}//end of initial_condition2()

void initial_condition3()
{// 2D square particle in the system
//double *cr=&c_r[0][0][0];
//double *etar=&eta_r[0][0][0];
int x_i,x_j,x_k,xijk;
double C=(N-1)/2.0;
 init_homo();
for_xijk	
		if( fabs(x_j-C)<=N/4.0 && fabs(x_k-C)<=N/4.0 )
		{
		c_r[x_i][x_j][x_k]=8.2;
		eta_r[0][x_i][x_j][x_k]=1;
		}
//		else 
//		{
//		cr[xijk]=0;
//		etar[xijk]=0;
//		}
efor_xijk
}//end of initial_condition3()

void initial_condition4()
{//Flat interface across z-direction between (c=1,eta=-1) and (c=1,eta=1)
double *cr=&c_r[0][0][0];
//double *etar=&eta_r[0][0][0];
int x_i,x_j,x_k,xijk;
 init_homo();
for_xijk	
		if(x_k<N/4)
/*		{
		cr[xijk]=5;
		etar[xijk]=0;
		}
		else 
		*/
		{
		cr[xijk]=8.0;
		eta_r[1][x_i][x_j][x_k]=1;
		}

efor_xijk
}//end of initial_condition4()

void initial_condition5(double r)
{// Circular particle in 2D spherical particle in 3D
double *cr=&c_r[0][0][0];
//double *etar=&eta_r[0][0][0];
int x_i,x_j,x_k,xijk,i;
double C[3];
 init_homo();
//Center of the system
C[0]=(L-1)/2.0;
C[1]=(M-1)/2.0;
C[2]=(N-1)/2.0;

if(L==2 && M!=2)//For2 2D case
{
for_xijk
		if( (x_k-C[2])*(x_k-C[2]) + (x_j-C[1])*(x_j-C[1]) <= r*r )
		{
		cr[xijk]=8.3;
		eta_r[0][x_i][x_j][x_k]=1;
		}
/*		else 
		{
		cr[xijk]=6.2;
		etar[xijk]=0;
		}
		*/
efor_xijk
}//end of 2D case if
else if(L!=2 && M!=2)//For 3D case
{
for_xijk
		if( (x_i-C[0])*(x_i-C[0]) + (x_j-C[1])*(x_j-C[1]) + (x_k-C[2])*(x_k-C[2]) <= r*r )
		{
		eta_r[0][x_i][x_j][x_k]=1.0;
	//	cr[xijk]=1;
/*		for(i=0;i<V;i++){
		eta_r[i][x_i][x_j][x_k]=0;
		}
		*/
		}
/*		else 
		{
		cr[xijk]=6.3;
		etar[xijk]=0;
		}
		*/
efor_xijk
}//end of 3D case if

}//end of initial_condition5()

void initial_condition6(double a, double b)
{//ellipsoidal particle in 2D
double *cr=&c_r[0][0][0];
//double *etar=&eta_r[0][0][0];
int x_i,x_j,x_k,xijk;
double C[3];
//Center of the system
 init_homo();
C[0]=(L-1)/2.0;
C[1]=(M-1)/2.0;
C[2]=(N-1)/2.0;

if(L==2 && M!=2)//For2 2D case
{
for_xijk

		if( (x_k-C[2])*(x_k-C[2])/b/b + (x_j-C[1])*(x_j-C[1])/a/a <= 1 )
/*		{
		cr[xijk]=0;
		etar[xijk]=0;
		}
		else 
		*/
		{
		cr[xijk]=8.3;
		eta_r[0][x_i][x_j][x_k]=1;
		}
efor_xijk
}

}//end of initial_condition6()

void initial_particle(int a, int b)//b is the num of the particle
{
int x_i,x_j,x_k,xijk;
double *cr=&c_r[0][0][0];
//double *etar=&eta_r[0][0][0];
int ip,iv,ix,iy,iz,ih,ix_old,iy_old,iz_old;
long int indx2;
float tp1,tp2,tp3,t1,t2,t3,m1,m2,m3;
float radius;
 init_homo();
        srand(100);
/*for_xijk
    cr[xijk]=6.5;
    etar[xijk]=0.0; 
efor_xijk
*/
//homogeneous configuration of particle

       // Init_HOMO(eta);
        printf("Generation of sphere particles with random distribution size.\n");
        for(ip=0;ip<b;ip++)
        {
        radius=rand()%a>8?(rand()%a):8;
     //   cout<<radius/2<<endl;
      //  iv=0;
            iv=rand()%V;
                ix=rand()%L;
                iy=rand()%M;
                iz=rand()%N;
for_xijk 
                                {

                                        tp1=(float)(x_i-ix);
                                        if(fabs(tp1)>L/2)
                                                tp1=tp1>0?(tp1-L):(tp1+L);
                                        tp2=(float)(x_j-iy);
                                        if(fabs(tp2)>M/2)
                                                tp2=tp2>0?(tp2-M):(tp2+M);
                                        tp3=(float)(x_k-iz);
                                        if(fabs(tp3)>N/2)
                                                tp3=tp3>0?(tp3-N):(tp3+N);
                                                //place H phase at (1/4.1/4,1/4)
                                                                                  //      if(tp1*tp1+tp2*tp2+tp3*tp3>radius*radius/4.0)
                                                                                                                          //      continue;
                                                                                                                                                                                  if(tp1*tp1/radius/radius+tp2*tp2/radius/radius+tp3*tp3/radius/radius<=1.0/4.0)
        {      
                cr[xijk]=0.0;
                eta_r[iv][x_i][x_j][x_k]=1.0;//X_Ni
}

        }       //end of space loop
efor_xijk
     }//end of particle number loop
//        return 0;


}

void init_homo()
{
int x_i,x_j,x_k,xijk;
    for(int i=0;i<V;i++)
	{
for_xijk
   eta_r[i][x_i][x_j][x_k]=0;
   c_r[x_i][x_j][x_k]=0.0;   
efor_xijk
}
}
