#include <stdio.h>
#include "inhom_v3.h"
voigt A[5];
int main()
{
V6_loop{A[2][mi]=2.0;}		
voigt s0;
double *s=&A[2][0];
#pragma acc parallel loop  
V6_loop{s0[mi]=s[mi];}
printf("s0[3]=%lf\n",s0[3]);
return 0;
}
