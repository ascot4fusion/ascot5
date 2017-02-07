/**
 * @file test_math.c
 * @brief Test program for math functions
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "ascot5.h"
#include "math.h"

double myfunc(double x){
  return sin(10*x)*exp(-x);
}

int main(void) {
  int i;
  srand48(0);

  double accurate=(10.0-(sin(100.0)+10.0*cos(100.0))*exp(-10.0))/101.0;
  double approx=math_simpson(myfunc,0.0,10.0,1.e-10);
  
  printf("accurate result: %20.18g\n",accurate);
  printf("numerical value: %20.18g\n",approx);
  printf("relative error: %20.18g\n",fabs((accurate-approx)/accurate));
  
  
  for(i=0; i < 1000000; i++) {
    real r = math_normal_rand();
    printf("%lf\n", r);
  }
  

  printf("Testing matrix operations\n");
  real matA[6] = {1, 2, 
		  3, 4, 
		  5, 6};/* 3 x 2 matrix*/
  real matB[8] = { 7,  8,  9, 10, 
		  11, 12, 13, 14};/* 2 x 4 matrix*/
  real matC[12];/* 3 x 4 matrix*/
  math_matmul(matA, matB, 3, 2, 4, matC);
  printf("Matrix multiplication\n");
  printf("Calculated    Correct\n");
  printf(" %g %g %g %g   %g %g %g %g\n",matC[0],matC[1], matC[2], matC[3], 39.0, 49.0,  59.0,  69.0);
  printf(" %g %g %g %g   %g %g %g %g\n",matC[4],matC[5], matC[6], matC[7], 54.0, 68.0,  82.0,  96.0);
  printf(" %g %g %g %g   %g %g %g %g\n",matC[8],matC[9],matC[10],matC[11], 69.0, 87.0, 105.0, 123.0);

}
