/**
 * @file test_math.c
 * @brief Test program for math functions
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "../ascot5.h"
#include "../math.h"

/**
 * Evaluate some mathematical function to be used for testing purposes
 */
double myfunc(double x){
  return sin(10*x)*exp(-x);
}

/**
 * Main function for the test program.
 */
int main(void) {

  int i;

  double accurate=(10.0-(sin(100.0)+10.0*cos(100.0))*exp(-10.0))/101.0;
  double approx=math_simpson(myfunc,0.0,10.0,1.e-10);


  printf("Checking simpson rule.\n");
  printf("accurate result: %20.18g\n",accurate);
  printf("numerical value: %20.18g\n",approx);
  printf("relative error: %20.18g\n",fabs((accurate-approx)/accurate));
  if ( fabs((accurate-approx)/accurate) > 1.0e-14){
	  printf("Fail!\n");
	  return 1;
  }else{
	  printf("ok!\n");
  }

  /*
  srand48(0);
  for(i=0; i < 1000000; i++) {
    real r = math_normal_rand();
    printf("%lf\n", r);
  }
  */

  printf("Testing matrix operations\n");
  real matA[6] = {1, 2,
                  3, 4,
                  5, 6};/* 3 x 2 matrix*/
  real matB[8] = { 7, 8, 9, 10,
                   11, 12, 13, 14};/* 2 x 4 matrix*/
  real matC[12];/* 3 x 4 matrix*/
  math_matmul(matA, matB, 3, 2, 4, matC);
  printf("Matrix multiplication\n");
  printf("Calculated    Correct\n");
  printf(" %g %g %g %g   %g %g %g %g\n",matC[0],matC[1], matC[2], matC[3], 39.0, 49.0,  59.0,  69.0);
  printf(" %g %g %g %g   %g %g %g %g\n",matC[4],matC[5], matC[6], matC[7], 54.0, 68.0,  82.0,  96.0);
  printf(" %g %g %g %g   %g %g %g %g\n",matC[8],matC[9],matC[10],matC[11], 69.0, 87.0, 105.0, 123.0);

  if (    matC[ 0] != 39.0 ||  matC[ 1] != 49.0 || matC[ 2] != 59.0 || matC[ 3] != 69.0 ||
		  matC[ 4] != 54.0 ||  matC[ 5] != 68.0 || matC[ 6] != 82.0 || matC[ 7] != 96.0 ||
		  matC[ 8] != 69.0 ||  matC[ 9] != 87.0 || matC[10] !=105.0 || matC[11] != 123.0 ){
	  printf("Fail!\n");
	  return 1;
  }else{
	  printf("ok!\n");
  }

  printf("Testing 3x3 matrix determinant\n");
  real matdeta[9],det;
  for(i=0;i<9;i++){
	  matdeta[i] = (float) i+1;
  }
  matdeta[8] = 10.0;

  det = math_determinant3x3(matdeta[0],matdeta[1],matdeta[2],matdeta[3],matdeta[4],matdeta[5],matdeta[6],matdeta[7],matdeta[8] );
  printf("det( [1 2 3; 4 5 6; 7 8 10] ) == %f  (correct: -3)\n",det);
  if ( fabs(det+3) > 1.0e-14){
	  return 1;
  }else{
	  printf("ok!\n");
  }
  return 0;
}
