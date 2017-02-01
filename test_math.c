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
}
