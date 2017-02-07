/**
 * @file math.c
 * @brief Mathematical utility functions
 *
 * @todo Implement as macros
 */

#include <math.h>
#include <stdlib.h>
#include "ascot5.h"
#include "math.h"

/**
 * @brief Return unit vector of a given 3D vector
 *
 * @param vec      3-length array
 * @param vec_unit corresponding unit vector
 */
void math_unit(real* vec, real* vec_unit) {
    real n = math_norm(vec);
    vec_unit[0] = vec[0]/n;
    vec_unit[1] = vec[1]/n;
    vec_unit[2] = vec[2]/n;
}

/**
 * @brief Convert cartesian to cylindrical coordinates 
 *
 * @param xyz input xyz coordinates in a 3-length array
 * @param rpz output rpz coordinates in a 3-length array
 */
void math_xyz2rpz(real* xyz, real* rpz) {
    rpz[0] = sqrt(xyz[0]*xyz[0]+xyz[1]*xyz[1]);
    rpz[1] = atan2(xyz[1], xyz[0]);
    rpz[2] = xyz[2];
}

/**
 * @brief Convert cylindrical to cartesian coordinates 
 *
 * @param rpz input rpz coordinates in a 3-length array
 * @param xyz output xyz coordinates in a 3-length array
 */
void math_rpz2xyz(real* rpz, real* xyz) {
    xyz[0] = rpz[0] * cos(rpz[1]);
    xyz[1] = rpz[0] * sin(rpz[1]);
    xyz[2] = rpz[2];
}

/**
 * @brief Convert a vector from cylindrical to cartesian coordinates
 *
 * This function converts a vector located at angle phi from cylindrical
 * to cartesian coordinates.
 *
 * @param rpz input rpz coordinates in a 3-length array
 * @param xyz output xyz coordinates in a 3-length array
 * @param phi phi coordinate of the vector
 *
 * @bug Does not work with SIMD (unknown reason)
 */
void math_vec_rpz2xyz(real* rpz, real* xyz, real phi) {
    xyz[0] = rpz[0] * cos(phi) - rpz[1] * sin(phi);
    xyz[1] = rpz[0] * sin(phi) + rpz[1] * cos(phi);
    xyz[2] = rpz[2];
}

/**
 * @brief Convert a vector from cartesian to cylindrical coordinates
 *
 * This function converts a vector located at angle phi from cartesian
 * to cylindrical coordinates.
 *
 * @param xyz input xyz coordinates in a 3-length array
 * @param rpz output rpz coordinates in a 3-length array
 * @param phi phi coordinate of the vector
 *
 * @bug Probably does not work with SIMD (unknown reason) (see vec_rpz2xyz) 
 */
void math_vec_xyz2rpz(real* xyz, real* rpz, real phi) {
    rpz[0] = xyz[0] * cos(phi) + xyz[1] * sin(phi);
    rpz[1] = -xyz[0] * sin(phi) + xyz[1] * cos(phi);
    rpz[2] = xyz[2];
}

/**
 * @brief Matrix multiplication
 *
 * This function performs matrix multiplication. The matrices are given as arrays,
 * the indexing being row-wise.
 *
 * @param matA input array representing a d1 x d2 matrix
 * @param matB input array representing a d2 x d3 matrix
 * @param d1   input scalar dimension (columns in matA and matC)
 * @param d2   input scalar dimension (rows in matA and columns in matB)
 * @param d3   input scalar dimension (rows in matB and matC)
 * @param matC output array representing a d1 x d3 matrix
 */
void math_matmul(real* matA, real* matB, int d1, int d2, int d3, real* matC) {
  int i;
  int j;
  int k;
  real sum;
  for (i = 0; i < d1; i=i+1) {
    for (j = 0; j < d3; j=j+1) {
      sum = 0.0;
      for (k = 0; k < d2; k=k+1){
	sum = sum + matA[k * d1 + i]*matB[j * d2 + k];
      }
      matC[i * d3 + j] = sum;
    }
  }
}

/**
 * @brief Generate normally distributed random numbers
 *
 * This function generates normally distributed random numbers with expected
 * value of 0 and variance of 1 using the polar method [1]. The rand48 random
 * number generator should be initialized with srand48() before calling this 
 * function.
 *
 * [1] G. Marsaglia, T. A. Bray. A convenient method for generating normal  
 *     variables. Siam Review 6.3 (1964): 260-264.
 *     http://epubs.siam.org/doi/abs/10.1137/1006063:w
 *
 * @todo Currently only one of the generated numbers is returned; a second
 *       independent variate is Y = v2 * sqrt(-2*log(s) / s)
 * @todo Try other random number generators such as those in Intel MKL
 */
real math_normal_rand(void) {
    real X;
    real v1, v2, s;
    do {
        v1 = drand48() * 2 - 1;
        v2 = drand48() * 2 - 1;
        s = v1*v1 + v2*v2;
    } while(s >= 1);

    X = v1 * sqrt(-2*log(s) / s);
        
    return X;
}

int math_ipow(int a, int p) {
    int i;
    if(p == 0) {
        return 0;
    } else {
        int pow = a;
        for(i = 0; i < p-1; i++) {
            pow *= a;
        }
        return pow;
    }
}

/**
 * @brief Helper routine for "math_simpson"
 */
double math_simpson_helper(double (*f)(double), double a, double b, double eps, double S, double fa, double fb, double fc, int bottom) {
  double c = (a + b)/2, h = b - a;
  double d = (a + c)/2, e = (c + b)/2;
  double fd = f(d), fe = f(e);                                               
  double Sleft = (h/12)*(fa + 4*fd + fc);
  double Sright = (h/12)*(fc + 4*fe + fb);
  double S2 = Sleft + Sright;
  //printf("[a,b]=[%lf,%lf], (S2-S)/S=%lf,eps=%lf\n",a,b,fabs((S2-S)/S),eps);
  if (bottom <= 0 || fabs(S2 - S) <= eps*fabs(S)) {
    return  S2 + (S2 - S)/15;      
  }
  return math_simpson_helper(f, a, c, eps, Sleft,  fa, fc, fd, bottom-1)
    +math_simpson_helper(f, c, b, eps, Sright, fc, fb, fe, bottom-1);                     
}         

/**
 * @brief Adaptive Simpsons rule for integral
 *
 * This function uses recursive splitting of the interval
 * until either maximum number of intervals is reached or
 * given relative error tolerance is obtained
 *
 * @param f: pointer to the integrand function (of one variable)
 * @param a: lower limit for the interval
 * @param b: upper limit for the interval
 * @param eps: absolute tolerance 
 * @param maxDepth: maximum number of splits for the interval
 *
 * @bug Probably does not work with SIMD 
 */
double math_simpson(double (*f)(double), double a, double b, double eps) {
  double c = (a + b)/2, h = b - a;                                             
  double fa = f(a), fb = f(b), fc = f(c);                                       
  double S = (h/6)*(fa + 4*fc + fb); 
  return math_simpson_helper(f, a, b, eps, S, fa, fb, fc, math_maxSimpsonDepth);        
}                                                                               

void math_linspace(real* vec, real a, real b, int n) {
    if(n == 1) {
        vec[0] = b;
    } else {
        real d = (b-a)/(n-1);
        int i;
        for(i = 0; i < n; i++) {
            vec[i] = a+i*d;
        }
    }
}
