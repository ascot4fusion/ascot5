#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../math.h"

int main() {
    int n = 30;
    real vec_min = (real) -n;
    real vec_max = (real) n;
    real vec[n];
    math_linspace(vec, vec_min, vec_max, n);

    int success;

    int n_samples = n * 1000;
    printf("Checking %d random values...", n_samples);
    success = 1;
    for(int i = 0; i < n_samples; i++) {
        real sval = vec_min + (vec_max - vec_min)*(real)rand()/RAND_MAX;
        real* res = math_rsearch(sval, vec, n);
        if(!res) {
            printf("\nERROR: Number %.4f not in vector.\n", sval);
            success = 0;
            break;
        }
        else if(*res <= sval && *(res + 1) > sval) {
            /* Do nothing */
        }
        else {
            printf("\nERROR: Number %.4f found between "
                   "%.4f and %.4f = vec[%d] and vec[%d].\n",
                   sval, *res, *(res + 1),
                   (int)(res - vec),
                   (int)(res + 1 - vec));
            success = 0;
            break;
        }
    }
    if (success) printf("\t\t\t\tOK!\n");

    printf("Checking %d first grid points...", n - 1);
    success = 1;
    for(int i = 0; i < n - 1; i++) {
        real* res = math_rsearch(vec[i], vec, n);
        if(!res) {
            printf("\nERROR: Grid point %.4f not found.\n", vec[i]);
            success = 0;
            break;
        }
        else if(res == vec + i) {
            /* Do nothing */
        }
        else {
            printf("\nERROR: Grid point %.4f at wrong place %.4f = vec[%d]\n",
                   vec[i], *res, i);
            success = 0;
            break;
        }
    }
    if(success) printf("\t\t\tOK!\n");

    printf("Checking last grid point...");
    real* res = math_rsearch(vec[n-1], vec, n);
    if(res) {
        printf("\nERROR: Found last grid point %.4f at %.4f = vec[%d].\n",
               vec[n-1], *res, (int)(res - vec));
        success = 0;
    }
    else {
        printf("\t\t\t\tOK!\n");
    }

    printf("Checking two points just outside grid...");
    real tol = 0.0001;
    real sval1 = vec[0] - tol * (vec[1] - vec[0]);
    real sval2 = vec[n-1] + tol * (vec[n-1] - vec[n-2]);
    if(math_rsearch(sval1, vec, n)) {
        printf("\nERROR: Number %.4f found (less than %.4f = vec[0]).\n",
               sval1, vec[0]);
        success = 0;
    }
    if(math_rsearch(sval2, vec, n)) {
        printf("\nERROR: Number %.4f found (more than %.4f = vec[n-1]).\n",
               sval2, vec[n-1]);
        success = 0;
    }
    else {
        printf("\t\tOK!\n");
    }

    return 0;
}
