/**
 * @file test_spline.c
 * @brief Unit test program for spline interpolation convergence
 *
 * Make (compile) and run from ascot5/ folder by:
 *     >> make test_spline
 *     >> ./test_spline
 *
 * @todo Solve the problem that CPU time measurement is sometimes giving one
 *       or a few absurdly large values, e.g. 8.887821e+252
 */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../math.h"
#include "../consts.h"
#include "../spline/spline.h"
#include "../spline/interp.h"
#include <time.h>

/**
 * @brief Identifiers for alternative spline representations
 */
enum representation {
    EXPL = 0, /**< Explicit representation of cubic splines */
    COMP = 1  /**< Compact representation of cubic splines  */
};

/**
 * @brief Identifiers for alternative analytical functions
 */
enum analyticalFunction {
    TRIG = 0, /**< Trigonometric analytical function */
    EXP  = 1  /**< Exponential analytical function */
};

int test_interp1D();
int test_interp2D();
int test_interp3D();

/**
 * Main function for testing spline interpolation convergence
 */
int main() {

    printf("\nStarting spline interpolation convergence unit test.\n");

    /* Initialize parameters */
    int anl_func = TRIG;
    int rep = COMP;
    int n_rnd = 1000000; /* Number of random points for error calculation */
    /* Common number of grid points and interval min and max for all
       coordinates, since we want h_x = h_y = h_z for simple monitoring of
       error scaling. Also, min and max are chosen to give trigonometric
       functions perfect BC matching. */
    int n[5] = {8, 16, 32, 64, 128};
    real min = 0.0;
    real max = 2*CONST_PI;

    /* Error arrays */
    real err1D_df[5*3] = {0.0};
    real err2D_df[5*6] = {0.0};
    real err3D_df[5*10] = {0.0};

    /* Count CPU time used */
    clock_t start = clock();
    clock_t end = clock();
    /* CPU times used for different operations, for each dimensionality and
       resolution. The order within one resolution: [init, eval] */
    double cpu_time_1D[5*2];
    double cpu_time_2D[5*2];
    double cpu_time_3D[5*2];

    /* Loop through resolutions and test all dimensions. Currently available
       boundary conditions: NATURALBC and PERIODICBC. Note that the h for
       PERIODICBCs will be slightly smaller, because of the one extra
       subinterval between min and max. However, the relative significance of
       this fades quickly with increasing n. We cannot tweak the NATURALBC
       max to eliminate this, since max = 2*pi is needed for perfect BC match
       for trigonometric functions also for NATURALBC. */
    for(int i_n = 0; i_n < 5; i_n++) {
        test_interp1D(&err1D_df[i_n*3], anl_func, rep, n_rnd,
                      n[i_n], NATURALBC,
                      min, max,
                      start, end, &cpu_time_1D[i_n*2]);
        test_interp2D(&err2D_df[i_n*6], anl_func, rep, n_rnd,
                      n[i_n], n[i_n], NATURALBC, NATURALBC,
                      min, max,
                      min, max,
                      start, end, &cpu_time_2D[i_n*2]);
        test_interp3D(&err3D_df[i_n*10], anl_func, rep, n_rnd,
                      n[i_n], n[i_n], n[i_n], NATURALBC, PERIODICBC, NATURALBC,
                      min, max,
                      min, max,
                      min, max,
                      start, end, &cpu_time_3D[i_n*2]);
    }

    /* Print results, in Matlab compatible format. If you copy-paste the output
       to Matlab, nice vectors and matrices are written. */
    /* Grid interval h. Be careful when mixing other BCs with PERIODICBC,
       which has one more subinterval. Now we are using the h for NATURALBC,
       even knowing that it does not exactly equal the h for PERIODICBC. */
    printf("\nh = [");
    for(int i = 0; i < 5; i++) {
        printf(" %le", (max-min)/(n[i]-1));
    }
    printf("];\n");
    /* 1D error */
    printf("\nerr1D_df = ...\n");
    for(int i = 0; i < 3; i++) {
        for(int j = 0; j < 5; j++) {
            if(i == 0 && j == 0) {
                printf("[%le", err1D_df[j*3+i]);
            } else {
                printf(" %le", err1D_df[j*3+i]);
            }
        }
        if(i == 2) {
            printf("];\n");
        } else {
            printf("\n");
        }
    }
    /* 2D error */
    printf("\nerr2D_df = ...\n");
    for(int i = 0; i < 6; i++) {
        for(int j = 0; j < 5; j++) {
            if(i == 0 && j == 0) {
                printf("[%le", err2D_df[j*6+i]);
            } else {
                printf(" %le", err2D_df[j*6+i]);
            }
        }
        if(i == 5) {
            printf("];\n");
        } else {
            printf("\n");
        }
    }
    /* 3D error */
    printf("\nerr3D_df = ...\n");
    for(int i = 0; i < 10; i++) {
        for(int j = 0; j < 5; j++) {
            if(i == 0 && j == 0) {
                printf("[%le", err3D_df[j*10+i]);
            } else {
                printf(" %le", err3D_df[j*10+i]);
            }
        }
        if(i == 9) {
            printf("];\n");
        } else {
            printf("\n");
        }
    }
    /* Number of grid points */
    printf("\nn = [");
    for(int i = 0; i < 5; i++) {
        printf(" %d", n[i]);
    }
    printf("];\n");
    /* CPU time */
    printf("\ncpu_time = ...\n");
    for(int i = 0; i < 5; i++) {
        if(i == 0) {
            printf("[%le", cpu_time_1D[i*2]);
        } else {
            printf(" %le", cpu_time_1D[i*2]);
        }
        printf(" %le", cpu_time_1D[i*2+1]);
        printf(" %le", cpu_time_2D[i*2]);
        printf(" %le", cpu_time_2D[i*2+1]);
        printf(" %le", cpu_time_3D[i*2]);
        printf(" %le", cpu_time_3D[i*2+1]);
        if(i == 4) {
            printf("];\n");
        } else {
            printf("\n");
        }
    }

    printf("\nThe above print-out can be copy-pasted to Matlab.\n");
    printf("Ending spline interpolation convergence unit test.\n");
}

/**
 * Function that tests convergence of 1D spline interpolation
 */
int test_interp1D(real* err_df, int anl_func, int rep, int n_rnd,
                  int n_x, int bc_x,
                  real x_min, real x_max,
                  clock_t start, clock_t end, double* cpu_time) {

    /* Initialize variable and test data */
    real* x = (real*) malloc((n_x+1*(bc_x==PERIODICBC))*sizeof(real));
    math_linspace(x, x_min, x_max, n_x+1*(bc_x==PERIODICBC));
    real* f_inp = (real*) malloc(n_x*sizeof(real));
    for(int i = 0; i < n_x; i++) {
        if(anl_func == TRIG) {
            f_inp[i] = sin(x[i]);
        } else if(anl_func == EXP) {
            f_inp[i] = exp(-(pow(x[i] - CONST_PI, 2)));
        }
    }

    /* Construct spline, and time it */
    start = clock();
    real* c;
    interp1D_data str;
    if(rep == EXPL) {
        // TODO: Is this -1 for natural only possible with 1D, bc e.g. in 2D
        // temporarily need also the last one, when finding second step coefs?
        c = (real*) malloc((n_x-1*(bc_x==NATURALBC))*NSIZE_EXPL1D*sizeof(real));
        interp1Dexpl_init_coeff(c, f_inp, n_x, bc_x, x_min, x_max);
        interp1Dexpl_init_spline(&str, c, n_x, bc_x, x_min, x_max);
    } else if(rep == COMP) {
        c = (real*) malloc(n_x*NSIZE_COMP1D*sizeof(real));
        interp1Dcomp_init_coeff(c, f_inp, n_x, bc_x, x_min, x_max);
        interp1Dcomp_init_spline(&str, c, n_x, bc_x, x_min, x_max);
    } else {
        /* Dummy clause to avoid compiler warning */
        c = (real*) malloc(sizeof(real));
    }
    end = clock();
    cpu_time[0] = ((double) (end-start))/CLOCKS_PER_SEC;

    /* Evaluate and compare values at random points to determine error */
    real x_rnd;
    real f_spl; /* Function value */
    real err_f = 0.0;
    real df_anl[3]; /* Function value and derivatives */
    real df_spl[3];
    for(int i = 0; i < n_rnd; i++) {
        /* Draw random point */
        x_rnd = x_min + (real)rand()/(real)RAND_MAX*(x_max-x_min);
        /* Calculate the analytical value and derivatives */
        if(anl_func == TRIG) {
            df_anl[0] = sin(x_rnd);
            df_anl[1] = cos(x_rnd);
            df_anl[2] = -sin(x_rnd);
        } else if(anl_func == EXP) {
            df_anl[0] = exp(-(pow(x_rnd - CONST_PI, 2)));
            df_anl[1] = -2*exp(-(pow(x_rnd - CONST_PI, 2)))*(x_rnd - CONST_PI);
            df_anl[2] = exp(-(pow(x_rnd - CONST_PI, 2)))*
                (4*pow(x_rnd, 2) - 8*CONST_PI*x_rnd + 4*pow(CONST_PI, 2) - 2);
        }
        /* Evaluate spline interpolant, and time it cumulatively */
        if(rep == EXPL) {
            start = clock();
            interp1Dexpl_eval_f(&f_spl, &str, x_rnd);
            interp1Dexpl_eval_df(df_spl, &str, x_rnd);
            end = clock();
        }else if(rep == COMP) {
            start = clock();
            interp1Dcomp_eval_f(&f_spl, &str, x_rnd);
            interp1Dcomp_eval_df(df_spl, &str, x_rnd);
            end = clock();
        }
        cpu_time[1] = cpu_time[1] + ((double) (end-start))/CLOCKS_PER_SEC;
        /* Cumulate error */
        err_f += fabs(df_anl[0]-f_spl);
        for(int j = 0; j < 3; j++) {
            err_df[j] += fabs(df_anl[j]-df_spl[j]);
        }
    }
    /* Average error */
    err_f /= n_rnd; // TODO: Check that f_spl == df_spl[0]!
    for(int j = 0; j < 3; j++) {
        err_df[j] /= n_rnd;
    }

    /* Free allocated memory */
    free(x);
    free(f_inp);
    free(c);

    return 0;
}

/**
 * Function that tests convergence of 2D spline interpolation
 */
int test_interp2D(real* err_df, int anl_func, int rep, int n_rnd,
                  int n_x, int n_y, int bc_x, int bc_y,
                  real x_min, real x_max,
                  real y_min, real y_max,
                  clock_t start, clock_t end, double* cpu_time) {

    /* Initialize variables and test data */
    real* x = (real*) malloc((n_x+1*(bc_x==PERIODICBC))*sizeof(real));
    math_linspace(x, x_min, x_max, n_x+1*(bc_x==PERIODICBC));
    real* y = (real*) malloc((n_y+1*(bc_y==PERIODICBC))*sizeof(real));
    math_linspace(y, y_min, y_max, n_y+1*(bc_y==PERIODICBC));
    real* f_inp = (real*) malloc(n_y*n_x*sizeof(real));
    for(int i_y = 0; i_y < n_y; i_y++) {
        for(int i_x = 0; i_x < n_x; i_x++) {
            if(anl_func == TRIG) {
                f_inp[i_y*n_x+i_x] = sin(x[i_x])*sin(y[i_y]);
            } else if(anl_func == EXP) {
                f_inp[i_y*n_x+i_x] = exp(-(pow(x[i_x]-CONST_PI, 2)))
                                    *exp(-(pow(y[i_y]-CONST_PI, 2)));
            }
        }
    }

    /* Construct spline, and time it */
    start = clock();
    real* c;
    interp2D_data str;
    if(rep == EXPL) {
        c = (real*) malloc(n_y*n_x*NSIZE_EXPL2D*sizeof(real));
        interp2Dexpl_init_coeff(c, f_inp, n_x, n_y, bc_x, bc_y,
                                x_min, x_max, y_min, y_max);
        interp2Dexpl_init_spline(&str, c, n_x, n_y, bc_x, bc_y,
                                 x_min, x_max, y_min, y_max);
    } else if(rep == COMP) {
        c = (real*) malloc(n_y*n_x*NSIZE_COMP2D*sizeof(real));
        interp2Dcomp_init_coeff(c, f_inp, n_x, n_y, bc_x, bc_y,
                                x_min, x_max, y_min, y_max);
        interp2Dcomp_init_spline(&str, c, n_x, n_y, bc_x, bc_y,
                                 x_min, x_max, y_min, y_max);
    } else {
        /* Dummy clause to avoid compiler warning */
        c = (real*) malloc(sizeof(real));
    }
    end = clock();
    cpu_time[0] = ((double) (end-start))/CLOCKS_PER_SEC;

    /* Evaluate and compare values at random points to determine error */
    real x_rnd;
    real y_rnd;
    real f_spl; /* Function value */
    real err_f = 0.0;
    real df_anl[6]; /* Function value and derivatives */
    real df_spl[6];
    for(int i = 0; i < n_rnd; i++) {
        /* Draw random point */
        x_rnd = x_min + (real)rand()/(real)RAND_MAX*(x_max-x_min);
        y_rnd = y_min + (real)rand()/(real)RAND_MAX*(y_max-y_min);
        /* Calculate the analytical value and derivatives */
        if(anl_func == TRIG) {
            df_anl[0] = sin(x_rnd)*sin(y_rnd);
            df_anl[1] = cos(x_rnd)*sin(y_rnd);
            df_anl[2] = sin(x_rnd)*cos(y_rnd);
            df_anl[3] = -sin(x_rnd)*sin(y_rnd);
            df_anl[4] = sin(x_rnd)*(-sin(y_rnd));
            df_anl[5] = cos(x_rnd)*cos(y_rnd);
        } else if(anl_func == EXP) {
            df_anl[0] = exp(-pow(x_rnd-CONST_PI, 2)
                            -pow(y_rnd-CONST_PI, 2));
            df_anl[1] = -2*(x_rnd - CONST_PI)*exp(-pow(x_rnd-CONST_PI, 2)
                                              -pow(y_rnd-CONST_PI, 2));
            df_anl[2] = -2*(y_rnd-CONST_PI)*exp(-pow(x_rnd-CONST_PI, 2)
                                            -pow(y_rnd-CONST_PI, 2));
            df_anl[3] = (4*pow(x_rnd, 2) - 8*CONST_PI*x_rnd
                         + 4*pow(CONST_PI, 2) - 2)
                *exp(-pow(x_rnd - CONST_PI, 2)
                     -pow(y_rnd - CONST_PI, 2));
            df_anl[4] = (4*pow(y_rnd, 2) - 8*CONST_PI*y_rnd
                         + 4*pow(CONST_PI, 2) - 2)
                *exp(-pow(x_rnd-CONST_PI, 2)
                     -pow(y_rnd-CONST_PI, 2));
            df_anl[5] = 4*(x_rnd-CONST_PI)*(y_rnd-CONST_PI)*
                exp(-pow(x_rnd-CONST_PI, 2)
                    -pow(y_rnd-CONST_PI, 2));
        }
        /* Evaluate spline interpolant, and time it cumulatively */
        if(rep == EXPL) {
            start = clock();
            interp2Dexpl_eval_f(&f_spl, &str, x_rnd, y_rnd);
            interp2Dexpl_eval_df(df_spl, &str, x_rnd, y_rnd);
            end = clock();
        }else if(rep == COMP) {
            start = clock();
            interp2Dcomp_eval_f(&f_spl, &str, x_rnd, y_rnd);
            interp2Dcomp_eval_df(df_spl, &str, x_rnd, y_rnd);
            end = clock();
        }
        cpu_time[1] = cpu_time[1] + ((double) (end-start))/CLOCKS_PER_SEC;
        /* Cumulate error */
        err_f += fabs(df_anl[0]-f_spl);
        for(int j = 0; j < 6; j++) {
            err_df[j] += fabs(df_anl[j]-df_spl[j]);
        }
    }
    /* Average and print error */
    err_f /= n_rnd;
    for(int j = 0; j < 6; j++) {
        err_df[j] /= n_rnd;
    }

    /* Free allocated memory */
    free(x);
    free(y);
    free(f_inp);
    free(c);

    return 0;
}

/**
 * Function that tests convergence of 3D spline interpolation
 */
int test_interp3D(real* err_df, int anl_func, int rep, int n_rnd,
                  int n_x, int n_y, int n_z, int bc_x, int bc_y, int bc_z,
                  real x_min, real x_max,
                  real y_min, real y_max,
                  real z_min, real z_max,
                  clock_t start, clock_t end, double* cpu_time) {

    /* Initialize variables and test data */
    real* x = (real*) malloc((n_x+1*(bc_x==PERIODICBC))*sizeof(real));
    math_linspace(x, x_min, x_max, n_x+1*(bc_x==PERIODICBC));
    real* y = (real*) malloc((n_y+1*(bc_y==PERIODICBC))*sizeof(real));
    math_linspace(y, y_min, y_max, n_y+1*(bc_y==PERIODICBC));
    real* z = (real*) malloc((n_z+1*(bc_z==PERIODICBC))*sizeof(real));
    math_linspace(z, z_min, z_max, n_z+1*(bc_z==PERIODICBC));
    real* f_inp = (real*) malloc(n_z*n_y*n_x*sizeof(real));
    for(int i_z = 0; i_z < n_z; i_z++) {
        for(int i_y = 0; i_y < n_y; i_y++) {
            for(int i_x = 0; i_x < n_x; i_x++) {
                if(anl_func == TRIG) {
                    f_inp[i_z*n_y*n_x+i_y*n_x+i_x] = sin(x[i_x])
                                                    *sin(y[i_y])
                                                    *sin(z[i_z]);
                } else if(anl_func == EXP) {
                    f_inp[i_z*n_y*n_x+i_y*n_x+i_x] =
                        exp(-pow(x[i_x]-CONST_PI, 2)
                            -pow(y[i_y]-CONST_PI, 2)
                            -pow(z[i_z]-CONST_PI, 2));
                }
            }
        }
    }

    /* Construct spline, and time it */
    start = clock();
    real* c;
    interp3D_data str;
    if(rep == EXPL) {
        c = (real*) malloc(n_z*n_y*n_x*NSIZE_EXPL3D*sizeof(real));
        interp3Dexpl_init_coeff(c, f_inp, n_x, n_y, n_z, bc_x, bc_y, bc_z,
                                x_min, x_max, y_min, y_max, z_min, z_max);
        interp3Dexpl_init_spline(&str, c, n_x, n_y, n_z, bc_x, bc_y, bc_z,
                                 x_min, x_max, y_min, y_max, z_min, z_max);
    } else if(rep == COMP) {
        c = (real*) malloc(n_z*n_y*n_x*NSIZE_COMP3D*sizeof(real));
        interp3Dcomp_init_coeff(c, f_inp, n_x, n_y, n_z, bc_x, bc_y, bc_z,
                                x_min, x_max, y_min, y_max, z_min, z_max);
        interp3Dcomp_init_spline(&str, c, n_x, n_y, n_z, bc_x, bc_y, bc_z,
                                 x_min, x_max, y_min, y_max, z_min, z_max);
    } else {
        /* Dummy clause to avoid compiler warning */
        c = (real*) malloc(sizeof(real));
    }
    end = clock();
    cpu_time[0] = ((double) (end-start))/CLOCKS_PER_SEC;

    /* Evaluate and compare values at random points to determine error */
    real x_rnd;
    real y_rnd;
    real z_rnd;
    real f_spl; /* Function value */
    real err_f = 0.0;
    real df_anl[10]; /* Function value and derivatives */
    real df_spl[10];
    for(int i = 0; i < n_rnd; i++) {
        /* Draw random point */
        x_rnd = x_min + (real)rand()/(real)RAND_MAX*(x_max-x_min);
        y_rnd = y_min + (real)rand()/(real)RAND_MAX*(y_max-y_min);
        z_rnd = z_min +  (real)rand()/(real)RAND_MAX*(z_max-z_min);
        /* Calculate the analytical value and derivatives */
        if(anl_func == TRIG) {
            df_anl[0] = sin(x_rnd)*sin(y_rnd)*sin(z_rnd);
            df_anl[1] = cos(x_rnd)*sin(y_rnd)*sin(z_rnd);
            df_anl[2] = sin(x_rnd)*cos(y_rnd)*sin(z_rnd);
            df_anl[3] = sin(x_rnd)*sin(y_rnd)*cos(z_rnd);
            df_anl[4] = -sin(x_rnd)*sin(y_rnd)*sin(z_rnd);
            df_anl[5] = -sin(x_rnd)*sin(y_rnd)*sin(z_rnd);
            df_anl[6] = -sin(x_rnd)*sin(y_rnd)*sin(z_rnd);
            df_anl[7] = cos(x_rnd)*cos(y_rnd)*sin(z_rnd);
            df_anl[8] = cos(x_rnd)*sin(y_rnd)*cos(z_rnd);
            df_anl[9] = sin(x_rnd)*cos(y_rnd)*cos(z_rnd);
        } else if(anl_func == EXP) {
            df_anl[0] = exp(-pow(x_rnd-CONST_PI, 2)
                            -pow(y_rnd-CONST_PI, 2)
                            -pow(z_rnd-CONST_PI, 2));
            df_anl[1] = -2*exp(-pow(x_rnd-CONST_PI, 2)
                               -pow(y_rnd-CONST_PI, 2)
                               -pow(z_rnd-CONST_PI, 2))*(x_rnd-CONST_PI);
            df_anl[2] = -2*exp(-pow(x_rnd-CONST_PI, 2)
                               -pow(y_rnd-CONST_PI, 2)
                               -pow(z_rnd-CONST_PI, 2))*(y_rnd-CONST_PI);
            df_anl[3] = -2*exp(-pow(x_rnd-CONST_PI, 2)
                               -pow(y_rnd-CONST_PI, 2)
                               -pow(z_rnd-CONST_PI, 2))*(z_rnd-CONST_PI);
            df_anl[4] = -2*exp(-pow(x_rnd-CONST_PI, 2)
                               -pow(y_rnd-CONST_PI, 2)
                               -pow(z_rnd-CONST_PI, 2))
                +4*exp(-pow(x_rnd-CONST_PI, 2)
                       -pow(y_rnd-CONST_PI, 2)
                       -pow(z_rnd-CONST_PI, 2))*pow(x_rnd-CONST_PI, 2);
            df_anl[5] = -2*exp(-pow(x_rnd - CONST_PI, 2)
                               -pow(y_rnd-CONST_PI, 2)
                               -pow(z_rnd-CONST_PI, 2))
                +4*exp(-pow(x_rnd-CONST_PI, 2)
                       -pow(y_rnd-CONST_PI, 2)
                       -pow(z_rnd-CONST_PI, 2))*pow(y_rnd-CONST_PI, 2);
            df_anl[6] = -2*exp(-pow(x_rnd-CONST_PI, 2)
                               -pow(y_rnd-CONST_PI, 2)
                               -pow(z_rnd-CONST_PI, 2))
                +4*exp(-pow(x_rnd-CONST_PI, 2)
                       -pow(y_rnd-CONST_PI, 2)
                       -pow(z_rnd-CONST_PI, 2))*pow(z_rnd-CONST_PI, 2);
            df_anl[7] = 4*exp(-pow(x_rnd-CONST_PI, 2)
                              -pow(y_rnd-CONST_PI, 2)
                              -pow(z_rnd-CONST_PI, 2))
                *(x_rnd-CONST_PI)*(y_rnd-CONST_PI);
            df_anl[8] = 4*exp(-pow(x_rnd-CONST_PI, 2)
                              -pow(y_rnd-CONST_PI, 2)
                              -pow(z_rnd-CONST_PI, 2))
                *(x_rnd-CONST_PI)*(z_rnd-CONST_PI);
            df_anl[9] = 4*exp(-pow(x_rnd-CONST_PI, 2)
                              -pow(y_rnd-CONST_PI, 2)
                              -pow(z_rnd-CONST_PI, 2))
                *(y_rnd-CONST_PI)*(z_rnd-CONST_PI);
        }
        /* Evaluate spline interpolant, and time it cumulatively */
        if(rep == EXPL) {
            start = clock();
            interp3Dexpl_eval_f(&f_spl, &str, x_rnd, y_rnd, z_rnd);
            interp3Dexpl_eval_df(df_spl, &str, x_rnd, y_rnd, z_rnd);
            end = clock();
        }else if(rep == COMP) {
            start = clock();
            interp3Dcomp_eval_f(&f_spl, &str, x_rnd, y_rnd, z_rnd);
            interp3Dcomp_eval_df(df_spl, &str, x_rnd, y_rnd, z_rnd);
            end = clock();
        }
        cpu_time[1] = cpu_time[1] + ((double) (end-start))/CLOCKS_PER_SEC;
        /* Cumulate error */
        err_f += fabs(df_anl[0]-f_spl);
        for(int j = 0; j < 10; j++) {
            err_df[j] += fabs(df_anl[j]-df_spl[j]);
        }
    }
    /* Average and print error */
    err_f /= n_rnd;
    for(int j = 0; j < 10; j++) {
        err_df[j] /= n_rnd;
    }

    /* Free allocated memory */
    free(x);
    free(y);
    free(z);
    free(f_inp);
    free(c);

    return 0;
}
