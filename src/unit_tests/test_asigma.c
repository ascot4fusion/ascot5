/**
 * @file test_asigma.c
 * @brief Test program for asigma module functions
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../asigma.h"
#include "../math.h"

int main() {
    printf("Starting atomic sigma unit test\n");
    int err = 0;

    /* Initialize atomic sigma */
    printf("Initializing atomic sigma\n");
    int N_reac = 1;
    int* z_1 = malloc(N_reac*sizeof(int));
    int* a_1 = malloc(N_reac*sizeof(int));
    real* m_1 = malloc(N_reac*sizeof(real));
    int* z_2 = malloc(N_reac*sizeof(int));
    int* a_2 = malloc(N_reac*sizeof(int));
    real* m_2 = malloc(N_reac*sizeof(real));
    int* reac_type = malloc(N_reac*sizeof(int));
    z_1[0] = 1;
    a_1[0] = 1;
    m_1[0] = 1.6726219*1e-27;
    z_2[0] = 1;
    a_2[0] = 1;
    m_2[0] = 1.6726219*1e-27;
    reac_type[0] = 6;
    asigma_data data;
    err = asigma_init(N_reac, z_1, a_1, m_1,
                      z_2, a_2, m_2,
                      reac_type, &data);

    /* Initialize abscissae for test grid */
    printf("Calculating sigmav test grids\n");
    real E_min = 1.0*1e+2;
    real E_max = 1.5*1e+5;
    int N_E = 150;
    real* E = malloc(N_E*sizeof(real));
    if(E_min >= data.sigmav[0].x_min && E_max <= data.sigmav[0].x_max) {
        math_linspace(E, E_min, E_max, N_E);
    } else {
        printf("Error: Energy span goes outside spline domain\n");
    }
    real T_min = 1.0*1e+2;
    real T_max = 1.0*1e+4;
    int N_T = 3;
    real T[N_T];
    if(T_min >= data.sigmav[0].y_min && T_max <= data.sigmav[0].y_max) {
        T[0] = T_min;
        T[1] = 1.0*1e+3;
        T[N_T-1] = T_max;
    } else {
        printf("Error: Energy span goes outside spline domain\n");
    }
    real* sigmav = malloc(N_reac*N_T*N_E*sizeof(real));
    real n_e = 0.0; /* Dummy value for electron density */
    real n_i[1];
    n_i[0] = n_e;
    /* Evaluate rate coefficients at test grid points for each reaction */
    for(int i_reac = 0; i_reac < N_reac; i_reac++) {
        for(int i_T = 0; i_T < N_T; i_T++) {
            for(int i_E = 0; i_E < N_E; i_E++) {
                sigmav[i_reac*N_T*N_E+i_T*N_E+i_E] =
                    asigma_eval_sigmav(z_1[i_reac], a_1[i_reac],
                                       z_2[i_reac], a_2[i_reac],
                                       reac_type[i_reac], &data,
                                       E[i_E], T[i_T], T[i_T], T[i_T],
                                       n_e, n_i);
            }
        }
    }

    /* Write test grid to file */
    char* file_name = "/u/97/ollusp1/unix/Documents/git/ascot5/"
        "test_asigma_output.dat";
    FILE* f_ptr;
    f_ptr = fopen(file_name, "w");
    if(f_ptr == NULL) {
        printf("Error: Error while opening file\n");
        exit(EXIT_FAILURE);
    }
    err = fprintf(f_ptr, "%d %d %d\n", N_E, N_T, N_reac);
    for(int i_E = 0; i_E < N_E; i_E++) {
        err = fprintf(f_ptr, "%le ", E[i_E]);
    }
    err = fprintf(f_ptr, "\n");
    for(int i_T = 0; i_T < N_T; i_T++) {
        err = fprintf(f_ptr, "%le ", T[i_T]);
    }
    err = fprintf(f_ptr, "\n");
    for(int i_reac = 0; i_reac < N_reac; i_reac++) {
        for(int i_T = 0; i_T < N_T; i_T++) {
            for(int i_E = 0; i_E < N_E; i_E++) {
                err = fprintf(f_ptr,"%le ",sigmav[i_reac*N_T*N_E+i_T*N_E+i_E]);
            }
        }
    }
    fclose(f_ptr);

    /* Free memory */
    free(z_1);
    free(a_1);
    free(z_2);
    free(a_2);
    free(reac_type);

    printf("Finishing atomic sigma unit test\n");
    return err;
}
