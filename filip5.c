/**
 * @file filip5.c
 * @brief FIeld LIne Poincareplotter
 *
 * This program traces field lines and outputs data for creating poincare plots.
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include "ascot5.h"
#include "wall.h"
#include "B_field.h"
#include "simulate.h"
#include "filip5.h"
#include "math.h"

#pragma omp declare target
void simulate_fieldline(int id, int n, fieldline* l, sim_offload_data sim,
                        real* B_offload_array, real* wall_offload_array);
void simulate_fieldline_rk4(int id, int n, fieldline* l, sim_offload_data sim,
                        real* B_offload_array, real* wall_offload_array);
#pragma omp end declare target

int main(int argc, char** argv) {
    int mpi_rank, mpi_size;

    #ifdef MPI
    MPI_Status status;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
    #else
    mpi_rank = 0;
    mpi_size = 1;
    #endif

    #if VERBOSE >= 1
    printf("Initialized MPI, rank %d, size %d.\n", mpi_rank, mpi_size);
    #endif

    /* Prepare simulation parameters and data for offload */
    sim_offload_data sim;
    real* B_offload_array;
    real* wall_offload_array;
    B_field_init_offload(&sim.B_offload_data, &B_offload_array);
    wall_init_offload(&sim.wall_offload_data, &wall_offload_array);

    int n_r = 1;
    int n_phi = 1;
    int n = 2*n_r*n_phi;
    fieldline* lines = (fieldline*) malloc(n*sizeof(fieldline));
    real rmin = 7;
    real rmax = 7.5;
    real z0 = 0.6;
    int i, j;
    for(i = 0; i < n_r; i++) {
        for(j = 0; j < n_phi; j++) {
            lines[i*n_r+j].r = rmin + (rmax-rmin) * ((real) i) / ((real) n_r);
            lines[i*n_r+j].phi = ((real) j) / ((real) n_phi) * 2 * math_pi;
            lines[i*n_r+j].z = z0;
            lines[i*n_r+j].dir = 1;
            lines[i*n_r+j].running = 1;
            lines[i*n_r+j+n_phi*n_r].r = lines[i*n_r+j].r;
            lines[i*n_r+j+n_phi*n_r].phi = lines[i*n_r+j].phi;
            lines[i*n_r+j+n_phi*n_r].z = lines[i*n_r+j].z;
            lines[i*n_r+j+n_phi*n_r].dir = -1;
            lines[i*n_r+j+n_phi*n_r].running = 1;
        }
    }

    int start_index = mpi_rank * (n / mpi_size);
    lines += start_index;

    if(mpi_rank == mpi_size-1) {
        n = n - mpi_rank * (n / mpi_size);
    }
    else {
        n = n / mpi_size;
    }

    #ifndef NOTARGET
    int n_mic = n / 2;
    int n_host = 0;
    #else
    int n_mic = 0;
    int n_host = n;
    #endif

    double mic0_start, mic0_end, mic1_start, mic1_end, host_start, host_end;

    fflush(stdout);
    omp_set_nested(1);
    #pragma omp parallel sections num_threads(3)
    {
    #ifndef NOTARGET
    #pragma omp section
    {
        mic0_start = omp_get_wtime();
        #pragma omp target device(0) map( \
            p[0:n_mic], \
            B_offload_array[0:sim.B_offload_data.offload_array_length], \
            wall_offload_array[0:sim.wall_offload_data.offload_array_length] \
         )
        simulate_fieldline(1, n_mic, lines, sim, B_offload_array,
                           wall_offload_array);
        mic0_end = omp_get_wtime();
    }

    #pragma omp section
    {
        mic1_start = omp_get_wtime();
        #pragma omp target device(1) map( \
            p[n_mic:n_mic], \
            B_offload_array[0:sim.B_offload_data.offload_array_length], \
            wall_offload_array[0:sim.wall_offload_data.offload_array_length] \
         )
        simulate_fieldline(2, n_mic, lines+n_mic, sim, B_offload_array,
                           wall_offload_array);
        mic1_end = omp_get_wtime();
    } 
    #else
    #pragma omp section
    {
        host_start = omp_get_wtime();
        simulate_fieldline(0, n_host, lines+2*n_mic, sim, B_offload_array,
                           wall_offload_array);
        host_end = omp_get_wtime();
    }
    #endif
    }
    /* Code excution returns to host. */

    #if VERBOSE >= 1
        printf("mic0 %lf s, mic1 %lf s, host %lf s\n", mic0_end-mic0_start,
               mic1_end-mic1_start, host_end-host_start);
    #endif

    B_field_free_offload(&sim.B_offload_data, &B_offload_array);
    wall_free_offload(&sim.wall_offload_data, &wall_offload_array);
    free(lines);

    return 0;
}

void simulate_fieldline(int id, int n, fieldline* l, sim_offload_data sim,
                        real* B_offload_array, real* wall_offload_array) {
    wall_data wall_data;
    wall_init(&wall_data, &sim.wall_offload_data, wall_offload_array);
    B_field_data B_data;
    B_field_init(&B_data, &sim.B_offload_data, B_offload_array);

    FILE* f_tor = fopen("filip5.poincareTor","w");
    FILE* f_pol = fopen("filip5.poincarePol","w");
    FILE* f_sum = fopen("filip5.summary","w");
    FILE* f_debug = fopen("filip5.debug","w");

    real tol = 1e-8;
    real h0 = 1e-2;
    real axis_r = B_field_get_axis_r(&B_data);
    real axis_z = B_field_get_axis_z(&B_data);
    double tt_start = omp_get_wtime();
    int i;
//    #pragma omp parallel for
    for(i = 0; i < n; i++) {
        l[i].ini_r = l[i].r;
        l[i].ini_phi = l[i].phi;
        l[i].ini_z = l[i].z;

        real h = h0 * l[i].dir;

        int orbits_pol = 0;
        int orbits_tor = 0;

        while(l[i].running && l[i].length < 200000 && orbits_pol < 5000
              && orbits_tor < 5000) {

            if(yerr[0][0] > tol || yerr[1][0] > tol || yerr[2][0] > tol) {
                h *= 0.5;
            } else {
                h *= 1.1;
                l[i].r = y[0][0];
                l[i].phi = y[1][0];
                l[i].z = y[2][0];
                l[i].length += sqrt((y[0][0]-yprev[0][0])*(y[0][0]-yprev[0][0])
                   + y[0][0]*(y[1][0]-yprev[1][0])*y[0][0]*(y[1][0]-yprev[1][0])
                   + (y[2][0]-yprev[2][0])*(y[2][0]-yprev[2][0]));

                if(wall_hit_wall(yprev[0][0], yprev[1][0], yprev[2][0],
                                 l[i].r, l[i].phi, l[i].z, &wall_data)) {
                    l[i].running = 0;
                }

                if(fabs(fmod(yprev[1][0],2*math_pi) - fmod(y[1][0],2*math_pi))
                   > math_pi) {
                    orbits_tor++;
                    fprintf(f_pol, "%d %lf %lf %lf\n", i, y[0][0], y[2][0], l[i].length);
                }

                if((yprev[2][0]-axis_z)*(y[2][0]-axis_z)<0 && y[0][0]>axis_r) {
                    orbits_pol++;
                    real psi[1][NSIMD], rho[1][NSIMD];
                    B_field_eval_psi(0,psi, y[0][0], y[1][0], y[2][0], &B_data);
                    B_field_eval_rho(0, rho, psi[0][0], &B_data);
                    fprintf(f_tor, "%d %lf %lf %lf\n", i, fmod(y[1][0],2*math_pi), rho[0][0], l[i].length);
                }
            }
        }
        real ini_psi[1][NSIMD], ini_rho[1][NSIMD];
        B_field_eval_psi(0, ini_psi, l[i].ini_r, l[i].ini_phi, l[i].ini_z,
                         &B_data);
        B_field_eval_rho(0, ini_rho, ini_psi[0][0], &B_data);
        fprintf(f_sum, "%d %lf %lf %lf %lf %lf %ld %d %d\n", i, l[i].ini_r, l[i].ini_phi, l[i].ini_z, ini_rho[0][0], l[i].length, !l[i].running, orbits_pol, orbits_tor);
    }
    double tt_end = omp_get_wtime();

    #if VERBOSE >= 1
    printf("%d: %d lines done in %lf s.\n", id, n,
           tt_end-tt_start);
    #endif

    fclose(f_tor);
    fclose(f_pol);
    fclose(f_sum);
    fclose(f_debug);
}

void step_fieldline(
            real b21 = 0.2;
            real b31 = 0.075;
            real b32 = 0.225;
            real b41 = 0.3;
            real b42 = -0.9;
            real b43 = 1.2;
            real b51 = -0.203703703703704; 
            real b52 = 2.5;
            real b53 = -2.592592592592593;
            real b54 = 1.296296296296296;
            real b61 = 0.029495804398148;
            real b62 = 0.341796875;
            real b63 = 0.041594328703704;
            real b64 = 0.400345413773148;
            real b65 = 0.061767578125;
            real c1 = 0.097883597883598;
            real c3 = 0.402576489533011;
            real c4 = 0.210437710437710;
            real c6 = 0.289102202145680;
            real dc1 = c1-0.102177372685185;
            real dc3 = c3-0.383907903439153;
            real dc4 = c4-0.244592737268519;
            real dc5 = -0.019321986607143;
            real dc6 = c6-0.25;

            real k1[3];
            real k2[3];
            real k3[3];
            real k4[3];
            real k5[3];
            real k6[3];
            real yprev[3];
            real y[3];
            real yerr[3];

            yprev[0] = l[i].r;
            yprev[1] = l[i].phi;
            yprev[2] = l[i].z;

            B_field_eval_B(k1, yprev[0], yprev[1], yprev[2],
                           &B_data);
            k1[1] /= yprev[0];
            real normB = math_normc(k1[0], k1[1], k1[2]);
            k1[0] /= normB;
            k1[1] /= normB;
            k1[2] /= normB;
            int j;
            for(j = 0; j < 3; j++) {
                y[j] = yprev[j] + b21*h*k1[j];
            }

            B_field_eval_B(0, k2, y[0], y[1], y[2], &B_data);
            k2[1] /= y[0];
            normB = math_normc(k2[0], k2[1], k2[2]);
            k2[0] /= normB;
            k2[1] /= normB;
            k2[2] /= normB;
            for(j = 0; j < 3; j++) {
                y[j] = yprev[j] + h*(b31*k1[j] + b32*k2[j]);
            }

            B_field_eval_B(0, k3, y[0], y[1], y[2], &B_data);
            k3[1] /= y[0];
            normB = math_normc(k3[0], k3[1], k3[2]);
            k3[0] /= normB;
            k3[1] /= normB;
            k3[2] /= normB;
            for(j = 0; j < 3; j++) {
                y[j] = yprev[j] + h*(b41*k1[j] + b42*k2[j]
                                           + b43*k3[j]);
            }

            B_field_eval_B(0, k4, y[0], y[1], y[2], &B_data);
            k4[1] /= y[0];
            normB = math_normc(k4[0], k4[1], k4[2]);
            k4[0] /= normB;
            k4[1] /= normB;
            k4[2] /= normB;
            for(j = 0; j < 3; j++) {
                y[j] = yprev[j] + h*(b51*k1[j] + b52*k2[j]
                                           + b53*k3[j] + b54*k4[j]);
            }

            B_field_eval_B(0, k5, y[0], y[1], y[2], &B_data);
            k5[1] /= y[0];
            normB = math_normc(k5[0], k5[1], k5[2]);
            k5[0] /= normB;
            k5[1] /= normB;
            k5[2] /= normB;
            for(j = 0; j < 3; j++) {
                y[j] = yprev[j] + h*(b61*k1[j] + b62*k2[j]
                                           + b63*k3[j] + b64*k4[j]
                                           + b65*k5[j]);
            }

            B_field_eval_B(0, k6, y[0], y[1], y[2], &B_data);
            k6[1] /= y[0];
            normB = math_normc(k6[0], k6[1], k6[2]);
            k6[0] /= normB;
            k6[1] /= normB;
            k6[2] /= normB;
            for(j = 0; j < 3; j++) {
                y[j] = yprev[j] + h * (c1*k1[j] + c3*k3[j]
                                             + c4*k4[j] + c6*k6[j]);
                yerr[j] = h * (dc1*k1[j] + dc3*k3[j] + dc4*k4[j]
                                  + dc5*k5[j] + dc6*k6[j]);
            }

void simulate_fieldline_rk4(int id, int n, fieldline* l, sim_offload_data sim,
                        real* B_offload_array, real* wall_offload_array) {
    wall_data wall_data;
    wall_init(&wall_data, &sim.wall_offload_data, wall_offload_array);
    B_field_data B_data;
    B_field_init(&B_data, &sim.B_offload_data, B_offload_array);

    FILE* f_tor = fopen("filip5.poincareTor","w");
    FILE* f_pol = fopen("filip5.poincarePol","w");
    FILE* f_sum = fopen("filip5.summary","w");
    FILE* f_debug = fopen("filip5.debug","w");

    real h = 1e-0;
    real axis_r = B_field_get_axis_r(&B_data);
    real axis_z = B_field_get_axis_z(&B_data);
    double tt_start = omp_get_wtime();
    int i;
    #pragma omp parallel for
    for(i = 0; i < n; i++) {
        l[i].ini_r = l[i].r;
        l[i].ini_phi = l[i].phi;
        l[i].ini_z = l[i].z;

        h *= l[i].dir;

        int orbits_pol = 0;
        int orbits_tor = 0;

        while(l[i].running && l[i].length < 200000 && orbits_pol < 5000
              && orbits_tor < 5000) {
            real k1[3][NSIMD];
            real k2[3][NSIMD];
            real k3[3][NSIMD];
            real k4[3][NSIMD];
            real yprev[3][NSIMD];
            real y[3][NSIMD];

            yprev[0][0] = l[i].r;
            yprev[1][0] = l[i].phi;
            yprev[2][0] = l[i].z;

            B_field_eval_B(0, k1, yprev[0][0], yprev[1][0], yprev[2][0],
                           &B_data);
            real normB = math_normc(k1[0][0], k1[1][0], k1[2][0]);
            k1[0][0] /= normB;
            k1[1][0] /= normB;
            k1[2][0] /= normB;
            k1[1][0] /= yprev[0][0];
            int j;
            for(j = 0; j < 3; j++) {
                y[j][0] = yprev[j][0] + h/2.0*k1[j][0];
            }

            B_field_eval_B(0, k2, y[0][0], y[1][0], y[2][0], &B_data);
            normB = math_normc(k2[0][0], k2[1][0], k2[2][0]);
            k2[0][0] /= normB;
            k2[1][0] /= normB;
            k2[2][0] /= normB;
            k2[1][0] /= y[0][0];
            for(j = 0; j < 3; j++) {
                y[j][0] = yprev[j][0] + h/2.0*k2[j][0];
            }

            B_field_eval_B(0, k3, y[0][0], y[1][0], y[2][0], &B_data);
            normB = math_normc(k3[0][0], k3[1][0], k3[2][0]);
            k3[0][0] /= normB;
            k3[1][0] /= normB;
            k3[2][0] /= normB;
            k3[1][0] /= y[0][0];
            for(j = 0; j < 3; j++) {
                y[j][0] = yprev[j][0] + h*k3[j][0];
            }

            B_field_eval_B(0, k4, y[0][0], y[1][0], y[2][0], &B_data);
            normB = math_normc(k4[0][0], k4[1][0], k4[2][0]);
            k4[0][0] /= normB;
            k4[1][0] /= normB;
            k4[2][0] /= normB;
            k4[1][0] /= y[0][0];
            for(j = 0; j < 3; j++) {
                y[j][0] = yprev[j][0]
                    + h/6.0 * (k1[j][0] + 2*k2[j][0] + 2*k3[j][0] + k4[j][0]);
            }

            /*fprintf(f_debug, "%d %lf %lf %lf\n", i, y[0][0], y[1][0], y[2][0]);*/
            l[i].r = y[0][0];
            l[i].phi = y[1][0];
            l[i].z = y[2][0];
            l[i].length += sqrt((y[0][0]-yprev[0][0])*(y[0][0]-yprev[0][0])
                   + y[0][0]*(y[1][0]-yprev[1][0])*y[0][0]*(y[1][0]-yprev[1][0])
                   + (y[2][0]-yprev[2][0])*(y[2][0]-yprev[2][0]));

            if(wall_hit_wall(yprev[0][0], yprev[1][0], yprev[2][0],
                             l[i].r, l[i].phi, l[i].z, &wall_data)) {
                l[i].running = 0;
            }

            if(fabs(fmod(yprev[1][0],2*math_pi) - fmod(y[1][0],2*math_pi))
               > math_pi) {
                orbits_tor++;
                real s = (0.0 - fmod(yprev[1][0],2*math_pi)) / (y[1][0] - yprev[1][0]);
                fprintf(f_pol, "%d %lf %lf %lf\n", i, yprev[0][0]+s*(y[0][0]-yprev[0][0]), yprev[2][0]+s*(y[2][0]-yprev[2][0]), l[i].length);
            }

            if((yprev[2][0]-axis_z)*(y[2][0]-axis_z) < 0 && y[0][0] > axis_r) {
                orbits_pol++;
                real s = (axis_z - yprev[2][0]) / (y[2][0] - yprev[2][0]);
                real psi[1][NSIMD], rho[1][NSIMD];
                B_field_eval_psi(0, psi, yprev[0][0]+s*(y[0][0]-yprev[0][0]), yprev[1][0]+s*(y[1][0]-yprev[1][0]), axis_z, &B_data);
                B_field_eval_rho(0, rho, psi[0][0], &B_data);
                fprintf(f_tor, "%d %lf %lf %lf\n", i, fmod(yprev[1][0],2*math_pi)+s*(y[1][0] - yprev[1][0]), rho[0][0], l[i].length);
            }
        }
        real psi[1][NSIMD], rho[1][NSIMD];
        B_field_eval_psi(0, psi, l[i].ini_r, l[i].ini_phi, l[i].ini_z, &B_data);
        B_field_eval_rho(0, rho, psi[0][0], &B_data);
        fprintf(f_sum, "%d %lf %lf %lf %lf %lf %ld %d %d\n", i, l[i].ini_r, l[i].ini_phi, l[i].ini_z, rho[0][0], l[i].length, !l[i].running, orbits_pol, orbits_tor);
    }
    double tt_end = omp_get_wtime();

    #if VERBOSE >= 1
    printf("%d: %d lines done in %lf s.\n", id, n,
           tt_end-tt_start);
    #endif

    fclose(f_tor);
    fclose(f_pol);
    fclose(f_sum);
    fclose(f_debug);
}
