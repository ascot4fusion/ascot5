/**
 * @file test_interact.c
 * @brief Test program for Coulomb collision interactions
 *
 * This test program calculates some coulomb collision parameters using the
 * functions in interact.c. This can be compared to the values given by
 * ASCOT4's test_interact.
 */
#include <stdio.h>
#include <math.h>
#include "ascot5.h"
#include "interact.h"

int main(int argc, char** argv) {
    real dens[2] = { 1e20, 1e20 };
    real temp[2] = { 6.96e7, 4.64e7 };
    real charge[2] = { -1.6e-19, 1.6e-19 };
    real mass[2] = { 9.1e-31, 3.3e-27 };

    real v = 1.3e7;
    real m = 6.6e-27;
    real q = 3.2e-19;

    real debye_length2 = interact_debye_length2(temp, dens, charge, 2);
    real coulomb_log_e = interact_coulomb_log(v, m, q, temp[0], dens[0], charge[0], mass[0], debye_length2);
    real coulomb_log_i = interact_coulomb_log(v, m, q, temp[1], dens[1], charge[1], mass[1], debye_length2);
    real nu = interact_nu(v, m, q, temp, dens, charge, mass, 2);
    real D_par = interact_D_par(v, m, q, temp, dens, charge, mass, 2);
    real D_perp = interact_D_perp(v, m, q, temp, dens, charge, mass, 2);
 
    printf("Debye length: %le\n", debye_length2);
    printf("Coulomb log e: %le\n", coulomb_log_e);
    printf("Coulomb log i: %le\n", coulomb_log_i);
    printf("nu_s: %le\n", nu);
    printf("D_par: %le\n", D_par);
    printf("D_perp: %le\n", D_perp);

    return 0;
}
