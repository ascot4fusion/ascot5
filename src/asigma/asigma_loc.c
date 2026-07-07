/**
 * @file asigma_loc.c
 * @brief Atomic reaction data from local files
 *
 * Atomic reaction data (sigmas) originating from local files and
 * interpolated using splines. If the data for a reaction is missing,
 * a corresponding analytical model implemented in ASCOT5 might be used.
 */
#include <stdio.h>
#include <stdlib.h>
#include "../ascot5.h"
#include "../print.h"
#include "../error.h"
#include "../spline/interp.h"
#include "../consts.h"
#include "../math.h"
#include "../physlib.h"
#include "../suzuki.h"
#include "../asigma.h"
#include "asigma_loc.h"

/**
 * @brief Initialize local file atomic data and check inputs
 *
 * @param data pointer to the data struct
 *
 * @return zero if initialization success
 */
int asigma_loc_init(asigma_loc_data* data, int nreac,
                    int* z1, int* a1, int* z2, int* a2, int* reactype,
                    int* ne, real* emin, real* emax,
                    int* nn, real* nmin, real* nmax,
                    int* nT, real* Tmin, real* Tmax, real* sigma) {

    int err = 0;
    data->N_reac = nreac;
    data->z_1 = (int*) malloc( nreac * sizeof(int) );
    data->a_1 = (int*) malloc( nreac * sizeof(int) );
    data->z_2 = (int*) malloc( nreac * sizeof(int) );
    data->a_2 = (int*) malloc( nreac * sizeof(int) );
    data->reac_type = (int*) malloc( nreac * sizeof(int) );
    data->sigma = (interp1D_data*) malloc( nreac * sizeof(interp1D_data) );
    data->sigmav = (interp2D_data*) malloc( nreac * sizeof(interp2D_data) );
    data->BMSsigmav = (interp3D_data*) malloc( nreac * sizeof(interp3D_data) );
    real* pos = sigma;
    for(int i_reac = 0; i_reac < nreac; i_reac++) {
        data->z_1[i_reac] = z1[i_reac];
        data->a_1[i_reac] = a1[i_reac];
        data->z_2[i_reac] = z2[i_reac];
        data->a_2[i_reac] = a2[i_reac];
        data->reac_type[i_reac] = reactype[i_reac];

        /* Initialize spline struct according to dimensionality of
           reaction data (and mark reaction availability) */
        int dim = (ne[i_reac] > 1) + (nn[i_reac] > 1) + (nT[i_reac] > 1);
        switch(dim) {
            case 1:
                err = interp1Dcomp_setup(&data->sigma[i_reac], pos,
                                         ne[i_reac], NATURALBC,
                                         emin[i_reac], emax[i_reac]);
                pos += ne[i_reac];
                break;
            case 2:
                err = interp2Dcomp_setup(&data->sigmav[i_reac], pos,
                                         ne[i_reac], nT[i_reac],
                                         NATURALBC, NATURALBC,
                                         emin[i_reac], emax[i_reac],
                                         Tmin[i_reac], Tmax[i_reac]);
                pos += ne[i_reac] * nT[i_reac];
                break;
            case 3:
                err = interp3Dcomp_setup(&data->BMSsigmav[i_reac], pos,
                                         ne[i_reac], nn[i_reac], nT[i_reac],
                                         NATURALBC, NATURALBC, NATURALBC,
                                         emin[i_reac], emax[i_reac],
                                         nmin[i_reac], nmax[i_reac],
                                         Tmin[i_reac], Tmax[i_reac]);
                pos += ne[i_reac] * nn[i_reac] * nT[i_reac];
                break;
            default:
                /* Unrecognized dimensionality. Produce error. */
                print_err("Error: Unrecognized abscissa dimensionality\n");
                err = 1;
                break;
        }
    }

    print_out(VERBOSE_IO, "\nFound data for %d atomic reactions:\n",
              data->N_reac);
    for(int i_reac = 0; i_reac < data->N_reac; i_reac++) {
        print_out(VERBOSE_IO,
              "Reaction number / Total number of reactions = %d / %d\n"
              "  Reactant species Z_1 / A_1, Z_2 / A_2     = %d / %d, %d / %d\n"
              "  Min/Max energy                            = %1.2le / %1.2le\n"
              "  Min/Max density                           = %1.2le / %1.2le\n"
              "  Min/Max temperature                       = %1.2le / %1.2le\n"
              "  Number of energy grid points              = %d\n"
              "  Number of density grid points             = %d\n"
              "  Number of temperature grid points         = %d\n",
              i_reac+1, data->N_reac, data->z_1[i_reac], data->a_1[i_reac],
              data->z_2[i_reac], data->a_2[i_reac], emin[i_reac], emax[i_reac],
              nmin[i_reac], nmax[i_reac], Tmin[i_reac], Tmax[i_reac],
              ne[i_reac], nn[i_reac], nT[i_reac]);
    }

    return err;
}

/**
 * @brief Free allocated resources
 *
 * @param offload_data pointer to the data struct
*/
void asigma_loc_free(asigma_loc_data* data) {
    for(int i_reac = 0; i_reac < data->N_reac; i_reac++) {
        if(data->reac_type[i_reac] == sigmav_CX) {
            free(data->sigmav[i_reac].c);
        }
        else if(data->reac_type[i_reac] == sigmav_EII) {
            free(data->sigmav[i_reac].c);
        }
        else if(data->reac_type[i_reac] == sigmav_BMS) {
            free(data->BMSsigmav[i_reac].c);
        }
    }
    free(data->z_1);
    free(data->a_1);
    free(data->z_2);
    free(data->a_2);
    free(data->reac_type);
    free(data->sigma);
    free(data->sigmav);
    free(data->BMSsigmav);
}

/**
 * @brief Offload data to the accelerator.
 *
 * @param data pointer to the data struct
 */
void asigma_loc_offload(asigma_loc_data* data) {
    // TODO: Implement
}

/**
 * @brief Evaluate charge-exchange rate coefficient
 *
 * This function evaluates the rate coefficient (<sigma*v>) for the atomic
 * reaction corresponding to the reaction identifiers given as parameters
 * at the given fast particle energy and bulk plasma conditions.
 *
 * This is a SIMD function.
 *
 * @param ratecoeff pointer to evaluated rate coefficient
 * @param z_1 atomic number of fast particle
 * @param a_1 atomic mass number of fast particle
 * @param E energy of fast particle
 * @param mass mass of fast particle
 * @param znum atomic numbers of bulk neutrals
 * @param anum atomic mass numbers of bulk neutrals
 * @param T_0 temperature of bulk neutrals
 * @param n_0 neutral densities
 * @param extrapolate don't raise error but set values outside abscissae to zero
 * @param asigma_data pointer to atomic data struct
 *
 * @return zero if evaluation succeeded
 */
a5err asigma_loc_eval_cx(
    real* ratecoeff, int z_1, int a_1, real E, real mass, int nspec,
    const int* znum, const int* anum, real T_0, real* n_0, int extrapolate,
    asigma_loc_data* asigma_data) {
    a5err err = 0;

    /* Convert Joule to eV */
    E   /= CONST_E;
    T_0 /= CONST_E;
    *ratecoeff = 0;
    for(int i_spec = 0; i_spec < nspec; i_spec++) {

        /* Find the matching reaction */
        int reac_found = -1, i_reac;
        for(i_reac = 0; i_reac < asigma_data->N_reac; i_reac++) {
            if(asigma_data->z_1[i_reac] == z_1 &&
               asigma_data->a_1[i_reac] == a_1 &&
               asigma_data->z_2[i_reac] == znum[i_spec] &&
               asigma_data->a_2[i_reac] == anum[i_spec] &&
               asigma_data->reac_type[i_reac] == sigmav_CX) {
                reac_found = i_reac;
            }
        }
        i_reac = reac_found;

        if(reac_found < 0) {
            /* Reaction not found. Raise error. */
            err = error_raise(ERR_INPUT_EVALUATION, __LINE__, EF_ASIGMA_LOC);
        } else {
            real sigmav;
            int interperr = interp2Dcomp_eval_f(
                                &sigmav, &asigma_data->sigmav[i_reac], E, T_0);

            /* Interpolation error means the data has to be extrapolated */
            if(interperr) {
                if(extrapolate) {
                    sigmav = 0.0;
                } else {
                    err = error_raise( ERR_INPUT_EVALUATION, __LINE__,
                                        EF_ASIGMA_LOC );
                }
            }
            *ratecoeff += sigmav*n_0[i_spec];
        }
    }

    return err;
}

/**
 * @brief Evaluate beam stopping rate coefficient
 *
 * This function first tries to evaluate BMS with ADAS data. If not present,
 * the Suzuki model is used instead.
 *
 * This is a SIMD function.
 *
 * @param ratecoeff pointer to evaluated rate coefficient
 * @param z_1 atomic number of fast particle
 * @param a_1 atomic mass number of fast particle
 * @param E energy of fast particle
 * @param mass mass of fast particle
 * @param nion number of bulk ion species
 * @param znum atomic numbers of bulk particles
 * @param anum atomic mass numbers of bulk particles
 * @param T_e electron temperature of bulk plasma
 * @param n_i densities of bulk ions
 * @param extrapolate don't raise error but set values outside abscissae to zero
 * @param asigma_data pointer to atomic data struct
 *
 * @return zero if evaluation succeeded
 */
a5err asigma_loc_eval_bms(
    real* ratecoeff, int z_1, int a_1, real E, real mass, int nion,
    const int* znum, const int* anum, real T_e, real* n_i, int extrapolate,
    asigma_loc_data* asigma_data) {
    a5err err = 0;

    /* Convert Joule to eV */
    real E_eV = E / CONST_E;
    T_e /= CONST_E;

    /* Find the matching reaction. Note that BMS data is same for all
     * isotopes, so we don't compare anums */
    int reac_found = -1; real n_e = 0; *ratecoeff = 0;
    for(int i_spec = 0; i_spec < nion; i_spec++) {
        n_e += znum[i_spec] * n_i[i_spec];
        for(int i_reac = 0; i_reac < asigma_data->N_reac; i_reac++) {
            if(asigma_data->z_1[i_reac]       == z_1 &&
               asigma_data->z_2[i_reac]       == znum[i_spec] &&
               asigma_data->reac_type[i_reac] == sigmav_BMS) {
                reac_found = i_reac;
                real sigmav;
                int interperr = \
                        interp3Dcomp_eval_f(
                            &sigmav, &asigma_data->BMSsigmav[i_reac],
                            E_eV/a_1, znum[i_spec] * n_i[i_spec], T_e);

                /* Interpolation error means the data has to be extrapolated */
                if(interperr) {
                    if(extrapolate) {
                        sigmav = 0.0;
                    } else {
                        err = error_raise( ERR_INPUT_EVALUATION, __LINE__,
                                           EF_ASIGMA_LOC );
                    }
                }
                *ratecoeff += sigmav * ( znum[i_spec] * n_i[i_spec]);
            }
        }
    }
    *ratecoeff /= n_e;

    if(reac_found < 0) {
        /* Reaction not found. Try Suzuki before throwing error. */
        T_e *= CONST_E;
        real gamma = physlib_gamma_Ekin(mass, E);
        real vnorm = physlib_vnorm_gamma(gamma);
        if(suzuki_sigmav(ratecoeff, E/a_1, vnorm, n_e, T_e, nion, n_i, anum,
                         znum)) {
            err = error_raise(ERR_INPUT_EVALUATION, __LINE__, EF_ASIGMA_LOC);
        }
    }

    return err;
}

/**
 * @brief Evaluate electron impact ionization rate coefficient
 *
 * This function evaluates the rate coefficient (<sigma*v>) for the atomic
 * reaction corresponding to the reaction identifiers given as parameters
 * at the given fast particle energy and bulk plasma conditions.
 *
 * This is a SIMD function.
 *
 * @param ratecoeff pointer to evaluated rate coefficient
 * @param z_1 atomic number of fast particle
 * @param a_1 atomic mass number of fast particle
 * @param E energy of fast particle
 * @param mass mass of fast particle
 * @param Te electron temperature
 * @param ne electron density
 * @param extrapolate don't raise error but set values outside abscissae to zero
 * @param asigma_data pointer to atomic data struct
 *
 * @return zero if evaluation succeeded
 */
a5err asigma_loc_eval_eii(
    real* ratecoeff, int z_1, int a_1, real E, real mass, real Te, real ne,
    int extrapolate, asigma_loc_data* asigma_data) {

    a5err err = 0;

    /* Convert Joule to eV */
    E   /= CONST_E;
    Te /= CONST_E;
    *ratecoeff = 0;

    /* Find the matching reaction */
    int reac_found = -1, i_reac;
    for(i_reac = 0; i_reac < asigma_data->N_reac; i_reac++) {
        if(asigma_data->z_1[i_reac] == z_1 &&
            asigma_data->a_1[i_reac] == a_1 &&
            asigma_data->z_2[i_reac] == 0 &&
            asigma_data->a_2[i_reac] == 0 &&
            asigma_data->reac_type[i_reac] == sigmav_EII) {
            reac_found = i_reac;
        }
    }
    i_reac = reac_found;

    if(reac_found < 0) {
        /* Reaction not found. Raise error. */
        err = error_raise(ERR_INPUT_EVALUATION, __LINE__, EF_ASIGMA_LOC);
    } else {
        real sigmav;
        int interperr = interp2Dcomp_eval_f(
                            &sigmav, &asigma_data->sigmav[i_reac], E, Te);

        /* Interpolation error means the data has to be extrapolated */
        if(interperr) {
            if(extrapolate) {
                sigmav = 0.0;
            } else {
                err = error_raise( ERR_INPUT_EVALUATION, __LINE__,
                                    EF_ASIGMA_LOC );
            }
        }
        *ratecoeff += sigmav*ne;
    }

    return err;

}
