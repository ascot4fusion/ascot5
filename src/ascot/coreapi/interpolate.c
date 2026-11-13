/**
 * Routines for interpolating input data and derived quantities for external
 * use.
 *
 * Functions in this file allows to evaluate input data and quantities using
 * the same methods as is used in the actual simulation.
 */
#include "ascot.h"
#include "consts.h"
#include "data/atomic.h"
#include "data/bfield.h"
#include "data/boozer.h"
#include "data/efield.h"
#include "data/mhd.h"
#include "data/neutral.h"
#include "data/plasma.h"
#include "data/rfof.h"
#include "data/wall.h"
#include "datatypes.h"
#include "defines.h"
#include "parallel.h"
#include "utils/mathlib.h"
#include "utils/physlib.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/**
 * Store value on the output array.
 *
 * The output arrays have the format arr[icomp][ipnt], where icomp is the
 * component index and ipnt is the point index. No value is stored if the output
 * array is NULL.
 *
 * @param icomp Component index for output arrays that have multiple components
 *        (e.g. B which has components Br, Bphi, and Bz).
 * @param ipnt Point index corresponding to the point where the value is
 *        evaluated.
 * @param val Value to be stored.
 * @param arr Pointer to the output array or NULL to skip storing this
 *        quantiity.
 */
#define STORE(ipnt, icomp, val, arr)                                           \
    do                                                                         \
    {                                                                          \
        if (arr)                                                               \
            (arr)[(ipnt)][(icomp)] = (val);                                    \
    } while (0)

void ascot_interpolate(
    Bfield *bfield, Efield *efield, Plasma *plasma, Neutral *neutral,
    Boozer *boozer, Mhd *mhd, size_t npnt, int modenumber, real R[npnt],
    real phi[npnt], real z[npnt], real t[npnt], real B[3][npnt],
    real Bjac[9][npnt], real psi[4][npnt], real rho[2][npnt], real E[3][npnt],
    real n[][npnt], real T[2][npnt], real n0[][npnt], real T0[][npnt],
    real theta[4][npnt], real zeta[4][npnt], real alpha[5][npnt],
    real Phi[5][npnt], real mhd_b[3][npnt], real mhd_e[3][npnt],
    real mhd_phi[][npnt])
{
    int ONLY_PERTURBATIONS = 1;
#pragma omp parallel for
    for (size_t k = 0; k < npnt; k++)
    {
        real Bq[15], psival[4], rhoval[2], Eq[3], ns[MAX_SPECIES],
            Ts[MAX_SPECIES], psithetazeta[12], epert[3], bpert[3], Phipert[1],
            alpha_dalpha[5], Phi_dPhi[5];
        int psi_valid = 0, rho_valid = 0;
        int n_species, isinside;
        if (bfield && !Bfield_eval_b_db(Bq, R[k], phi[k], z[k], t[k], bfield))
        {
            STORE(0, k, Bq[0], B);
            STORE(1, k, Bq[1], B);
            STORE(2, k, Bq[2], B);
            STORE(0, k, Bq[3], Bjac);
            STORE(1, k, Bq[4], Bjac);
            STORE(2, k, Bq[5], Bjac);
            STORE(3, k, Bq[6], Bjac);
            STORE(4, k, Bq[7], Bjac);
            STORE(5, k, Bq[8], Bjac);
            STORE(6, k, Bq[9], Bjac);
            STORE(7, k, Bq[10], Bjac);
            STORE(8, k, Bq[11], Bjac);
        }
        if (bfield &&
            !Bfield_eval_psi_dpsi(psival, R[k], phi[k], z[k], t[k], bfield))
        {
            psi_valid = 1;
            STORE(0, k, psival[0], psi);
            STORE(1, k, psival[1], psi);
            STORE(2, k, psival[2], psi);
            STORE(3, k, psival[3], psi);
        }
        if (bfield && psi_valid && !Bfield_eval_rho(rhoval, psival[0], bfield))
        {
            rho_valid = 1;
            STORE(0, k, rhoval[0], rho);
            STORE(1, k, rhoval[1], rho);
        }
        if (efield && bfield &&
            !Efield_eval_e(Eq, R[k], phi[k], z[k], t[k], efield, bfield))
        {
            STORE(0, k, Eq[0], E);
            STORE(1, k, Eq[1], E);
            STORE(2, k, Eq[2], E);
        }
        if (plasma)
        {
            n_species = Plasma_get_n_species(plasma);
        }
        if (plasma && rho_valid &&
            !Plasma_eval_nT(
                ns, Ts, rhoval[0], R[k], phi[k], z[k], t[k], plasma))
        {
            STORE(0, k, Ts[0] / CONST_E, T);
            STORE(1, k, Ts[1] / CONST_E, T);
            for (int i = 0; i < n_species; i++)
            {
                STORE(i, k, ns[i], n);
            }
        }
        if (neutral)
        {
            n_species = Neutral_get_n_species(neutral);
        }
        if (neutral && rho_valid &&
            !Neutral_eval_nT(
                ns, Ts, rhoval[0], R[k], phi[k], z[k], t[k], neutral))
        {
            for (int i = 0; i < n_species; i++)
            {
                STORE(i, k, ns[i], n0);
                STORE(i, k, Ts[i], T0);
            }
        }
        if (boozer && bfield &&
            !Boozer_map_coordinates(
                psithetazeta, &isinside, R[k], phi[k], z[k], t[k], boozer,
                bfield))
        {
            if (isinside)
            {
                STORE(0, k, psithetazeta[4], theta);
                STORE(0, k, psithetazeta[8], zeta);
                STORE(1, k, psithetazeta[5], theta);
                STORE(2, k, psithetazeta[6], theta);
                STORE(3, k, psithetazeta[7], theta);
                STORE(1, k, psithetazeta[9], zeta);
                STORE(2, k, psithetazeta[10], zeta);
                STORE(3, k, psithetazeta[11], zeta);
            }
        }
        if (mhd && boozer && bfield &&
            !Mhd_eval_alpha_Phi(
                alpha_dalpha, Phi_dPhi, R[k], phi[k], z[k], t[k], modenumber,
                mhd, bfield, boozer))
        {
            STORE(0, k, alpha_dalpha[0], alpha);
            STORE(1, k, alpha_dalpha[2], alpha);
            STORE(2, k, alpha_dalpha[3], alpha);
            STORE(3, k, alpha_dalpha[4], alpha);
            STORE(4, k, alpha_dalpha[1], alpha);
            STORE(0, k, Phi_dPhi[0], Phi);
            STORE(1, k, Phi_dPhi[2], Phi);
            STORE(2, k, Phi_dPhi[3], Phi);
            STORE(3, k, Phi_dPhi[4], Phi);
            STORE(4, k, Phi_dPhi[2], Phi);
        }
        if (mhd && boozer && bfield &&
            !Mhd_eval_perturbation(
                bpert, epert, Phipert, R[k], phi[k], z[k], t[k],
                ONLY_PERTURBATIONS, modenumber, mhd, bfield, boozer))
        {
            STORE(0, k, bpert[0], mhd_b);
            STORE(1, k, bpert[1], mhd_b);
            STORE(2, k, bpert[2], mhd_b);
            STORE(0, k, epert[0], mhd_e);
            STORE(1, k, epert[1], mhd_e);
            STORE(2, k, epert[2], mhd_e);
            STORE(0, k, Phipert[0], mhd_phi);
        }
    }
}

void ascot_eval_collcoefs(
    size_t nspecies, size_t npnt, int qa, real ma, int qb[nspecies],
    real mb[nspecies], real nb[nspecies], real Tb[nspecies], real va[npnt],
    real clog[nspecies][npnt], real mu0[nspecies][npnt],
    real dmu0[nspecies][npnt], real mu1[nspecies][npnt], real F[nspecies][npnt],
    real Dpara[nspecies][npnt], real Dperp[nspecies][npnt],
    real K[nspecies][npnt], real nu[nspecies][npnt], real Q[nspecies][npnt],
    real dQ[nspecies][npnt], real dDpara[nspecies][npnt])
{

    for (size_t i = 0; i < nspecies; i++)
    {
        CPU_PARALLELIZE_LOOP
        for (size_t k = 0; k < npnt; k++)
        {

            real clogab[MAX_SPECIES];
            mccc_coefs_clog(clogab, ma, qa, va[k], nspecies, mb, qb, nb, Tb);

            real vb = sqrt(2 * Tb[i] / mb[i]);
            real x = va[k] / vb;
            real mufun[3] = {0., 0., 0.};
            mccc_coefs_mufun(mufun, x);

            real Fb = mccc_coefs_F(
                ma, qa, mb[i], qb[i], nb[i], vb, clogab[i], mufun[0]);
            real Qb = mccc_coefs_Q(
                ma, qa, mb[i], qb[i], nb[i], vb, clogab[i], mufun[0]);
            real dQb = mccc_coefs_dQ(
                ma, qa, mb[i], qb[i], nb[i], vb, clogab[i], mufun[2]);
            real Dparab = mccc_coefs_Dpara(
                ma, qa, va[k], qb[i], nb[i], vb, clogab[i], mufun[0]);
            real Dperpb = mccc_coefs_Dperp(
                ma, qa, va[k], qb[i], nb[i], vb, clogab[i], mufun[1]);
            real dDparab = mccc_coefs_dDpara(
                ma, qa, va[k], qb[i], nb[i], vb, clogab[i], mufun[0], mufun[2]);
            real Kb = mccc_coefs_K(va[k], Dparab, dDparab, Qb);
            real nub = mccc_coefs_nu(va[k], Dperpb);

            STORE(i, k, mufun[0], mu0);
            STORE(i, k, mufun[1], mu1);
            STORE(i, k, mufun[2], dmu0);
            STORE(i, k, clogab[i], clog);
            STORE(i, k, Fb, F);
            STORE(i, k, Qb, Q);
            STORE(i, k, dQb, dQ);
            STORE(i, k, Dparab, Dpara);
            STORE(i, k, Dperpb, Dperp);
            STORE(i, k, Kb, K);
            STORE(i, k, nub, nu);
            STORE(i, k, dDparab, dDpara);
        }
    }
}

/**
 * Evaluate ICRH electric field and the resonance condition.
 *
 * The evaluated electric field consists of left-hand (+) and right-hand (-)
 * circularly polarized components. The circularly polarised components are
 * notorious of their negatice effects on the mental health of their consumers.
 * In the context of this function, they have the following definitions:
 * E_plus_real = Re{E_LH ( cos( phase(E_LH) ) + i sin( phase(E_LH) ) } and
 * E_minus_real = Re{E_RH ( cos( phase(E_RH) ) - i sin( phase(E_RH) ) }.
 * E_LH, E_RH and their phases are all real. Furthermore, E_LH is called E_+ and
 * R_RH is called E_- in the context of RFOF and also more generally. Sometimes
 * the definition E+- = E_x +- iE_y. It is seen that |E+| = sqrt(E+  E+*) where
 * E+* is the complex conjugate. But because E+*=E-, one has
 * |E+| = sqrt(E+  E-) = |E-|. This is obviously not the case in this function,
 * as the magnitudes are different but this definition is also used sometimes.
 *
 * The resonance condition is given by
 *
 * omega_wave - n * omega_gyro - k_parallel * v_parallel - k_perp dot v_drift
 * = 0. Whether the Doppler shift (k_par v_par) or/and the drift term
 * (k_per v_drif) is used, is specified in the rfof_codeparam.xml. The drift
 * has not yet been implemented (Jan 2025).
 *
 * @param bfield magnetic field data.
 * @param rfof RFOF data.
 * @param Neval number of evaluation points.
 * @param R R coordinates of the evaluation points [m].
 * @param phi phi coordinates of the evaluation points [rad].
 * @param z z coordinates of the evaluation points [m].
 * @param t time coordinates of the evaluation points [s].
 * @param mass test particle mass (for computing resonance) [kg].
 * @param q test particle charge (for computing resonance) [C].
 * @param vpar test particle parallel velocity (for computing resonance) [m/s].
 * @param Eplus_real Real part of the right-handed electric field component of
 *        the wave [V/m]. (See comment above!)
 * @param Eminus_real Real part of the left-handed electric field component of
 *        the wave [V/m]. (See comment above!)
 * @param Eplus_imag Imaginary part of the right-handed electric field component
 *        of the wave [V/m]. (See comment above!)
 * @param Eminus_imag Imaginary part of the left-handed electric field
 *        component of the wave [V/m]. (See comment above!)
 * @param res_cond value of the resonance condition where zero is the resonance
 *        [1].
 */
void ascot_eval_rfof(
    Bfield *bfield, Rfof *rfof, int Neval, real *R, real *phi, real *z, real *t,
    real mass, real q, real vpar, real *Eplus_real, real *Eminus_real,
    real *Eplus_imag, real *Eminus_imag, real *res_cond)
{

#pragma omp parallel
    {
        /* The function that evaluates resonance condition takes an RFOF marker
         * as an input. However, only the R and vpar values are actually used.
         * Therefore, we initialize a dummy marker and adjust only the values of
         * R and vpar. */
        rfof_marker rfof_mrk;
        int dummy_int = 1;
        real dummy_real = -999.0; /*-999.0 to be on the safe side */
        rfof_set_up(&rfof_mrk, rfof);

#pragma omp for
        for (int k = 0; k < Neval; k++)
        {
            real B[3];
            if (Bfield_eval_b(B, R[k], phi[k], z[k], t[k], bfield))
            {
                continue;
            }
            real B_magn = sqrt(B[0] * B[0] + B[1] * B[1] + B[2] * B[2]);
            real gyrofreq = q * B_magn / mass;
            rfof_set_marker_manually(
                &rfof_mrk, &dummy_int, &dummy_real, &(R[k]), &dummy_real,
                &dummy_real, &dummy_real, &dummy_real, &dummy_real, &dummy_real,
                &dummy_real, &dummy_real, &dummy_real, &vpar, &dummy_real,
                &gyrofreq, &dummy_real, &dummy_real, &dummy_int, &dummy_int);

            int nharm; /* For storing return value which is not used */
            rfof_eval_resonance_function(
                &(res_cond[k]), &nharm, &rfof_mrk, rfof);

            // TODO: this should return a non-zero value for failed evaluations
            rfof_eval_rf_wave(
                &(Eplus_real[k]), &(Eminus_real[k]), &(Eplus_imag[k]),
                &(Eminus_imag[k]), R[k], z[k], rfof);
        }
        rfof_tear_down(&rfof_mrk);
    }
}

#undef STORE
