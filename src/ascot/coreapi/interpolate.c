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
 * Store val at index idx in output array arr if the output array is not NULL.
 */
#define STORE(idx, val, arr)                                                   \
    do                                                                         \
    {                                                                          \
        if (arr)                                                               \
            (arr)[(idx)] = (val);                                              \
    } while (0)

void ascot_interpolate(
    Bfield *bfield, Efield *efield, Plasma *plasma, Neutral *neutral,
    Boozer *boozer, Mhd *mhd,
    // Atomic* atomic,
    size_t npnt, int modenumber, real R[npnt], real phi[npnt], real z[npnt],
    real t[npnt], real B[3][npnt], real Bjac[9][npnt], real psi[4][npnt],
    real rho[2][npnt], real E[3][npnt], real n[][npnt], real T[2][npnt],
    real n0[][npnt], real T0[][npnt], real theta[4][npnt], real zeta[4][npnt],
    real alpha[5][npnt], real Phi[5][npnt], real mhd_b[3][npnt],
    real mhd_e[3][npnt], real mhd_phi[npnt])
{
    int ONLY_PERTURBATIONS = 1;
    OMP_PARALLEL_CPU_ONLY
    for (size_t k = 0; k < npnt; k++)
    {
        real Bq[15], psival[4], rhoval[2], Eq[3], ns[MAX_SPECIES],
            Ts[MAX_SPECIES], psithetazeta[12], epert[3], bpert[3], Phipert[1],
            alpha_dalpha[5], Phi_dPhi[5];
        int psi_valid = 0, rho_valid = 0;
        int n_species, isinside;

        if (bfield && !Bfield_eval_b_db(Bq, R[k], phi[k], z[k], t[k], bfield))
        {
            STORE(0 * npnt + k, Bq[0], *B);
            STORE(1 * npnt + k, Bq[4], *B);
            STORE(2 * npnt + k, Bq[8], *B);
            STORE(0 * npnt + k, Bq[1], *Bjac);
            STORE(1 * npnt + k, Bq[2], *Bjac);
            STORE(2 * npnt + k, Bq[3], *Bjac);
            STORE(3 * npnt + k, Bq[5], *Bjac);
            STORE(4 * npnt + k, Bq[6], *Bjac);
            STORE(5 * npnt + k, Bq[7], *Bjac);
            STORE(6 * npnt + k, Bq[9], *Bjac);
            STORE(7 * npnt + k, Bq[10], *Bjac);
            STORE(8 * npnt + k, Bq[11], *Bjac);
        }
        if (bfield &&
            !Bfield_eval_psi_dpsi(psival, R[k], phi[k], z[k], t[k], bfield))
        {
            psi_valid = 1;
            STORE(0 * npnt + k, psival[0], *psi);
            STORE(1 * npnt + k, psival[1], *psi);
            STORE(2 * npnt + k, psival[2], *psi);
            STORE(3 * npnt + k, psival[3], *psi);
        }
        if (bfield && psi_valid && !Bfield_eval_rho(rhoval, psival[0], bfield))
        {
            rho_valid = 1;
            STORE(0 * npnt + k, rhoval[0], *rho);
            STORE(1 * npnt + k, rhoval[1], *rho);
        }
        if (efield && bfield &&
            !Efield_eval_e(Eq, R[k], phi[k], z[k], t[k], efield, bfield))
        {
            STORE(0 * npnt + k, Eq[0], *E);
            STORE(1 * npnt + k, Eq[1], *E);
            STORE(2 * npnt + k, Eq[2], *E);
        }
        if (plasma)
        {
            n_species = Plasma_get_n_species(plasma);
        }
        if (plasma && rho_valid &&
            !Plasma_eval_nT(
                ns, Ts, rhoval[0], R[k], phi[k], z[k], t[k], plasma))
        {
            STORE(0 * npnt + k, Ts[0] / CONST_E, *T);
            STORE(1 * npnt + k, Ts[1] / CONST_E, *T);
            for (int i = 0; i < n_species; i++)
            {
                STORE(i * npnt + k, ns[i], *n);
            }
        }
        if (neutral)
        {
            n_species = Neutral_get_n_species(neutral);
        }
        if (neutral && rho_valid &&
            !Neutral_eval_density(ns, rhoval[0], R[k], phi[k], z[k], t[k], neutral))
        {
            for (int i = 0; i < n_species; i++)
            {
                STORE(i * npnt + k, ns[i], *n0);
            }
        }
        if (neutral && rho_valid &&
            !Neutral_eval_temperature(Ts, rhoval[0], R[k], phi[k], z[k], t[k], neutral))
        {
            for (int i = 0; i < n_species; i++)
            {
                STORE(i * npnt + k, Ts[i], *T0);
            }
        }
        if (boozer && bfield &&
            !Boozer_map_coordinates(
                psithetazeta, &isinside, R[k], phi[k], z[k], t[k], boozer, bfield))
        {
            if (isinside)
            {
                STORE(0 * npnt + k, psithetazeta[4], *theta);
                STORE(0 * npnt + k, psithetazeta[8], *zeta);
                STORE(1 * npnt + k, psithetazeta[5], *theta);
                STORE(2 * npnt + k, psithetazeta[6], *theta);
                STORE(3 * npnt + k, psithetazeta[7], *theta);
                STORE(1 * npnt + k, psithetazeta[9], *zeta);
                STORE(2 * npnt + k, psithetazeta[10], *zeta);
                STORE(3 * npnt + k, psithetazeta[11], *zeta);
            }
        }
        if (mhd && boozer && bfield &&
            !Mhd_eval_alpha_Phi(
                alpha_dalpha, Phi_dPhi, R[k], phi[k], z[k], t[k], modenumber, mhd,
                bfield, boozer))
        {
            STORE(0 * npnt + k, alpha_dalpha[0], *alpha);
            STORE(1 * npnt + k, alpha_dalpha[2], *alpha);
            STORE(2 * npnt + k, alpha_dalpha[3], *alpha);
            STORE(3 * npnt + k, alpha_dalpha[4], *alpha);
            STORE(4 * npnt + k, alpha_dalpha[1], *alpha);
            STORE(0 * npnt + k, Phi_dPhi[0], *Phi);
            STORE(1 * npnt + k, Phi_dPhi[2], *Phi);
            STORE(2 * npnt + k, Phi_dPhi[3], *Phi);
            STORE(3 * npnt + k, Phi_dPhi[4], *Phi);
            STORE(4 * npnt + k, Phi_dPhi[2], *Phi);
        }
        if (mhd && boozer && bfield &&
            !Mhd_eval_perturbation(
                bpert, epert, Phipert, R[k], phi[k], z[k], t[k], ONLY_PERTURBATIONS,
                modenumber, mhd, bfield, boozer))
        {
            STORE(0 * npnt + k, bpert[0], *mhd_b);
            STORE(1 * npnt + k, bpert[1], *mhd_b);
            STORE(2 * npnt + k, bpert[2], *mhd_b);
            STORE(0 * npnt + k, epert[0], *mhd_e);
            STORE(1 * npnt + k, epert[1], *mhd_e);
            STORE(2 * npnt + k, epert[2], *mhd_e);
            STORE(0 * npnt + k, Phipert[0], mhd_phi);
        }
    }
}

void ascot_eval_collcoefs(
    Bfield *bfield, Plasma *plasma, size_t npnt, real ma, real qa, real R[npnt],
    real phi[npnt], real z[npnt], real t[npnt], real va[npnt],
    real coefficients[6][npnt], real clog[][npnt], real mu[3][npnt])
{
    int n_species = Plasma_get_n_species(plasma);
    const real *qb = Plasma_get_species_charge(plasma);
    const real *mb = Plasma_get_species_mass(plasma);
    mccc_data mccc;

    OMP_PARALLEL_CPU_ONLY
    for (size_t k = 0; k < npnt; k++)
    {
        real mufun[3] = {0., 0., 0.};

        /* Evaluate rho as it is needed to evaluate plasma parameters */
        real psi, rho[2];
        if (Bfield_eval_psi(&psi, R[k], phi[k], z[k], t[k], bfield))
        {
            continue;
        }
        if (Bfield_eval_rho(rho, psi, bfield))
        {
            continue;
        }

        real nb[MAX_SPECIES], Tb[MAX_SPECIES];
        if (Plasma_eval_nT(
                nb, Tb, rho[0], R[k], phi[k], z[k], t[k], plasma))
        {
            continue;
        }

        /* Loop through all plasma species */
        for (int ib = 0; ib < n_species; ib++)
        {

            /* Coulomb logarithm */
            real clogab[MAX_SPECIES];
            mccc_coefs_clog(clogab, ma, qa, va[k], n_species, mb, qb, nb, Tb);

            /* Special functions */
            real vb = sqrt(2 * Tb[ib] / mb[ib]);
            real x = va[k] / vb;
            mccc_coefs_mufun(mufun, x, &mccc);

            (void)clog;
            (void)mu;
            (void)coefficients;

            /* Coefficients */
            /*
            real Fb = mccc_coefs_F(
                ma, qa, mb[ib], qb[ib], nb[ib], vb, clogab[ib], mufun[0]);
            real Qb = mccc_coefs_Q(
                ma, qa, mb[ib], qb[ib], nb[ib], vb, clogab[ib], mufun[0]);
            real dQb = mccc_coefs_dQ(
                ma, qa, mb[ib], qb[ib], nb[ib], vb, clogab[ib], mufun[2]);
            real Dparab = mccc_coefs_Dpara(
                ma, qa, va[k], qb[ib], nb[ib], vb, clogab[ib], mufun[0]);
            real Dperpb = mccc_coefs_Dperp(
                ma, qa, va[k], qb[ib], nb[ib], vb, clogab[ib], mufun[1]);
            real dDparab = mccc_coefs_dDpara(
                ma, qa, va[k], qb[ib], nb[ib], vb, clogab[ib], mufun[0],
                mufun[2]);
            real Kb = mccc_coefs_K(va[k], Dparab, dDparab, Qb);
            real nub = mccc_coefs_nu(va[k], Dperpb);
            */

            // int idx = ib * npnt + k;
            //  STORE(n_species*npnt*0 + ib*npnt + k, mufun[0], **mu);
            /**
            if(mu0 != NULL)    { mu0[idx]    = mufun[0];   }
            if(mu1 != NULL)    { mu1[idx]    = mufun[1];   }
            if(dmu0 != NULL)   { dmu0[idx]   = mufun[2];   }
            if(clog != NULL)   { clog[idx]   = clogab[ib]; }
            if(F != NULL)      { F[idx]      = Fb;         }
            if(Dpara != NULL)  { Dpara[idx]  = Dparab;     }
            if(Dperp != NULL)  { Dperp[idx]  = Dperpb;     }
            if(K != NULL)      { K[idx]      = Kb;         }
            if(nu != NULL)     { nu[idx]     = nub;        }
            if(Q != NULL)      { Q[idx]      = Qb;         }
            if(dQ != NULL)     { dQ[idx]     = dQb;        }
            if(dDpara != NULL) { dDpara[idx] = dDparab;    } */
        }
    }
}

/**
 * Evaluate atomic reaction rate coefficient.
 *
 * @param Simulation initialized simulation data struct
 * @param Neval number of evaluation points in (R, phi, z, t).
 * @param R R coordinates of the evaluation points [m].
 * @param phi phi coordinates of the evaluation points [rad].
 * @param z z coordinates of the evaluation points [m].
 * @param t time coordinates of the evaluation points [s].
 * @param Nv number of evaluation points in velocity space.
 * @param va test particle velocities [m/s].
 * @param Aa test particle atomic mass number.
 * @param Za test particle charge number.
 * @param ma test particle mass.
 * @param reac_type reaction type
 * @param ratecoeff output array where evaluated values are stored [1/m^2].
 */
void ascot_eval_ratecoeff(
    Bfield *bfield, Plasma *plasma, Neutral *neutral, Atomic *atomic, int Neval,
    real *R, real *phi, real *z, real *t, int Nv, real *va, int Aa, int Za,
    real ma, int reac_type, real *ratecoeff)
{

    const int *Zb = Plasma_get_species_znum(plasma);
    const int *Ab = Plasma_get_species_anum(plasma);
    int nion = Plasma_get_n_species(plasma) - 1;
    int nspec = Neutral_get_n_species(neutral);

    OMP_PARALLEL_CPU_ONLY
    for (int k = 0; k < Neval; k++)
    {
        real psi[1], rho[2], T0[1], n[MAX_SPECIES], T[MAX_SPECIES],
            n0[MAX_SPECIES];
        if (Bfield_eval_psi(psi, R[k], phi[k], z[k], t[k], bfield))
        {
            continue;
        }
        if (Bfield_eval_rho(rho, psi[0], bfield))
        {
            continue;
        }
        if (Plasma_eval_nT(
                n, T, rho[0], R[k], phi[k], z[k], t[k], plasma))
        {
            continue;
        }
        if (Neutral_eval_temperature(T0, rho[0], R[k], phi[k], z[k], t[k], neutral))
        {
            continue;
        }
        if (Neutral_eval_density(n0, rho[0], R[k], phi[k], z[k], t[k], neutral))
        {
            continue;
        }
        for (int j = 0; j < Nv; j++)
        {
            real E =
                (physlib_gamma_vnorm(va[j]) - 1.0) * ma * CONST_C * CONST_C;
            real val;
            switch (reac_type)
            {
            case sigmav_CX:
                if (Atomic_eval_cx(
                        &val, Za, Aa, E, ma, nspec, Zb, Ab, T0[0], n0, atomic))
                {
                    continue;
                }
                ratecoeff[Nv * k + j] = val;
                break;
            case sigmav_BMS:
                if (Atomic_eval_bms(
                        &val, Za, Aa, E, ma, nion, Zb, Ab, T[0], n, atomic))
                {
                    continue;
                }
                ratecoeff[Nv * k + j] = val * n[0];
                break;
            default:
                break;
            }
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
