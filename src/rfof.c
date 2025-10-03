/**
 * Implements rfof.h.
 */
#include "rfof.h"
#include "consts.h"
#include "particle.h"
#include "physlib.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef RFOF
void __ascot5_icrh_routines_MOD_call_initev_excl_marker_stuff(
    const char *xml_filename, int **xml_filename_len, void **cptr_rfglobal,
    void **cptr_rfof_input_params);
void __ascot5_icrh_routines_MOD_call_initialise_res_mem(
    void **cptr_mem, int *cptr_mem_shape_i, int *cptr_mem_shape_j,
    void **cptr_rfglobal, void **cptr_rfof_input_param);
void __ascot5_icrh_routines_MOD_call_initialise_diagnostics(
    void **cptr_RFglobal, void **cptr_diagno);

void __ascot5_icrh_routines_MOD_call_set_marker_pointers(
    void **cptr_marker, int **id, real **weight, real **R, real **phi, real **z,
    real **psi, real **charge, real **mass, real **Ekin, real **velocity,
    real **mu, real **pphicanonical, real **vpar, real **vperp, real **gyrof,
    real **tauB, real **vdriftRho, real **acc, int *isOrbitTimeAccelerated,
    int *is_already_allocated);

void __ascot5_icrh_routines_MOD_call_rf_kick(
    double *time, double *dtin, int *mpi_rank, void **cptr_marker,
    void **cptr_mem, void **cptr_rfglobal, void **cptr_rfdiagno,
    void **cptr_rfof_input, int *mem_shape_i, int *mem_shape_j, int *err,
    rfof_output *out);

void __ascot5_icrh_routines_MOD_call_reset_res_mem(
    void **rfof_mem_pointer, int *mem_shape_i, int *mem_shape_j);

void __ascot5_icrh_routines_MOD_call_deallocate_rfof_input_param(
    void **cptr_rfof_input_param);
void __ascot5_icrh_routines_MOD_call_deallocate_rfglobal(void **cptr_rfglobal);
void __ascot5_icrh_routines_MOD_call_deallocate_res_mem(
    void **cptr_res_mem, int *cptr_mem_shape_i, int *cptr_mem_shape_j);
void __ascot5_icrh_routines_MOD_call_deallocate_diagnostics(void **cptr_diagno);
void __ascot5_icrh_routines_MOD_deallocate_marker(void **cptr_rfof_marker);

void __ascot5_icrh_routines_MOD_get_rf_wave_local_v2(
    real *R, real *z, real *rho_tor, real *theta, void **cptr_rfglobal,
    real *e_plus_real, real *e_minus_real, real *e_plus_imag,
    real *e_minus_imag);

void __ascot5_icrh_routines_MOD_eval_resonance_function(
    void **cptr_marker, void **cptr_rfglobal, real *omega_res, int *nharm);

void __ascot5_icrh_routines_MOD_print_marker_stuff(void **marker_pointer);

void __ascot5_icrh_routines_MOD_print_mem_stuff(void **mem_pointer);
#endif

void rfof_init(rfof_data *rfof_data)
{
#ifdef RFOF
    const char xml_filename[128] = RFOF_CODEPARAM_XML;
    int xml_filename_len = strlen(xml_filename);
    int *xml_filename_len_ptr = &xml_filename_len;
    __ascot5_icrh_routines_MOD_call_initev_excl_marker_stuff(
        xml_filename, &xml_filename_len_ptr, &(rfof_data->rfglobal),
        &(rfof_data->rfof_input_params));
#endif
}

void rfof_free(rfof_data *rfof)
{
#ifdef RFOF
    __ascot5_icrh_routines_MOD_call_deallocate_rfof_input_param(
        &rfof->rfof_input_params);
    __ascot5_icrh_routines_MOD_call_deallocate_rfglobal(&rfof->rfglobal);
#endif
}

void rfof_set_up(rfof_marker *rfof_mrk, rfof_data *rfof)
{
#ifdef RFOF
    for (int i = 0; i < NSIMD; i++)
    {
        /* Initialize marker data with dummy values */
        real dummy = -999.0, *ptr = &dummy;  /* -999.0 should be used for real
                                                dummy when dealing with RFOF */
        int dummyint = 0, *iptr = &dummyint; /* -999 should be used for int
                                                dummy with RFOF but in this
                                                case it is irrelevant as iptr
                                                only goes to RFOF marker id  */
        int is_accelerated = 0;
        int is_already_allocated = 0;
        __ascot5_icrh_routines_MOD_call_set_marker_pointers(
            &(rfof_mrk->p[i]), &iptr, &ptr, &ptr, &ptr, &ptr, &ptr, &ptr, &ptr,
            &ptr, &ptr, &ptr, &ptr, &ptr, &ptr, &ptr, &ptr, &ptr, &ptr,
            &is_accelerated, &is_already_allocated);

        /* Initialize resonance history */
        __ascot5_icrh_routines_MOD_call_initialise_res_mem(
            &rfof_mrk->history_array[i], &rfof_mrk->nrow[i], &rfof_mrk->ncol[i],
            &rfof->rfglobal, &rfof->rfof_input_params);

        /* Initialize diagnostics */
        __ascot5_icrh_routines_MOD_call_initialise_diagnostics(
            &rfof->rfglobal, &(rfof_mrk->diag_array[i]));
    }
#endif
}

void rfof_tear_down(rfof_marker *rfof_mrk)
{
#ifdef RFOF
    for (int i = 0; i < NSIMD; i++)
    {
        __ascot5_icrh_routines_MOD_call_deallocate_res_mem(
            &rfof_mrk->history_array[i], &rfof_mrk->nrow[i],
            &rfof_mrk->ncol[i]);
        __ascot5_icrh_routines_MOD_call_deallocate_diagnostics(
            &rfof_mrk->diag_array[i]);
        __ascot5_icrh_routines_MOD_deallocate_marker(&rfof_mrk->p[i]);
    }
#endif
}

void rfof_clear_history(rfof_marker *rfof_mrk, int i)
{
#ifdef RFOF
    __ascot5_icrh_routines_MOD_call_reset_res_mem(
        &rfof_mrk->history_array[i], &rfof_mrk->nrow[i], &rfof_mrk->ncol[i]);
#endif
}

void rfof_resonance_check_and_kick_gc(
    particle_simd_gc *p, real *hin, real *hout, rfof_marker *rfof_mrk,
    rfof_data *rfof_data, B_field_data *Bdata)
{
#ifdef RFOF
    for (int i = 0; i < NSIMD; i++)
    {
        if (p->id[i] > 0 && p->running[i])
        {
            /* Evaluate derived quantities needed by librfof */
            /* NOTE: It is possible that the use of (multiple) physlib
               functions to evalutate some quantity introduces some error which,
               at least when cumulated, grows intolerable.                    */
            real psi, B, Ekin, vnorm, P_phi, v_par, v_perp, gyrof, tauB;
            B_field_eval_psi(
                &psi, p->r[i], p->phi[i], p->z[i], p->time[i], Bdata);
            psi *= CONST_2PI; // librfof is COCOS 13
            B = math_normc(p->B_r[i], p->B_phi[i], p->B_z[i]);
            v_par = p->ppar[i] / p->mass[i];
            Ekin = p->ppar[i] * p->ppar[i] / (2 * p->mass[i]) + p->mu[i] * B;
            vnorm = sqrt(
                (p->ppar[i] / p->mass[i]) * (p->ppar[i] / p->mass[i]) +
                2 * p->mu[i] * B / p->mass[i]);
            v_perp = sqrt(2 * p->mu[i] * B / p->mass[i]);

            P_phi = phys_ptoroid_gc(
                p->charge[i], p->r[i], p->ppar[i], psi, B, p->B_phi[i]);
            gyrof = phys_gyrofreq_ppar(
                p->mass[i], p->charge[i], p->mu[i], p->ppar[i], B);
            real q_safe = 1.0;
            real majR = 1.65; // For now AUG
            real minR = 0.6;  // AUG
            tauB = CONST_2PI * q_safe * p->r[i] / fabs(v_par) *
                   sqrt(2 * majR / minR);
            real vdriftRho = 0; // Assuming this is not needed in librfof
            real acceleration = 1.0;
            int is_accelerated = 0;
            int is_preallocated = 1;

            int dummy_id = 1;
            int *dummy_Id_ptr = &dummy_id;

            real *weight_ptr = &(p->weight[i]);
            real *r_ptr = &(p->r[i]);
            real *phi_ptr = &(p->phi[i]);
            real *z_ptr = &(p->z[i]);
            real *charge_ptr = &(p->charge[i]);
            real *mass_ptr = &(p->mass[i]);
            real *mu_ptr = &(p->mu[i]);
            real *Ekin_ptr = &Ekin;
            real *psi_ptr = &psi;
            real *speed_ptr = &vnorm;
            real *P_phi_ptr = &P_phi;
            real *v_par_ptr = &v_par;
            real *v_perp_ptr = &v_perp;
            real *gyrof_ptr = &gyrof;
            real *tauB_ptr = &tauB;
            real *vdriftRho_ptr = &vdriftRho;
            real *acc_ptr = &acceleration;

            /* Update the fields of RFOF marker */
            __ascot5_icrh_routines_MOD_call_set_marker_pointers(
                &(rfof_mrk->p[i]), &(dummy_Id_ptr), &(weight_ptr), &(r_ptr),
                &(phi_ptr), &(z_ptr), &(psi_ptr), &(charge_ptr), &(mass_ptr),
                &(Ekin_ptr), &(speed_ptr), &(mu_ptr), &(P_phi_ptr),
                &(v_par_ptr), &(v_perp_ptr), &(gyrof_ptr), &(tauB_ptr),
                &(vdriftRho_ptr), &(acc_ptr), &is_accelerated,
                &is_preallocated);

            /* Contains the change in physical quantities and the suggested
             * time-step after the kick */
            rfof_output rfof_data_pack = {
                .dmu = 0.0,
                .dvpar = 0.0,
                .de = 0.0,
                .deCumulative = 0.0,
                .dpitch = 0.0,
                .maxAcc = 0.0,
                .RFdt = 0.0,
            };

            int err = 0;
            int mpi_rank = 0; // RFOF does not work with MPI yet

            /* Ready to kick some ash (if in resonance) */
            __ascot5_icrh_routines_MOD_call_rf_kick(
                &(p->time[i]), &(hin[i]), &mpi_rank, &(rfof_mrk->p[i]),
                &(rfof_mrk->history_array[i]), &(rfof_data->rfglobal),
                &(rfof_mrk->diag_array[i]), &(rfof_data->rfof_input_params),
                &(rfof_mrk->nrow[i]), &(rfof_mrk->ncol[i]), &err,
                &rfof_data_pack);

            /* Most marker phase-space coordinates are updated automatically
             * via the pointers in rfof_mrk except ppar which we update here */
            p->ppar[i] = p->ppar[i] + p->mass[i] * (rfof_data_pack.dvpar);

            if (err == 7)
            {
                /* Overshot the resonance. Mark the suggested time-step as
                 * negative to inform ASCOT to retry the time step */
                hout[i] = -rfof_data_pack.RFdt;
            }
            else
            {
                /* Interaction was successful. The suggested time-step is
                 * a guess how long till the marker enters the resonance */
                hout[i] = rfof_data_pack.RFdt;
            }
        }
    }
#endif
}

void rfof_set_marker_manually(
    rfof_marker *rfof_mrk, int *id, real *weight, real *R, real *phi, real *z,
    real *psi, real *charge, real *mass, real *Ekin, real *vnorm, real *mu,
    real *Pphi, real *vpar, real *vperp, real *gyrof, real *vdriftRho,
    real *acc, int *is_accelerated, int *is_already_allocated)
{
#ifdef RFOF
    // TODO: Check if this is needed and implement properly.
    real bigR = 2.0;
    real tauB_target = CONST_2PI * (*charge) * (bigR) / (*vpar);
    real *tauB = &tauB_target;
    __ascot5_icrh_routines_MOD_call_set_marker_pointers(
        &rfof_mrk->p[0], &id, &weight, &R, &phi, &z, &psi, &charge, &mass,
        &Ekin, &vnorm, &mu, &Pphi, &vpar, &vperp, &gyrof, &tauB, &vdriftRho,
        &acc, is_accelerated, is_already_allocated);
#endif
}

void rfof_eval_rf_wave(
    real *e_plus_real, real *e_minus_real, real *e_plus_imag,
    real *e_minus_imag, real R, real z, rfof_data *rfof)
{
#ifdef RFOF
    /* The current implementation assumes that the ICRH wave field is given in
        R,z co-ordinates. It would be possible to give it in rho_tor,theta as
        well. Information about the choice of co-ordinates is saved in the
        rfglobal fortran struct in the fields rfglobal%coord_name_x1 and
        rfglobal_coord_name_x2. */
    real rho_tor = -999.0; // currently dummy (-999.0 is standard dummy in RFOF)
    real theta = -999.0;   // currently dummy
    __ascot5_icrh_routines_MOD_get_rf_wave_local_v2(
        &R, &z, &rho_tor, &theta, &rfof->rfglobal, e_plus_real, e_minus_real,
        e_plus_imag, e_minus_imag);
#endif
}

void rfof_eval_resonance_function(
    real *omega_res, int *nharm, rfof_marker *rfof_mrk, rfof_data *rfof)
{
#ifdef RFOF
    __ascot5_icrh_routines_MOD_eval_resonance_function(
        &rfof_mrk->p[0], &rfof->rfglobal, omega_res, nharm);
#endif
}
