/**
 * @file rfof_interface.c
 * @brief Contains the function to be called during the simulation when using
 * ICRH. Requires librfof.so library which contains the Fortran routines.
**/
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include "physlib.h"
#include "consts.h"
#include "particle.h"
#include "rfof.h"

/**
 * @brief
 */
typedef struct rfof_output {
    double dmu;          /**< Change in magnetic moment due to ICRH kick.     */
    double dvpar;        /**< Change in parallel velocity component due to ICRH
                              kick.                                           */
    double de;           /**< Change in energy due to a single ICRH kick [J]. */
    double deCumulative; /**< Change in energy due to possibly several ICRH
                              kicks during an orbit time step [J]             */
    double dpitch;       /**< Change in pitch due to ICRH kick                */
    double maxAcc;       /**< Maximum acceleration allowed by RFOF            */
    double RFdt;         /**< time step suggested by RFOF                     */
} rfof_output;

#ifdef RFOF
void __ascot5_icrh_routines_MOD_call_initev_excl_marker_stuff(
    const char* xml_filename, int **xml_filename_len, void** cptr_rfglobal,
    void** cptr_rfof_input_params, int* n_waves, int* n_modes);
void __ascot5_icrh_routines_MOD_call_initialise_res_mem(void** cptr_mem,
    int* cptr_mem_shape_i, int* cptr_mem_shape_j, void** cptr_rfglobal,
    void** cptr_rfof_input_param);
void __ascot5_icrh_routines_MOD_call_initialise_diagnostics(
    void** cptr_RFglobal, void** cptr_diagno);

void __ascot5_icrh_routines_MOD_call_set_marker_pointers(void** cptr_marker,
    int** id, real** weight, real** R, real** phi, real** z, real** psi,
    real** charge, real** mass, real** Ekin, real** velocity, real** mu,
    real** pphicanonical, real** vpar, real** vperp, real** gyrof, real** tauB,
    real** vdriftRho, real** acc, int* isOrbitTimeAccelerated,
    int* is_already_allocated);

void __ascot5_icrh_routines_MOD_call_rf_kick(double*time, double*dtin,
    int* mpi_rank, void** cptr_marker, void** cptr_mem, void** cptr_rfglobal,
    void** cptr_rfdiagno, void** cptr_rfof_input, int* mem_shape_i,
    int* mem_shape_j, int *rfof_err, rfof_output* out,
    real* cptr_de_rfof_during_step);

void __ascot5_icrh_routines_MOD_call_reset_res_mem(void** rfof_mem_pointer,
    int* mem_shape_i, int* mem_shape_j);

void __ascot5_icrh_routines_MOD_call_deallocate_rfof_input_param(
    void** cptr_rfof_input_param);
void __ascot5_icrh_routines_MOD_call_deallocate_rfglobal(void** cptr_rfglobal);
void __ascot5_icrh_routines_MOD_call_deallocate_res_mem(void** cptr_res_mem,
    int* cptr_mem_shape_i, int* cptr_mem_shape_j);
void __ascot5_icrh_routines_MOD_call_deallocate_diagnostics(void** cptr_diagno);
void __ascot5_icrh_routines_MOD_deallocate_marker(void** cptr_rfof_marker);

void __ascot5_icrh_routines_MOD_get_rf_wave_local_v2(real* R, real* z,
    real* rho_tor, real* theta, void** cptr_rfglobal, real* e_plus_real,
    real* e_minus_real, real* e_plus_imag, real* e_minus_imag);

void __ascot5_icrh_routines_MOD_eval_resonance_function(void** cptr_marker,
    void** cptr_rfglobal, real* omega_res, int* nharm);

void __ascot5_icrh_routines_MOD_print_marker_stuff(void** marker_pointer);

void __ascot5_icrh_routines_MOD_print_mem_stuff(void** mem_pointer);
#endif

/**
 * @brief Initialise input data.
 *
 * Reads the ICRH (RFOF) inputs (xml, xsd, ASCII) and initialises the wave
 * field.
 *
 * @param rfof_data pointer to the RFOF data structure
 */
void rfof_init(rfof_data* rfof_data) {
#ifdef RFOF
    const char xml_filename[128] = RFOF_CODEPARAM_XML;
    int xml_filename_len = strlen(xml_filename);
    int*xml_filename_len_ptr = &xml_filename_len;
    int n_waves = -999;
    int n_modes = -999;
    __ascot5_icrh_routines_MOD_call_initev_excl_marker_stuff(xml_filename,
        &xml_filename_len_ptr, &(rfof_data->rfglobal),
        &(rfof_data->rfof_input_params), &n_waves, &n_modes);
    rfof_data->n_waves = n_waves;
    rfof_data->n_modes = n_modes;

    rfof_data->summed_timesteps = 0.0;

    /* Allocate memory for the 1d array that contains the accumulated energy
    changes for each wave and each mode in each wave */
    rfof_data->dE_RFOF_modes_and_waves = (real*)calloc(n_waves * n_modes, sizeof(real));
#endif
}

/**
 * @brief Deallocates the rfof_input_param struct on the fortran side.
 *
 * There exists only one copy of this struct and therefore it is to be
 * deallocated in the simulate.c after the loop is completed.
 */
void rfof_free(rfof_data* rfof) {
#ifdef RFOF
    __ascot5_icrh_routines_MOD_call_deallocate_rfof_input_param(
        &rfof->rfof_input_params);
    __ascot5_icrh_routines_MOD_call_deallocate_rfglobal(&rfof->rfglobal);
    free(rfof->dE_RFOF_modes_and_waves);
#endif
}

/**
 * @brief Initialises resonance history, diagnostics, and the marker struct
 *
 * This function is to be called before the simulation loop.
 *
 * @param rfof_mrk pointer to the local RFOF marker data
 * @param rfof pointer to the shared RFOF data
 */
void rfof_set_up(rfof_marker* rfof_mrk, rfof_data* rfof) {
#ifdef RFOF
    for(int i=0; i< NSIMD; i++) {
        /* Initialize marker data with dummy values */
        real dummy = -999.0, *ptr = &dummy;   /* -999.0 should be used for real
                                                 dummy when dealing with RFOF */
        int dummyint = 0, *iptr = &dummyint;  /* -999 should be used for int
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
            &rfof_mrk->history_array[i], &rfof_mrk->nrow[i],
            &rfof_mrk->ncol[i], &rfof->rfglobal,
            &rfof->rfof_input_params);

        /* Initialize diagnostics */
        __ascot5_icrh_routines_MOD_call_initialise_diagnostics(
            &rfof->rfglobal, &(rfof_mrk->diag_array[i]));
    }
#endif
}

/**
 * @brief Deallocates the data structs used by the RFOF marker simulation data
 *
 * This function is to be called after the simulation loop.
 *
 * @param rfof_mrk pointer to the RFOF marker simulation data
 */
void rfof_tear_down(rfof_marker* rfof_mrk) {
#ifdef RFOF
    for(int i=0; i< NSIMD; i++) {
        __ascot5_icrh_routines_MOD_call_deallocate_res_mem(
            &rfof_mrk->history_array[i], &rfof_mrk->nrow[i],
            &rfof_mrk->ncol[i]);
        __ascot5_icrh_routines_MOD_call_deallocate_diagnostics(
            &rfof_mrk->diag_array[i]);
        __ascot5_icrh_routines_MOD_deallocate_marker(&rfof_mrk->p[i]);
    }
#endif
}

/**
 * @brief Clears resonance history of an RFOF marker
 *
 * History should be cleared whenever a marker finishes simulation. The new
 * marker cannot receive ICRH kicks during the first two time steps as its
 * resonance history must have at least two data points stored to estimate
 * the resonance location.
 *
 * @param rfof_mrk pointer to the RFOF marker simulation data
 * @param i index in the NSIMD array for the marker whose history is cleared
 */
void rfof_clear_history(rfof_marker* rfof_mrk, int i) {
#ifdef RFOF
    __ascot5_icrh_routines_MOD_call_reset_res_mem(
        &rfof_mrk->history_array[i], &rfof_mrk->nrow[i], &rfof_mrk->ncol[i]);
#endif
}


void rfof_update_energy_array_of_the_process(rfof_data* rfof_data,
    real** energy_arrays_for_NSIMD_markers,
    real* accumulated_time_for_NSIMD_markers, int* cycle_array) {
#ifdef RFOF
    /* Update the RFOF energy and time data in the rfof_data struct for those
     markers that reached the end condition just not */
        for (int i = 0; i < NSIMD; i++){
            if (cycle_array[i] != 0) {
                // Add the particle's summed dt to the MPI process sum
                #pragma omp critical
                {
                    rfof_data->summed_timesteps += accumulated_time_for_NSIMD_markers[i];
                    accumulated_time_for_NSIMD_markers[i] = 0.0;
                }

                // Add the particle's summed E's to the MPI process sum
                #pragma omp critical
                {
                    for (int j = 0; j < rfof_data->n_waves*rfof_data->n_modes; j++) {
                        // As a safety measure, NaNs are excluded here.
                        if (!isnan(energy_arrays_for_NSIMD_markers[i][j]-energy_arrays_for_NSIMD_markers[i][j])) {
                            rfof_data->dE_RFOF_modes_and_waves[j] += energy_arrays_for_NSIMD_markers[i][j];
                        } else {
                            printf("Trying to update accumulated RF energy with a nan\n");
                        }
                        energy_arrays_for_NSIMD_markers[i][j] = 0.0;
                    }
                }
            }
        }
#endif
}





/**
 * @brief Check if the marker is in resonance and apply kick
 *
 * 1. Updates the fields of the rfof_marker based on the given input
 * ascot_marker.
 *
 * 2. Calls the "kick" function, which
 *      a) Checks resonance condition and
 *      b) if in resonance, kicks marker and updates velocity of the rfof marker
 *          and consequently also ascot marker (as only pointers are passed when
 *          creating the rfof marker).
 *
 * The time step can fail if the marker overshoots the resonance.
 *
 * @param p pointer to marker simulation struct
 * @param hin current time step
 * @param hout suggestion for the next time step with negative sign
 * indicating a failed step
 * @param rfof_mrk pointer to the rfof marker simulation struct
 * @param rfof_data pointer to the shared rfof data
 * @param Bdata pointer to the magnetic field data needed to evaluate psi
 */
void rfof_resonance_check_and_kick_gc(
    particle_simd_gc* p, real* hin, real* hout, rfof_marker* rfof_mrk,
    rfof_data* rfof_data, B_field_data* Bdata, real** de_rfof_during_step) {
#ifdef RFOF
    for(int i=0; i<NSIMD; i++) {
        if(p->id[i] > 0 && p->running[i]) {
            a5err errflag = 0;

            /* Evaluate derived quantities needed by librfof */
            /* NOTE: It is possible that the use of (multiple) physlib
               functions to evalutate some quantity introduces some error which,
               at least when cumulated, grows intolerable.                    */
            real psi, B, Ekin, vnorm, P_phi, v_par, v_perp, gyrof,
            tauB;
            errflag = B_field_eval_psi(&psi, p->r[i], p->phi[i], p->z[i],
                p->time[i], Bdata);
            psi *= CONST_2PI; // librfof is COCOS 13
            B = math_normc(p->B_r[i], p->B_phi[i], p->B_z[i]);
            //gamma  = physlib_gamma_ppar(p->mass[i], p->mu[i], p->ppar[i], B);
            //Ekin   = physlib_Ekin_gamma(p->mass[i], gamma);
            //vnorm  = physlib_vnorm_gamma(gamma);

	        /* The relativisti formulas above did not work too well for low
            speeds but started introducing numerical errors. Hence, assume
            non-relativistic speeds. */
	        v_par = p->ppar[i]/p->mass[i];
	        Ekin = p->ppar[i]*p->ppar[i]/(2*p->mass[i]) + p->mu[i]*B;
            vnorm = sqrt( (p->ppar[i]/p->mass[i])*(p->ppar[i]/p->mass[i]) + 2*p->mu[i]*B/p->mass[i] );
            v_perp = sqrt( 2 * p->mu[i]*B/p->mass[i] );


	        P_phi  = phys_ptoroid_gc(p->charge[i], p->r[i], p->ppar[i], psi, B,
                p->B_phi[i]);
            gyrof  = phys_gyrofreq_ppar(p->mass[i], p->charge[i], p->mu[i],
                                        p->ppar[i], B);

            // Check for negative values and nans
            if (Ekin < 0 || isnan(Ekin)){
                printf("----------------------------------\nNONPHYSICAL Ekin = %.5e at t=%.3e (id=%ld)\n (r,z,phi)=(%.3e, %.3e, %.3e)\n-----------------------------\n",Ekin,p->mileage[i], p->id[i], p->r[i], p->z[i], p->phi[i]);
            }

            if (p->mu[i] < 0 || isnan(p->mu[i])){
                printf("----------------------------------\nNONPHYSICAL mu = %.5e at t=%.3e (id=%ld)\n (r,z,phi)=(%.3e, %.3e, %.3e)\n-----------------------------\n",p->mu[i],p->mileage[i], p->id[i], p->r[i], p->z[i], p->phi[i]);
            }

            if (B < 0 || isnan(B)){
                printf("----------------------------------\nNONPHYSICAL B = %.5e at t=%.3e (id=%ld)\n (r,z,phi)=(%.3e, %.3e, %.3e)\n-----------------------------\n", B, p->mileage[i], p->id[i], p->r[i], p->z[i], p->phi[i]);
            }

            if (vnorm < 0 || isnan(vnorm)){
                printf("----------------------------------\nNONPHYSICAL vnorm = %.5e at t=%.3e (id=%ld)\n (r,z,phi)=(%.3e, %.3e, %.3e)\n-----------------------------\n", vnorm, p->mileage[i], p->id[i], p->r[i], p->z[i], p->phi[i]);
            }

            if (v_perp < 0 || isnan(v_perp)){
                printf("----------------------------------\nNONPHYSICAL v_perp = %.5e at t=%.3e (id=%ld)\n (r,z,phi)=(%.3e, %.3e, %.3e)\n-----------------------------\n", v_perp, p->mileage[i], p->id[i], p->r[i], p->z[i], p->phi[i]);
            }

            if (gyrof < 0 || isnan(gyrof)){
                printf("----------------------------------\nNONPHYSICAL gyrof = %.5e at t=%.3e (id=%ld)\n (r,z,phi)=(%.3e, %.3e, %.3e)\n-----------------------------\n", gyrof, p->mileage[i], p->id[i], p->r[i], p->z[i], p->phi[i]);
            }

            /* This tauB formula does not work near the bounce points of a
            banana */
	        //real q_safe = 1.0;
            //real majR = 1.65;   // For now AUG
            //real minR = 0.6;    // AUG
	        //tauB = CONST_2PI*q_safe*p->r[i]/v_par*sqrt(2*majR/minR);
	        /* => Disable tauB for now. (Note: tauB not in use in RFOF-ASCOT4
            either) */
	        tauB = -999.0;

            real vdriftRho      = 0; /* As per Seppo, there should be a bug in
                                        RFOF related to this drift velocity;
                                        thus, this drift velocity should be zero
                                        for now. */
            real acceleration   = 1.0;
            int is_accelerated  = 0;
            int is_preallocated = 1;

            int dummy_id = 1;
            int* dummy_Id_ptr = &dummy_id;

            real* weight_ptr    = &(p->weight[i]);
            real* r_ptr         = &(p->r[i]);
            real* phi_ptr       = &(p->phi[i]);
            real* z_ptr         = &(p->z[i]);
            real* charge_ptr    = &(p->charge[i]);
            real* mass_ptr      = &(p->mass[i]);
            real* mu_ptr        = &(p->mu[i]);
            real* Ekin_ptr      = &Ekin;
            real* psi_ptr       = &psi;
            real* speed_ptr     = &vnorm;
            real* P_phi_ptr     = &P_phi;
            real* v_par_ptr     = &v_par;
            real* v_perp_ptr    = &v_perp;
            real* gyrof_ptr     = &gyrof;
            real* tauB_ptr      = &tauB;
            real* vdriftRho_ptr = &vdriftRho;
            real* acc_ptr       = &acceleration;

            /* Update the fields of RFOF marker */
            __ascot5_icrh_routines_MOD_call_set_marker_pointers(
                &(rfof_mrk->p[i]),
                &(dummy_Id_ptr),
                &(weight_ptr),
                &(r_ptr),
                &(phi_ptr),
                &(z_ptr),
                &(psi_ptr),
                &(charge_ptr),
                &(mass_ptr),
                &(Ekin_ptr),
                &(speed_ptr),
                &(mu_ptr),
                &(P_phi_ptr),
                &(v_par_ptr),
                &(v_perp_ptr),
                &(gyrof_ptr),
                &(tauB_ptr),
                &(vdriftRho_ptr),
                &(acc_ptr),
                &is_accelerated,
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

            int rfof_err = 0;
            int mpi_rank = 0; // RFOF does not do anything MPI-specific for now

            real old_mu = p->mu[i];

            /* Ready to kick some ash (if in resonance) */
            __ascot5_icrh_routines_MOD_call_rf_kick(
                &(p->time[i]), &(hin[i]), &mpi_rank,
                &(rfof_mrk->p[i]), &(rfof_mrk->history_array[i]),
                &(rfof_data->rfglobal), &(rfof_mrk->diag_array[i]),
                &(rfof_data->rfof_input_params), &(rfof_mrk->nrow[i]),
                &(rfof_mrk->ncol[i]), &rfof_err, &rfof_data_pack,
                (de_rfof_during_step[i]));

            if (p->mu[i] <= 0.0 || isnan(p->mu[i])) {
                printf("\nAt time = %.3e, mrk_id = %ld\n",p->time[i], p->id[i]);
                printf("after kick, mu = %.3e\n", p->mu[i]);
                printf("before kick, mu was %.3e\n", old_mu);
            }

            if (Ekin < 0 || isnan(Ekin)){
                printf("\nAfter kick\n");
                printf("----------------------------------\nNONPHYSICAL Ekin = %.5e at t=%.3e (id=%ld)\n (r,z,phi)=(%.3e, %.3e, %.3e)\n-----------------------------\n",Ekin,p->mileage[i], p->id[i], p->r[i], p->z[i], p->phi[i]);
            }

            if (p->mu[i] < 0 || isnan(p->mu[i])){
                printf("\nAfter kick\n");
                printf("----------------------------------\nNONPHYSICAL mu = %.5e at t=%.3e (id=%ld)\n (r,z,phi)=(%.3e, %.3e, %.3e)\n-----------------------------\n",p->mu[i],p->mileage[i], p->id[i], p->r[i], p->z[i], p->phi[i]);
            }

            if (B < 0 || isnan(B)){
                printf("\nAfter kick\n");
                printf("----------------------------------\nNONPHYSICAL B = %.5e at t=%.3e (id=%ld)\n (r,z,phi)=(%.3e, %.3e, %.3e)\n-----------------------------\n", B, p->mileage[i], p->id[i], p->r[i], p->z[i], p->phi[i]);
            }

            if (vnorm < 0 || isnan(vnorm)){
                printf("\nAfter kick\n");
                printf("----------------------------------\nNONPHYSICAL vnorm = %.5e at t=%.3e (id=%ld)\n (r,z,phi)=(%.3e, %.3e, %.3e)\n-----------------------------\n", vnorm, p->mileage[i], p->id[i], p->r[i], p->z[i], p->phi[i]);
            }

            if (v_perp < 0 || isnan(v_perp)){
                printf("\nAfter kick\n");
                printf("----------------------------------\nNONPHYSICAL v_perp = %.5e at t=%.3e (id=%ld)\n (r,z,phi)=(%.3e, %.3e, %.3e)\n-----------------------------\n", v_perp, p->mileage[i], p->id[i], p->r[i], p->z[i], p->phi[i]);
            }

            if (gyrof < 0 || isnan(gyrof)){
                printf("\nAfter kick\n");
                printf("----------------------------------\nNONPHYSICAL gyrof = %.5e at t=%.3e (id=%ld)\n (r,z,phi)=(%.3e, %.3e, %.3e)\n-----------------------------\n", gyrof, p->mileage[i], p->id[i], p->r[i], p->z[i], p->phi[i]);
            }


            /* Most marker phase-space coordinates are updated automatically
             * via the pointers in rfof_mrk except ppar which we update here */
            p->ppar[i] = p->ppar[i] + p->mass[i]*(rfof_data_pack.dvpar);


            if (rfof_err == 7) {
                /* Overshot the resonance. Mark the suggested time-step as
                 * negative to inform ASCOT to retry the time step */
                hout[i] = -rfof_data_pack.RFdt;
            } else {
                /* Interaction was successful. The suggested time-step is
                 * a guess how long till the marker enters the resonance */
                hout[i] = rfof_data_pack.RFdt;
            }

            /* Error handling */
            if(errflag) {
                p->err[i]     = errflag;
                p->running[i] = 0;
            }
        }
   }
#endif
}

/**
 * @brief Explicitly set the coordinates in a marker struct
 *
 * This routine is only used for testing and in libascot. Only the first marker
 * in the NSIMD array of rfof_marker is manipulated.
 *
 * @param rfof_mrk pointer to the RFOF marker data
 * @param id marker identifier
 * @param weight marker weight [prt/s]
 * @param R major radius [m]
 * @param phi toroidal angle [rad]
 * @param z z coordinate [m]
 * @param psi poloidal flux [Wb/rad]
 * @param charge particle charge [C]
 * @param mass particle mass [kg]
 * @param Ekin kinetic energy [J]
 * @param vnorm velocity [m/s]
 * @param mu magnetic moment [J/T]
 * @param Pphi canonical toroidal angular momentum [kg*m/s]
 * @param vpar perpendicular velocity [m/s]
 * @param vperp parallel velocity [m/s]
 * @param gyrof gyrofrequency [rad/s]
 * @param vdriftRho drift velocity in radial direction
 * @param acc accleration factor, should be 1.0 always
 * @param is_accelerated flag indicating whether acceleration is used, false
 * always
 * @param is_already_allocated flag whether to allocate a new marker struct,
 * should be false always
*/
void rfof_set_marker_manually(rfof_marker* rfof_mrk, int* id,
    real* weight, real* R, real* phi, real* z, real* psi, real* charge,
    real* mass, real* Ekin, real* vnorm, real* mu, real* Pphi,
    real* vpar, real* vperp, real* gyrof, real* vdriftRho, real* acc,
    int* is_accelerated, int* is_already_allocated) {
#ifdef RFOF
    // TODO: Check if this is needed and implement properly.
    real bigR = 2.0;
    real tauB_target = CONST_2PI*(*charge)*(bigR)/(*vpar);
    real* tauB = &tauB_target;
    __ascot5_icrh_routines_MOD_call_set_marker_pointers(
        &rfof_mrk->p[0], &id, &weight, &R, &phi, &z, &psi, &charge, &mass,
        &Ekin, &vnorm, &mu, &Pphi, &vpar, &vperp, &gyrof, &tauB, &vdriftRho, &acc,
        is_accelerated, is_already_allocated);
#endif
}

/**
 * @brief Calculate the local E+ and E- values of the ICRH field
 *
 * The definitions of E+ and E- sometimes differ. In this context,
 * E+ = E_LH * (cos(phi) + i sin(phi)),
 * where E_LH is the magnitude of the left-hand polarised (rotating) component
 * and phi is its phase. That is, E_LH and phi are real. Often, however, E_LH is
 * called E+ which creates confusion. Afterall, in this function, E+ is a
 * complex number.
 *
 * @param e_plus_real Re("E+"") component of the local wave field
 * @param e_minus_real Re("E-"") component of the local wave field
 * @param e_plus_imag Im("E+"") component of the local wave field
 * @param e_minus_imag Im("E-"") component of the local wave field
 * @param R major radius coordinate [m]
 * @param z z-coordinate [m]
 * @param rfof pointer to the RFOF data structure
 */
void rfof_eval_rf_wave(
    real* e_plus_real, real* e_minus_real, real* e_plus_imag,
    real* e_minus_imag, real R, real z, rfof_data* rfof) {
#ifdef RFOF
    /* The current implementation assumes that the ICRH wave field is given in
        R,z co-ordinates. It would be possible to give it in rho_tor,theta as
        well. Information about the choice of co-ordinates is saved in the
        rfglobal fortran struct in the fields rfglobal%coord_name_x1 and
        rfglobal_coord_name_x2. */
    real rho_tor = -999.0; //currently dummy (-999.0 is standard dummy in RFOF)
    real theta = -999.0; //currently dummy
    __ascot5_icrh_routines_MOD_get_rf_wave_local_v2(
        &R, &z, &rho_tor, &theta, &rfof->rfglobal, e_plus_real, e_minus_real,
        e_plus_imag, e_minus_imag);
#endif
}

/**
 * @brief Evaluate the value of resonance function (zero at the resonance)
 *
 * This function finds the closest resonance and in addition to evaluating the
 * resonance function it also returns the corresponding harmonic value.
 *
 * The resonance is evaluated for the first marker in the NSIMD array in
 * rfof_marker fields.
 *
 * @param cptr_marker void pointer to the rfof_marker
 * @param rfof_data pointer to the RFOF data structure
 * @param omega_res evaluated value of the resonance function
 * @param nharm the number of the closest harmonic found
 */
void rfof_eval_resonance_function(
    real* omega_res, int* nharm, rfof_marker* rfof_mrk, rfof_data* rfof) {
#ifdef RFOF
    __ascot5_icrh_routines_MOD_eval_resonance_function(
        &rfof_mrk->p[0], &rfof->rfglobal, omega_res, nharm);
#endif
}
