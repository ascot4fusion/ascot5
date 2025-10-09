#include "diag_energy_exchange.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "../B_field.h"
#include "../consts.h"
#include "../print.h"
#include "../plasma.h"
#include "../mhd/mhd_stat.h"

#define PSI_POINTS 1024 /**< Number of psi points for interpolation. */
#ifndef M_PI
#define M_PI 3.14159265358979323846 /**< Definition of pi if not defined. */
#endif

real simpson(real* f, int n, real dx) {
    real sum = f[0] + f[n-1];
    for (int i = 1; i < n-1; i++) {
        sum += (i % 2 == 0 ? 2 : 4) * f[i];
    }
    return sum * dx / 3.0;
}

/**
 * @brief Initialize the energy exchange diagnostic data structure.
 *
 * This function initializes the diagnostic data structure for energy exchange
 * calculations, allocating memory for source terms and energy change.
 *
 * @param data Pointer to the diagnostic data structure to be initialized.
 * @param nprt Number of particles.
 * @param boozerdata Pointer to Boozer coordinate data.
 * @param B_data Pointer to magnetic field data.
 * @param mhd_data Pointer to MHD mode data.
 */
void diag_energy_exchange_init(diag_energy_exchange_data* data, 
                               int nprt, boozer_data* boozerdata,
                               B_field_data* B_data,
                               mhd_stat_data* mhd_data, 
                               plasma_data* plasma_data) {

    if (data == NULL || boozerdata == NULL || \
        B_data == NULL || mhd_data == NULL  || plasma_data == NULL) {
        print_err("Error: Null pointer passed to diag_energy_exchange_init.\n");
        return;
    }

    if(data->enabled == 1) diag_energy_exchange_free(data);

    data->n_modes = mhd_data->n_modes;
    data->nprt = nprt;
    data->mhd_data = mhd_data;
    data->boozerdata = boozerdata;
    data->B_data = B_data;
    data->plasma_data = plasma_data;

    // The energy exchange diagnostic, that can be used for understanding
    // mode evolution requires that the MHD has a non-zero frequency. 
    // We will set here a flag indicating which modes require this mode evolution.
    if (data->n_modes <= 0) {
        print_err("Error: No MHD modes defined.\n");
        return;
    }
    data->is_mode_evol = (int*)malloc(data->n_modes * sizeof(int));

    for (int i = 0; i < data->n_modes; i++) {
        if (mhd_data->omega_nm[i] == 0.0) {
            data->is_mode_evol[i] = 0; // No evolution for this mode
        } else {
            data->is_mode_evol[i] = 1; // Evolution is required
        }
    }

    // Allocate memory for source terms and energy change
    data->S1 = (real*)malloc(data->n_modes * nprt * sizeof(real));
    data->S2 = (real*)malloc(data->n_modes * nprt * sizeof(real));
    data->dEnergy = (real*)malloc(data->n_modes * nprt * sizeof(real));
    memset(data->S1, 0, data->n_modes * nprt * sizeof(real));
    memset(data->S2, 0, data->n_modes * nprt * sizeof(real));
    memset(data->dEnergy, 0, data->n_modes * nprt * sizeof(real));

    // Finally we get from the plasma data the mass of the main ion species.
    data->thrmass = plasma_get_species_mass(plasma_data)[1];

    data->enabled = 1; // Set enabled flag
}

/**
 * @brief Free allocated resources for the energy exchange diagnostic data.
 *
 * This function frees the memory allocated for the diagnostic data structure.
 *
 * @param data Pointer to the diagnostic data structure to be freed.
 */
void diag_energy_exchange_free(diag_energy_exchange_data* data) {
    if(data == NULL) return;
    if(data->enabled == 0) return;
    if (data->S1) free(data->S1);
    if (data->S2) free(data->S2);
    if (data->dEnergy) free(data->dEnergy);
    data->mhd_data = NULL;
    data->boozerdata = NULL;
    data->n_modes = 0;
    data->nprt = 0;

    data->enabled = 0; // Reset enabled flag
}

/**
 * @brief Offload data to the accelerator.
 *
 * This function prepares the diagnostic data for offloading to an accelerator.
 * It is a placeholder for future implementation.
 *
 * @param data Pointer to the diagnostic data structure.
 */
void diag_energy_exchange_offload(diag_energy_exchange_data* data) {
    GPU_MAP_TO_DEVICE(data->S1[0: data->n_modes * data->nprt], 
                    data->S2[0: data->n_modes * data->nprt], 
                    data->dEnergy[0: data->n_modes * data->nprt]);
}


void diag_energy_exchange_onload(diag_energy_exchange_data* data) {
    GPU_MAP_FROM_DEVICE(data->S1[0: data->n_modes * data->nprt], 
                        data->S2[0: data->n_modes * data->nprt], 
                        data->dEnergy[0: data->n_modes * data->nprt]);

}



/**
 * @brief Energy exchange calculation for full-orbit particles.
 * 
 * This function updates the energy exchange diagnostic data based on the
 * interaction between a full-orbit particle and the MHD modes, based on the
 * end and start conditions of the particle for a given time step.
 * 
 * Since the total energy difference in the p_i and p_f may also include other
 * factors, such as drag and diffusion, RF, etc and it is the total energy 
 * difference, not discriminating by mode numbers (n, m), this function will
 * give instead a comprehensive diagnostic on the energy exchange.
 * 
 * To compute the energy exchange and the current coupling terms, we will use
 * the trapezoidal rule to integrate over the time step of the particle.
 * 
 * @todo This should be done more efficiently, there is no need to evaluate again
 * the initial and final electric fields: just the final one is enough.
 * 
 * @param data Pointer to the diagnostic data structure.
 * @param p_f Pointer to the full-orbit particle at the end of the time step.
 * @param p_i Pointer to the full-orbit particle at the start of the time step.
 */
void diag_energy_exchange_update_fo(diag_energy_exchange_data* data, 
                                    particle_simd_fo* p_f, 
                                    particle_simd_fo* p_i){

    // Loop over all particles
    GPU_PARALLEL_LOOP_ALL_LEVELS
    for(int i = 0; i < p_f->n_mrk; i++) {
        int index = p_f->index[i]; // Get the global index of the particle
        if(p_i->running[i] == 0) continue; // Skip if the particle is not running
        // Electric field at the particle position.
        real Er_i, Ez_i, Ephi_i;
        real Er_f, Ez_f, Ephi_f;

        // Change in electric field at the particle position.
        real dEr_i, dEz_i, dEphi_i;
        real dEr_f, dEz_f, dEphi_f;

        // Auxiliary variables for the particle.
        real pert_field[14];

        // Common values for the particle.
        real dt = p_f->time[i] - p_i->time[i];
        real qm = p_i->charge[i] / p_i->mass[i]; // Charge/mass ratio

        // Computing the Alfven velocity at the particle position.
        real iondens_f, iondens_i;
        plasma_eval_dens(&iondens_f, p_f->rho[i], 
            p_f->r[i], p_f->phi[i], p_f->z[i], p_f->time[i], 1, data->plasma_data);
        plasma_eval_dens(&iondens_i, p_i->rho[i], 
            p_i->r[i], p_i->phi[i], p_i->z[i], p_i->time[i], 1, data->plasma_data);

        // Evaluating the magnetic field strenght at the particle position.
        real B[3];
        real Babs2_f, Babs2_i;
        B_field_eval_B(B, p_f->r[i], p_f->phi[i], p_f->z[i], p_f->time[i], data->B_data);
        Babs2_f = B[0] * B[0] + B[1] * B[1] + B[2] * B[2];
        B_field_eval_B(B, p_i->r[i], p_i->phi[i], p_i->z[i], p_f->time[i], data->B_data);
        Babs2_i = B[0] * B[0] + B[1] * B[1] + B[2] * B[2];

        // Getting the poloidal field strenght times the radius.
        real R2Bpol2_f = p_f->r[i] * p_f->r[i] * (B[0] * B[0] + B[2] * B[2]);
        real R2Bpol2_i = p_i->r[i] * p_i->r[i] * (B[0] * B[0] + B[2] * B[2]);

        // Evaluating the Alfven velocity at the particle position.
        real vA2_f = Babs2_f / (iondens_f * data->thrmass) * R2Bpol2_f;
        real vA2_i = Babs2_i / (iondens_i * data->thrmass) * R2Bpol2_i;
        
        for (int j = 0; i < data->n_modes; i++) {
            if( data->is_mode_evol[j] == 0 ) continue; // Skip if no evolution for this mode
            // 1. Getting the current electric field at the 
            //    particle position.
            mhd_stat_eval_perturbations_dt(pert_field, p_f->r[i], p_f->phi[i], 
                p_f->z[i], p_f->time[i], 1, j, data->boozerdata, 
                data->mhd_data, data->B_data);
            Er_f = pert_field[3];
            Ephi_f = pert_field[4];
            Ez_f = pert_field[5];
            dEr_f = pert_field[10];
            dEphi_f = pert_field[11];
            dEz_f = pert_field[12];

            mhd_stat_eval_perturbations_dt(pert_field, p_i->r[i], p_i->phi[i], 
                p_i->z[i], p_i->time[i], 1, j, data->boozerdata, 
                data->mhd_data, data->B_data);
            Er_i = pert_field[3];
            Ephi_i = pert_field[4];
            Ez_i = pert_field[5];
            dEr_i = pert_field[10];
            dEphi_i = pert_field[11];
            dEz_i = pert_field[12];
            
            real Edotv_f = Er_f * p_f->p_r[i] + Ephi_f * p_f->p_phi[i] + Ez_f * p_f->p_z[i];
            real Edotv_i = Er_i * p_i->p_r[i] + Ephi_i * p_i->p_phi[i] + Ez_i * p_i->p_z[i];
            real dEdotv_f = dEr_f * p_f->p_r[i] + dEphi_f * p_f->p_phi[i] + dEz_f * p_f->p_z[i];
            real dEdotv_i = dEr_i * p_i->p_r[i] + dEphi_i * p_i->p_phi[i] + dEz_i * p_i->p_z[i];
            
            // 1. Calculating the energy change for the particle.
            real deltaE = 0.5 * Edotv_f + 0.5 * Edotv_i;
            deltaE *= (p_f->time[i] - p_i->time[i]) * qm;
            data->dEnergy[j * data->nprt + i] += deltaE;

            // 2.3. Computing the source term S1.
            real omega = data->mhd_data->omega_nm[j];
            real factor1 = dt * qm / (2.0 * omega * omega);
            data->S1[j * data->nprt + index] += 0.5 * (vA2_f * Edotv_f * p_f->weight[i] + \
                                                       vA2_i * Edotv_i * p_i->weight[i] ) * factor1;
            // 2.4. Computing the source term S2.
            real factor2 = factor1 / omega;
            data->S2[j * data->nprt + index] += 0.5 * (vA2_f * dEdotv_f * p_f->weight[i] + \
                                                   vA2_i * dEdotv_i * p_i->weight[i] ) * factor2;
        }
    }

}


/**
* @brief Energy exchange calculation for guiding center particles.
*
* This function updates the energy exchange diagnostic data based on the
* interaction between a guiding center particle and the MHD modes, based on the
* end and start conditions of the particle for a given time step.
*
* Since the total energy difference in the p_i and p_f may also include other
* factors, such as drag and diffusion, RF, etc and it is the total energy
* difference, not discriminating by mode numbers (n, m), this function will
* give instead a comprehensive diagnostic on the energy exchange.
*
* To compute the energy exchange and the current coupling terms, we will use
* the trapezoidal rule to integrate over the time step of the particle.
*
* The version of the energy exchange here deployed is the one from the book
* R. B. White, "The theory of toroidally confined plasmas", 3rd edition,
* Oxford University Press, 2014, section 6.9.
*
* @todo This should be done more efficiently, there is no need to evaluate again
* the initial and final electric fields: just the final one is enough.
* @param data Pointer to the diagnostic data structure.
* @param p_f Pointer to the guiding center particle at the end of the time step.
* @param p_i Pointer to the guiding center particle at the start of the time step.
*/      
void diag_energy_exchange_update_gc(diag_energy_exchange_data* data, 
                                     particle_simd_gc* p_f, 
                                     particle_simd_gc* p_i){
    // Loop over all particles
    GPU_PARALLEL_LOOP_ALL_LEVELS
    for(int i = 0; i < NSIMD; i++) {
        int index = p_f->index[i]; // Get the global index of the particle
        if(p_i->running[i] == 0) continue; // Skip if the particle is not running

        // Common values for the particle.
        real dt = p_f->time[i] - p_i->time[i];
        real qm = p_i->charge[i] / p_i->mass[i]; // Charge/mass ratio
        real w0 = p_i->weight[i];
        real w1 = p_f->weight[i];

        // Computing the Alfven velocity at the particle position.
        real iondens_f, iondens_i;
        plasma_eval_dens(&iondens_f, p_f->rho[i], 
            p_f->r[i], p_f->phi[i], p_f->z[i], p_f->time[i], 1, data->plasma_data);
        plasma_eval_dens(&iondens_i, p_i->rho[i], 
            p_i->r[i], p_i->phi[i], p_i->z[i], p_i->time[i], 1, data->plasma_data);

        // Evaluating the magnetic field strenght at the particle position.
        real B[3];
        real Babs2_f, Babs2_i;
        B_field_eval_B(B, p_f->r[i], p_f->phi[i], p_f->z[i], p_f->time[i], data->B_data);
        Babs2_f = B[0] * B[0] + B[1] * B[1] + B[2] * B[2];
        B_field_eval_B(B, p_i->r[i], p_i->phi[i], p_i->z[i], p_i->time[i], data->B_data);
        Babs2_i = B[0] * B[0] + B[1] * B[1] + B[2] * B[2];

        // Computing the poloidal field strenght times the radius.
        real R2Bpol2_f = p_f->r[i] * p_f->r[i] * (B[0] * B[0] + B[2] * B[2]);
        real R2Bpol2_i = p_i->r[i] * p_i->r[i] * (B[0] * B[0] + B[2] * B[2]);

        // Evaluating the Alfven velocity at the particle position.
        // Technically we are evaluating the square of the Alfven velocity time mu_0.
        real vA2_f = Babs2_f / (iondens_f * data->thrmass);
        real vA2_i = Babs2_i / (iondens_i * data->thrmass);
        real Babs_f = sqrt(Babs2_f);
        real Babs_i = sqrt(Babs2_i);

        // Evaluating the parallel gyroradius at the particle position.
        real vpar_f = p_f->ppar[i] / p_f->mass[i];
        real vpar_i = p_i->ppar[i] / p_i->mass[i];

        real mhd_dmhd[10];
        for (int j = 0; j < data->n_modes; j++) {
            if( data->is_mode_evol[j] == 0 ) continue; // Skip if no evolution for this mode
            // Getting the mode data.
            mhd_stat_eval(mhd_dmhd, p_f->r[i], p_f->phi[i], 
                          p_f->z[i], p_f->time[i], j, data->boozerdata, 
                          data->mhd_data, data->B_data);

            real alpha = mhd_dmhd[0]; // Magnetic potential
            real phipot = mhd_dmhd[5]; // Electric potential
            real alpha_dot = mhd_dmhd[1]; // Time derivative of magnetic potential
            real phipot_dot = mhd_dmhd[6]; // Time derivative of electric potential
            real omega = data->mhd_data->omega_nm[j]; // Toroidal frequency
            real mass = p_i->mass[i]; // Particle mass
            real charge = p_i->charge[i]; // Particle charge

            // 1. Evaluating the mode energy exchange.
            real dE_f = charge * (vpar_f * alpha * Babs_f - phipot) * omega;
            real dE_i = charge * (vpar_i * alpha * Babs_i - phipot) * omega;
            data->dEnergy[j * data->nprt + index] += 0.5 * (dE_f + dE_i) * dt;

            // 2. Computing the source term S1.
            data->S1[j * data->nprt + index] += 0.5 * (vA2_f * dE_f * w1 * R2Bpol2_f + \
                                                       vA2_i * dE_i * w0 * R2Bpol2_i) * dt;

            // 3. Computing the source term S2.
            real dEdotv_f = charge * (vpar_f * alpha_dot * Babs_f - phipot_dot) * omega*omega;
            real dEdotv_i = charge * (vpar_i * alpha_dot * Babs_i - phipot_dot) * omega*omega;
            data->S2[j * data->nprt + index] += 0.5 * (vA2_f * dEdotv_f * w1 * R2Bpol2_f + \
                                                       vA2_i * dEdotv_i * w0 * R2Bpol2_i) * dt;
        }
    }
}

/**
 * @brief Compute the total source terms for the energy exchange diagnostic, 
 * combining all the particle contributions.
 * 
 * This function computes the total source terms for the energy exchange
 * diagnostic by summing the contributions from all particles.
 * 
 * @param data Pointer to the diagnostic data structure.
 * @param S1 Output array for the source term S1.
 * @param S2 Output array for the source term S2.
 * @param clear Whether to clear the source terms after the summation.
 */

void diag_energy_exchange_compact(diag_energy_exchange_data* data, 
                                  real* S1, real* S2, int clear){
    // Check if the diagnostic is enabled
    if (data->enabled == 0) return;

    // Check if the data and source term pointers are valid
    if (data == NULL || S1 == NULL || S2 == NULL) return;

    // Summing all the source terms for all the particles.
    GPU_PARALLEL_LOOP_ALL_LEVELS_REDUCTION(S1, S2)
    for (int i = 0; i < data->nprt; i++) {
        for (int j = 0; j < data->n_modes; j++) {
            real ismodeevol = (real) data->is_mode_evol[j];
            S1[j] += data->S1[j * data->nprt + i] * ismodeevol;
            S2[j] += data->S2[j * data->nprt + i] * ismodeevol;
        }

    }
    
    // If clear is true, we reset the source terms to zero.
    if (clear) {
        memset(data->S1, 0, data->n_modes * data->nprt * sizeof(real));
        memset(data->S2, 0, data->n_modes * data->nprt * sizeof(real));
        memset(data->dEnergy, 0, data->n_modes * data->nprt * sizeof(real));
    }

}

/** 
 * @brief Update the number of particles in the energy exchange diagnostic.
 * This function updates the number of particles in the diagnostic data structure
 * and reallocates memory for the source terms and energy change arrays.
 * @param data Pointer to the diagnostic data structure.
 * @param nprt New number of particles.
*/
void diag_energy_exchange_update_nprt(diag_energy_exchange_data* data, 
                                        int nprt) {
    // Check if the diagnostic is enabled
    if (data->enabled == 0) return;

    // Removing the data from the GPU.
    GPU_MAP_DELETE_DEVICE(data->S1[0: data->n_modes * data->nprt], 
                          data->S2[0: data->n_modes * data->nprt], 
                          data->dEnergy[0: data->n_modes * data->nprt]);

    // Free the previous source terms and reallocate for the new number of particles
    free(data->S1);
    free(data->S2);
    free(data->dEnergy);

    data->nprt = nprt;
    data->S1 = (real*)malloc(data->n_modes * nprt * sizeof(real));
    data->S2 = (real*)malloc(data->n_modes * nprt * sizeof(real));
    data->dEnergy = (real*)malloc(data->n_modes * nprt * sizeof(real));
    
    memset(data->S1, 0, data->n_modes * nprt * sizeof(real));
    memset(data->S2, 0, data->n_modes * nprt * sizeof(real));
    memset(data->dEnergy, 0, data->n_modes * nprt * sizeof(real));

    // Offloading the data to the GPU again.
    GPU_MAP_TO_DEVICE(data->S1[0: data->n_modes * data->nprt], 
                      data->S2[0: data->n_modes * data->nprt], 
                      data->dEnergy[0: data->n_modes * data->nprt]);
}