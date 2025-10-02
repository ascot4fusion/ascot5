#include "RF_stix_particle_history.h"
#include "../particle.h"
#include <stdlib.h>
#include <string.h>


/**
 * @brief Initializes the RF particle history structure with the current particle state.
 * 
 * @param hist Pointer to the RF_particle_history structure to be initialized.
 * @param p Pointer to the particle_simd_gc structure containing the particle state.
 * @param nwaves Number of waves in the Stix diffusion operator.
 * @param imrk Index of the marker in the particle array.
 * @param omega Array of wave frequencies.
 * @param ntor Array of toroidal mode numbers.
 * @param lhigh Flag indicating whether high-frequency data is used.
 */
void RF_particle_history_init(RF_particle_history* hist, particle_simd_gc* p, int nwaves, real h,
                              int imrk, real* omega, int* ntor, int lhigh) {
    hist->resn = (real*) malloc(nwaves * sizeof(real) * lhigh);
    hist->resp = (real*) malloc(nwaves * sizeof(real) * lhigh);
    if(!p->running[imrk]) return; // Particle is not running, skip.

    hist->qm = p->charge[imrk] / p->mass[imrk]; // Charge/mass ratio
    
    // We need to add to the data the current values.
    for(int j = 0; j < RF_N_HISTORY; j++) {
        hist->dt[j] = h; // Current time
        hist->bnorm[j] = sqrt(p->B_r[imrk]   * p->B_r[imrk]   +
                              p->B_phi[imrk] * p->B_phi[imrk] +
                              p->B_z[imrk]   * p->B_z[imrk]);
        hist->R[j] = p->r[imrk];
        hist->rhopara[j] = p->ppar[imrk] / hist->bnorm[j];
    }

    for (int i = 0; i < nwaves; i++) {
        for(int j = 0; j < lhigh; j++) {
            real kpara = ntor[i] / p->r[imrk];
            real doppler = kpara * p->ppar[imrk] / p->mass[imrk];
            real Babs = hist->bnorm[0]; // Use the first history point for Babs
            int lres = j + 1; // Resonance level starts from 1
            hist->resn[i * lhigh + j] = omega[i] - lres * hist->qm * Babs - doppler;
            hist->resp[i * lhigh + j] = omega[i] - lres * hist->qm * Babs + doppler;
        }
    }

    hist->nwaves = nwaves;
    hist->lhigh = lhigh;

    hist->omega = (real*) malloc(nwaves * sizeof(real));
    memcpy(hist->omega, omega, nwaves * sizeof(real));
    hist->ntor = (int*) malloc(nwaves * sizeof(int));
    memcpy(hist->ntor, ntor, nwaves * sizeof(int));
}

/**
 * @brief Evaluates the number of kicks and resonance levels for a particle based on its history.
 * 
 * Evaluates the crossings for a given wave and marker whether it has crossed the resonance condition for
 * both positive and negative Doppler shifted resonance conditions.
 * 
 * @param hist Pointer to the RF_particle_history structure containing the particle's history.
 * @param p Pointer to the particle_simd_gc structure containing the particle state.
 * @param imrk Index of the marker in the particle array.
 * @param iwave Index of the wave for which the kicks are evaluated.
 * @param nkicks Pointer to an integer where the number of kicks will be stored.
 * @param lres Pointer to an integer where the resonance level will be stored.
 */
void RF_particle_eval_nkicks(RF_particle_history* hist, particle_simd_gc* p, 
                             int imrk, int iwave, int* nkicks,
                             int *lres) {

    *nkicks = 0; // Initialize the number of kicks
    *lres = -1; // Initialize l to -1
    real qm = hist->qm;
    real Babs = sqrt(p->B_r[imrk] * p->B_r[imrk] +
                     p->B_phi[imrk] * p->B_phi[imrk] +
                     p->B_z[imrk] * p->B_z[imrk]);
    
    
    // Checking if the particle has advanced...
    if((hist->bnorm[0] == hist->bnorm[1]) && \
       (hist->bnorm[1] == hist->bnorm[2])){
        *nkicks = 0;
        *lres = -1;
        return; // No change in B, no resonance crossing possible.
    }

    for(int j = 0; j < hist->lhigh; j++){
        real resn_prv = hist->resn[iwave * hist->lhigh + j];
        real resp_prv = hist->resp[iwave * hist->lhigh + j];

        real kpara = hist->ntor[iwave] / p->r[imrk];
        real doppler = kpara * p->ppar[imrk] / p->mass[imrk];

        // Evaluating current resonance locations
        int l = j + 1; // Resonance level starts from 1
        hist->resn[iwave * hist->lhigh + j] = hist->omega[iwave] - l * qm * Babs - doppler;
        hist->resp[iwave * hist->lhigh + j] = hist->omega[iwave] - l * qm * Babs + doppler;

        // Check if the resonance condition has changed sign
        if ((resn_prv * hist->resn[iwave * hist->lhigh + j] < 0.0)){
            (*nkicks)++;
            *lres = l; // Update the resonance level

            // Checking the positive resonance condition
            if(resp_prv * hist->resp[iwave * hist->lhigh + j] < 0.0) (*nkicks)++;
            break;
        }
    }
}

/**
 * @brief Frees the memory allocated for the RF particle history structure.
 * 
 * @param hist Pointer to the RF_particle_history structure to be freed.
 */
void RF_particle_history_free(RF_particle_history* hist) {
    if (hist->resn) free(hist->resn);
    if (hist->resp) free(hist->resp);
    if (hist->omega) free(hist->omega);
    if (hist->ntor) free(hist->ntor);

    hist->resn = NULL;
    hist->resp = NULL;
    hist->omega = NULL;
    hist->ntor = NULL;
}

/**
 * @brief Updates the RF particle history with the current state of the particle.
 * 
 * @param hist 
 * @param p 
 * @param imrk 
 */
void RF_particle_history_update(RF_particle_history* hist, particle_simd_gc* p, int imrk, real h){
    // Shift history data to make room for the new entry
    for(int j = RF_N_HISTORY - 1; j > 0; j--) {
        hist->dt[j] = hist->dt[j - 1];
        hist->bnorm[j] = hist->bnorm[j - 1];
        hist->R[j] = hist->R[j - 1];
        hist->rhopara[j] = hist->rhopara[j - 1];
    }

    // Add the current state as the newest entry
    hist->dt[0] = h; // Current time step.
    hist->bnorm[0] = sqrt(p->B_r[imrk] * p->B_r[imrk] +
                          p->B_phi[imrk] * p->B_phi[imrk] +
                          p->B_z[imrk] * p->B_z[imrk]);
    hist->R[0] = p->r[imrk];
    hist->rhopara[0] = p->ppar[imrk] / hist->bnorm[0];
}