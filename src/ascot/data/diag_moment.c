/**
 * Implements "diag_moment.h".
 */
#include "diag_moment.h"
#include "consts.h"
#include "data/bfield.h"
#include "defines.h"
#include "marker.h"
#include "parallel.h"
#include "utils/mathlib.h"
#include "utils/physlib.h"
#include <stdlib.h>

void DiagMoment_offload(DiagMoment *moment)
{
    SUPPRESS_UNUSED_WARNING(moment);
    GPU_MAP_TO_DEVICE(moment->values[0:moment->axis.n*MOMENT_NUM])
}

void DiagMoment_onload(DiagMoment *moment)
{
    SUPPRESS_UNUSED_WARNING(moment);
    GPU_MAP_FROM_DEVICE(moment->values[0:moment->axis.n*MOMENT_NUM])
}

void DiagMoment_update_go(
    DiagMoment *moment, Bfield *bfield, MarkerGyroOrbit *mrk_f,
    MarkerGyroOrbit *mrk_i)
{
    GPU_PARALLEL_LOOP_ALL_LEVELS
    for (size_t i = 0; i < mrk_f->size; i++)
    {
        int valid = 1;
        HistAxis *axis = &(moment->axis);
        real computed_moments[MOMENT_NUM];

        /*Norm of the magnetic field |B|*/
        real bnorm_f = math_normc(mrk_f->br[i], mrk_f->bphi[i], mrk_f->bz[i]);

        /*Norm of the magnetic field |B| at the start of the timestep*/
        real bnorm_i = math_normc(mrk_i->br[i], mrk_i->bphi[i], mrk_i->bz[i]);

        /*Norm of the momentum p |p|*/
        real pnorm_f = math_normc(mrk_f->p_r[i], mrk_f->p_phi[i], mrk_f->p_z[i]);

        /*Norm of the momentum p |p| at start of timestep*/
        real pnorm_i = math_normc(mrk_i->p_r[i], mrk_i->p_phi[i], mrk_i->p_z[i]);

        /*Momentum p along the magnetic field B. p \dot B / |B|*/
        real ppar_f =
            (
                mrk_f->p_r[i] * mrk_f->br[i] + 
                mrk_f->p_phi[i] * mrk_f->bphi[i] +
                mrk_f->p_z[i] * mrk_f->bz[i]
            ) / bnorm_f;
        /*Momentum p along the magnetic field B at the start of the timestep. p \dot B / |B|*/
        real ppar_i =
            (
                mrk_i->p_r[i] * mrk_i->br[i] + 
                mrk_i->p_phi[i] * mrk_i->bphi[i] +
                mrk_i->p_z[i] * mrk_i->bz[i]
            ) / bnorm_i;

        
        /*Lorentz factor gamma*/
        real gamma_f = physlib_gamma_pnorm(mrk_f->mass, pnorm_f);

        /*Lorentz factor gamma at the start of the timestep*/
        real gamma_i = physlib_gamma_pnorm(mrk_i->mass, pnorm_i);

        /*Parallel velocity*/
        real vpar_f = ppar_f / (gamma_f * mrk_f->mass);

        /*Toroidal velocity*/
        real vtor_f = mrk_f->p_phi[i] / (gamma_f * mrk_f->mass);

        /*Pitch xi*/
        real xi_f = pnorm_f > 0 ? ppar_f / pnorm_f : 0; 

        /*Total kinetic energy from mass and gamma*/
        real ekin_f = physlib_Ekin_gamma(mrk_f->mass, gamma_f);

        /*Total kinetic energy from mass and gamma at the start of the timestep*/
        real ekin_i = physlib_Ekin_gamma(mrk_i->mass, gamma_i);

        /*Change in total kinetic energy*/
        real delta_e = ekin_f - ekin_i; 

        /*Length of the timestep*/
        real delta_t = (mrk_f->time[i] - mrk_i->time[i]);

        /*Charge e of the marker*/
        real charge_f = CONST_E * mrk_f->charge[i];

        /*Charge e of the marker at the start of the timestep*/
        real charge_i = CONST_E * mrk_i->charge[i];

        /*Psi at the start and end of the timestep*/
        real psi_f, psi_i;
        valid *= !Bfield_eval_psi(&psi_f, mrk_f->r[i], mrk_f->phi[i], mrk_f->z[i], mrk_f->time[i], bfield); 
        valid *= !Bfield_eval_psi(&psi_i, mrk_i->r[i], mrk_i->phi[i], mrk_i->z[i], mrk_i->time[i], bfield); 

        /*Change in psi during the timestep*/
        real delta_psi = psi_f - psi_i;

        /*Change in parallel momentum during the timestep*/
        real delta_ppar = ppar_f - ppar_i; 


        /**Change in momentum conjugate to phi*/
        real delta_P_phi = 
             phys_ptoroid_fo(charge_f, mrk_f->r[i], mrk_f->p_phi[i], psi_f)
            -phys_ptoroid_fo(charge_i, mrk_i->r[i], mrk_i->p_phi[i], psi_i);

        /*Weight of the marker at the end of the timestep w*/
        real marker_weight_w = mrk_f->weight[i];

        computed_moments[MOMENT_DENSITY] = marker_weight_w * delta_t; /**Particle density W_n = w*\delta t */
        computed_moments[MOMENT_ENERGYDENSITY] = ekin_f * marker_weight_w * delta_t; /**Energy density W_e = E*w*\delta t */
        computed_moments[MOMENT_ENERGYDENSITYPAR] = ekin_f * xi_f * xi_f * marker_weight_w * delta_t; /**Parallel energy density W_e_|| = E*xi^2*w*\delta t */
        computed_moments[MOMENT_JPAR] = charge_f * vpar_f * marker_weight_w * delta_t; /**Parallel current W_j_|| = e*v_||*w*\delta t */
        computed_moments[MOMENT_JTOR] = charge_f * vtor_f * marker_weight_w * delta_t; /**Toroidal current W_j = e*v*w*\delta t*/
        computed_moments[MOMENT_CPDEP] = -delta_e * marker_weight_w; /**Collisional power deposition to plasma W_P,c=-\delta E_d*w */
        computed_moments[MOMENT_TORTORQ] = -charge_f * delta_psi * marker_weight_w; /**Toroidal torque j x B, W_\tau,jxB=-e*\delta\psi_p*w */
        computed_moments[MOMENT_TORCTORQ] = -mrk_f->r[i]*(delta_ppar*mrk_f->bphi[i]/bnorm_f) * marker_weight_w; /**Toroidal collision torque W_p,c=-R*(\delta\p_||*b_\phi/|B|) */
        computed_moments[MOMENT_TORCTORQPPHI] = -delta_P_phi * marker_weight_w; /**Toroidal torque from changes in toroidal canonical momentum W_\tau,p_\phi = -e*\delta P_\phi*w*/

        real rho_f = mrk_f->rho[i];
        real rho_i = mrk_i->rho[i];

        valid = valid && (mrk_i->running[i]
            && rho_f >= axis->min && rho_f <= axis->max
            && rho_i >= axis->min && rho_i <= axis->max);

        if(valid){
            size_t bin_index_f = math_bin_index(rho_f, axis->n, axis->min, axis->max);
            size_t bin_index_i = math_bin_index(rho_i, axis->n, axis->min, axis->max);

            size_t start = bin_index_f < bin_index_i ? bin_index_f: bin_index_i;
            size_t end = bin_index_f > bin_index_i ? bin_index_f: bin_index_i;
            real bin_weight = 1 / (real)(end-start+1);
            for(size_t final_bin_index = start; final_bin_index <= end; final_bin_index++){
                for(size_t moment_index = 0; moment_index < MOMENT_NUM; moment_index++)
                {
                    GPU_ATOMIC
                    moment->values[final_bin_index*MOMENT_NUM + moment_index] += computed_moments[moment_index] * bin_weight;
                }
            }
        }
    }
}

void DiagMoment_update_gc(
    DiagMoment *moment, Bfield *bfield, MarkerGuidingCenter *mrk_f,
    MarkerGuidingCenter *mrk_i)
{
    GPU_PARALLEL_LOOP_ALL_LEVELS
    for (size_t i = 0; i < mrk_f->size; i++)
    {
        int valid = 1;
        HistAxis *axis = &(moment->axis);
        real computed_moments[MOMENT_NUM];

        /*Norm of the magnetic field |B|*/
        real bnorm_f = math_normc(mrk_f->br[i], mrk_f->bphi[i], mrk_f->bz[i]);

        /*Norm of the magnetic field |B|*/
        real bnorm_i = math_normc(mrk_i->br[i], mrk_i->bphi[i], mrk_i->bz[i]);

        /*Norm of the momentum p |p|*/
        real pnorm_f = physlib_gc_p(mrk_f->mass, mrk_f->mu[i], mrk_f->ppar[i], bnorm_f);

        /*Norm of the momentum p |p| at the start of the timestep*/
        real pnorm_i = physlib_gc_p(mrk_i->mass, mrk_i->mu[i], mrk_i->ppar[i], bnorm_i);

        /*Momentum p along the magnetic field B. p \dot B / |B|*/
        real ppar_f = mrk_f->ppar[i];

        /*Momentum p along the magnetic field B. p \dot B / |B|*/
        real ppar_i = mrk_i->ppar[i];

        /*Lorentz factor from mass*/
        real gamma_f = physlib_gamma_pnorm(mrk_f->mass, pnorm_f);

        /*Lorentz factor from mass at the start of the timestep*/
        real gamma_i = physlib_gamma_pnorm(mrk_i->mass, pnorm_i);
        
        /*Parallel velocity*/
        real vpar_f = ppar_f / (gamma_f * mrk_f->mass);

        /*Toroidal velocity*/
        real vtor_f = vpar_f * mrk_f->bphi[i] / bnorm_f;

        /*Pitch xi*/
        real xi_f = pnorm_f > 0 ? ppar_f / pnorm_f : 0; 

        /*Total kinetic energy from mass and gamma*/
        real ekin_f = physlib_Ekin_gamma(mrk_f->mass, gamma_f);

        /*Total kinetic energy from mass and gamma at the start of timestep*/
        real ekin_i = physlib_Ekin_gamma(mrk_i->mass, gamma_i);

        /*Change in total kinetic energy*/
        real delta_e = ekin_f - ekin_i;

        /*Length of the timestep*/
        real delta_t = (mrk_f->time[i] - mrk_i->time[i]);

        /*Psi at the start and end of the timestep*/
        real psi_f, psi_i;
        valid *= !Bfield_eval_psi(&psi_f, mrk_f->r[i], mrk_f->phi[i], mrk_f->z[i], mrk_f->time[i], bfield); 
        valid *= !Bfield_eval_psi(&psi_i, mrk_i->r[i], mrk_i->phi[i], mrk_i->z[i], mrk_i->time[i], bfield); 
        /*Change in psi during the timestep*/
        real delta_psi = psi_f - psi_i;

        /*Change in parallel momentum during the timestep*/
        real delta_ppar = ppar_f - ppar_i; 

        /*Charge of marker e*/
        real charge_f = CONST_E * mrk_f->charge[i];

        /*Charge of marker e at start of timestep*/
        real charge_i = CONST_E * mrk_i->charge[i];

        real delta_P_phi = 
             phys_ptoroid_gc(charge_f, mrk_f->r[i], ppar_f, psi_f, bnorm_f, mrk_f->bphi[i])
            -phys_ptoroid_gc(charge_i, mrk_i->r[i], ppar_i, psi_i, bnorm_i, mrk_i->bphi[i]);


        /*Weight of the marker at the end of the timestep w*/
        real marker_weight_w = mrk_f->weight[i];

        computed_moments[MOMENT_DENSITY] = marker_weight_w * delta_t; /**Particle density W_n = w*\delta t */
        computed_moments[MOMENT_ENERGYDENSITY] = ekin_f * marker_weight_w * delta_t; /**Energy density W_e = E*w*\delta t */
        computed_moments[MOMENT_ENERGYDENSITYPAR] = ekin_f * xi_f * xi_f * marker_weight_w * delta_t; /**Parallel energy density W_e_|| = E*xi^2*w*\delta t */
        computed_moments[MOMENT_JPAR] = charge_f * vpar_f * marker_weight_w * delta_t; /**Parallel current W_j_|| = e*v_||*w*\delta t */
        computed_moments[MOMENT_JTOR] = charge_f * vtor_f * marker_weight_w * delta_t; /**Toroidal current W_j = e*v*w*\delta t*/
        computed_moments[MOMENT_CPDEP] = -delta_e * marker_weight_w; /**Collisional power deposition to plasma W_P,c=-\delta E_d*w */
        computed_moments[MOMENT_TORTORQ] = -charge_f * delta_psi * marker_weight_w; /**Toroidal torque j x B, W_\tau,jxB=-e*\delta\psi_p*w */
        computed_moments[MOMENT_TORCTORQ] = -mrk_f->r[i]*(delta_ppar*mrk_f->bphi[i]/bnorm_f) * marker_weight_w; /**Toroidal collision torque W_p,c=-R*(\delta\p_||*b_\phi/|B|) */
        computed_moments[MOMENT_TORCTORQPPHI] = -delta_P_phi * marker_weight_w; /**Toroidal torque from changes in toroidal canonical momentum W_\tau,p_\phi = -e*\delta P_\phi*w*/

        real rho_f = mrk_f->rho[i];
        real rho_i = mrk_i->rho[i];

        valid = valid && (mrk_i->running[i]
            && rho_f >= axis->min && rho_f <= axis->max
            && rho_i >= axis->min && rho_i <= axis->max);


        if(valid){
            size_t bin_index_f = math_bin_index(rho_f, axis->n, axis->min, axis->max);
            size_t bin_index_i = math_bin_index(rho_i, axis->n, axis->min, axis->max);

            size_t start = bin_index_f < bin_index_i ? bin_index_f: bin_index_i;
            size_t end = bin_index_f > bin_index_i ? bin_index_f: bin_index_i;
            real bin_weight = 1 / (real)(end-start+1);
            for(size_t final_bin_index = start; final_bin_index <= end; final_bin_index++){
                for(size_t moment_index = 0; moment_index < MOMENT_NUM; moment_index++)
                {
                    GPU_ATOMIC
                    moment->values[final_bin_index*MOMENT_NUM + moment_index] += computed_moments[moment_index] * bin_weight;
                }
            }
        }
    }
}
