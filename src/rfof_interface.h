/**
 * @file rfof_interface.h
 * @brief Contains the functions to be called from the simulation loop when
 * using ICRH.
**/

#ifndef RFOF_INTERFACE_H
#define RFOF_INTERFACE_H

#include <stdlib.h>

/**
 * @brief A one-to-one copy of the same struct from ASCOT4. Passed to RFOF kick
 * routine which stores some of the results in these fields.
*/
typedef struct prt_rfof {
    double dmu;          /**< Change in magnetic moment due to ICRH kick.     */
    double dvpar;        /**< Change in parallel velocity component due to ICRH
                              kick.                                           */
    double de;           /**< Change in energy due to a single ICRH kick [J]. */
    double deCumulative; /**< Change in energy due to possibly several ICRH 
                              kicks during an orbit time step [J]             */
    double dpitch;       /**< Change in pitch due to ICRH kick                */
    double maxAcc;       /**< Maximum acceleration allowed by RFOF            */
    double RFdt;         /**< time step suggested by RFOF                     */
} prt_rfof;

/** @brief Struct containing void pointers to RFOF structs */
typedef struct {
    void* cptr_rfof_input_params; /**< Pointer to rfof_input_param struct on the
                                       fortran side. This is needed when
                                       allocating resonance memory.           */
    void* cptr_rfglobal;          /**< Wave field; same for all markers.      */
    int icrh_initialised;         /**< 0 if not initialised, 1 if initialised.
                                      Currently used only in gui when plotting*/
} rfof_data;




/********************************* FUNCTIONS **********************************/


/* INITIALISATION */

void rfof_interface_initev_excl_marker_stuff(char* xml_filename,
    int **xml_filename_len,void** cptr_rfglobal, void** cptr_rfof_input_params);

void rfof_interface_initialise_res_mem(void** cptr_mem, int* cptr_mem_shape_i,
    int* cptr_mem_shape_j, void** cptr_rfglobal, void** cptr_rfof_input_param);

void rfof_interface_initialise_diagnostics(void** cptr_RFglobal,
    void** cptr_diagno);

void rfof_interface_allocate_rfof_marker(void** rfof_marker_pointer);


/* STUFF TO DO BETWEEN KICKS */

void rfof_interface_set_marker_pointers(void** cptr_marker, int* id, 
    real* weight, real* R, real* phi, real* z, real* psi, real* charge,
    real* mass, real* Ekin, real* velocity, real* mu, real* pphicanonical,
    real* vpar, real* vperp, real* gyrof, real* vdriftRho, real* acc,
    int* isOrbitTimeAccelerated, int* is_already_allocated);


/* KICKS */

void rfof_interface_do_rfof_stuff_gc(particle_simd_gc* ascot_marker, real* hin,
    real* hout_rfof, rfof_data rfof_data, B_field_data* Bdata,
    void** rfof_marker_pointer_array, void** rfof_mem_pointer_array,
    void** rfof_diag_pointer_array, int* mem_shape_i, int* mem_shape_j);

// To be implemented shortly
/*
void rfof_interface_do_rfof_stuff_fo(particle_simd_fo* ascot_marker, real* hin,
    real* hout_rfof, rfof_data rfof_data, B_field_data* Bdata,
    void** rfof_marker_pointer_array, void** rfof_mem_pointer_array,
    void** rfof_diag_pointer_array, int* mem_shape_i, int* mem_shape_j);
*/


/* RESET RESONANCE MEMORY */

void rfof_interface_reset_icrh_mem(void** rfof_marker_pointer, int* mem_shape_i,
    int* mem_shape_j);


/* DEALLOCATION ROUTINES */

void rfof_interface_deallocate_rfof_input_param(void** cptr_rfof_input_param);

void rfof_interface_deallocate_rfglobal(void** cptr_rfglobal);

void rfof_interface_deallocate_res_mem(void** cptr_res_mem,
    int* cptr_mem_shape_i, int* cptr_mem_shape_j);

void rfof_interface_deallocate_diagnostics(void** cptr_diagno);

void rfof_interface_deallocate_marker(void** cptr_rfof_marker);


/* FOR VISUALISING ICRH WAVE FIELD AND RESONANCE */

void rfof_interface_get_rf_wave_local(real* R, real* z, real* rho_tor, real* theta,
    void** cptr_wi, real* e_plus_out, real* e_minus_out);

void rfof_interface_eval_resonance_function(void** cptr_marker, void** cptr_rfglobal, real* omega_res, int* nharm);
#endif
