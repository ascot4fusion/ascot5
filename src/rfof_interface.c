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
#include "particle.h"
#include "rfof_interface.h"

#ifdef RFOF


/* INITIALISATION */

void __ascot5_icrh_routines_MOD_call_initev_excl_marker_stuff(char* xml_filename,
    int **xml_filename_len,void** cptr_rfglobal, void** cptr_rfof_input_params);
void __ascot5_icrh_routines_MOD_call_initialise_res_mem(void** cptr_mem,
    int* cptr_mem_shape_i, int* cptr_mem_shape_j, void** cptr_rfglobal,
    void** cptr_rfof_input_param);
void __ascot5_icrh_routines_MOD_call_initialise_diagnostics(void** cptr_RFglobal,
    void** cptr_diagno);

/* NOTE: There is no separate routine for initialising the markers; it is done
   using the call_set_marker_pointers function with the argument
   is_already_addlocated=0                                                    */


/* STUFF TO DO BETWEEN KICKS */

void __ascot5_icrh_routines_MOD_call_set_marker_pointers(void** cptr_marker,
    int** id, real** weight, real** R, real** phi, real** z, real** psi,
    real** charge, real** mass, real** Ekin, real** velocity, real** mu,
    real** pphicanonical, real** vpar, real** vperp, real** gyrof,
    real** vdriftRho, real** acc, int* isOrbitTimeAccelerated,
    int* is_already_allocated);


/* KICK */

void __ascot5_icrh_routines_MOD_call_rf_kick(double*time, double*dtin, int* myMPIprocID,
    prt_rfof* rfof_data_pack, void** cptr_marker, void**cptr_mem,
    void** cptr_rfglobal, void** cptr_rfdiagno, int *err, int*mem_shape_i,
    int*mem_shape_j);


/* RESET RESONANCE MEMORY */

void __ascot5_icrh_routines_MOD_call_reset_res_mem(void** rfof_mem_pointer,
    int* mem_shape_i, int* mem_shape_j);


/* DEALLOCATIONS*/

void __ascot5_icrh_routines_MOD_call_deallocate_rfof_input_param(
    void** cptr_rfof_input_param);
void __ascot5_icrh_routines_MOD_call_deallocate_rfglobal(void** cptr_rfglobal);
void __ascot5_icrh_routines_MOD_call_deallocate_res_mem(void** cptr_res_mem,
    int* cptr_mem_shape_i, int* cptr_mem_shape_j);
void __ascot5_icrh_routines_MOD_call_deallocate_diagnostics(void** cptr_diagno);
void __ascot5_icrh_routines_MOD_deallocate_marker(void** cptr_rfof_marker);


/* FOR VISUALISING ICRH WAVE FIELD  */

void __ascot5_icrh_routines_MOD_get_rf_wave_local_v2(real* R, real* z,
    real* rho_tor, real* theta, void** cptr_wi, real* e_plus_out,
    real* e_minus_out);

void __ascot5_icrh_routines_MOD_eval_resonance_function(void** cptr_marker,
    void** cptr_rfglobal, real* omega_res, int* nharm);


/* FOR DEBUGGING (n.b. TOTALLY OBSOLETE, OBVIOUSLY) */

void __ascot5_icrh_routines_MOD_print_marker_stuff(void** marker_pointer);
#endif




/********************************* FUNCTIONS **********************************/


/* INITIALISATION */

/**
 * @brief Initialise everyting  excluding marker stuff. Reads the ICRH (RFOF)
 * inputs (xml, xsd, ASCII) and initialises the wave field (variable name:
 * RFglobal).
 *
 * @param xml_filename Name of the xml file (less than 124 char)
 * @param xml_filename_len Length of the xml_filename (excluding the '\0' at the
 * end)
 * @param cptr_rfglobal void pointer to the constructed wave field. Cannot be
 * used to access the wave field but acts as a reference.
 * @param cptr_rfof_input_params void pointer to an RFOF struct containing the
 * input parameters. Only relevant when constructing the resonance memorys later
 * on.
*/
const char* xml_filename = "rfof_codeparam.xml";
void rfof_interface_initev_excl_marker_stuff(rfof_data* rfof_data) {

#ifdef RFOF
    //const char* xml_filename = RFOF_CODEPARAM_XML;
    int xml_filename_len = strlen(xml_filename);
    int*xml_filename_len_ptr = &xml_filename_len;
    printf("xml_filename = %s\n", xml_filename);
    printf("len(eml_file_name) = %d\n", xml_filename_len);
    __ascot5_icrh_routines_MOD_call_initev_excl_marker_stuff(xml_filename,
        &xml_filename_len_ptr, &(rfof_data->cptr_rfglobal), &(rfof_data->cptr_rfof_input_params));
#endif
};

/**
 * @brief Initialize RFOF data on target.
 *
 * Offloading for the RFOF data is not yet implemented, so this function just
 * copies pointers from the offload data to the target data for now.
 *
 * @param rfof_data rfof data struct on target
 * @param rfof_offload_data rfof data on host
 */
void rfof_init(rfof_data* rfof, rfof_data* rfof_offload_data) {
    rfof->cptr_rfglobal = rfof_offload_data->cptr_rfglobal;
    rfof->cptr_rfof_input_params = rfof_offload_data->cptr_rfof_input_params;
}

/**
 * @brief Initialises resonance memory for rfof markers. To be called before the
 * time step loop. (Each RFOF marker has its ows resonance memory matrix which
 * in turn has elements corresponding to each wave and its mode.)
 * @param cptr_mem Handle to corresponding Fortran resonance memory matrix pointer
 * @param cptr_mem_shape_i Size of the first dimension of the resonance memory
 * matrix
 * @param cptr_mem_shape_j Size of the second dimension of the resonance memory
 * matrix
 * @param cptr_rfglobal Handle to the Fortan wave field struct
 * @param cptr_rfof_input_param Handle to the Fortran input parameter struct
 * */
void rfof_interface_initialise_res_mem(void** cptr_mem, int* cptr_mem_shape_i,
    int* cptr_mem_shape_j, void** cptr_rfglobal, void** cptr_rfof_input_param) {
#ifdef RFOF
    __ascot5_icrh_routines_MOD_call_initialise_res_mem(cptr_mem, cptr_mem_shape_i,
        cptr_mem_shape_j, cptr_rfglobal, cptr_rfof_input_param);
#endif
};

/**
 * @brief Initialises rfof diagnostics. These are not used but are given as
 * dummy inputs to the kick routine. To be called before the time step loop.
 * @param cptr_rfglobal Handle to the Fortan wave field struct
 * @param cptr_diagno Handle to the Fortran diagnostics struct
*/
void rfof_interface_initialise_diagnostics(void** cptr_rfglobal,
    void** cptr_diagno) {
#ifdef RFOF
    __ascot5_icrh_routines_MOD_call_initialise_diagnostics(cptr_rfglobal, cptr_diagno);
#endif
};

/**
 * @brief Allocates memory for an rfof marker on the fortran side. To be called
 *  before the time step loop. The memory should be deallocated at the end by
 * calling the corresponding deallocation routine at the end of this file.
 * @param rfof_marker_pointer Handle to rfof marker struct.
*/
void rfof_interface_allocate_rfof_marker(void** rfof_marker_pointer) {
#ifdef RFOF

    /* These are actually dummies and could be whatever you want              */
    real dummy_real_rfof = 41.99;
    real* dummy_real_rfof_ptr = &dummy_real_rfof;
    int dummy_int_rfof = 42;
    int *dummy_int_rfof_ptr = &dummy_int_rfof;

    /** @brief Needs to be zero. This is not stored in the marker; this input
     * parameter is only used when calling the call_set_marker_pointers to
     * determine whether to allocate new memory or use an existing marker     */
    int is_already_allocated = 0;

    /* Allocates memory and sets some dummy values to the marker fields       */
    __ascot5_icrh_routines_MOD_call_set_marker_pointers(
                rfof_marker_pointer,
                &dummy_int_rfof_ptr,
                &dummy_real_rfof_ptr,
                &dummy_real_rfof_ptr,
                &dummy_real_rfof_ptr,
                &dummy_real_rfof_ptr,
                &dummy_real_rfof_ptr,
                &dummy_real_rfof_ptr,
                &dummy_real_rfof_ptr,
                &dummy_real_rfof_ptr,
                &dummy_real_rfof_ptr,
                &dummy_real_rfof_ptr,
                &dummy_real_rfof_ptr,
                &dummy_real_rfof_ptr,
                &dummy_real_rfof_ptr,
                &dummy_real_rfof_ptr,
                &dummy_real_rfof_ptr,
                &dummy_real_rfof_ptr,
                &dummy_int_rfof, /* isOrbitTimeAccelerated: Currently a dummy
                                    (The fortran routine does nothing with this
                                    and sets orbit_time_accelerated to false
                                    regardless)                               */
                &is_already_allocated);
#endif
};


/* STUFF TO DO BETWEEN KICKS */

/**
 * @brief Sets the RFOF marker's fields (all Fortran pointers). If the RFOF
 * marker has not been allocated previously, this routine can be used to reate a
 * marker for the RFOF routines on the fortran side. This allocation of a new
 * RFOF marker should be done using the function
 * rfof_interface_allocate_rfof_marker defined in this file.
 * @param cptr_marker Handle to the RFOF marker struct in Fortran
 * @param id Id of the marker (unclear whether truly needed)
 * @param weight Marker weight
 * @param R Major radius co-ordinate
 * @param phi Toroidal angle co-ordinate
 * @param z Vertical co-ordinate
 * @param psi Poloidal flux function (following the ITM and ITER conventions;
 * COCOS 13/11; not divided by \f$ 2\pi \f$)
 * @param charge Charge of the particles represented by the marker.
 * @param mass Mass of the particles represented by the marker
 * @param Ekin Kinetic energy
 * @param Velocity Speed of the particle. Obviously.
 * @param mu Magnetic moment
 * @param pphicanonical Canonical momentum conjugate to phi
 * @param vpar Component of velocity parallel to B
 * @param vperp Magnitude of velocity perpendicular to B
 * @param gyrof Gyrofrequency
 * @param vdriftRho The component of the drift velocity w.r.t. the radial
 * direction (unclear if truly needed). This can used at least when evaluating
 * the Doppler shift in the resonance condition.
 * @param acc Time acceleration factor (not in use anymore, should be 1; NOT 0)
 * @param isOrbitTimeAcclerated False as implied above in acc. Logical in
 * fortran; represented by int (4 bytes) in C. Currently the Fortran routine
 * does nothing with this parameter but instead sets it to false for all
 * markers.
 * @param is_already_allocated If true, the routine tries to fetch the existing
 * particle corresponding to cptr_marker. If false, the routine allocates
 * memory for an RFOF marker and returns the c_location of that fortran struct
 * in cptr_marker.
*/
void rfof_interface_set_marker_pointers(void** cptr_marker, int* id,
    real* weight, real* R, real* phi, real* z, real* psi, real* charge,
    real* mass, real* Ekin, real* velocity, real* mu, real* pphicanonical,
    real* vpar, real* vperp, real* gyrof, real* vdriftRho, real* acc,
    int* isOrbitTimeAccelerated, int* is_already_allocated) {
#ifdef RFOF
    __ascot5_icrh_routines_MOD_call_set_marker_pointers(cptr_marker, &id,
        &weight, &R, &phi, &z, &psi, &charge, &mass, &Ekin, &velocity, &mu,
        &pphicanonical, &vpar, &vperp, &gyrof, &vdriftRho, &acc, isOrbitTimeAccelerated, is_already_allocated);
#endif
};


/**
 * @brief Function to be called in the main simulation loop during each step
 * when following the guiding centre.
 *
 * 1. Updates the fields of the rfof_marker based on the given input
 * ascot_marker.
 * 2. Calls the "kick" function, which
 *      a) Checks resonance condition and
 *      b) if in resonance, kicks marker and updates velocity of the rfof marker
 *          and consequently also ascot marker (as only pointers are passed when
 *          creating the rfof marker).
 * If the proposed time step, hin, was too large, hout_rfof will be negative
 * indicating that the whole time step has to be redone with a smaller time
 * step.
 * @param ascot_marker Ascot marker (no points for those who guessed it)
 * @param hin Time step proposed by the simulation loop
 * @param hout_rfof If negative, ICRH kick failed because hin was too large. In
 * that case the absolute value of hout_rfof should be used when redoing the
 * simulation loop. For successfull steps, this is the estimate to the next
 * resonance.
 * @param rfof_data A "package" of two ICRH related void pointers to Fortran
 * routines. Out of the two, only the pointer to the global wave field is used.
 * This variable is not to be confused with the local variable rfof_data_pack
 * defined inside this function.
 * @param Bdata The magnetic field. Needed for evaluating psi.
 * @param rfof_marker_pointer_array Contains void pointers which are handles to
 * rfof markers on the Fortran side.
 * @param rfof_mem_pointer_array Contains void pointers which are handles to
 * rfof markers' resonance memory matrices.
 * @param rfof_diag_pointer_array Contains void pointers which are handles to
 * rfof diagnostics structs (not in actual use but must be included).
 * @param mem_shape_i Array of rfof resonance memory matrix's first dimensions.
 * @param mem_shape_j Array of rfof resonance memory matrix's second dimensions.
*/
void rfof_interface_do_rfof_stuff_gc(particle_simd_gc* ascot_marker, real* hin,
    real* hout_rfof, rfof_data rfof_data, B_field_data* Bdata,
    void** rfof_marker_pointer_array, void** rfof_mem_pointer_array,
    void** rfof_diag_pointer_array, int* mem_shape_i, int* mem_shape_j) {
#ifdef RFOF
    for (int i=0; i<NSIMD;i++) {
        if (ascot_marker->id[i] > 0 && ascot_marker->running[i]) {
            /*Some of the fields needed by the RFOF_marker struct are not
            present in the ascot_marker struct and thus need to be evaluated
            first. */

            /** @brief Poloidal flux function (following the ITM and ITER
             * conventions; COCOS 13/11; not divided by \f$ 2\pi \f$)         */
            real psi;

            B_field_eval_psi(&psi, ascot_marker->r[i], ascot_marker->phi[i],
                ascot_marker->z[i], ascot_marker->time[i], Bdata);

            /** @brief Norm of B field                                        */
            real B = sqrt(pow(ascot_marker->B_r[i],2) +
                pow(ascot_marker->B_phi[i],2) + pow(ascot_marker->B_z[i],2));

            /** @brief Lorentz factor.*/
            real gamma = physlib_gamma_ppar(ascot_marker->mass[i],
                ascot_marker->mu[i], ascot_marker->ppar[i], B);

            /** @brief Energy [J]*/
            real Ekin = physlib_Ekin_gamma(ascot_marker->mass[i], gamma);

            /** @brief Speed */
            real speed = physlib_vnorm_gamma(gamma);

            /** @brief Canonical momentum conjugate to phi (toroidal momenutum).
             * Should be in SI units                                          */
            real p_phi = phys_ptoroid_gc(ascot_marker->charge[i],
                ascot_marker->r[i], ascot_marker->ppar[i], psi, B,
                ascot_marker->B_phi[i]);

            /** @brief Momentum norm */
            real p = physlib_pnorm_vnorm(ascot_marker->mass[i], speed);

            /** @brief pitch */
            real xi = physlib_gc_xi(ascot_marker->mass[i], ascot_marker->mu[i],
                ascot_marker->ppar[i], B);

            /** @brief Parallel momentum */
            real v_par = speed*xi;

            /** @brief Perpendicular speed */
            real v_perp = phys_vperp_gc(speed, v_par);

            /** @brief Gyrofrequency */
            real gyrof = phys_gyrofreq_ppar(ascot_marker->mass[i],
                ascot_marker->charge[i], ascot_marker->mu[i],
                ascot_marker->ppar[i], B);

            /** @brief Velocity drift in the direction of rho. Used when
             * evaluating the Doppler shift term in the resonance condition.
             * APPARENTLY NOT USED?                                           */
            real vdriftRho = 0;

            /** @brief Time acceleration factor. EVEN IF TIME ACCELERATION IS
             * OFF, THIS IS NOT A DUMMY BUT MUST BE SET TO 1.0                */
            real acc = 1.0;

            /** @brief Time accelecation on/off. NOTE: Fortran logical
             * corresponds to C int (4 bytes). NOTE2: Obsolete at the moment, as
             *  the Fortran routine currently forces this to be false (0)
             * regardless.                                                    */
            int isOrbitTimeAccelerated = 0;

            /** @brief Must be 1 (true) at this point as the marker is already
             * allocated.*/
            int is_already_allocated = 1;


            /* Now when we have all the needed parameters we call the rfof
            routine which sets the pointers of the rfof_marker (on the Fortran
            side. The situation is somewhat comparable to ordering pizza/kebab
            from Kontula: it's convenient that the place exists but you'd rather
            not go there yourself. We're in luck, though, as we happen to have a
            courier.) */

            int dummy_Id = 3141592;
            int* dummy_Id_ptr = &dummy_Id;

            /* These fields needed by RFOF marker we can pass straight away   */

            /** @brief The weight is only used for RFOF's own diagnostics which
             * are not of interest to us. Thusly, the weight given to the RFOF
             * marker could be arbitrary; from a certain point of view it's a
             * dummy. In other words, these are not the weights you are looking
             * for.*/
            real* weight_ptr = &(ascot_marker->weight[i]);
            real* r_ptr = &(ascot_marker->r[i]);
            real* phi_ptr = &(ascot_marker->phi[i]);
            real* z_ptr = &(ascot_marker->z[i]);
            real* charge_ptr = &(ascot_marker->charge[i]);
            real* mass_ptr = &(ascot_marker->mass[i]);
            real* mu_ptr = &(ascot_marker->mu[i]);

            /* These fields needed to be evaluated inside this function. They
            are not independent variables.                                    */
            real* Ekin_ptr = &Ekin;
            real* psi_ptr = &psi;
            real* speed_ptr = &speed;
            real* p_phi_ptr = &p_phi;
            real* v_par_ptr = &v_par;
            real* v_perp_ptr = &v_perp;
            real* gyrof_ptr = &gyrof;
            real* vdriftRho_ptr = &vdriftRho;
            real* acc_ptr = &acc;

            //For debugging store the old values
            real Ekin_old = Ekin;
            real p_phi_old = p_phi;
            real v_par_old = v_par;
            real v_perp_old = v_perp;
            real speed_old = speed;
            real ppar_old = ascot_marker->ppar[i];

            /* Update the fields of RFOF marker */
            __ascot5_icrh_routines_MOD_call_set_marker_pointers(
                &rfof_marker_pointer_array[i], /* Note that the pointer to the
                                                  pointer is passed.          */
                &(dummy_Id_ptr),
                &(weight_ptr),                 /* Number of real particles
                                                  represented by the marker.
                                                  Acts as a dummy.            */
                &(r_ptr),
                &(phi_ptr),
                &(z_ptr),
                &(psi_ptr),
                &(charge_ptr),
                &(mass_ptr),
                &(Ekin_ptr),
                &(speed_ptr),
                &(mu_ptr),
                &(p_phi_ptr),            /* pphicanonical                     */
                &(v_par_ptr),
                &(v_perp_ptr),
                &(gyrof_ptr),
                &(vdriftRho_ptr),        /* vdriftRho; EVALUATE IF ACTUALLY
                                            NEEDED                            */
                &(acc_ptr),

                &isOrbitTimeAccelerated, /**< Int in C, logical in Fortran    */
                &is_already_allocated);  /**< Int in C, logical in Fortran    */


            /** @brief Used for storing the "results" of calling RF kick, only
             * the RFdt field is utilised as it returns the time step
             * recommended by RFOF                                            */
            prt_rfof rfof_data_pack = {
                .dmu = 0.0,
                .dvpar = 0.0,
                .de = 0.0,
                .deCumulative = 0.0,
                .dpitch = 0.0,
                .maxAcc = 0.0,
                .RFdt = 0.0,
            };

            int err = 0;       /**< "empty" input to kick                     */

            /** @brief Number used to identify MPI nods during parallel
             * execution. This process is not to be paralleliz (?)            */
            int mpiprocid = 0;

            /* Ready to kick some ash */
            __ascot5_icrh_routines_MOD_call_rf_kick(&(ascot_marker->time[i]),
                &(hin[i]), &mpiprocid, &rfof_data_pack,
                &(rfof_marker_pointer_array[i]), &(rfof_mem_pointer_array[i]),
                &(rfof_data.cptr_rfglobal), &(rfof_diag_pointer_array[i]), &err,
                &(mem_shape_i[i]), &(mem_shape_j[i]));

            /* Some of the rfof marker's pointers where pointing to the fields
            of the ascot marker but some are pointing to the local variables
            inside this function (e.g. Ekin). We now need to update the rest of
            the ASCOT marker fields (= ppar) accordingly.*/


            //OLD VERSION

            //int sign_v_par_old = (v_par_old > 0) - (v_par_old < 0);
            //if(v_par_old*v_par < 0){
            //    /* The parallel velocity has flipped during the icrh kick, give opposite sign to v_par_old */
            //    ascot_marker->ppar[i] = -sign_v_par_old*phys_ppar_Ekin(ascot_marker->mass[i], Ekin, ascot_marker->mu[i], B);
            //}else{
            //    /* The parallel velocity has not flipped, return ppar with same sign as v_par_old */
            //    ascot_marker->ppar[i] = sign_v_par_old*phys_ppar_Ekin(ascot_marker->mass[i], Ekin, ascot_marker->mu[i], B);
            //}


            /* Current version for updating ppar of ASCOT marker based on the
            RFOF kick. Below there are several other methods listed for doing
            this. This method was found to conserve the ppar in the case that
            the marker would not receive a kick where as some of the methods
            below would slightly alter the markers ppar even when no kick is
            applied. This implementation corresponds to number 5 in the list
            below. In fear that this effect could accumulate, this method of
            evaluating ppar was deemed best. */
            ascot_marker->ppar[i] = phys_ppar_pphi(B, ascot_marker->r[i],
                                                   ascot_marker->B_phi[i],
                                                   p_phi,
                                                   ascot_marker->charge[i],
                                                   psi);


            // The other proposed methods for evaluating ppar:
            /*
            printf("\nCompute ppar using different input fields of RFOF marker to check consistency\n");

            // 0. ppar_old just for comparison
            printf("For reference, ppar_old = %e\n", ppar_old);

            // 1. Updata pitch and use pitch and v to get vpar. Then use ppar =
            //m*vpar
            xi = xi + rfof_data_pack.dpitch;
            real vpar1 = xi*speed;
            printf("1. ppar from v and pitch + dpitch = %e\n", ascot_marker->mass[i]*vpar1);

            // 2. Get ppar from v_par and m
            printf("2. ppar = m*vpar = %e\n", ascot_marker->mass[i]*v_par);

            // 3. Get vpar from speed and vperp
            real vpar2 = sqrt(speed*speed - v_perp*v_perp);
            printf("3. abs(ppar) from speed and vperp = %e\n", vpar2*ascot_marker->mass[i]);

            // 4. From Ekin get speed, from mu get vperp, then get vpar and ppar
            real speed2 = sqrt(2.0*Ekin/ascot_marker->mass[i]);
            real vperp2 = sqrt(2.0*ascot_marker->mu[i]*B/ascot_marker->mass[i]);
            real vpar3 = sqrt(speed2*speed2 - vperp2*vperp2);
            printf("4. abs(ppar) from Ekin and mu = %e\n",ascot_marker->mass[i]*vpar3);

            // 5. Get ppar from p_phi
            real ppar5 = B/(ascot_marker->r[i]*ascot_marker->B_phi[i])*(p_phi - ascot_marker->charge[i]*psi);
            printf("5. ppar from p_phi = %e\n", ppar5);

            // 6. Get ppar from gyrof and mu
            real ppar6 = ascot_marker->mass[i]*CONST_C*sqrt((ascot_marker->charge[i]*B/(ascot_marker->mass[i]*gyrof))*(ascot_marker->charge[i]*B/(ascot_marker->mass[i]*gyrof)) - 2*ascot_marker->mu[i]*B/(ascot_marker->mass[i]*CONST_C2) - 1.0);
            printf("6. abs(ppar) from gyrof and mu = %e\n", ppar6);

            // 7. Get ppar from vpar_old and dvpar in rfof_data_pack
            printf("7. ppar from vpar_old and dvpar in rfof_data_pack = %e\n", ascot_marker->mass[i]*(v_par_old + rfof_data_pack.dvpar));
            */


            /* Check if dt was sufficiently small and assign hout_rfof
            accordingly. */
            if (err == 7) {
                /*Interaction failed, particle overshot the resonance. Make
                hout_rfof negative. The absolute value of rfof_hout in this case
                is the estimate of what hin should have been so that the
                resonance would not have been missed (obviously less than hin).
                */
                hout_rfof[i] = -rfof_data_pack.RFdt;
            } else {
                /* Interaction was successful. dt returned by rfof is now the
                estimate to the next resonance. */
                if(rfof_data_pack.RFdt == 0) {
                    /* This if is only for debugginf purposes when it is
                    possible that the kick is not called. */
                    printf("rfof_datapack.RFdt = 0, "
                           "setting rf return dt to 1\n");
                    hout_rfof[i] = 1;
                } else {
                    /* This is where you normally go when you call kick and
                    there is no error. */
                    hout_rfof[i] = rfof_data_pack.RFdt;
                }
            }

        };
   };
#endif
};


/**
 * @brief Function to be called in the main simulation loop during each step for
 * full orbit.
 *
 * 1. Updates the fields of the rfof_marker based on the given input
 * ascot_marker.
 * 2. Calls the "kick" function, which
 *      a) Checks resonance condition and
 *      b) if in resonance, kicks marker and updates velocity of the rfof marker
 *          and consequently also ascot marker (as only pointers are passed when
 *          creating the rfof marker).
 * If the proposed time step, hin, was too large, hout_rfof will be negative
 * indicating that the whole time step has to be redone with a smaller time
 * step.
  * @param ascot_marker Ascot marker (no points for those who guessed it)
 * @param hin Time step proposed by the simulation loop
 * @param hout_rfof If negative, ICRH kick failed because hin was too large. In
 * that case the absolute value of hout_rfof should be used when redoing the
 * simulation loop. For successfull steps, this is the estimate to the next
 * resonance.
 * @param rfof_data A "package" of two ICRH related void pointers to Fortran
 * routines. Out of the two, only the pointer to the global wave field is used.
 * This variable is not to be confused with the local variable rfof_data_pack
 * defined inside this function.
 * @param Bdata The magnetic field. Needed for evaluating psi.
 * @param rfof_marker_pointer_array Contains void pointers which are handles to
 * rfof markers on the Fortran side.
 * @param rfof_mem_pointer_array Contains void pointers which are handles to
 * rfof markers' resonance memory matrices.
 * @param rfof_diag_pointer_array Contains void pointers which are handles to
 * rfof diagnostics structs (not in actual use but must be included).
 * @param mem_shape_i Array of rfof resonance memory matrix's first dimensions.
 * @param mem_shape_j Array of rfof resonance memory matrix's second dimensions.
*/
//void rfof_interface_do_rfof_stuff_fo(particle_simd_fo* ascot_marker, //real* hin, real* hout_rfof, rfof_data rfof_data) {  /* B_field_data* //Bdata ? */
//#ifdef RFOF
//    for (int i=0; i<NSIMD;i++) {
//        if (ascot_marker->id[i] > 0 && ascot_marker->running[i]) {
//            /*Some of the fields needed by the RFOF_marker struct are not
//            present in the ascot_marker struct and thus these need to be
//            evaluated first. */
//
//            real B = sqrt(pow(ascot_marker->B_r[i],2) + pow//(ascot_marker->B_phi[i],2) + pow(ascot_marker->B_z[i],2)); ///**< Norm of B field  */
//            real B_unit_vec[3] = {(ascot_marker->B_r[i])/B, //(ascot_marker->B_phi[i])/B, (ascot_marker->B_z[i])/B}; /**< //Unit vector parallel to B field.*/
//            real v_r = ascot_marker->p_r[i]/ascot_marker->mass[i]; /**< r
//            component of velocity. */
//            real v_phi = ascot_marker->p_phi[i]/ascot_marker->mass[i]; ///**< phi
//            component of velocity. */
//            real v_z = ascot_marker->p_z[i]/ascot_marker->mass[i]; /**< z
//            component of velocity. */
//            real v_vec[3] = {v_r, v_phi,v_z}; /**< Velocity vector. */
//            real v_par = 0; /**< Parallel component of velocity. */
//            real speed_squared = 0; /**< Speed -- who would have thought. //*/
//            for (int xi = 0; xi < 3; xi++) {
//                v_par += v_vec[xi] * B_unit_vec[xi];
//                speed_squared += pow(v_vec[xi],2);
//            }
//            real speed = sqrt(speed_squared); /**< Speed. */
//            real v_perp_squared = speed_squared - pow(v_par,2); /**<
//            Perpendicular speed squared.*/
//            real v_perp = sqrt(v_perp_squared); /**< Perpendicular speed. //*/
//            real Ekin = 1/2*speed_squared*ascot_marker->mass[i]; /**< //Kinetic
//            energy */
//            real mu; /**< Magnetic moment. */
//            if(v_par > 0) {
//                mu = ascot_marker->mass[i]*v_perp_squared/(2*B);
//            } else{
//                /* If vpar opposite to B field, add minus in front. */
//                mu = -ascot_marker->mass[i]*v_perp_squared/(2*B);
//            }
//            real gyrof = ascot_marker->charge[i]*B/ascot_marker->mass[i];
//
//            real vdriftRho = 0; /**< Velocity drift in the direction of //rho.
//            NOT NEEDED? */
//            real acc = 0; /**< Time acceleration factor. NOT NEEDED?*/
//            int isOrbitTimeAccelerated = 0; /**< Time accelecation on///off. NOTE:
//            fortran logical corresponds to c int (4 bytes).*/
//
//            /*Now when we have all the needed parameters we call the rfof
//            routine which sets the pointers of the rfof_marker (on the //Fortran
//            side. The situation is somewhat comparable to ordering pizza///kebab
//            from Kontula: it's convenient that the place exists but you'd //rather
//            not go there yourself. We're in luck, though, as we happen to //have
//            a courier.) */
//
//            __ascot5_icrh_routines_MOD_call_set_marker_pointers(
//                ascot_marker->id[i],
//                ascot_marker->weight[i],  /**<Number of real particles
//                represented
//                by the marker */
//                ascot_marker->r[i],
//                ascot_marker->phi[i],
//                ascot_marker->z[i],
//                ascot_marker->,    /* psi here */
//                ascot_marker->charge[i],
//                ascot_marker->mass[i],
//                Ekin,
//                speed,
//                mu,
//                ascot_marker->p_phi[i],   /* pphicanonical */
//                v_par,
//                v_perp,
//                gyrof,
//                vdriftRho, /* vdriftRho; EVALUATE IF ACTUALLY NEEDED */
//                acc,
//                isOrbitTimeAccelerated);
//
//            };
//    /* call rfof kick */
//
//
//    /* check if dt sufficiently small*/
//        /* update particle_simd_fo (p_r. p_phi, p_z) using Energy or //velocities
//        (check if consistent) */
//
//   };
//#endif
//};


/* RESET RESONANCE MEMORY */

/**
 * @brief Resets resonance memory of ICRH (RFOF) markers. Should be done when
 * the marker dies and a new one is born. Note that a marker with a newly
 * allocated or reseted memory cannot receive ICRH kicks during the first two
 * time steps as its resonance memory must have at least two data points stored
 * for it to be kicked.
 * @param rfof_mem_pointer Handle to rfof resonance memory matrix.
 * @param mem_shape_i RFOF resonance memory matrix's first dimension.
 * @param mem_shape_j RFOF resonance memory matrix's second dimension.
 */
void rfof_interface_reset_icrh_mem(void** rfof_mem_pointer, int* mem_shape_i,
    int* mem_shape_j) {
#ifdef RFOF
    __ascot5_icrh_routines_MOD_call_reset_res_mem(rfof_mem_pointer, mem_shape_i,
        mem_shape_j);
#endif
};


/* DEALLOCATION ROUTINES */

/**
 * @brief Deallocates the rfof_input_param struct on the fortran side. There
 * exists only one copy of this struct and therefore it is to be deallocated in
 * the simulate.c after the loop is completed.
 * @param cptr_rfof_input_param Handle to rfof input param struct on the Fortran
 * side.
 */
void rfof_interface_deallocate_rfof_input_param(void** cptr_rfof_input_param) {
#ifdef RFOF
    __ascot5_icrh_routines_MOD_call_deallocate_rfof_input_param(cptr_rfof_input_param);
#endif
};

/**
 * @brief Deallocates the rfglobal struct (wave field) on the fortran side.
 * There exists only one copy of this struct and therefore it is to be
 * deallocated in the simulate.c after the loop is completed -- much like the
 * input_param struct.
 * @param cptr_rfglobal Handle to rfof wave field on the Fortran side.
 */
void rfof_interface_deallocate_rfglobal(void** cptr_rfglobal) {
#ifdef RFOF
    __ascot5_icrh_routines_MOD_call_deallocate_rfglobal(cptr_rfglobal);
#endif
};

/**
 * @brief Deallocates the resonance memory matrix of a particle. To be done when
 * the simulation has finished. If a marker dies but there are still new ones in
 * the queue, the resonance memory matrix should only be resetted
 * (see:rfof_interface_reset_icrh_mem), not deallocated.
 * @param rfof_mem_pointer Handle to rfof resonance memory.
 * @param mem_shape_i RFOF resonance memory matrix's first dimension.
 * @param mem_shape_j RFOF resonance memory matrix's second dimension.
*/
void rfof_interface_deallocate_res_mem(void** cptr_res_mem,
    int* cptr_mem_shape_i, int* cptr_mem_shape_j) {
#ifdef RFOF
    __ascot5_icrh_routines_MOD_call_deallocate_res_mem(cptr_res_mem, cptr_mem_shape_i,
        cptr_mem_shape_j);
#endif
};

/**
 * @brief Deallocates the (dummy) dianostics of rfof markers.
 * @param cptr_diagno Handle to rfof diagnostics struct.
 */
void rfof_interface_deallocate_diagnostics(void** cptr_diagno) {
#ifdef RFOF
    __ascot5_icrh_routines_MOD_call_deallocate_diagnostics(cptr_diagno);
#endif
};

/**
 * @brief Deallocates the rfof marker.
 * @param cptr_rfof_marker Handle to rfof marker.
 */
void rfof_interface_deallocate_marker(void** cptr_rfof_marker) {
#ifdef RFOF
    __ascot5_icrh_routines_MOD_deallocate_marker(cptr_rfof_marker);
#endif
};


/* FOR VISUALISING ICRH WAVE FIELD  */

/**
 * @brief Return the local E+ and E- values of the ICRH field, given the
 * coordinates.
 * @param R Major radius
 * @param z Vertical co-ordinate
 * @param rho_tor
 * @param theta
 * @param cptr_rfglobal Void pointer to the RFglobal global wave field in
 * fortran
 * @param e_plus E+ component of the local wave field
 * @param e_minus E- component of the local wave field
 */
void rfof_interface_get_rf_wave_local(real* R, real* z, real* rho_tor,
    real* theta, void** cptr_wi, real* e_plus_out, real* e_minus_out) {
#ifdef RFOF
    __ascot5_icrh_routines_MOD_get_rf_wave_local_v2(R, z, rho_tor, theta ,
        cptr_wi, e_plus_out,  e_minus_out);
#endif
};

/**
 * @brief Function for evaluating the value of resonance function (0 = resonance)
 * @param cptr_marker void pointer to the rfof_marker
 * @param cptr_rfglobal void pointer to the rfof wave field
 * @param omega_res value of the resonance function
 * @param nharm harmonic index
 */
void rfof_interface_eval_resonance_function(void** cptr_marker,
    void** cptr_rfglobal, real* omega_res, int* nharm){
#ifdef RFOF
    __ascot5_icrh_routines_MOD_eval_resonance_function(cptr_marker,
        cptr_rfglobal, omega_res, nharm);
#endif
};