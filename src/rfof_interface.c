/**
 * @file rfof_interface.c
 * @brief Contains the function to be called from the simulation loop when using
 * ICRH.
**/

#include <stdio.h>
#include <string.h>
#include <math.h>
#include "physlib.h"
#include "particle.h"
#include "rfof_interface.h"
//#define RFOF 1

#ifdef RFOF


/* INITIALISATION */

void __valipalikka_MOD_call_initev_excl_marker_stuff(char* xml_filename,
    int **xml_filename_len,void** cptr_rfglobal, void** cptr_rfof_input_params);
void __valipalikka_MOD_call_initialise_res_mem(void** cptr_mem,
    int* cptr_mem_shape_i, int* cptr_mem_shape_j, void** cptr_rfglobal,
    void** cptr_rfof_input_param);
void __valipalikka_MOD_call_initialise_diagnostics(void** cptr_RFglobal,
    void** cptr_diagno);


/* STUFF TO DO BETWEEN KICKS */

void __valipalikka_MOD_call_set_marker_pointers(void** cptr_marker, int* id,
    real* weight, real* R, real* phi, real* z, real* psi, real* charge,
    real* mass, real* Ekin, real* velocity, real* mu, real* pphicanonical,
    real* vpar, real* vperp, real* gyrof, real* vdriftRho, real* acc,
    int* isOrbitTimeAccelerated, int* is_already_allocated);


/* KICK */

void __valipalikka_MOD_call_rf_kick(double*time, double*dtin, int* myMPIprocID,
    prt_rfof* rfof_data_pack, void** cptr_marker, void**cptr_mem,
    void** cptr_rfglobal, void** cptr_rfdiagno, int *err, int*mem_shape_i,
    int*mem_shape_j);


/* RESET RESONANCE MEMORY */

void __valipalikka_MOD_call_reset_res_mem(void** rfof_mem_pointer,
    int* mem_shape_i, int* mem_shape_j);


/* DEALLOCATIONS*/

void __valipalikka_MOD_call_deallocate_rfof_input_param(
    void** cptr_rfof_input_param);
void __valipalikka_MOD_call_deallocate_rfglobal(void** cptr_rfglobal);
void __valipalikka_MOD_call_deallocate_res_mem(void** cptr_res_mem,
    int* cptr_mem_shape_i, int* cptr_mem_shape_j);
void __valipalikka_MOD_call_deallocate_diagnostics(void** cptr_diagno);
void __valipalikka_MOD_deallocate_marker(void** cptr_rfof_marker);
#endif


/* TODO: Tänne voisi lisätä niitä for looppeja tuolta simulate_gc_adaptiven 
alusta, jossa allokoidaan markerit, muisti ja diagnostiikka niin sitten siellä 
näyttää siistimmältä.
*/

/* TODO: samaan tapaan voisi lisätä jonkun funktion joka käy for loopissa 
deallokoimassa muistit, markerit ja diagnostiikan. aaltokentän deallokointiin 
riittää että kutsutaan valipalikassa olevaa rutiinia joka kutsuu sisso_wrapperin
rutiinia, sillä aaltokentän deallokointi tehdään vain kerran. 
*/


/* FUNCTIONS */


/* INITIALISATION */

/**
 * @brief Initialise everyting  excluding marker stuff. Reads the ICRH (RFOF) 
 * inputs (xml, xsd, ASCII) and initialises the wave field.
 * @param xml_filename Name of the xml file (less than 124 char)
 * @param xml_filename_len Length of the xml_filename (excluding the '\0')
 * @param cptr_rfglobal void pointer to the constructed wave field. Cannot be 
 * used to access the wave field but acts as a reference.
 * @param cptr_rfof_input_param void pointer to an RFOF struct containing the 
 * input parameters. Only relevant when constructing the resonance memorys later
 * on.
*/
void rfof_interface_initev_excl_marker_stuff(char* xml_filename,
    int **xml_filename_len,void** cptr_rfglobal,void** cptr_rfof_input_params) {
#ifdef RFOF
    __valipalikka_MOD_call_initev_excl_marker_stuff(xml_filename,
        xml_filename_len,cptr_rfglobal, cptr_rfof_input_params);
#endif
};

/**
 * @brief Initialises resonance memory for rfof markers. To be called before the
 * time step loop.
 * @param cptr_mem Handle to corresponding Fortran resonance memory pointer
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
    __valipalikka_MOD_call_initialise_res_mem(cptr_mem, cptr_mem_shape_i,
        cptr_mem_shape_j, cptr_rfglobal, cptr_rfof_input_param);
#endif
};

/**
 * @brief Initialises rfof diagnostics. These are not used but are as of now 
 * given as dummy inputs to the kick routine. To be called before the time step 
 * loop.
 * @param cptr_rfglobal Handle to the Fortan wave field struct
 * @param cptr_diagno Handle to the Fortran diagnostics struct
*/
void rfof_interface_initialise_diagnostics(void** cptr_rfglobal,
    void** cptr_diagno) {
#ifdef RFOF  
    __valipalikka_MOD_call_initialise_diagnostics(cptr_rfglobal, cptr_diagno);
#endif
};

/**
 * @brief Allocates memory for an rfof marker on the fortran side. To be called
 *  before the time step loop. The memory should be deallocated at the end by 
 * calling the corresponding deallocation routine at the end of this file.
 * @param rfof_marker_pointer Handle to rfof marker struct.
*/
void rfof_interface_allocate_rfof_marker(void* rfof_marker_pointer) {
#ifdef RFOF
    real dummy_real_rfof = 0.0;
    int dummy_int_rfof = 0;
    int is_already_allocated = 0; /**< Needs to be zero. */
    __valipalikka_MOD_call_set_marker_pointers(
                &rfof_marker_pointer, 
                &dummy_int_rfof,
                &dummy_real_rfof,
                &dummy_real_rfof,
                &dummy_real_rfof,
                &dummy_real_rfof,
                &dummy_real_rfof, 
                &dummy_real_rfof,
                &dummy_real_rfof,
                &dummy_real_rfof,
                &dummy_real_rfof,
                &dummy_real_rfof,
                &dummy_real_rfof, 
                &dummy_real_rfof, 
                &dummy_real_rfof, 
                &dummy_real_rfof, 
                &dummy_real_rfof, 
                &dummy_real_rfof,
                &dummy_int_rfof,
                &is_already_allocated);
#endif
};


/* STUFF TO DO BETWEEN KICKS */

/**
 * @brief Creates a marker for the RFOF routines on the fortran side.
 * 
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
 * direction (unclear if truly needed)
 * @param acc Time acceleration factor (not in use anymore)
 * @param isOrbitTimeAcclerated False as implied above in acc. Logical in 
 * fortran; represented by int (4 bytes) in C.
*/
void rfof_interface_set_marker_pointers(void** cptr_marker, int* id,
    real* weight, real* R, real* phi, real* z, real* psi, real* charge,
    real* mass, real* Ekin, real* velocity, real* mu, real* pphicanonical,
    real* vpar, real* vperp, real* gyrof, real* vdriftRho, real* acc,
    int* isOrbitTimeAccelerated, int* is_already_allocated) {
#ifdef RFOF
    __valipalikka_MOD_call_set_marker_pointers(cptr_marker, id, weight, R, phi,
        z, psi, charge, mass, Ekin, velocity, mu, pphicanonical, vpar, vperp,
        gyrof, vdriftRho, acc, isOrbitTimeAccelerated, is_already_allocated);
#endif
};


/**
 * @brief Function to be called in the main simulation loop during each step 
 * when following the guiding centre.  
 * 
 * 1. Creates the rfof_marker based on the given input ascot_marker.
 * 2. Calls the "kick" function, which 
 *      a) Checks resonance condition and 
 *      b) if in resonance, kicks marker and updates velocity of the rfof marker 
 *          and consequently also ascot marker (as only pointers are passed when 
 *          creating the rfof marker).
 * @param ascot_marker Ascot marker (no points for those who guessed it)
 * @param hin Time step proposed by the simulation loop
 * @param hout_rfof If negative, ICRH kick failed because hin was too large. In 
 * that case the absolute value of hout_rfof should be used when redoing the 
 * simulation loop. For successfull steps, this is the estimate to the next 
 * resonance.
 * @param rfof_data A "package" of ICRH realted void pointers to Fortran 
 * routines
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
            
            real psi; /**< Poloidal flux function (following the ITM and ITER 
            conventions; COCOS 13/11; not divided by \f$ 2\pi \f$) */
            B_field_eval_psi(&psi, ascot_marker->r[i], ascot_marker->phi[i],
                ascot_marker->z[i], ascot_marker->time[i], Bdata);

            real B = sqrt(pow(ascot_marker->B_r[i],2) + 
                pow(ascot_marker->B_phi[i],2) + pow(ascot_marker->B_z[i],2)); 
                /**< Norm of B field  */

            real gamma = physlib_gamma_ppar(ascot_marker->mass[i],
                ascot_marker->mu[i], ascot_marker->ppar[i], B); /**< Lorentz 
                factor.*/

            real Ekin = physlib_Ekin_gamma(ascot_marker->mass[i], gamma); /**< 
                Kinetic energy*/

            real speed = physlib_vnorm_gamma(gamma); /**< Speed */

            real p_phi = phys_ptoroid_gc(ascot_marker->charge[i],
                ascot_marker->r[i], ascot_marker->ppar[i], psi, B,
                ascot_marker->B_phi[i]); /* canonical momentum conjugate to phi 
                (toroidal momenutum) */

            real p = physlib_pnorm_vnorm(ascot_marker->mass[i], speed); /**< 
                Momentum norm */

            real xi = physlib_gc_xi(ascot_marker->mass[i], ascot_marker->mu[i],
                ascot_marker->ppar[i], B); /**< pitch */

            real v_par = physlib_gc_ppar(p, xi); /**< Parallel momentum */

            real v_perp = phys_vperp_gc(speed, v_par); /**<Perpendicular speed*/

            real gyrof = phys_gyrofreq_ppar(ascot_marker->mass[i],
                ascot_marker->charge[i], ascot_marker->mu[i],
                ascot_marker->ppar[i], B); /**< Gyrofrequency */

            real vdriftRho = 0; /**< Velocity drift in the direction of rho. 
            NOT NEEDED?*/
            real acc = 0; /**< Time acceleration factor. NOT NEEDED? */
            int isOrbitTimeAccelerated = 0; /**< Time accelecation on/off. NOTE:
            fortran logical corresponds to c int (4 bytes). */
            int is_already_allocated = 1;

            /* Now when we have all the needed parameters we call the rfof 
            routine which sets the pointers of the rfof_marker (on the Fortran 
            side. The situation is somewhat comparable to ordering pizza/kebab 
            from Kontula: it's convenient that the place exists but you'd rather
            not go there yourself. We're in luck, though, as we happen to have a
            courier.) */
            
            /* TODO: Remove memory leak */
            int * dummy = (int*)malloc(sizeof(int)); 

            __valipalikka_MOD_call_set_marker_pointers(
                &rfof_marker_pointer_array[i], /* Note that the pointer to the 
                pointer is passed. */
                /*&(ascot_marker->id[i]),*/
                dummy,
                &(ascot_marker->weight[i]), /* Number of real particles 
                                                represented by the marker     */
                &(ascot_marker->r[i]),
                &(ascot_marker->phi[i]),
                &(ascot_marker->z[i]),
                &psi, 
                &(ascot_marker->charge[i]),
                &(ascot_marker->mass[i]),
                &Ekin,
                &speed,
                &(ascot_marker->mu[i]),
                &p_phi,                     /* pphicanonical                  */
                &v_par, 
                &v_perp, 
                &gyrof, 
                &vdriftRho,                 /* vdriftRho; EVALUATE IF ACTUALLY 
                                               NEEDED                         */
                &acc,
                &isOrbitTimeAccelerated,
                &is_already_allocated);
            
            /*Used for storing the "results" of calling RF kick, appears to be 
            redundant (?) */
            prt_rfof rfof_data_pack = {
                .dmu = 0.0,
                .dvpar = 0.0,
                .de = 0.0,
                .deCumulative = 0.0,
                .dpitch = 0.0,
                .maxAcc = 0.0,
                .RFdt = 0.0,
            };
            
            int err = 0;       /* "empty" input to kick                       */

            int mpiprocid = 1; /* Number used to identify MPI nods during 
                                  parallel execution. This process is not to be 
                                  parallelized (?)                            */

            /* Ready to kick some ash */
            __valipalikka_MOD_call_rf_kick(&(ascot_marker->time[i]), &(hin[i]),
                &mpiprocid, &rfof_data_pack, &rfof_marker_pointer_array[i],
                &(rfof_mem_pointer_array[i]), &(rfof_data.cptr_rfglobal),
                &(rfof_diag_pointer_array[i]), &err, &(mem_shape_i[i]),
                &(mem_shape_j[i]));

            /* Check if dt was sufficiently small and assign hout_rfof 
            accrodingly. */
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
                hout_rfof[i] = rfof_data_pack.RFdt;
            }

            /*The momenta of the ascot marker should now have been updated since
            pointers to those quantities were passed to RFOF. Thus, no need to
            manually update them here. */
        };  
   };
#endif
};


/**
 * @brief Function to be called in the main simulation loop during each step for 
 * full orbit.
 * 
 * 1. Creates the rfof_marker based on the given input ascot_marker.
 * 2. Calls the "kick" function, which 
 *      a) Checks resonance condition and 
 *      b) if in resonance, kicks marker and updates velocity of the rfof marker
 *          and consequently also ascot marker (as only pointers are passed when
 *          creating the rfof marker).
 * @param ascot_marker Ascot marker (no points for those who guessed it)
 * @param hin Time step proposed by the simulation loop
 * @param hout_rfof If negative, ICRH kick failed because hin was too large. In 
 * that case the absolute value of hout_rfof should be used when redoing the 
 * simulation loop. For successfull steps, this is the estimate to the next 
 * resonance.
 * @param rfof_data A "package" of ICRH realted void pointers to Fortran 
 * routines
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
//            __valipalikka_MOD_call_set_marker_pointers(
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
 * @brief Resets resonance memory of ICRH (RFOF) markers. 
 * @param rfof_mem_pointer Handle to rfof resonance memory.
 * @param mem_shape_i Array of rfof resonance memory matrix's first dimensions.
 * @param mem_shape_j Array of rfof resonance memory matrix's second dimensions.
 */
void rfof_interface_reset_icrh_mem(void** rfof_mem_pointer, int* mem_shape_i,
    int* mem_shape_j) {
#ifdef RFOF
    __valipalikka_MOD_call_reset_res_mem(rfof_mem_pointer, mem_shape_i,
        mem_shape_j);
#endif
};


/* DEALLOCATION ROUTINES */

/** 
 * @brief Deallocates the rfof_input_param struct on the fortran side. 
 * @brief cptr_rfof_input_param Handle to rfof input param struct on the Fortran
 * side.
 */
void rfof_interface_deallocate_rfof_input_param(void** cptr_rfof_input_param) {
#ifdef RFOF
    __valipalikka_MOD_call_deallocate_rfof_input_param(cptr_rfof_input_param);
#endif
};

/** 
 * @brief Deallocates the rfglobal struct (wave field) on the fortran side. *
 * @brief cptr_rfglobal Handle to rfof wave field on the Fortran side.
 */
void rfof_interface_deallocate_rfglobal(void** cptr_rfglobal) {
#ifdef RFOF
    __valipalikka_MOD_call_deallocate_rfglobal(cptr_rfglobal);
#endif
};

/** 
 * @brief Deallocates the resonance memory matrix of a particle.
 * @param rfof_mem_pointer Handle to rfof resonance memory. 
 * @param mem_shape_i Array of rfof resonance memory matrix's first dimensions.
 * @param mem_shape_j Array of rfof resonance memory matrix's second dimensions.
*/
void rfof_interface_deallocate_res_mem(void** cptr_res_mem,
    int* cptr_mem_shape_i, int* cptr_mem_shape_j) {
#ifdef RFOF
    __valipalikka_MOD_call_deallocate_res_mem(cptr_res_mem, cptr_mem_shape_i,
        cptr_mem_shape_j);
#endif
};

/** 
 * @brief Deallocates the (dummy) dianostics of rfof markers.
 * @param cptr_diagno Handle to rfof diagnostics struct.
 */
void rfof_interface_deallocate_diagnostics(void** cptr_diagno) {
#ifdef RFOF
    __valipalikka_MOD_call_deallocate_diagnostics(cptr_diagno);
#endif
};

/** 
 * @brief Deallocates the rfof marker.
 * @param cptr_rfof_marker Handle to rfof marker.
 */
void rfof_interface_deallocate_marker(void** cptr_rfof_marker) {
#ifdef RFOF
    __valipalikka_MOD_deallocate_marker(cptr_rfof_marker);
#endif
};