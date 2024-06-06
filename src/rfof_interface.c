/*
 * @brief Contains the function to be called from the simulation loop and all the functions called within it.
 *
 *
**/
#include <stdio.h>
#include <string.h>
#include "particle.h"
#include "rfof_interface.h"
#define RFOF 1
#ifdef RFOF
    


// init everything funktio



/**
 * @brief Function to be called in the main simulation loop. 
 * 
 * 1. Sets the rfof_marker based on the given input ascot_marker.
 * 2. Call the "kick" function, which 
 *      a) Checks resonance condition and 
 *      b) if in resonance, kicks marker and updates velocity.
*/

void do_rfof_stuff(particle_simd_fo* ascot_marker, real hin, real hout_rfof, rfof_data rfof_data) {
    for (int i=0; i<NSIMD;i++) {
    //If id > 0 and running
        if (ascot_marker->id(i) > 0 && ascot_marker->running(i)) {
            //Some of the parameters needed by the RFOF_marker are not members of the ascot_marker struct and thus these need to be evaluated first.

            
            //Non when we have all the needed parameters we call the rfof routine which sets the pointers of the rfof_marker

            set_marker_pointers(marker = marker   , &
                id                = id          , &
                weight            = weight      , &
                r                 = R           , &  !used to be rpz(1)
                phi               = phi         , &  !used to be rpz(2)
                z                 = z           , &  !used to be rpz(3)
                psi               = psi         , &
                rho_tor           = dummy       , &
                theta             = dummy       , &
                charge            = charge      , &
                mass              = mass        , &
                energy            = Ekin        , &
                energy_kinetic    = Ekin        , &
                velocity          = velocity    , &
                magneticMoment    = mu          , &
                pphi              = pphicanonical , &
                vpar              = vpar        , &
                vperp             = vperp       , &
                omega_gyro        = gyrof       , &
                taubounce         = dummy       , &
                vdrift            = dummy       , &
                vdriftRho         = vdriftRho   , &
                vdriftDia         = dummy       , &
                d_vpar_d_rho      = dummy       , &
                d_vpar_d_dia      = dummy       , &
                d_vperp_d_rho     = dummy       , &
                d_vperp_d_dia     = dummy       , &
                d_vDriftRho_d_rho = dummy       , &
                d_vDriftRho_d_dia = dummy       , &
                d_vDriftDia_d_rho = dummy       , &
                d_vDriftDia_d_dia = dummy       , &
                time_acceleration = acc           , &
                isOrbitTimeAccelerated = isOrbitTimeAccelerated);
            
            };
    //evaluate energy, omega_gyro,...
    //call_set_marker_ponter_c
    //call rfof kick
    //check if dt sufficiently small
        //update particle_simd_fo (p_r. p_phi, p_z) using Energy or velocities (check if consistent)
  
   };
};
#endif