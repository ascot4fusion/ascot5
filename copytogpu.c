/**
 * @file copytogpu.c
 * @brief Transfer data from and to the GPU
 */
#include "simulate.h"
#include "particle.h"
#include "B_field.h"
#include "E_field.h"
#include "copytogpu.h"

/**
 * @brief Copy data from CPU to GPU
*/
void simulate_fo_fixed_copy_to_gpu(sim_data* sim, particle_simd_fo *p_ptr, particle_simd_fo *p0_ptr, B_field_data* Bdata, E_field_data* Edata, particle_loc*  p_loc, real* hin, real* rnd, int n_queue_size) {

  GPU_MAP_TO_DEVICE(
		      p_loc[0:1],\
		      p_loc->r_arr1[0:n_queue_size],\
		      p_loc->r_arr2[0:n_queue_size],\
		      p_loc->r_arr3[0:n_queue_size],\
		      p_loc->r_arr4[0:n_queue_size],\
		      p_loc->r_arr5[0:n_queue_size],\
		      p_loc->i_arr1[0:n_queue_size],\
		      p_loc->i_arr2[0:n_queue_size],\
		      p_loc->i_arr3[0:n_queue_size],\
		      p_loc->i_arr4[0:n_queue_size],\
		      p_loc->i_arr5[0:n_queue_size],\
		      p_loc->i_arr6[0:n_queue_size],\
		      p_loc->i_arr7[0:n_queue_size],\
		      p_loc->i_arr8[0:n_queue_size],\
		      p_loc->i_arr9[0:n_queue_size],\
  		      p_ptr[0:1],\
		      p_ptr->running        [0:n_queue_size],\
		      p_ptr->r              [0:n_queue_size],\
		      p_ptr->phi            [0:n_queue_size],\
		      p_ptr->p_r            [0:n_queue_size],\
		      p_ptr->p_phi          [0:n_queue_size],\
		      p_ptr->p_z            [0:n_queue_size],\
		      p_ptr->mileage        [0:n_queue_size],\
		      p_ptr->z              [0:n_queue_size],\
		      p_ptr->charge         [0:n_queue_size],\
		      p_ptr->mass           [0:n_queue_size],\
		      p_ptr->B_r            [0:n_queue_size],\
		      p_ptr->B_r_dr         [0:n_queue_size],\
		      p_ptr->B_r_dphi       [0:n_queue_size],\
		      p_ptr->B_r_dz         [0:n_queue_size],\
		      p_ptr->B_phi          [0:n_queue_size],\
		      p_ptr->B_phi_dr       [0:n_queue_size],\
		      p_ptr->B_phi_dphi     [0:n_queue_size],\
		      p_ptr->B_phi_dz       [0:n_queue_size],\
		      p_ptr->B_z            [0:n_queue_size],\
		      p_ptr->B_z_dr         [0:n_queue_size],\
		      p_ptr->B_z_dphi       [0:n_queue_size],\
		      p_ptr->B_z_dz         [0:n_queue_size],\
		      p_ptr->rho            [0:n_queue_size],\
		      p_ptr->theta          [0:n_queue_size],\
		      p_ptr->err            [0:n_queue_size],\
		      p_ptr->time           [0:n_queue_size],\
		      p_ptr->weight         [0:n_queue_size],\
		      p_ptr->cputime        [0:n_queue_size],\
		      p_ptr->id             [0:n_queue_size],\
		      p_ptr->endcond        [0:n_queue_size],\
		      p_ptr->walltile       [0:n_queue_size],\
		      p_ptr->index          [0:n_queue_size],\
		      p_ptr->znum           [0:n_queue_size],\
		      p_ptr->anum           [0:n_queue_size],\
		      p_ptr->bounces        [0:n_queue_size],\
  		      p0_ptr[0:1],\
		      p0_ptr->running       [0:n_queue_size],\
		      p0_ptr->r             [0:n_queue_size],\
		      p0_ptr->phi           [0:n_queue_size],\
		      p0_ptr->p_r           [0:n_queue_size],\
		      p0_ptr->p_phi         [0:n_queue_size],\
		      p0_ptr->p_z           [0:n_queue_size],\
		      p0_ptr->mileage       [0:n_queue_size],\
		      p0_ptr->z             [0:n_queue_size],\
		      p0_ptr->charge        [0:n_queue_size],\
		      p0_ptr->mass          [0:n_queue_size],\
		      p0_ptr->B_r           [0:n_queue_size],\
		      p0_ptr->B_r_dr        [0:n_queue_size],\
		      p0_ptr->B_r_dphi      [0:n_queue_size],\
		      p0_ptr->B_r_dz        [0:n_queue_size],\
		      p0_ptr->B_phi         [0:n_queue_size],\
		      p0_ptr->B_phi_dr      [0:n_queue_size],\
		      p0_ptr->B_phi_dphi    [0:n_queue_size],\
		      p0_ptr->B_phi_dz      [0:n_queue_size],\
		      p0_ptr->B_z           [0:n_queue_size],\
		      p0_ptr->B_z_dr        [0:n_queue_size],\
		      p0_ptr->B_z_dphi      [0:n_queue_size],\
		      p0_ptr->B_z_dz        [0:n_queue_size],\
		      p0_ptr->rho           [0:n_queue_size],\
		      p0_ptr->theta         [0:n_queue_size],\
		      p0_ptr->err           [0:n_queue_size],\
		      p0_ptr->time          [0:n_queue_size],\
		      p0_ptr->weight        [0:n_queue_size],\
		      p0_ptr->cputime       [0:n_queue_size],\
		      p0_ptr->id            [0:n_queue_size],\
		      p0_ptr->endcond       [0:n_queue_size],\
		      p0_ptr->walltile      [0:n_queue_size],\
		      p0_ptr->index         [0:n_queue_size],\
		      p0_ptr->znum          [0:n_queue_size],\
		      p0_ptr->anum          [0:n_queue_size],\
		      p0_ptr->bounces       [0:n_queue_size],\
		      hin[0:n_queue_size],\
       		      sim[0:1],		\
		      Bdata[0:1],\
		      Edata[0:1],Edata->type, \
		      rnd[0:3*n_queue_size] \
			)

    switch(sim->wall_data.type) {
    case wall_type_2D:
      GPU_MAP_TO_DEVICE(
			sim->wall_data.w2d.wall_r[0:sim->wall_data.w2d.n],sim->wall_data.w2d.wall_z[0:sim->wall_data.w2d.n] )
      break;
    case wall_type_3D:
      GPU_MAP_TO_DEVICE(
			sim->wall_data.w3d.wall_tris[0:sim->wall_data.w3d.n*9+9],sim->wall_data.w3d.tree_array[0:sim->wall_data.w3d.tree_array_size] )
      break;
    default:
      break;
    }

    if(sim->diag_data.dist5D_collect) {
      GPU_MAP_TO_DEVICE(
			sim->diag_data.dist5D.histogram[0:sim->diag_data.dist5D.n_r * sim->diag_data.dist5D.n_phi * sim->diag_data.dist5D.n_z * sim->diag_data.dist5D.n_ppara * sim->diag_data.dist5D.n_pperp * sim->diag_data.dist5D.n_time * sim->diag_data.dist5D.n_q] )
    }

    if(sim->diag_data.dist6D_collect) {
      GPU_MAP_TO_DEVICE(
			sim->diag_data.dist6D.histogram[0:sim->diag_data.dist6D.n_r * sim->diag_data.dist6D.n_phi * sim->diag_data.dist6D.n_z * sim->diag_data.dist6D.n_pr * sim->diag_data.dist6D.n_pphi * sim->diag_data.dist6D.n_pz * sim->diag_data.dist6D.n_time * sim->diag_data.dist6D.n_q] )
    }

    if(sim->diag_data.distrho5D_collect) {
      GPU_MAP_TO_DEVICE(
			sim->diag_data.distrho5D.histogram[0:sim->diag_data.distrho5D.n_rho * sim->diag_data.distrho5D.n_theta * sim->diag_data.distrho5D.n_phi * sim->diag_data.distrho5D.n_ppara * sim->diag_data.distrho5D.n_pperp * sim->diag_data.distrho5D.n_time * sim->diag_data.distrho5D.n_q] )
    }

    if(sim->diag_data.distrho6D_collect) {
      GPU_MAP_TO_DEVICE(
			sim->diag_data.distrho6D.histogram[0:sim->diag_data.distrho6D.n_rho*sim->diag_data.distrho6D.n_theta*sim->diag_data.distrho6D.n_phi*sim->diag_data.distrho6D.n_pr*sim->diag_data.distrho6D.n_pphi*sim->diag_data.distrho6D.n_pz*sim->diag_data.distrho6D.n_time*sim->diag_data.distrho6D.n_q] )
    }

    if(sim->diag_data.distCOM_collect) {
      GPU_MAP_TO_DEVICE(
			sim->diag_data.distCOM.histogram[0:sim->diag_data.distCOM.n_mu*sim->diag_data.distCOM.n_Ekin*sim->diag_data.distCOM.n_Ptor] )
    }


  switch(Edata->type) {

    case E_field_type_1DS:
      GPU_MAP_TO_DEVICE(
			Edata->E1DS,Edata->E1DS.dV,Edata->E1DS.dV.c[0:Edata->E1DS.dV.n_x*NSIZE_COMP1D] )
      break;
      
    case E_field_type_TC:
      GPU_MAP_TO_DEVICE(
			Edata->ETC,Edata->ETC.Exyz[0:1] )
      break;
      
    default:
      break;
    }

    switch(sim->plasma_data.type) {

    case plasma_type_1D:
      GPU_MAP_TO_DEVICE(
		      sim->plasma_data.plasma_1D.mass       [0:MAX_SPECIES],\
		      sim->plasma_data.plasma_1D.charge     [0:MAX_SPECIES],\
		      sim->plasma_data.plasma_1D.anum       [0:MAX_SPECIES],\
		      sim->plasma_data.plasma_1D.znum       [0:MAX_SPECIES],\
		      sim->plasma_data.plasma_1D.rho        [0:sim->plasma_data.plasma_1D.n_rho],\
		      sim->plasma_data.plasma_1D.temp       [0:sim->plasma_data.plasma_1D.n_rho*sim->plasma_data.plasma_1D.n_species], \
  		      sim->plasma_data.plasma_1D.dens       [0:sim->plasma_data.plasma_1D.n_rho*sim->plasma_data.plasma_1D.n_species] )
      break;
			
      case plasma_type_1Dt:
	GPU_MAP_TO_DEVICE(
		      sim->plasma_data.plasma_1Dt.mass      [0:MAX_SPECIES],\
		      sim->plasma_data.plasma_1Dt.charge    [0:MAX_SPECIES],\
		      sim->plasma_data.plasma_1Dt.anum      [0:MAX_SPECIES],\
		      sim->plasma_data.plasma_1Dt.znum      [0:MAX_SPECIES],\
		      sim->plasma_data.plasma_1Dt.rho       [0:sim->plasma_data.plasma_1Dt.n_rho],\
		      sim->plasma_data.plasma_1Dt.temp      [0:sim->plasma_data.plasma_1Dt.n_time*sim->plasma_data.plasma_1Dt.n_rho*sim->plasma_data.plasma_1Dt.n_species],\
		      sim->plasma_data.plasma_1Dt.dens      [0:sim->plasma_data.plasma_1Dt.n_rho*sim->plasma_data.plasma_1Dt.n_species*sim->plasma_data.plasma_1Dt.n_time],\
		      sim->plasma_data.plasma_1Dt.time      [0:sim->plasma_data.plasma_1Dt.n_time] )
      break;

      case plasma_type_1DS:
	GPU_MAP_TO_DEVICE(
		      sim->plasma_data.plasma_1DS.mass      [0:MAX_SPECIES],\
		      sim->plasma_data.plasma_1DS.charge    [0:MAX_SPECIES],\
		      sim->plasma_data.plasma_1DS.anum      [0:MAX_SPECIES],\
		      sim->plasma_data.plasma_1DS.znum      [0:MAX_SPECIES],\
		      sim->plasma_data.plasma_1DS.temp      [0:2],\
		      sim->plasma_data.plasma_1DS.dens      [0:MAX_SPECIES],\
		      sim->plasma_data.plasma_1DS.temp[0].c [0:sim->plasma_data.plasma_1DS.temp[0].n_x*NSIZE_COMP1D],\
		      sim->plasma_data.plasma_1DS.temp[1].c [0:sim->plasma_data.plasma_1DS.temp[1].n_x*NSIZE_COMP1D] )
	for (int i=0;i<MAX_SPECIES;i++) {
	  GPU_MAP_TO_DEVICE(
			    sim->plasma_data.plasma_1DS.dens[i].c[0:sim->plasma_data.plasma_1DS.dens[i].n_x*NSIZE_COMP1D] )
    }
      break;

      default:
      break;
    }


    switch(Bdata->type) {

    case B_field_type_GS:
      GPU_MAP_TO_DEVICE(
			Bdata->BGS.psi_coeff[0:13] )
      break;
      
    case B_field_type_2DS:
      GPU_MAP_TO_DEVICE(
		      Bdata->B2DS.psi,    Bdata->B2DS.psi.c    [0:Bdata->B2DS.psi.n_x   *Bdata->B2DS.psi.n_y                          *NSIZE_COMP2D],\
  		      Bdata->B2DS.B_r,    Bdata->B2DS.B_r.c    [0:Bdata->B2DS.B_r.n_x   *Bdata->B2DS.B_r.n_y                          *NSIZE_COMP2D],\
		      Bdata->B2DS.B_phi,  Bdata->B2DS.B_phi.c  [0:Bdata->B2DS.B_phi.n_x *Bdata->B2DS.B_phi.n_y                        *NSIZE_COMP2D],\
		      Bdata->B2DS.B_z,    Bdata->B2DS.B_z.c    [0:Bdata->B2DS.B_z.n_x   *Bdata->B2DS.B_z.n_y                          *NSIZE_COMP2D] )
      break;
      
    case B_field_type_3DS:
      GPU_MAP_TO_DEVICE(
			Bdata->B3DS.psi,    Bdata->B3DS.psi.c    [0:Bdata->B3DS.psi.n_x   *Bdata->B3DS.psi.n_y                          *NSIZE_COMP2D],	\
			Bdata->B3DS.B_r,    Bdata->B3DS.B_r.c    [0:Bdata->B3DS.B_r.n_x   *Bdata->B3DS.B_r.n_y   *Bdata->B3DS.B_r.n_z   *NSIZE_COMP3D],	\
			Bdata->B3DS.B_phi,  Bdata->B3DS.B_phi.c  [0:Bdata->B3DS.B_phi.n_x *Bdata->B3DS.B_phi.n_y *Bdata->B3DS.B_phi.n_z *NSIZE_COMP3D],	\
			Bdata->B3DS.B_z,    Bdata->B3DS.B_z.c    [0:Bdata->B3DS.B_z.n_x   *Bdata->B3DS.B_z.n_y   *Bdata->B3DS.B_z.n_z   *NSIZE_COMP3D] )

      break;
      
    case B_field_type_STS:
      GPU_MAP_TO_DEVICE(

			Bdata->BSTS.axis_r, Bdata->BSTS.axis_r.c [0:Bdata->BSTS.axis_r.n_x                                                           ], \
			Bdata->BSTS.axis_z, Bdata->BSTS.axis_z.c [0:Bdata->BSTS.axis_z.n_x                                                           ],	\
			Bdata->BSTS.psi,    Bdata->BSTS.psi.c    [0:Bdata->BSTS.psi.n_x   *Bdata->BSTS.psi.n_y   *Bdata->BSTS.psi.n_z   *NSIZE_COMP3D],	\
			Bdata->BSTS.B_r,    Bdata->BSTS.B_r.c    [0:Bdata->BSTS.B_r.n_x   *Bdata->BSTS.B_r.n_y   *Bdata->BSTS.B_r.n_z   *NSIZE_COMP3D],	\
			Bdata->BSTS.B_z,    Bdata->BSTS.B_z.c    [0:Bdata->BSTS.B_z.n_x   *Bdata->BSTS.B_z.n_y   *Bdata->BSTS.B_z.n_z   *NSIZE_COMP3D],	\
			Bdata->BSTS.B_phi,  Bdata->BSTS.B_phi.c  [0:Bdata->BSTS.B_phi.n_x *Bdata->BSTS.B_phi.n_y *Bdata->BSTS.B_phi.n_z *NSIZE_COMP3D] )
      break;
      
    case B_field_type_TC:
      GPU_MAP_TO_DEVICE(
			Bdata->BTC.B[0:3],Bdata->BTC.dB[0:9] )
      break;
      
    default:
      break;
    }

}

/**
 * @brief Copy data from GPU to CPU
*/

void simulate_fo_fixed_copy_from_gpu(sim_data* sim, particle_simd_fo *p_ptr, int n_queue_size){

  GPU_UPDATE_FROM_DEVICE(
      p_ptr->running[0:n_queue_size],p_ptr->r[0:n_queue_size],p_ptr->phi[0:n_queue_size],p_ptr->p_r[0:n_queue_size],p_ptr->p_phi[0:n_queue_size],p_ptr->p_z[0:n_queue_size],p_ptr->mileage[0:n_queue_size], \
  p_ptr->z[0:n_queue_size],p_ptr->charge[0:n_queue_size],p_ptr->mass[0:n_queue_size],p_ptr->B_r[0:n_queue_size],p_ptr->B_r_dr[0:n_queue_size],p_ptr->B_r_dphi[0:n_queue_size],p_ptr->B_r_dz[0:n_queue_size], \
  p_ptr->B_phi[0:n_queue_size],p_ptr->B_phi_dr[0:n_queue_size],p_ptr->B_phi_dphi[0:n_queue_size],p_ptr->B_phi_dz[0:n_queue_size],p_ptr->B_z[0:n_queue_size],p_ptr->B_z_dr[0:n_queue_size],p_ptr->B_z_dphi[0:n_queue_size], \
  p_ptr->B_z_dz[0:n_queue_size],p_ptr->rho[0:n_queue_size],p_ptr->theta[0:n_queue_size],p_ptr->err[0:n_queue_size],p_ptr->time[0:n_queue_size],p_ptr->weight[0:n_queue_size],p_ptr->cputime[0:n_queue_size], \
      p_ptr->id[0:n_queue_size],p_ptr->endcond[0:n_queue_size],p_ptr->walltile[0:n_queue_size],p_ptr->index[0:n_queue_size],p_ptr->znum[0:n_queue_size],p_ptr->anum[0:n_queue_size],p_ptr->bounces[0:n_queue_size] )

    GPU_MAP_FROM_DEVICE(
			      sim[0:1]  )
}
