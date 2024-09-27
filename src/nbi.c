/**
 * @file nbi.c
 * @brief Functions for NBI simulation and particle generation
 */
#include <stdlib.h>
#include <math.h>
#include "print.h"
#include "ascot5.h"
#include "consts.h"
#include "physlib.h"
#include "math.h"
#include "random.h"
#include "nbi.h"

/**
 * @brief Initialize NBI data struct on target
 *
 * @param nbi pointer to data struct on target
 * @param offload_data pointer to offload data struct
 * @param offload_array pointer to offload array
 */
void nbi_init(nbi_data* data, int ninj, int* id, int* anum, int* znum,
              real* mass, real* power, real* efrac, real* energy,
              real* div_h, real* div_v, real* div_halo_v, real* div_halo_h,
              real* div_halo_frac, int* nbeamlet, real* beamlet_xyz) {
    int idx = 0;
    data->ninj =ninj;
    data->inj = (nbi_injector*) malloc( ninj*sizeof(nbi_injector) );
    for(int i=0; i<data->ninj; i++) {
        data->inj[i].anum          = anum[i];
        data->inj[i].znum          = znum[i];
        data->inj[i].mass          = mass[i];
        data->inj[i].power         = power[i];
        data->inj[i].energy        = energy[i];
        data->inj[i].efrac[0]      = efrac[3*i+0];
        data->inj[i].efrac[1]      = efrac[3*i+1];
        data->inj[i].efrac[2]      = efrac[3*i+2];
        data->inj[i].div_h         = div_h[i];
        data->inj[i].div_v         = div_v[i];
        data->inj[i].div_halo_frac = div_halo_frac[i];
        data->inj[i].div_halo_h    = div_halo_h[i];
        data->inj[i].div_halo_v    = div_halo_v[i];
        data->inj[i].id            = id[i];
        data->inj[i].n_beamlet     = nbeamlet[i];

        int n_beamlet = data->inj[i].n_beamlet;
        data->inj[i].beamlet_x = (real*) malloc( n_beamlet*sizeof(real) );
        data->inj[i].beamlet_y = (real*) malloc( n_beamlet*sizeof(real) );
        data->inj[i].beamlet_z = (real*) malloc( n_beamlet*sizeof(real) );
        data->inj[i].beamlet_dx = (real*) malloc( n_beamlet*sizeof(real) );
        data->inj[i].beamlet_dy = (real*) malloc( n_beamlet*sizeof(real) );
        data->inj[i].beamlet_dz = (real*) malloc( n_beamlet*sizeof(real) );
        for(int j = 0; j < n_beamlet; j++) {
            data->inj[i].beamlet_x[j]  = beamlet_xyz[idx + 0*n_beamlet + j];
            data->inj[i].beamlet_y[j]  = beamlet_xyz[idx + 1*n_beamlet + j];
            data->inj[i].beamlet_z[j]  = beamlet_xyz[idx + 2*n_beamlet + j];
            data->inj[i].beamlet_dx[j] = beamlet_xyz[idx + 3*n_beamlet + j];
            data->inj[i].beamlet_dy[j] = beamlet_xyz[idx + 4*n_beamlet + j];
            data->inj[i].beamlet_dz[j] = beamlet_xyz[idx + 5*n_beamlet + j];
        }
        idx += 6 * n_beamlet;
    }

    int err = 0;
    print_out(VERBOSE_IO, "\nNBI input\n");
    print_out(VERBOSE_IO, "Number of injectors %d:\n", data->ninj);
    for(int i=0; i < data->ninj; i++) {
        print_out(VERBOSE_IO, "\n  Injector ID %d (%d beamlets) Power: %1.1e\n",
                  data->inj[i].id, data->inj[i].n_beamlet,
                  data->inj[i].power);
        print_out(VERBOSE_IO,
                  "    Anum %d Znum %d mass %1.1e amu energy %1.1e eV\n",
                  data->inj[i].anum, data->inj[i].znum,
                  data->inj[i].mass / CONST_U, data->inj[i].energy / CONST_E);
        print_out(VERBOSE_IO,
                  "    Energy fractions: %1.1e (Full) %1.1e (1/2) %1.1e (1/3)\n",
                  data->inj[i].efrac[0], data->inj[i].efrac[1],
                  data->inj[i].efrac[2]);

        /* Even if halo fraction is zero, the divergences should be nonzero
           to avoid division by zero during evaluation. Do this after the
           input has been printed as to not confuse the user */
        if(data->inj[i].div_halo_frac == 0) {
            data->inj[i].div_halo_h = 1e-10;
            data->inj[i].div_halo_v = 1e-10;
        }
    }
    return err;
}

/**
 * @brief Free allocated resources
 *
 * @param data pointer to the data struct
 */
void nbi_free(nbi_data* data) {
    for(int i=0; i<data->ninj; i++) {
        data->inj[i].n_beamlet = 0;
        free(data->inj[i].beamlet_x);
        free(data->inj[i].beamlet_y);
        free(data->inj[i].beamlet_z);
        free(data->inj[i].beamlet_dx);
        free(data->inj[i].beamlet_dy);
        free(data->inj[i].beamlet_dz);
    }
    data->ninj = 0;
    free(data->inj);
}

/**
 * @brief Sample injected marker's coordinates.
 *
 * @param inj pointer to injector data
 * @param xyz initialized marker's position in cartesian coordinates [m]
 * @param vxyz initialized marker's velocity in cartesian coordinates [m/s]
 * @param rng pointer to random number generator data
 */
void nbi_inject(real* xyz, real* vxyz, nbi_injector* inj, random_data* rng) {
    /* Pick a random beamlet and initialize marker there */
    int i_beamlet = floor(random_uniform(rng) * inj->n_beamlet);
    xyz[0] = inj->beamlet_x[i_beamlet];
    xyz[1] = inj->beamlet_y[i_beamlet];
    xyz[2] = inj->beamlet_z[i_beamlet];

    /* Pick marker energy based on the energy fractions */
    real energy;
    real r = random_uniform(rng);
    if(r < inj->efrac[0]) {
        energy = inj->energy;
    } else if(r < inj->efrac[0] + inj->efrac[1]) {
        energy = inj->energy / 2;
    } else {
        energy = inj->energy / 3;
    }

    /* Calculate vertical and horizontal normals for beam divergence */
    real dir[3], normalv[3], normalh[3], tmp[3];
    dir[0] = inj->beamlet_dx[i_beamlet];
    dir[1] = inj->beamlet_dy[i_beamlet];
    dir[2] = inj->beamlet_dz[i_beamlet];

    real phi   = atan2(dir[1], dir[0]);
    real theta = acos(dir[2]);

    normalv[0] = sin(theta+CONST_PI/2) * cos(phi);
    normalv[1] = sin(theta+CONST_PI/2) * sin(phi);
    normalv[2] = cos(theta+CONST_PI/2);

    math_cross(dir, normalv, tmp);
    math_unit(tmp, normalh);

    /* Generate a random number on interval [0,1]. If the number is smaller
     * than the halo fraction, this marker will be considered as part
     * of the halo. The divergences are different for the core and halo.
     * Additionally, two normally distributed random variables are generated:
     * the divergences in horizontal and the vertical directions. */
    real div_h, div_v;
    if(random_uniform(rng) < inj->div_halo_frac) {
        // Use halo divergences instead
        div_h = inj->div_halo_h * random_normal(rng) / sqrt(2.0);
        div_v = inj->div_halo_v * random_normal(rng) / sqrt(2.0);
    } else {
        div_h = inj->div_h * random_normal(rng) / sqrt(2.0);
        div_v = inj->div_v * random_normal(rng) / sqrt(2.0);
    }

    /* Convert the divergence angle to an unit vector. The marker velocity
     * vector points in the direction of v_beam + v_divergence. */
    tmp[0] = cos(div_h) * ( cos(div_v) * dir[0] + sin(div_v) * normalv[0] )
        + sin(div_h) * normalh[0];
    tmp[1] = cos(div_h) * ( cos(div_v) * dir[1] + sin(div_v) * normalv[1] )
        + sin(div_h) * normalh[1];
    tmp[2] = cos(div_h) * ( cos(div_v) * dir[2] + sin(div_v) * normalv[2] )
        + sin(div_h) * normalh[2];
    math_unit(tmp, dir);

    real gamma = physlib_gamma_Ekin(inj->mass, energy);
    real absv  = sqrt( 1.0 - 1.0 / (gamma * gamma) ) * CONST_C;
    vxyz[0] = absv * dir[0];
    vxyz[1] = absv * dir[1];
    vxyz[2] = absv * dir[2];
}
