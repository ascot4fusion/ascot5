/**
 * Implements "nbi_source.h".
 */
#include "nbi_source.h"
#include "ascot.h"
#include "consts.h"
#include "data/atomic.h"
#include "data/bfield.h"
#include "data/diag.h"
#include "data/marker.h"
#include "data/neutral.h"
#include "data/plasma.h"
#include "data/wall.h"
#include "datatypes.h"
#include "defines.h"
#include "endcond.h"
#include "parallel.h"
#include "simulate.h"
#include "utils/mathlib.h"
#include "utils/physlib.h"
#include "utils/random.h"
#include <getopt.h>
#include <math.h>
#include <omp.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

void nbi_source_sample_beam(real xyz[3], real vxyz[3], real rnd[5], Nbi *nbi)
{

    size_t random_beamlet = floor(rnd[0] * nbi->n_beamlet);
    xyz[0] = nbi->beamlet_x[random_beamlet];
    xyz[1] = nbi->beamlet_y[random_beamlet];
    xyz[2] = nbi->beamlet_z[random_beamlet];

    real ekin;
    if (rnd[1] < nbi->efrac[0])
        ekin = nbi->energy;
    else if (rnd[1] < nbi->efrac[0] + nbi->efrac[1])
        ekin = nbi->energy / 2;
    else
        ekin = nbi->energy / 3;
    real gamma = physlib_gamma_Ekin(nbi->mass, ekin);
    real vnorm = sqrt(1.0 - 1.0 / (gamma * gamma)) * CONST_C;

    real dir[3] = {
        nbi->beamlet_dx[random_beamlet],
        nbi->beamlet_dy[random_beamlet],
        nbi->beamlet_dz[random_beamlet],
    };

    real up[3] = {0.0, 0.0, 1.0};

    real n_h[3];
    math_cross(up, dir, n_h);
    math_unit(n_h, n_h);

    real n_v[3];
    math_cross(dir, n_h, n_v);
    math_unit(n_v, n_v);

    real div[2];
    div[0] = (rnd[2] < nbi->halo_fraction ? nbi->div_halo[0] : nbi->div[0]) *
             rnd[3] / sqrt(2.);
    div[1] = (rnd[2] < nbi->halo_fraction ? nbi->div_halo[1] : nbi->div[1]) *
             rnd[4] / sqrt(2.);

    real tmp[3] = {
        dir[0] + div[0] * n_v[0] + div[1] * n_h[0],
        dir[1] + div[0] * n_v[1] + div[1] * n_h[1],
        dir[2] + div[0] * n_v[2] + div[1] * n_h[2],
    };
    math_unit(tmp, dir);

    vxyz[0] = vnorm * dir[0];
    vxyz[1] = vnorm * dir[1];
    vxyz[2] = vnorm * dir[2];
}

void nbi_source_inject_markers(
    Simulation *sim, Nbi *nbi, real time_step, real max_dist, real tlim[2],
    size_t nmrk, State mrk[nmrk])
{

    /* Set marker weights assuming a large number is created so that the energy
     * fractions of generated markers are close to the injector values */
    real f = 1.0 * nbi->efrac[0] + (1.0 / 2) * nbi->efrac[1] +
             (1.0 / 3) * nbi->efrac[2];
    real weight = (nbi->power / nbi->energy) / (f * nmrk);

    CPU_PARALLELIZE_LOOP
    for (size_t i = 0; i < nmrk; i++)
    {
        real rnd[6];
        size_t rng_count = 0;
        random_uniform(rnd, 4, &i, &rng_count, 0);
        random_norm(&rnd[3], 2, &i, &rng_count, 0);
        real time = tlim[0] + rnd[0] * (tlim[1] - tlim[0]);

        real xyz[3], vxyz[3], rpz[3], vhat[3];
        nbi_source_sample_beam(xyz, vxyz, &rnd[1], nbi);
        math_xyz2rpz(xyz, rpz);
        math_unit(vxyz, vhat);

        real psi;
        err_t err = 0;
        err_t outside_field =
            Bfield_eval_psi(&psi, rpz[0], rpz[1], rpz[2], time, &sim->bfield);
        while (outside_field)
        {
            xyz[0] += time_step * vhat[0];
            xyz[1] += time_step * vhat[1];
            xyz[2] += time_step * vhat[2];
            math_xyz2rpz(xyz, rpz);
            outside_field = Bfield_eval_psi(
                &psi, rpz[0], rpz[1], rpz[2], time, &sim->bfield);
            if (sqrt(xyz[0] * xyz[0] + xyz[1] * xyz[1] + xyz[2] * xyz[2]) >
                max_dist)
            {
                err = ERROR_RAISE(ERR_NBI_MAX_DIST, SIMULATE_NBI_SOURCE_C);
                break;
            }
        }

        real vrpz[3];
        math_vec_xyz2rpz(vxyz, vrpz, rpz[1]);
        real gamma = physlib_gamma_vnorm(math_norm(vrpz));

        mrk[i].err = err;
        mrk[i].rprt = rpz[0];
        mrk[i].zprt = rpz[2];
        mrk[i].phiprt = rpz[1];
        mrk[i].pr = vrpz[0] * gamma * nbi->mass;
        mrk[i].pz = vrpz[2] * gamma * nbi->mass;
        mrk[i].pphi = vrpz[1] * gamma * nbi->mass;
        mrk[i].mass = nbi->mass;
        mrk[i].anum = nbi->anum;
        mrk[i].znum = nbi->znum;
        mrk[i].time = time;
        mrk[i].weight = weight;
        mrk[i].charge = 0;
        mrk[i].mileage = 0.;
        mrk[i].cputime = 0.;
        mrk[i].endcond = 0;
        mrk[i].walltile = 0;
    }
}
