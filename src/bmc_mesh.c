/**
 * @file bmc_mesh.c
 * @brief Routines for initializing and processing BMC mesh
 */
#include <stdlib.h>
#include <string.h>
#include "ascot5.h"
#include "mpi_interface.h"
#include "simulate/simulate_bmc.h"
#include "bmc_mesh.h"

/**
 * @brief Initialize BMC mesh
 *
 * Right now this function uses cartesian grid to initialize the mesh.
 *
 * @param rmin minimum  value
 * @param rmax maximum  value
 * @param nr number of radial bins
 * @param phimin minimum  value
 * @param phimax maximum  value
 * @param nphi number of toroidal bins
 * @param zmin minimum  value
 * @param zmax maximum  value
 * @param nz number of z bins
 * @param pparamin minimum  value
 * @param pparamax maximum  value
 * @param nppara number of parallel momentum bins
 * @param pperpmin minimum  value
 * @param pperpmax maximum  value
 * @param npperp number of perpendicular momentum bins
 * @param mesh pointer to the mesh object to be initialized
 *
 * @return zero if initialization succeeded
 */
int bmc_mesh_init(real min_r,     real max_r,     int n_r,
                  real min_phi,   real max_phi,   int n_phi,
                  real min_z,     real max_z,     int n_z,
                  real min_ppara, real max_ppara, int n_ppara,
                  real min_pperp, real max_pperp, int n_pperp,
                  bmc_mesh* mesh) {
    mesh->n_r     = n_r + 1;
    mesh->n_phi   = n_phi;
    mesh->n_z     = n_z + 1;
    mesh->n_ppara = n_ppara + 1;
    mesh->n_pperp = n_pperp + 1;
    mesh->size = ((size_t) mesh->n_r) * ((size_t) mesh->n_phi) * ((size_t) mesh->n_z)
               * ((size_t) mesh->n_ppara) * ((size_t) mesh->n_pperp);
    mesh->r        = (real*) malloc(mesh->n_r     * sizeof(real));
    mesh->phi      = (real*) malloc(mesh->n_phi   * sizeof(real));
    mesh->z        = (real*) malloc(mesh->n_z     * sizeof(real));
    mesh->ppara    = (real*) malloc(mesh->n_ppara * sizeof(real));
    mesh->pperp    = (real*) malloc(mesh->n_pperp * sizeof(real));
    mesh->val_next = (real*) malloc(mesh->size    * sizeof(real));
    mesh->val_prev = (real*) malloc(mesh->size    * sizeof(real));

    /* Store abscissae */
    for(int i=0; i<=n_r; i++)
        mesh->r[i] = min_r + i * (max_r - min_r) / n_r ;
    for(int i=0; i<=n_z; i++)
        mesh->z[i] = min_z + i * (max_z - min_z) / n_z;
    for(int i=0; i<n_phi; i++)
        mesh->phi[i] = min_phi + i * (max_phi - min_phi) / (n_phi+1);
    for(int i=0; i<=n_ppara; i++)
        mesh->ppara[i] = min_ppara + i * (max_ppara - min_ppara) / n_ppara;
    for(int i=0; i<=n_pperp; i++)
        mesh->pperp[i] = min_pperp + i * (max_pperp - min_pperp) / n_pperp;

    /* Set initial probability to zero */
    memset(mesh->val_prev, 0.0, mesh->size * sizeof(real));
    memset(mesh->val_next, 0.0, mesh->size * sizeof(real));

    return 0;
}

/**
 * @brief Free resources used by a mesh object
 *
 * @param mesh mesh object to be deallocated
 */
void bmc_mesh_free(bmc_mesh* mesh) {
    free(mesh->r);
    free(mesh->phi);
    free(mesh->z);
    free(mesh->ppara);
    free(mesh->pperp);
    free(mesh->val_next);
    free(mesh->val_prev);

    mesh->size = 0;
}

/**
 * @brief Return phase-space coordinates corresponding to given mesh element
 *
 * @param mesh pointer to the mesh object
 * @param idx index of the mesh element
 * @param coords array where coordinates (R,phi,z,ppara,pperp) will be stored
 */
void bmc_mesh_index2pos(bmc_mesh* mesh, size_t idx, real coords[5]) {
    int pperp = round( idx /
        ( mesh->n_r * mesh->n_z * mesh->n_phi * mesh->n_ppara ) );
    idx -= pperp * mesh->n_r * mesh->n_z * mesh->n_phi * mesh->n_ppara;
    int ppara = round( idx / ( mesh->n_r * mesh->n_z * mesh->n_phi ) );
    idx -= ppara * mesh->n_r * mesh->n_z * mesh->n_phi;
    int phi   = round( idx / ( mesh->n_r * mesh->n_z ) );
    idx -= phi * mesh->n_r * mesh->n_z;
    int z     = round( idx / mesh->n_r );
    idx -=     z * mesh->n_r;
    int r     = idx;

    coords[0] = mesh->r[r];
    coords[1] = mesh->phi[phi];
    coords[2] = mesh->z[z];
    coords[3] = mesh->ppara[ppara];
    coords[4] = mesh->pperp[pperp];
}

real bmc_mesh_interpolate(bmc_mesh* mesh, real r, real phi, real z, real ppara,
                          real pperp) {
    /* Find the matrix indices where this pseudo-particle belongs to */
    size_t i_r = ((r - mesh->r[0])
        / (mesh->r[1] - mesh->r[0]));
    size_t i_z = ((z - mesh->z[0])
        / (mesh->z[1] - mesh->z[0]));
    size_t i_ppara = ((ppara - mesh->ppara[0])
    / (mesh->ppara[1] - mesh->ppara[0]));
    size_t i_pperp = ((pperp - mesh->pperp[0])
        / (mesh->pperp[1] - mesh->pperp[0]));

    /* Periodic variable */
    size_t i_phi, i_phi1;
    if(mesh->n_phi == 1) {
        i_phi = 0;
    }
    else {
        i_phi = ((phi - mesh->phi[0]) / (mesh->phi[1] - mesh->phi[0]));
    }
    i_phi1 = i_phi + 1;
    if(i_phi == mesh->n_phi-1) i_phi1 = 0;

    /* Zero outside */
    real val = 0;
    if(i_ppara >= 0 && i_ppara < mesh->n_ppara-2 &&
       i_pperp >= 0 && i_pperp < mesh->n_pperp-2 &&
       i_r >= 0 && i_r < mesh->n_r-2 && i_z >= 0 && i_z < mesh->n_z-2) {

        /* Interpolate the probability value at the particle location
        *  using the nearby nodes (for linear interpolation) */
        real dr[2]   = {mesh->r[i_r+1] - r, r - mesh->r[i_r]};
        real dphi[2] = {mesh->phi[i_phi1] - phi, phi - mesh->phi[i_phi]};
        real dz[2]   = {mesh->z[i_z+1] - z, z - mesh->z[i_z]};
        real dppa[2] = {mesh->ppara[i_ppara+1] - ppara,
                        ppara - mesh->ppara[i_ppara]};
        real dppe[2] = {mesh->pperp[i_pperp+1] - pperp,
                        pperp - mesh->pperp[i_pperp]};
        real vol = (mesh->r[i_r+1]         - mesh->r[i_r])
                * fmax(fabs(mesh->phi[i_phi1] - mesh->phi[i_phi]), 1.0)
                * (mesh->z[i_z+1]         - mesh->z[i_z])
                * (mesh->ppara[i_ppara+1] - mesh->ppara[i_ppara])
                * (mesh->pperp[i_pperp+1] - mesh->pperp[i_pperp]);

        if(mesh->n_phi == 1) {
            dphi[0] = 0.5;
            dphi[1] = 0.5;
        }

        /* 5D linear interpolation */
        for(int i1=0; i1<2; i1++)
        for(int i2=0; i2<2; i2++)
        for(int i3=0; i3<2; i3++)
        for(int i4=0; i4<2; i4++)
        for(int i5=0; i5<2; i5++) {
            size_t idx0 =
                  (i_pperp + i5) * (mesh->n_r) * (mesh->n_z) * (mesh->n_phi)
                                 * (mesh->n_ppara)
                + (i_ppara + i4) * (mesh->n_r) * (mesh->n_z) * (mesh->n_phi)
                + (i_phi   + i3 * (i_phi1 - i_phi) ) * (mesh->n_r) * (mesh->n_z)
                + (i_z     + i2) * (mesh->n_r)
                + i_r      + i1;
            val += mesh->val_prev[idx0] * dr[i1] * dz[i2] * dphi[i3]
                * dppa[i4] * dppe[i5];
        }
        val /= vol;
    }
    return val;
}

/**
 * @brief Prepare mesh for the next time step
 *
 * This function should be called after the calculation for the current step
 * has finished.
 *
 * @param mesh pointer to the mesh data
 */
void bmc_mesh_finishstep(bmc_mesh* mesh) {
#ifdef MPI
    MPI_Allreduce(mesh->val_next, mesh->val_prev, mesh->size,
                  MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#else
    memcpy(mesh->val_prev, mesh->val_next, mesh->size * sizeof(real));
#endif
    memset(mesh->val_next, 0.0, mesh->size * sizeof(real));
}

/**
 * @brief Take a time step and update the probability
 *
 * For MPI parallelization it is possible to only update mesh elements
 * [start,stop) assuming that the push-matrix values have been calculated for
 * this same interval.
 *
 * The push matrix is assumed to have a format [imesh*HERMITE_KNOT + iknot].
 *
 * @param mesh pointer to the mesh data
 * @param start the first mesh index
 * @param stop the final mesh index
 * @param r push-matrix marker final R-coordinates [m]
 * @param phi push-matrix marker final (periodic, not cumulative)
 *        phi-coordinates [rad]
 * @param z push-matrix marker final z-coordinates [m]
 * @param ppara push-matrix marker final parallel momentum coordinate [kg*m/s]
 * @param pperp push-matrix marker final perpendicular momentum coordinate
 *        [kg*m/s]
 * @param fate flag indicating whether the push-matrix marker terminated with
 *        error [-1], hit wall [1], hit FILD [2], or finished normally [0]
 */
void bmc_mesh_update(bmc_mesh* mesh, size_t start, size_t stop,
                     real* r, real* phi, real* z, real* ppara, real* pperp,
                     int* fate) {

    real hermite_w[HERMITE_KNOTS] = HERMITE_W;
    #pragma omp parallel for \
        shared(start, stop, mesh, r, phi, z, ppara, pperp, fate)
    for(size_t iprt=start; iprt < stop; iprt++) {

        size_t idx = (iprt-start)*HERMITE_KNOTS;
        for(int i_knot=0; i_knot<HERMITE_KNOTS; i_knot++) {

            real ri = r[idx + i_knot];
            real zi = z[idx + i_knot];
            real phii = phi[idx + i_knot];
            real pparai = ppara[idx + i_knot];
            real pperpi = pperp[idx + i_knot];

            real val = 0;
            if(fate[idx + i_knot] == 2) {
                val = 1.0;
            }
            else if (fabs(fate[idx + i_knot]) == 1) {
                val = 0;
            }
            else {
                val = bmc_mesh_interpolate(mesh, ri, phii, zi, pparai, pperpi);
            }
            mesh->val_next[iprt] += val * hermite_w[i_knot];
        }
    }
}