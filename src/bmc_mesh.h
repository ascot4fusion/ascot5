/**
 * @file bmc_mesh.h
 * @brief Header file for bmc_mesh.c
 *
 * Along with the function declarations, this header contains the definition
 * of the mesh struct.
 */
#ifndef BMC_MESH_H
#define BMC_MESH_H

#include <stdlib.h>
#include "ascot5.h"

/**
 * @brief BMC mesh in 5D phase space
 */
typedef struct {
    int n_r;        /**< Number of R grid points                       */
    int n_phi;      /**< Number of phi grid points                     */
    int n_z;        /**< Number of z grid points                       */
    int n_ppara;    /**< Number of parallel momentum grid points       */
    int n_pperp;    /**< Number of perpendicular momentum grid points  */
    real* r;        /**< Mesh R-abscissa [m]                           */
    real* phi;      /**< Mesh phi-abscissa [rad]                       */
    real* z;        /**< Mesh z-abscissa [m]                           */
    real* ppara;    /**< Mesh parallel momentum abscissa [kg*m/s]      */
    real* pperp;    /**< Mesh perpendicular momentum abscissa [kg*m/s] */
    real* val_prev; /**< Probability evaluated in the previous step    */
    real* val_next; /**< Probability currently being evaluated         */
    size_t size;    /**< Size of the mesh (i.e. the two arrays above)  */
} bmc_mesh;

int bmc_mesh_init();
void bmc_mesh_free(bmc_mesh* mesh);
#pragma omp declare simd uniform(mesh)
void bmc_mesh_index2pos(bmc_mesh* mesh, size_t idx, real coords[5]);
real bmc_mesh_interpolate(bmc_mesh* mesh, real r, real phi, real z, real ppara,
                          real pperp);
void bmc_mesh_finishstep(bmc_mesh* mesh);
void bmc_mesh_update(bmc_mesh* mesh, size_t start, size_t stop,
                     real* r, real* phi, real* z, real* ppara, real* pperp,
                     int* fate);

#endif