/**
 * @file nbi_data.h
 * Defines neutral beam injector data.
 */
#ifndef NBI_DATA_H
#define NBI_DATA_H

#include <stddef.h>

/**
 * @brief Structure for describing an NBI injector
 */
typedef struct
{
    size_t id;            /**< Unique identifier for this injector.           */
    size_t n_beamlet;     /**< Number of beamlets in this injector.           */
    double *beamlet_x;    /**< Beamlet x coordinates [m].                     */
    double *beamlet_y;    /**< Beamlet y coordinates [m].                     */
    double *beamlet_z;    /**< Beamlet z coordinates [m].                     */
    double *beamlet_dx;   /**< Beamlet direction (unit) vector x component.   */
    double *beamlet_dy;   /**< Beamlet direction (unit) vector y component.   */
    double *beamlet_dz;   /**< Beamlet direction (unit) vector z component.   */
    double power;         /**< Injection power [W].                           */
    double energy;        /**< Full energy of injected particles [J].         */
    double efrac[3];      /**< Fractions of full, 1/2 and 1/3 energy.         */
    double div_h;         /**< Vertical divergence [rad].                     */
    double div_v;         /**< Horizontal divergence [rad].                   */
    double div_halo_frac; /**< Fraction of power in the halo.                 */
    double div_halo_h;    /**< Horizontal divergence of the halo [rad].       */
    double div_halo_v;    /**< Vertical divergence of the halo [rad].         */
    int anum;             /**< Mass number of injected species.               */
    int znum;             /**< Charge number of injected species.             */
    double mass;          /**< Mass of injected species [kg].                 */
} Nbi;

#endif
