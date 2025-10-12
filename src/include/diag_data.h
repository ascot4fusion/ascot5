/**
 * @file diag_data.h
 * Diagnostics data.
 */
#ifndef DIAG_DATA_H
#define DIAG_DATA_H

#include "defines.h"
#include <stddef.h>

/**
 * @brief Orbit diagnostics data struct.
 *
 * The pointers are asssigned to the offload array but only those pointers
 * are used which corresponds to the simulation mode (e.g. GC simulation).
 *
 * The marker orbit data are stored in maxpoints length chunks in the
 * offload array. The chuncks are in no particular order. Once chunk is
 * filled and the marker is still recording, new points replace the old
 * ones from the start.
 */
typedef struct
{
    real *id;      /**< Marker ID                                            */
    real *mileage; /**< Time marker has been simulated for [s]               */
    real *r;       /**< Marker R coordinate [m]                              */
    real *phi;     /**< Marker phi coordinate [rad]                          */
    real *z;       /**< Marker z coordiante [m]                              */
    real *p_r;     /**< Particle momentum R component [kg m/s]               */
    real *p_phi;   /**< Particle momentum phi component [kg m/s]             */
    real *p_z;     /**< Particle momentum z component [kg m/s]               */
    real *ppar;    /**< Guiding center parallel momentum [kg m/s]            */
    real *mu;      /**< Guiding center magnetic moment [J/T]                 */
    real *zeta;    /**< Guiding center gyroangle [rad]                       */
    real *weight;  /**< Marker weight [1]                                    */
    real *charge;  /**< Marker charge [C]                                    */
    real *rho;     /**< Normalized poloidal flux at marker position [1]      */
    real *theta;   /**< Marker poloidal angle [rad]                          */
    real *B_r;     /**< Magnetic field R component at marker position [T]    */
    real *B_phi;   /**< Magnetic field phi component at marker position [T]  */
    real *B_z;     /**< Magnetic field z component at marker position [T]    */
    real *simmode; /**< In what simulation mode data point was recorded      */
    real *pncrid;  /**< Id for the poincare plot a point corresponds to      */
    real *pncrdi;  /**< Direction in which Poincare plane was crossed        */

    size_t *mrk_pnt;   /**< Index of the last recorded point             */
    real *mrk_recorded; /**< Last time (in seconds) a marker was updated  */

} diag_orb_data;

/**
 * @brief Simple linked list link for storing data points.
 */
typedef struct diag_transcoef_link
{
    real rho;      /**< Current radial coordinate (R or rho)                */
    real time;     /**< Current time                                        */
    int pitchsign; /**< Sign of the pitch at this point                     */
    struct diag_transcoef_link *prevlink; /**< Pointer to the previous link */
} diag_transcoef_link;

/**
 * @brief Transport coefficient diagnostics offload data struct.
 */
typedef struct
{
    diag_transcoef_link **datapoints; /**< Temporary data storage             */

    int *id;     /**< Marker ID whose data is stored at this index            */
    real *Kcoef; /**< Calculated drift coefficients                           */
    real *Dcoef; /**< Calculated diffusion coefficients where negative value
                    means coefficients are/were not calculated                */

} diag_transcoef_data;

/**
 * @brief Histogram parameters
 */
typedef struct
{
    size_t step_1; /**< step for 2nd fastest running index   */
    size_t step_2; /**< step for 3rd fastest running index   */
    size_t step_3; /**< step for 4th fastest running index   */
    size_t step_4; /**< step for 5th fastest running index   */
    size_t step_5; /**< step for 6th fastest running index   */
    size_t step_6; /**< step for 7th fastest running index   */

    real *histogram; /**< pointer to start of histogram array */
} dist_5D_data;

/**
 * @brief Histogram parameters on target
 */
typedef struct
{
    int n_r;    /**< number of r bins           */
    real min_r; /**< value of lowest r bin      */
    real max_r; /**< value of highest r bin     */

    int n_phi;    /**< number of r bins           */
    real min_phi; /**< value of lowest r bin      */
    real max_phi; /**< value of highest r bin     */

    int n_z;    /**< number of z bins           */
    real min_z; /**< value of lowest z bin      */
    real max_z; /**< value of highest z bin     */

    int n_pr;    /**< number of p_r bins         */
    real min_pr; /**< value of lowest p_r bin    */
    real max_pr; /**< value of highest p_r bin   */

    int n_pphi;    /**< number of p_phi bins       */
    real min_pphi; /**< value of lowest p_phi bin  */
    real max_pphi; /**< value of highest p_phi bin */

    int n_pz;    /**< number of p_z bins         */
    real min_pz; /**< value of lowest p_z bin    */
    real max_pz; /**< value of highest p_z bin   */

    int n_time;    /**< number of r bins           */
    real min_time; /**< value of lowest r bin      */
    real max_time; /**< value of highest r bin     */

    int n_q;    /**< number of r bins           */
    real min_q; /**< value of lowest r bin      */
    real max_q; /**< value of highest r bin     */

    size_t step_1; /**< step for 2nd fastest running index   */
    size_t step_2; /**< step for 3rd fastest running index   */
    size_t step_3; /**< step for 4th fastest running index   */
    size_t step_4; /**< step for 5th fastest running index   */
    size_t step_5; /**< step for 6th fastest running index   */
    size_t step_6; /**< step for 7th fastest running index   */
    size_t step_7; /**< step for 8th fastest running index   */

    real *histogram; /**< pointer to start of histogram array */
} dist_6D_data;

/**
 * @brief Histogram parameters on target
 */
typedef struct
{
    int n_mu;    /**< number of mu bins                     */
    real min_mu; /**< value of lowest mu bin                */
    real max_mu; /**< value of highest ,u bin               */

    int n_Ekin;    /**< number of Ekin bins                   */
    real min_Ekin; /**< value of lowest Ekin bin              */
    real max_Ekin; /**< value of highest Ekin bin             */

    int n_Ptor;    /**< number of Ptor bins                   */
    real min_Ptor; /**< value of lowest Ptor bin              */
    real max_Ptor; /**< value of highest Ptor bin             */

    size_t step_1; /**< step for 2nd fastest running index    */
    size_t step_2; /**< step for 3rd fastest running index    */

    real *histogram; /**< pointer to start of histogram array */
} dist_COM_data;

/**
 * @brief Histogram parameters
 */
typedef struct
{
    int n_rho;    /**< number of rho bins                   */
    real min_rho; /**< value of lowest rho bin              */
    real max_rho; /**< value of highest rho bin             */

    int n_theta;    /**< number of poloidal angle bins        */
    real min_theta; /**< value of lowest pol bin              */
    real max_theta; /**< value of highest pol bin             */

    int n_phi;    /**< number of phi bins                   */
    real min_phi; /**< value of lowest phi bin              */
    real max_phi; /**< value of highest phi bin             */

    int n_ppara;    /**< number of p_parallel bins            */
    real min_ppara; /**< value of lowest p_parallel bin       */
    real max_ppara; /**< value of highest p_parallel bin      */

    int n_pperp;    /**< number of p_perpendicular bins       */
    real min_pperp; /**< value of lowest p_perpendicular bin  */
    real max_pperp; /**< value of highest p_perpendicular bin */

    int n_time;    /**< number of time bins                  */
    real min_time; /**< value of lowest time bin             */
    real max_time; /**< value of highest time bin            */

    int n_q;    /**< number of charge bins                */
    real min_q; /**< value of lowest charge bin           */
    real max_q; /**< value of highest charge bin          */

    size_t step_1; /**< step for 2nd fastest running index   */
    size_t step_2; /**< step for 3rd fastest running index   */
    size_t step_3; /**< step for 4th fastest running index   */
    size_t step_4; /**< step for 5th fastest running index   */
    size_t step_5; /**< step for 6th fastest running index   */
    size_t step_6; /**< step for 7th fastest running index   */

    real *histogram; /**< pointer to start of histogram array */
} dist_rho5D_data;

/**
 * @brief Histogram parameters on target
 */
typedef struct
{
    int n_rho;    /**< number of rho bins            */
    real min_rho; /**< value of lowest rho bin       */
    real max_rho; /**< value of highest rho bin      */

    int n_theta;    /**< number of poloidal angle bins */
    real min_theta; /**< value of lowest theta bin     */
    real max_theta; /**< value of highest theta bin    */

    int n_phi;    /**< number of phi bins            */
    real min_phi; /**< value of lowest phi bin       */
    real max_phi; /**< value of highest phi bin      */

    int n_pr;    /**< number of p_r bins            */
    real min_pr; /**< value of lowest p_r bin       */
    real max_pr; /**< value of highest p_r bin      */

    int n_pphi;    /**< number of p_phi bins          */
    real min_pphi; /**< value of lowest p_phi bin     */
    real max_pphi; /**< value of highest p_phi bin    */

    int n_pz;    /**< number of p_z bins            */
    real min_pz; /**< value of lowest p_z bin       */
    real max_pz; /**< value of highest p_z bin      */

    int n_time;    /**< number of time bins           */
    real min_time; /**< value of lowest time bin      */
    real max_time; /**< value of highest time bin     */

    int n_q;    /**< number of charge bins         */
    real min_q; /**< value of lowest charge bin    */
    real max_q; /**< value of highest charge bin   */

    size_t step_1; /**< step for 2nd fastest running index   */
    size_t step_2; /**< step for 3rd fastest running index   */
    size_t step_3; /**< step for 4th fastest running index   */
    size_t step_4; /**< step for 5th fastest running index   */
    size_t step_5; /**< step for 6th fastest running index   */
    size_t step_6; /**< step for 7th fastest running index   */
    size_t step_7; /**< step for 8th fastest running index   */

    real *histogram; /**< pointer to start of histogram array */
} dist_rho6D_data;

#define HIST_ALLDIM 16

/**
 * @brief Quantities that can be used as histogram axis coordinates.
 */
typedef enum
{
    R,      /**< The R coordinate in cylindrical basis [m].                   */
    PHI,    /**< The phi coordinate in cylindrical basis [rad].               */
    Z,      /**< The z coordinate in cylindrical basis [m].                   */
    RHO,    /**< Square root of normalized poloidal flux [1].                 */
    THETA,  /**< Poloidal angle [rad].                                        */
    PPAR,   /**< Momentum component parallel to the magnetic field [kg*m/s].  */
    PPERP,  /**< Momentum component orthogonal to the magnetic field [kg*m/s].*/
    PR,     /**< Momentum R-component [kg*m/s].                               */
    PPHI,   /**< Momentum phi-component [kg*m/s].                             */
    PZ,     /**< Momentum z-component [kg*m/s].                               */
    EKIN,   /**< Kinetic energy [J].                                          */
    XI,     /**< Pitch [1].                                                   */
    MU,     /**< Magnetic moment [J/T].                                       */
    PTOR,   /**< Canonical toroidal angular momentum [kg*m/s].                */
    TIME,   /**< Time instant (laboratory time) [s].                          */
    CHARGE, /**< Charge state [e].                                            */
} hist_coordinate;

/**
 * @brief Coordinate axis for the histogram.
 */
typedef struct
{
    hist_coordinate name; /**< Coordinate mapped to this axis                 */
    real min;             /**< Lower limit of the coordinate interval         */
    real max;             /**< Upper limit of the coordinate interval         */
    size_t n;             /**< Number of bins in this axis                    */
} hist_axis;

/**
 * @brief Histogram parameters.
 *
 * The bins are stored as a flattened array. The element (i0,i1,...,in)
 * corresponds to bins[i0*strides[0] + i1*strides[1] + ... + in].
 */
typedef struct
{
    hist_axis axes[HIST_ALLDIM]; /**< The coordinate axes.                  */
    size_t strides[HIST_ALLDIM - 1]; /**< Row length for each dimension. */
    size_t nbin; /**< Number of bins.                       */
    real *bins;  /**< The bin array.                        */
} histogram;

typedef struct
{
    diag_orb_data *orbit;
    diag_transcoef_data *transport_coefficient;
    dist_5D_data *dist5d;
    dist_6D_data *dist6d;
    dist_rho5D_data *dist5drho;
    dist_rho6D_data *dist6drho;
    dist_COM_data *distcom;
} Diagnostics;

#endif