/**
 * @file interp.h
 * @brief Spline interpolation library
 *
 * Spline interpolation fits cubic splines on data given on a uniform grid. Each
 * axis may have either natural or periodic boundary condition.
 *
 * There exists two representations for the splines: compact and explicit. Both
 * give identical results but the difference is that compact requires fewer
 * coefficients to be stored (and fetched from the memory at each evaluation)
 * than the explicit. On the other hand, explicit requires less computations
 * per evaluations, but compact is usually faster, and conserves memory,
 * especially in 3D. Therefore compact splines are preferred.
 *
 * To initialize splines, first call corresponding init_coeff() function which
 * evaluates coefficients to a pre-allocated array (i.e. to the offload array).
 * Then (after offloading is done) call init_spline() which assigns the
 * coefficients to a spline struct.
 *
 * In order to allocate the array for storing the coefficients, one needs to
 * know how many coefficients are stored per data grid point:
 *
 * - 1D compact  2, explicit 4
 * - 2D compact  4, explicit 16
 * - 3D compact  8, explicit 64
 */
#ifndef INTERP_H
#define INTERP_H
#include "../ascot5.h"
#include "../offload.h"
#include "../error.h"

/**
 * @brief Boundary conditions for the spline interpolation.
 */
enum boundaryCondition {
    NATURALBC  = 0, /**< Second derivative is zero at both ends               */
    PERIODICBC = 1  /**< Function has same value and derivatives at both ends */
};

/**
 * @brief Number of coefficients stored for each data point.
 */
enum splinesize {
    NSIZE_COMP1D =  2,
    NSIZE_COMP2D =  4,
    NSIZE_COMP3D =  8,
    NSIZE_EXPL1D =  4,
    NSIZE_EXPL2D = 16,
    NSIZE_EXPL3D = 64
};

/**
 * @brief Cubic interpolation struct.
 */
typedef struct {
    int n_x;     /**< number of x grid points                        */
    int bc_x;    /**< boundary condition for x coordinate            */
    real x_min;  /**< minimum x coordinate in the grid               */
    real x_max;  /**< maximum x coordinate in the grid               */
    real x_grid; /**< interval between two adjacent points in x grid */
    real* c;     /**< pointer to array with spline coefficients      */
} interp1D_data;

/**
 * @brief Bicubic interpolation struct.
 */
typedef struct {
    int n_x;     /**< number of x grid points                        */
    int n_y;     /**< number of y grid points                        */
    int bc_x;    /**< boundary condition for x coordinate            */
    int bc_y;    /**< boundary condition for y coordinate            */
    real x_min;  /**< minimum x coordinate in the grid               */
    real x_max;  /**< maximum x coordinate in the grid               */
    real x_grid; /**< interval between two adjacent points in x grid */
    real y_min;  /**< minimum y coordinate in the grid               */
    real y_max;  /**< maximum y coordinate in the grid               */
    real y_grid; /**< interval between two adjacent points in y grid */
    real* c;     /**< pointer to array with spline coefficients      */
} interp2D_data;

/**
 * @brief Tricubic interpolation struct.
 */
typedef struct {
    int n_x;     /**< number of x grid points                        */
    int n_y;     /**< number of y grid points                        */
    int n_z;     /**< number of z grid points                        */
    int bc_x;    /**< boundary condition for x coordinate            */
    int bc_y;    /**< boundary condition for y coordinate            */
    int bc_z;    /**< boundary condition for z coordinate            */
    real x_min;  /**< minimum x coordinate in the grid               */
    real x_max;  /**< maximum x coordinate in the grid               */
    real x_grid; /**< interval between two adjacent points in x grid */
    real y_min;  /**< minimum y coordinate in the grid               */
    real y_max;  /**< maximum y coordinate in the grid               */
    real y_grid; /**< interval between two adjacent points in y grid */
    real z_min;  /**< minimum z coordinate in the grid               */
    real z_max;  /**< maximum z coordinate in the grid               */
    real z_grid; /**< interval between two adjacent points in z grid */
    real* c;     /**< pointer to array with spline coefficients      */
    real x_inv_grid, y_inv_grid, z_inv_grid; /* 1/grid */
    real xg2_over6, yg2_over6, zg2_over6;    /* grid^2 / 6 */
    real xg_over6, yg_over6, zg_over6;    /* grid^2 / 6 */
    int y_stride8; /* n_x*8 */
    int z_stride8; /* n_y*n_x*8 */
} interp3D_data;

int interp1Dcomp_init_coeff(real* c, real* f,
                            int n_x, int bc_x,
                            real x_min, real x_max);

int interp2Dcomp_init_coeff(real* c, real* f,
                            int n_x, int n_y, int bc_x, int bc_y,
                            real x_min, real x_max,
                            real y_min, real y_max);

int interp3Dcomp_init_coeff(real* c, real* f,
                            int n_x, int n_y, int n_z,
                            int bc_x, int bc_y, int bc_z,
                            real x_min, real x_max,
                            real y_min, real y_max,
                            real z_min, real z_max);

int interp1Dexpl_init_coeff(real* c, real* f,
                            int n_x, int bc_x,
                            real x_min, real x_max);

int interp2Dexpl_init_coeff(real* c, real* f,
                            int n_x, int n_y, int bc_x, int bc_y,
                            real x_min, real x_max,
                            real y_min, real y_max);

int interp3Dexpl_init_coeff(real* c, real* f,
                            int n_x, int n_y, int n_z,
                            int bc_x, int bc_y, int bc_z,
                            real x_min, real x_max,
                            real y_min, real y_max,
                            real z_min, real z_max);

void interp1Dcomp_init_spline(interp1D_data* str, real* c,
                              int n_x, int bc_x,
                              real x_min, real x_max);

void interp2Dcomp_init_spline(interp2D_data* str, real* c,
                              int n_x, int n_y, int bc_x, int bc_y,
                              real x_min, real x_max,
                              real y_min, real y_max);

void interp3Dcomp_init_spline(interp3D_data* str, real* c,
                              int n_x, int n_y, int n_z,
                              int bc_x, int bc_y, int bc_z,
                              real x_min, real x_max,
                              real y_min, real y_max,
                              real z_min, real z_max);

void interp1Dexpl_init_spline(interp1D_data* str, real* c,
                              int n_x, int bc_x,
                              real x_min, real x_max);

void interp2Dexpl_init_spline(interp2D_data* str, real* c,
                              int n_x, int n_y, int bc_x, int bc_y,
                              real x_min, real x_max,
                              real y_min, real y_max);

void interp3Dexpl_init_spline(interp3D_data* str, real* c,
                              int n_x, int n_y, int n_z,
                              int bc_x, int bc_y, int bc_z,
                              real x_min, real x_max,
                              real y_min, real y_max,
                              real z_min, real z_max);

int interp1Dcomp_setup(interp1D_data* str, real* f, int n_x, int bc_x,
                       real x_min, real x_max);

int interp2Dcomp_setup(interp2D_data* str, real* f, int n_x, int n_y,
                       int bc_x, int bc_y, real x_min, real x_max,
                       real y_min, real y_max);

int interp3Dcomp_setup(interp3D_data* str, real* f,
                       int n_x, int n_y, int n_z, int bc_x, int bc_y, int bc_z,
                       real x_min, real x_max, real y_min, real y_max,
                       real z_min, real z_max);

GPU_DECLARE_TARGET_SIMD_UNIFORM(str)
a5err interp1Dcomp_eval_f(real* f, interp1D_data* str, real x);
DECLARE_TARGET_END
GPU_DECLARE_TARGET_SIMD_UNIFORM(str)
a5err interp2Dcomp_eval_f(real* f, interp2D_data* str, real x, real y);
DECLARE_TARGET_END
GPU_DECLARE_TARGET_SIMD_UNIFORM(str)
a5err interp3Dcomp_eval_f(real* f, interp3D_data* str,
                         real x, real y, real z);
DECLARE_TARGET_END

DECLARE_TARGET_SIMD_UNIFORM(str)
a5err interp1Dexpl_eval_f(real* f, interp1D_data* str, real x);
DECLARE_TARGET_SIMD_UNIFORM(str)
a5err interp2Dexpl_eval_f(real* f, interp2D_data* str, real x, real y);
DECLARE_TARGET_SIMD_UNIFORM(str)
a5err interp3Dexpl_eval_f(real* f, interp3D_data* str,
                          real x, real y, real z);

GPU_DECLARE_TARGET_SIMD_UNIFORM(str)
a5err interp1Dcomp_eval_df(real* f_df, interp1D_data* str, real x);
DECLARE_TARGET_END
GPU_DECLARE_TARGET_SIMD_UNIFORM(str)
a5err interp2Dcomp_eval_df(real* f_df, interp2D_data* str, real x, real y);
DECLARE_TARGET_END
GPU_DECLARE_TARGET_SIMD_UNIFORM(str)
a5err interp3Dcomp_eval_df(real* f_df, interp3D_data* str,
                           real x, real y, real z);
DECLARE_TARGET_END

DECLARE_TARGET_SIMD_UNIFORM(str)
a5err interp1Dexpl_eval_df(real* f_df, interp1D_data* str, real x);
DECLARE_TARGET_SIMD_UNIFORM(str)
a5err interp2Dexpl_eval_df(real* f_df, interp2D_data* str, real x, real y);
DECLARE_TARGET_SIMD_UNIFORM(str)
a5err interp3Dexpl_eval_df(real* f_df, interp3D_data* str,
                           real x, real y, real z);

static inline real wrap_periodic(real v, real vmin, real vmax)
{
    const real L = vmax - vmin;
    real t = (v - vmin) / L;   /* peut être négatif */
    t = t - floor(t);          /* frac dans [0,1) */
    return vmin + t * L;
}


DECLARE_TARGET_SIMD_UNIFORM(str)
static inline int interp3Dquad_eval_f_and_grad(
    real * __restrict__ f,
    real * __restrict__ g, // g[0]=fx, g[1]=fy, g[2]=fz
    const interp3D_data * __restrict__ str,
    const real * __restrict__ finp,
    real x, real y, real z)
{
    // 1) BC (reprend ton style)
    if (str->bc_x == PERIODICBC) x = wrap_periodic(x, str->x_min, str->x_max);
    if (str->bc_y == PERIODICBC) y = wrap_periodic(y, str->y_min, str->y_max);
    if (str->bc_z == PERIODICBC) z = wrap_periodic(z, str->z_min, str->z_max);

    if (str->bc_x == NATURALBC && (x < str->x_min || x > str->x_max)) return 1;
    if (str->bc_y == NATURALBC && (y < str->y_min || y > str->y_max)) return 1;
    if (str->bc_z == NATURALBC && (z < str->z_min || z > str->z_max)) return 1;

    // 2) indices cellule + t
    const real fx = (x - str->x_min) * str->x_inv_grid;
    const real fy = (y - str->y_min) * str->y_inv_grid;
    const real fz = (z - str->z_min) * str->z_inv_grid;

    int ix = (int)fx;  // cellule [ix, ix+1]
    int iy = (int)fy;
    int iz = (int)fz;

    ix -= (x >= str->x_max);
    iy -= (y >= str->y_max);
    iz -= (z >= str->z_max);

    const real tx = fx - (real)ix; // in [0,1)
    const real ty = fy - (real)iy;
    const real tz = fz - (real)iz;

    // 3) choisir stencil 3 points: i-1, i, i+1 autour de ix
    // Donc on a besoin des indices gx = ix-1, ix, ix+1.
    int x0 = ix - 1, x1 = ix, x2 = ix + 1;
    int y0 = iy - 1, y1 = iy, y2 = iy + 1;
    int z0 = iz - 1, z1 = iz, z2 = iz + 1;

    // gestion bord (simple):
    // - périodique: wrap modulo n
    // - naturel: clamp aux bornes (dégrade en bord, mais OK si tu acceptes moins précis)
    // NOTE: str->n_x etc existent chez toi plus haut dans d'autres structs; ici ton interp3D_data les a.
    // Dans ton typedef actuel, n_x/n_y/n_z sont là -> OK.

    #define WRAP(i,n) ((i) < 0 ? (i)+(n) : ((i) >= (n) ? (i)-(n) : (i)))
    #define CLAMP(i,n) ((i) < 0 ? 0 : ((i) >= (n) ? (n)-1 : (i)))

    if (str->bc_x == PERIODICBC) { x0=WRAP(x0,str->n_x); x1=WRAP(x1,str->n_x); x2=WRAP(x2,str->n_x); }
    else                         { x0=CLAMP(x0,str->n_x); x1=CLAMP(x1,str->n_x); x2=CLAMP(x2,str->n_x); }

    if (str->bc_y == PERIODICBC) { y0=WRAP(y0,str->n_y); y1=WRAP(y1,str->n_y); y2=WRAP(y2,str->n_y); }
    else                         { y0=CLAMP(y0,str->n_y); y1=CLAMP(y1,str->n_y); y2=CLAMP(y2,str->n_y); }

    if (str->bc_z == PERIODICBC) { z0=WRAP(z0,str->n_z); z1=WRAP(z1,str->n_z); z2=WRAP(z2,str->n_z); }
    else                         { z0=CLAMP(z0,str->n_z); z1=CLAMP(z1,str->n_z); z2=CLAMP(z2,str->n_z); }

    // 4) poids quadratiques Lagrange sur [-1,0,1] avec t in [0,1]
    const real tx2 = tx*tx;
    const real wx0 = (real)0.5 * tx * (tx - (real)1.0);
    const real wx1 = (real)1.0 - tx2;
    const real wx2 = (real)0.5 * tx * (tx + (real)1.0);

    const real dwx0 = tx - (real)0.5;
    const real dwx1 = (real)-2.0 * tx;
    const real dwx2 = tx + (real)0.5;

    const real ty2 = ty*ty;
    const real wy0 = (real)0.5 * ty * (ty - (real)1.0);
    const real wy1 = (real)1.0 - ty2;
    const real wy2 = (real)0.5 * ty * (ty + (real)1.0);

    const real dwy0 = ty - (real)0.5;
    const real dwy1 = (real)-2.0 * ty;
    const real dwy2 = ty + (real)0.5;

    const real tz2 = tz*tz;
    const real wz0 = (real)0.5 * tz * (tz - (real)1.0);
    const real wz1 = (real)1.0 - tz2;
    const real wz2 = (real)0.5 * tz * (tz + (real)1.0);

    const real dwz0 = tz - (real)0.5;
    const real dwz1 = (real)-2.0 * tz;
    const real dwz2 = tz + (real)0.5;

    // 5) accumulate 27 points
    // layout finp[z*n_y*n_x + y*n_x + x]
    const int nx = str->n_x;
    const int nxy = str->n_x * str->n_y;

    // helper lambdas (en C pur: macro)
    #define FVAL(xx,yy,zz) finp[(zz)*nxy + (yy)*nx + (xx)]

    // on fait somme separable: d'abord combiner x, puis y, puis z
    // mais ici version simple (27 MAC)
    real val = 0.0;
    real dtx = 0.0;
    real dty_acc = 0.0;
    real dtz_acc = 0.0;

    const int xs[3] = {x0,x1,x2};
    const int ys[3] = {y0,y1,y2};
    const int zs[3] = {z0,z1,z2};
    const real wx[3] = {wx0,wx1,wx2};
    const real wy[3] = {wy0,wy1,wy2};
    const real wz[3] = {wz0,wz1,wz2};
    const real dwx[3] = {dwx0,dwx1,dwx2};
    const real dwy[3] = {dwy0,dwy1,dwy2};
    const real dwz[3] = {dwz0,dwz1,dwz2};

    for (int kk=0; kk<3; ++kk) {
        for (int jj=0; jj<3; ++jj) {
            // pré-accum en x pour réduire ops (optionnel)
            real sx  = 0.0;
            real sdx = 0.0;

            const int yy = ys[jj];
            const int zz = zs[kk];

            // 3 points en x
            real f0 = FVAL(xs[0], yy, zz);
            real f1 = FVAL(xs[1], yy, zz);
            real f2 = FVAL(xs[2], yy, zz);

            sx  = wx[0]*f0 + wx[1]*f1 + wx[2]*f2;
            sdx = dwx[0]*f0 + dwx[1]*f1 + dwx[2]*f2;

            const real wywz   = wy[jj]*wz[kk];
            const real dwywz  = dwy[jj]*wz[kk];
            const real wydwz  = wy[jj]*dwz[kk];

            val     += wywz  * sx;
            dtx     += wywz  * sdx;     // d/dtx
            dty_acc += dwywz * sx;      // d/dty
            dtz_acc += wydwz * sx;      // d/dtz
        }
    }

    *f = val;
    g[0] = dtx     * str->x_inv_grid;
    g[1] = dty_acc * str->y_inv_grid;
    g[2] = dtz_acc * str->z_inv_grid;

    return 0;
}


#endif
