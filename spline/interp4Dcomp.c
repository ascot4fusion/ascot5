/**
 * @file interp4Dcomp.c
 * @brief 4 coordinates cubic spline interpolation in compact form
 */
#include <stdlib.h>
#include <math.h>
#include "../ascot5.h"
#include "../consts.h"
#include "interp.h"
#include "spline.h"

/**
 * @brief Calculate 4D cubic spline interpolation coefficients for 4D data
 *
 * This function calculates the 4D cubic spline interpolation coefficients for
 * the given data and stores them in an array. Compact cofficients are
 * calculated directly.
 *
 * @param c allocated array of length n_t*n_z*n_y*n_x*16 to store the coefficients
 * @param f 3D data to be interpolated
 * @param n_x number of data points in the x direction
 * @param n_y number of data points in the y direction
 * @param n_z number of data points in the z direction
 * @param n_t number of data points in the t direction
 * @param bc_x boundary condition for x axis
 * @param bc_y boundary condition for y axis
 * @param bc_z boundary condition for z axis
 * @param bc_t boundary condition for t axis
 * @param x_min minimum value of the x axis
 * @param x_max maximum value of the x axis
 * @param y_min minimum value of the y axis
 * @param y_max maximum value of the y axis
 * @param z_min minimum value of the z axis
 * @param z_max maximum value of the z axis
 * @param t_min minimum value of the t axis
 * @param t_max maximum value of the t axis
 *
 * @return zero if initialization succeeded
 */
int interp4Dcomp_init_coeff(real* c, real* f,
                            int n_x, int n_y, int n_z, int n_t,
                            int bc_x, int bc_y, int bc_z, int bc_t,
                            real x_min, real x_max,
                            real y_min, real y_max,
                            real z_min, real z_max,
                            real t_min, real t_max) {

    /* For periodic boundary condition, grid maximum value and the last data
       point are not the same. Take this into account in grid interval       */
    real x_grid, y_grid, z_grid, t_grid;
    if(bc_x == NATURALBC || bc_x == PERIODICBC) {
        x_grid = (x_max - x_min) / ( n_x - 1 * (bc_x == NATURALBC) );
    }
    else {
        return 1;
    }

    if(bc_y == NATURALBC || bc_y == PERIODICBC) {
        y_grid = (y_max - y_min) / ( n_y - 1 * (bc_y == NATURALBC) );
    }
    else {
        return 1;
    }

    if(bc_z == NATURALBC || bc_z == PERIODICBC) {
        z_grid = (z_max - z_min) / ( n_z - 1 * (bc_z == NATURALBC) );
    }
    else {
        return 1;
    }

    if(bc_t == NATURALBC || bc_t == PERIODICBC) {
        t_grid = (t_max - t_min) / ( n_t - 1 * (bc_t == NATURALBC) );
    }
    else {
        return 1;
    }

    /* Allocate helper quantities */
    real* f_x = malloc(n_x*sizeof(real));
    real* f_y = malloc(n_y*sizeof(real));
    real* f_z = malloc(n_z*sizeof(real));
    real* f_t = malloc(n_t*sizeof(real));
    real* c_x = malloc(n_x*2*sizeof(real));
    real* c_y = malloc(n_y*2*sizeof(real));
    real* c_z = malloc(n_z*2*sizeof(real));
    real* c_t = malloc(n_t*2*sizeof(real));


    if(f_x == NULL || f_y == NULL || f_z == NULL ||  f_t == NULL ||
       c_x == NULL || c_y == NULL || c_z == NULL || c_t == NULL) {
        return 1;
    }

    /* 4D cubic spline volume coefficients: For i_x, k_y, j_z, m_t the 16
       coefficients are:
       c[NSIZE_COMP4D*(m_t*n_x*n_z*n_y + k_y*n_x*n_z + j_z*n_x + i_x) + 0] = D1111 = f;
       c[NSIZE_COMP4D*(m_t*n_x*n_z*n_y + k_y*n_x*n_z + j_z*n_x + i_x) + 1] = D2111 = fxx;
       c[NSIZE_COMP4D*(m_t*n_x*n_z*n_y + k_y*n_x*n_z + j_z*n_x + i_x) + 2] = D1211 = fzz;
       c[NSIZE_COMP4D*(m_t*n_x*n_z*n_y + k_y*n_x*n_z + j_z*n_x + i_x) + 3] = D2211 = fxxzz;
       c[NSIZE_COMP4D*(m_t*n_x*n_z*n_y + k_y*n_x*n_z + j_z*n_x + i_x) + 4] = D1121 = fyy;
       c[NSIZE_COMP4D*(m_t*n_x*n_z*n_y + k_y*n_x*n_z + j_z*n_x + i_x) + 5] = D2121 = fxxyy;
       c[NSIZE_COMP4D*(m_t*n_x*n_z*n_y + k_y*n_x*n_z + j_z*n_x + i_x) + 6] = D1221 = fzzyy;
       c[NSIZE_COMP4D*(m_t*n_x*n_z*n_y + k_y*n_x*n_z + j_z*n_x + i_x) + 7] = D2221 = fxxzzyy;
       c[NSIZE_COMP4D*(m_t*n_x*n_z*n_y + k_y*n_x*n_z + j_z*n_x + i_x) + 8] = D1112 = ftt;
       c[NSIZE_COMP4D*(m_t*n_x*n_z*n_y + k_y*n_x*n_z + j_z*n_x + i_x) + 9] = D2112 = fxxtt;
       c[NSIZE_COMP4D*(m_t*n_x*n_z*n_y + k_y*n_x*n_z + j_z*n_x + i_x) +10] = D1212 = fzztt;
       c[NSIZE_COMP4D*(m_t*n_x*n_z*n_y + k_y*n_x*n_z + j_z*n_x + i_x) +11] = D2212 = fxxzztt;
       c[NSIZE_COMP4D*(m_t*n_x*n_z*n_y + k_y*n_x*n_z + j_z*n_x + i_x) +12] = D1122 = fyytt;
       c[NSIZE_COMP4D*(m_t*n_x*n_z*n_y + k_y*n_x*n_z + j_z*n_x + i_x) +13] = D2122 = fxxyytt;
       c[NSIZE_COMP4D*(m_t*n_x*n_z*n_y + k_y*n_x*n_z + j_z*n_x + i_x) +14] = D1222 = fzzyytt;
       c[NSIZE_COMP4D*(m_t*n_x*n_z*n_y + k_y*n_x*n_z + j_z*n_x + i_x) +15] = D2222 = fxxzzyytt;
       Note how we account for normalized grid. */

    for(int m_t=0; m_t<n_t; m_t++) {
        for(int k_y=0; k_y<n_z; k_y++) {
            for(int j_z=0; j_z<n_y; j_z) {
                /* Cubic spline of f along x
                   to get fxx (for each j_z,k_y,m_t) */
                for(int i_x=0; i_x<n_x; j_x) {
                    f_x[i_x] = f[m_t*n_x*n_z*n_y + k_y*n_x*n_z + j_z*n_x + i_x];
                }
                splinecomp(f_x, n_x, bc_x, c_x);
                for(int i_x=0; i_x<n_x; i_x++) {
                    c[NSIZE_COMP4D*(m_t*n_x*n_z*n_y + k_y*n_x*n_z + j_z*n_x + i_x) + 0] =
                        c_x[i_x*2];
                    c[NSIZE_COMP4D*(m_t*n_x*n_z*n_y + k_y*n_x*n_z + j_z*n_x + i_x) + 1] =
                        c_x[i_x*2 + 1]/(x_grid*x_grid);
                }
            }
            /* Cubic spline of f, fxx along z
               to get fzz, fxxzz (for each i_x,k_y,m_t) */
            for(int i_x=0; i_x<n_x; i_x++) {
                /* fzz */
                for(int j_z=0; j_z<n_z; j_z++) {
                    f_z[j_z] = f[m_t*n_x*n_z*n_y + k_y*n_x*n_z + j_z*n_x + i_x];
                }
                splinecomp(f_z, n_z, bc_z, c_z);
                for(int j_z=0; j_z<n_z; j_z++) {
                    c[NSIZE_COMP4D*(m_t*n_x*n_z*n_y + k_y*n_x*n_z + j_z*n_x + i_x) + 2] =
                        c_z[j_z*2 + 1]/(z_grid*z_grid);
                }
                /* fxxzz */
                for(int j_z=0; j_z<n_z; j_z++) {
                    f_z[j_z] =
                        c[NSIZE_COMP4D*(m_t*n_x*n_z*n_y + k_y*n_x*n_z + j_z*n_x + i_x) + 1];
                }
                splinecomp(f_z, n_z, bc_z, c_z);
                for(int j_z=0; j_z<n_z; j_z++) {
                    c[NSIZE_COMP4D*(m_t*n_x*n_z*n_y + k_y*n_x*n_z + j_z*n_x + i_x) + 3] =
                        c_z[j_z*2 + 1]/(z_grid*z_grid);
                }
            }
        }
        /* Cubic spline of f, fxx, fzz, fxxzz along y
           to get fyy, fxxyy, fzzyy, fxxzzyy (for each i_x,j_z,m_t) */
        for(int j_z=0; j_z<n_z; j_z++) {
            for(int i_x=0; i_x<n_x; i_x++) {
                /* fyy */
                for(int k_y=0; k_y<n_y; k_y++) {
                    f_y[k_y] = f[m_t*n_x*n_z*n_y + k_y*n_x*n_z + j_z*n_x + i_x];
                }
                splinecomp(f_y, n_y, bc_y, c_y);
                for(int k_y=0; k_y<n_y; k_y++) {
                    c[NSIZE_COMP4D*(m_t*n_x*n_z*n_y + k_y*n_x*n_z + j_z*n_x + i_x) + 4] =
                        c_y[k_y*2 + 1]/(y_grid*y_grid);
                }
                /* fxxyy */
                for(int k_y=0; k_y<n_y; k_y++) {
                    f_y[k_y] =
                        c[NSIZE_COMP4D*(m_t*n_x*n_z*n_y + k_y*n_x*n_z + j_z*n_x + i_x) + 1];
                }
                splinecomp(f_y, n_y, bc_y, c_y);
                for(int k_y=0; k_y<n_y; k_y++) {
                    c[NSIZE_COMP4D*(m_t*n_x*n_z*n_y + k_y*n_x*n_z + j_z*n_x + i_x) + 5] =
                        c_y[k_y*2 + 1]/(y_grid*y_grid);
                }
                /* fzzyy */
                for(int k_y=0; k_y<n_y; k_y++) {
                    f_y[k_y] =
                        c[NSIZE_COMP4D*(m_t*n_x*n_z*n_y + k_y*n_x*n_z + j_z*n_x + i_x) + 2];
                }
                splinecomp(f_y, n_y, bc_y, c_y);
                for(int k_y=0; k_y<n_y; k_y++) {
                    c[NSIZE_COMP4D*(m_t*n_x*n_z*n_y + k_y*n_x*n_z + j_z*n_x + i_x) + 6] =
                        c_y[k_y*2 + 1]/(y_grid*y_grid);
                }
                /* fxxzzyy */
                for(int k_y=0; k_y<n_y; k_y++) {
                    f_y[k_y] =
                        c[NSIZE_COMP4D*(m_t*n_x*n_z*n_y + k_y*n_x*n_z + j_z*n_x + i_x) + 3];
                }
                splinecomp(f_y, n_y, bc_y, c_y);
                for(int k_y=0; k_y<n_y; k_y++) {
                    c[NSIZE_COMP4D*(m_t*n_x*n_z*n_y + k_y*n_x*n_z + j_z*n_x + i_x) + 7] =
                        c_y[k_y*2 + 1]/(y_grid*y_grid);
                }
            }
        }
    }
    /* Cubic spline of f, fxx, fzz, fxxzz, fyy, fxxyy, fzzyy, fxxzzyy along t
       to get ftt, fxxtt, fzztt, fxxzztt, fyytt, fxxyytt, fzzyytt, fxxzzyytt
       (for each i_x,j_z,k_y) */

    /* Cubic spline along y for each xz-pair to find the compact coefficients
       of the tricubic spline volume */
    for(int k_y=0; k_y<n_y; k_y++) {
        for(int j_z=0; j_z<n_z; j_z++) {
            for(int i_x=0; i_x<n_x; i_x++) {
                /* ftt */
                for(int m_t=0; m_t<n_t; m_t++) {
                    f_t[m_t] = f[m_t*n_x*n_z*n_y + k_y*n_x*n_z + j_z*n_x + i_x];
                }
                splinecomp(f_t, n_t, bc_t, c_t);
                for(int m_t=0; m_t<n_t; m_t++) {
                    c[NSIZE_COMP4D*(m_t*n_x*n_z*n_y + k_y*n_x*n_z + j_z*n_x + i_x) + 8] =
                        c_t[m_t*2 + 1]/(t_grid*t_grid);
                }
                /* c[9:15] are calculated from c[1:7] in a loop */
                for(int i_c=1; i_c<8; i_c++) {
                    for(int m_t=0; m_t<n_t; m_t++) {
                        f_t[m_t] =
                            c[NSIZE_COMP4D*(m_t*n_x*n_z*n_y + k_y*n_x*n_z + j_z*n_x + i_x) + i_c];
                    }
                    splinecomp(f_t, n_t, bc_t, c_t);
                    for(int m_t=0; m_t<n_t; m_t++) {
                        c[NSIZE_COMP4D*(m_t*n_x*n_z*n_y + k_y*n_x*n_z + j_z*n_x + i_x) + i_c + 8] =
                            c_t[m_t*2 + 1]/(t_grid*t_grid);
                    }
                }
            }
        }
    }

    /* Free allocated memory */
    free(f_x);
    free(f_y);
    free(f_z);
    free(f_t);
    free(c_x);
    free(c_y);
    free(c_z);
    free(c_t);


    return 0;
}

/**
 * @brief Initialize a 4D cubic spline
 *
 * @param str pointer to spline to be initialized
 * @param c array where coefficients are stored
 * @param n_x number of data points in the x direction
 * @param n_y number of data points in the y direction
 * @param n_z number of data points in the z direction
 * @param n_t number of data points in the t direction
 * @param bc_x boundary condition for x axis
 * @param bc_y boundary condition for y axis
 * @param bc_z boundary condition for z axis
 * @param bc_t boundary condition for t axis
 * @param x_min minimum value of the x axis
 * @param x_max maximum value of the x axis
 * @param y_min minimum value of the y axis
 * @param y_max maximum value of the y axis
 * @param z_min minimum value of the z axis
 * @param z_max maximum value of the z axis
 * @param t_min minimum value of the t axis
 * @param t_max maximum value of the t axis
 */
void interp4Dcomp_init_spline(interp4D_data* str, real* c,
                              int n_x, int n_y, int n_z, int n_t,
                              int bc_x, int bc_y, int bc_z, int bc_t,
                              real x_min, real x_max,
                              real y_min, real y_max,
                              real z_min, real z_max,
                              real t_min, real t_max) {

    real x_grid = (x_max - x_min) / ( n_x - 1 * (bc_x == NATURALBC) );
    real y_grid = (y_max - y_min) / ( n_y - 1 * (bc_y == NATURALBC) );
    real z_grid = (z_max - z_min) / ( n_z - 1 * (bc_z == NATURALBC) );
    real t_grid = (t_max - t_min) / ( n_t - 1 * (bc_t == NATURALBC) );


    str->n_x    = n_x;
    str->n_y    = n_y;
    str->n_z    = n_z;
    str->n_t    = n_t;
    str->bc_x   = bc_x;
    str->bc_y   = bc_y;
    str->bc_z   = bc_z;
    str->bc_t   = bc_t;
    str->x_min  = x_min;
    str->x_max  = x_max;
    str->x_grid = x_grid;
    str->y_min  = y_min;
    str->y_max  = y_max;
    str->y_grid = y_grid;
    str->z_min  = z_min;
    str->z_max  = z_max;
    str->z_grid = z_grid;
    str->t_min  = t_min;
    str->t_max  = t_max;
    str->t_grid = t_grid;
    str->c      = c;
}

/**
 * @brief Evaluate interpolated value of 4D scalar field
 *
 * This function evaluates the interpolated value of a 4D scalar field using
 * 4D cubic spline interpolation coefficients of the compact form.
 *
 * @param f variable in which to place the evaluated value
 * @param str data struct for data interpolation
 * @param x x-coordinate
 * @param y y-coordinate
 * @param z z-coordinate
 * @param t t-coordinate 
 *
 * @return zero on success and one if (x,y,z,t) point is outside the grid.
 */
int interp4Dcomp_eval_f(real* f, interp4D_data* str, real x, real y, real z, real t) {

    /* Make sure periodic coordinates are within [max, min] region. */
    if(str->bc_x == PERIODICBC) {
        x = fmod(x - str->x_min, str->x_max - str->x_min) + str->x_min;
        x = x + (x < str->x_min) * (str->x_max - str->x_min);
    }
    if(str->bc_y == PERIODICBC) {
        y = fmod(y - str->y_min, str->y_max - str->y_min) + str->y_min;
        y = y + (y < str->y_min) * (str->y_max - str->y_min);
    }
    if(str->bc_z == PERIODICBC) {
        z = fmod(z - str->z_min, str->z_max - str->z_min) + str->z_min;
        z = z + (z < str->z_min) * (str->z_max - str->z_min);
    }
    if(str->bc_t == PERIODICBC) {
        t = fmod(t - str->t_min, str->t_max - str->t_min) + str->t_min;
        t = t + (t < str->t_min) * (str->t_max - str->t_min);
    }

    /* index for x variable */
    int i_x   = (x - str->x_min) / str->x_grid;
    /* Normalized x coordinate in current cell
       (the order facilitates the evaluation)*/
    real dx[4];
    dx[2] = (x - (str->x_min + i_x*str->x_grid)) / str->x_grid; /* p_x */
    dx[0] = 1.0 - dx[2];                                        /* q_x */
    dx[1] = dx[0]*(dx[0]*dx[0]-1)*str->x_grid*str_x_grid/6.0;   /* s_x */
    dx[3] = dx[2]*(dx[2]*dx[2]-1)*str->x_grid*str_x_grid/6.0;   /* r_x */

    /* index for y variable */
    int k_y   = (y - str->y_min) / str->y_grid;
   /* Normalized y coordinate in current cell
       (the order facilitates the evaluation)*/
    real dy[4];
    dy[2] = (y - (str->y_min + k_y*str->y_grid)) / str->y_grid; /* p_y */
    dy[0] = 1.0 - dy[3];                                        /* q_y */
    dy[1] = dy[0]*(dy[0]*dy[0]-1)*str->y_grid*str_y_grid/6.0;   /* s_y */
    dy[3] = dy[2]*(dy[2]*dy[2]-1)*str->y_grid*str_y_grid/6.0;   /* r_y */

    /* index for z variable */
    int j_z   = (z - str->z_min) / str->z_grid;
   /* Normalized z coordinate in current cell
       (the order facilitates the evaluation)*/
    real dz[4];
    dz[2] = (z - (str->z_min + j_z*str->z_grid)) / str->z_grid; /* p_z */
    dz[0] = 1.0 - dz[3];                                        /* q_z */
    dz[1] = dz[0]*(dz[0]*dz[0]-1)*str->z_grid*str_z_grid/6.0;   /* s_z */
    dz[3] = dz[2]*(dz[2]*dz[2]-1)*str->z_grid*str_z_grid/6.0;   /* r_z */

    /* index for t variable */
    int m_t   = (t - str->t_min) / str->t_grid;
   /* Normalized t coordinate in current cell
       (the order facilitates the evaluation)*/
    real dt[4];
    dt[2] = (t - (str->t_min + m_t*str->t_grid)) / str->t_grid; /* p_t */
    dt[0] = 1.0 - dt[3];                                        /* q_t */
    dt[1] = dt[0]*(dt[0]*dt[0]-1)*str->t_grid*str_t_grid/6.0;   /* s_t */
    dt[3] = dt[2]*(dt[2]*dt[2]-1)*str->t_grid*str_t_grid/6.0;   /* r_t */

    /**< Index jump to cell */
    int n  = NSIZE_COMP4D*
        (m_t*str->n_x*str->n_z*str->n_y + k_y*str->n_x*str->n_z + j_z*str->n_x + i_x);
    int x1 = NSIZE_COMP4D;                           /* Index jump one x forward */
    int y1 = NSIZE_COMP4D*str->n_z*str->n_x;         /* Index jump one y forward */
    int z1 = NSIZE_COMP4D*str->n_x;                  /* Index jump one z forward */
    int t1 = NSIZE_COMP4D*str->n_x*str->n_z*str->n_y;/* Index jump one z forward */

    int err = 0;

    /* Enforce periodic BC or check that the coordinate is within the grid. */
    if( str->bc_x == PERIODICBC && i_x == str->n_x-1 ) {
        x1 = -(str->n_x-1)*x1;
    }
    else if( str->bc_x == NATURALBC && (x < str->x_min || x > str->x_max) ) {
        err = 1;
    }
    if( str->bc_y == PERIODICBC && k_y == str->n_y-1 ) {
        y1 = -(str->n_y-1)*y1;
    }
    else if( str->bc_y == NATURALBC && (y < str->y_min || y > str->y_max) ) {
        err = 1;
    }
    if( str->bc_z == PERIODICBC && j_z == str->n_z-1 ) {
        z1 = -(str->n_z-1)*z1;
    }
    else if( str->bc_z == NATURALBC && (z < str->z_min || z > str->z_max) ) {
        err = 1;
    }
    if( str->bc_t == PERIODICBC && m_t == str->n_t-1 ) {
        t1 = -(str->n_t-1)*t1;
    }
    else if( str->bc_t == NATURALBC && (t < str->t_min || t > str->t_max) ) {
        err = 1;
    }

    if(!err) {

        /* Evaluate splines */
        real d_aux1, d_aux2; /**< Auxiliary normalized coordinates */
        *f = 0;
        /* loops to move through the nodes */
        for(int i_node_x=0; i_node_x<2; i_node_x++) {
            for(int i_node_z=0; i_node_z<2; i_node_z++) {
                for(int i_node_y=0; i_node_y<2; i_node_y++) {
                    for(int i_node_t=0; i_node_t<2; i_node_t++) {
                        /* loops to move through the coefficients */
                        for(int i_coef_t=0; i_coef_t<2; i_coef_t++) {
                            for(int i_coef_y=0; i_coef_y<2; i_coef_y++) {
                                d_aux1 = dt[i_node_t*2+i_coef_t]*dy[i_node_y*2+i_coef_y];
                                for(int i_coef_z=0; i_coef_z<2; i_coef_z++) {
                                    d_aux2 = d_aux1*dz[i_node_z*2+i_coef_z];
                                    for(int i_coef_x=0; i_coef_x<2; i_coef_x++) {
                                        *f +=
                                            str->c[n +
                                                   i_node_x*x1 + i_node_z*z1 +
                                                   i_node_y*y1 + i_node_t*t1 +
                                                   i_coef_t*8 + i_coef_y*4 +
                                                   i_coef_z*2 + i_coef_x]*
                                            d_aux2*dx[i_node_x*2+i_coef_x];
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    return err;
}

/**
 * @brief Evaluate interpolated value of 4D field and 1st and 2nd derivatives
 *
 * This function evaluates the interpolated value of a 4D scalar field and
 * its 1st and 2nd derivatives using bicubic spline interpolation coefficients
 * of the compact form.
 *
 * The evaluated  values are returned in an array with following elements:
 * - f_df[0] = f
 * - f_df[1] = f_x
 * - f_df[2] = f_y
 * - f_df[3] = f_z
 * - f_df[4] = f_xx
 * - f_df[5] = f_yy
 * - f_df[6] = f_zz
 * - f_df[7] = f_xy
 * - f_df[8] = f_xz
 * - f_df[9] = f_yz
 *
 * @param f_df array in which to place the evaluated values
 * @param str data struct for data interpolation
 * @param x x-coordinate
 * @param y y-coordinate
 * @param z z-coordinate
 * @param t t-coordinate
 *
 * @return zero on success and one if (x,y,z,t) point is outside the grid.
 */
int interp4Dcomp_eval_df(real* f_df, interp4D_data* str,
                         real x, real y, real z, real t) {
    /* Make sure periodic coordinates are within [max, min] region. */
    if(str->bc_x == PERIODICBC) {
        x = fmod(x - str->x_min, str->x_max - str->x_min) + str->x_min;
        x = x + (x < str->x_min) * (str->x_max - str->x_min);
    }
    if(str->bc_y == PERIODICBC) {
        y = fmod(y - str->y_min, str->y_max - str->y_min) + str->y_min;
        y = y + (y < str->y_min) * (str->y_max - str->y_min);
    }
    if(str->bc_z == PERIODICBC) {
        z = fmod(z - str->z_min, str->z_max - str->z_min) + str->z_min;
        z = z + (z < str->z_min) * (str->z_max - str->z_min);
    }
    if(str->bc_t == PERIODICBC) {
        t = fmod(t - str->t_min, str->t_max - str->t_min) + str->t_min;
        t = t + (t < str->t_min) * (str->t_max - str->t_min);
    }

    /* index for x variable */
    int i_x   = (x - str->x_min) / str->x_grid;
    /* Normalized x coordinate in current cell
       (the order facilitates the evaluation)*/
    real dx[4];
    dx[2] = (x - (str->x_min + i_x*str->x_grid)) / str->x_grid; /* p_x */
    dx[0] = 1.0 - dx[2];                                        /* q_x */
    dx[1] = dx[0]*(dx[0]*dx[0]-1)*str->x_grid*str_x_grid/6.0;   /* s_x */
    dx[3] = dx[2]*(dx[2]*dx[2]-1)*str->x_grid*str_x_grid/6.0;   /* r_x */
    /* Normalized x coordinate 1st derivates in current cell
       (the order facilitates the evaluation)*/
    real d_dx[4];
    d_dx[2] = 1 /str->x_grid;                             /* d(p_x)/dx */
    d_dx[0] =-1 /str->x_grid;                             /* d(q_x)/dx */
    d_dx[1] = str->x_grid/6.0*(1-3*dx[0]*dx[0]);          /* d(s_x)/dx */
    d_dx[3] = str->x_grid/6.0*(3*dx[2]*dx[2]-1);          /* d(r_x)/dx */
    /* Normalized x coordinate 2nd derivates in current cell
       (the order facilitates the evaluation)*/
    real dd_dx[4];
    dd_dx[2] = 0;                                       /* d2(p_x)/dx2 */
    dd_dx[0] = 0;                                       /* d2(q_x)/dx2 */
    dd_dx[1] = dx[0];                                   /* d2(s_x)/dx2 */
    dd_dx[3] = dx[2];                                   /* d2(r_x)/dx2 */

    /* index for y variable */
    int k_y   = (y - str->y_min) / str->y_grid;
   /* Normalized y coordinate in current cell
       (the order facilitates the evaluation)*/
    real dy[4];
    dy[2] = (y - (str->y_min + k_y*str->y_grid)) / str->y_grid; /* p_y */
    dy[0] = 1.0 - dy[3];                                        /* q_y */
    dy[1] = dy[0]*(dy[0]*dy[0]-1)*str->y_grid*str_y_grid/6.0;   /* s_y */
    dy[3] = dy[2]*(dy[2]*dy[2]-1)*str->y_grid*str_y_grid/6.0;   /* r_y */
    /* Normalized y coordinate 1st derivates in current cell
       (the order facilitates the evaluation)*/
    real d_dy[4];
    d_dy[2] = 1 /str->y_grid;                             /* d(p_y)/dy */
    d_dy[0] =-1 /str->y_grid;                             /* d(q_y)/dy */
    d_dy[1] = str->y_grid/6.0*(1-3*dy[0]*dy[0]);          /* d(s_y)/dy */
    d_dy[3] = str->y_grid/6.0*(3*dy[2]*dy[2]-1);          /* d(r_y)/dy */
    /* Normalized y coordinate 2nd derivates in current cell
       (the order facilitates the evaluation)*/
    real dd_dy[4];
    dd_dy[2] = 0;                                       /* d2(p_y)/dy2 */
    dd_dy[0] = 0;                                       /* d2(q_y)/dy2 */
    dd_dy[1] = dy[0];                                   /* d2(s_y)/dy2 */
    dd_dy[3] = dy[2];                                   /* d2(r_y)/dy2 */

    /* index for z variable */
    int j_z   = (z - str->z_min) / str->z_grid;
   /* Normalized z coordinate in current cell
       (the order facilitates the evaluation)*/
    real dz[4];
    dz[2] = (z - (str->z_min + j_z*str->z_grid)) / str->z_grid; /* p_z */
    dz[0] = 1.0 - dz[3];                                        /* q_z */
    dz[1] = dz[0]*(dz[0]*dz[0]-1)*str->z_grid*str_z_grid/6.0;   /* s_z */
    dz[3] = dz[2]*(dz[2]*dz[2]-1)*str->z_grid*str_z_grid/6.0;   /* r_z */
    /* Normalized z coordinate 1st derivates in current cell
       (the order facilitates the evaluation)*/
    real d_dz[4];
    d_dz[2] = 1 /str->z_grid;                             /* d(p_z)/dz */
    d_dz[0] =-1 /str->z_grid;                             /* d(q_z)/dz */
    d_dz[1] = str->z_grid/6.0*(1-3*dz[0]*dz[0]);          /* d(s_z)/dz */
    d_dz[3] = str->z_grid/6.0*(3*dz[2]*dz[2]-1);          /* d(r_z)/dz */
    /* Normalized z coordinate 2nd derivates in current cell
       (the order facilitates the evaluation)*/
    real dd_dz[4];
    dd_dz[2] = 0;                                       /* d2(p_z)/dz2 */
    dd_dz[0] = 0;                                       /* d2(q_z)/dz2 */
    dd_dz[1] = dz[0];                                   /* d2(s_z)/dz2 */
    dd_dz[3] = dz[2];                                   /* d2(r_z)/dz2 */

    /* index for t variable */
    int m_t   = (t - str->t_min) / str->t_grid;
   /* Normalized t coordinate in current cell
       (the order facilitates the evaluation)*/
    real dt[4];
    dt[2] = (t - (str->t_min + m_t*str->t_grid)) / str->t_grid; /* p_t */
    dt[0] = 1.0 - dt[3];                                        /* q_t */
    dt[1] = dt[0]*(dt[0]*dt[0]-1)*str->t_grid*str_t_grid/6.0;   /* s_t */
    dt[3] = dt[2]*(dt[2]*dt[2]-1)*str->t_grid*str_t_grid/6.0;   /* r_t */

    /**< Index jump to cell */
    int n  = NSIZE_COMP4D*
        (m_t*str->n_x*str->n_z*str->n_y + k_y*str->n_x*str->n_z + j_z*str->n_x + i_x);
    int x1 = NSIZE_COMP4D;                           /* Index jump one x forward */
    int y1 = NSIZE_COMP4D*str->n_z*str->n_x;         /* Index jump one y forward */
    int z1 = NSIZE_COMP4D*str->n_x;                  /* Index jump one z forward */
    int t1 = NSIZE_COMP4D*str->n_x*str->n_z*str->n_y;/* Index jump one z forward */

    int err = 0;

    /* Enforce periodic BC or check that the coordinate is within the grid. */
    if( str->bc_x == PERIODICBC && i_x == str->n_x-1 ) {
        x1 = -(str->n_x-1)*x1;
    }
    else if( str->bc_x == NATURALBC && (x < str->x_min || x > str->x_max) ) {
        err = 1;
    }
    if( str->bc_y == PERIODICBC && k_y == str->n_y-1 ) {
        y1 = -(str->n_y-1)*y1;
    }
    else if( str->bc_y == NATURALBC && (y < str->y_min || y > str->y_max) ) {
        err = 1;
    }
    if( str->bc_z == PERIODICBC && j_z == str->n_z-1 ) {
        z1 = -(str->n_z-1)*z1;
    }
    else if( str->bc_z == NATURALBC && (z < str->z_min || z > str->z_max) ) {
        err = 1;
    }
    if( str->bc_t == PERIODICBC && m_t == str->n_t-1 ) {
        t1 = -(str->n_t-1)*t1;
    }
    else if( str->bc_t == NATURALBC && (t < str->t_min || t > str->t_max) ) {
        err = 1;
    }

    if(!err) {

        /* Evaluate splines */
        real d_aux1, d_aux2; /**< Auxiliary normalized coordinates */

        /* f */
        f_df[0] = 0;
        /* loops to move through the nodes */
        for(int i_node_x=0; i_node_x<2; i_node_x++) {
            for(int i_node_z=0; i_node_z<2; i_node_z++) {
                for(int i_node_y=0; i_node_y<2; i_node_y++) {
                    for(int i_node_t=0; i_node_t<2; i_node_t++) {
                        /* loops to move through the coefficients */
                        for(int i_coef_t=0; i_coef_t<2; i_coef_t++) {
                            for(int i_coef_y=0; i_coef_y<2; i_coef_y++) {
                                d_aux1 = dt[i_node_t*2+i_coef_t]*dy[i_node_y*2+i_coef_y];
                                for(int i_coef_z=0; i_coef_z<2; i_coef_z++) {
                                    d_aux2 = d_aux1*dz[i_node_z*2+i_coef_z];
                                    for(int i_coef_x=0; i_coef_x<2; i_coef_x++) {
                                        f_df[0] +=
                                            str->c[n +
                                                   i_node_x*x1 + i_node_z*z1 +
                                                   i_node_y*y1 + i_node_t*t1 +
                                                   i_coef_t*8 + i_coef_y*4 +
                                                   i_coef_z*2 + i_coef_x]*
                                            d_aux2*dx[i_node_x*2+i_coef_x];
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        /* f_x */
        f_df[1] = 0;
        /* loops to move through the nodes */
        for(int i_node_x=0; i_node_x<2; i_node_x++) {
            for(int i_node_z=0; i_node_z<2; i_node_z++) {
                for(int i_node_y=0; i_node_y<2; i_node_y++) {
                    for(int i_node_t=0; i_node_t<2; i_node_t++) {
                        /* loops to move through the coefficients */
                        for(int i_coef_t=0; i_coef_t<2; i_coef_t++) {
                            for(int i_coef_y=0; i_coef_y<2; i_coef_y++) {
                                d_aux1 = dt[i_node_t*2+i_coef_t]*dy[i_node_y*2+i_coef_y];
                                for(int i_coef_z=0; i_coef_z<2; i_coef_z++) {
                                    d_aux2 = d_aux1*dz[i_node_z*2+i_coef_z];
                                    for(int i_coef_x=0; i_coef_x<2; i_coef_x++) {
                                        f_df[1] +=
                                            str->c[n +
                                                   i_node_x*x1 + i_node_z*z1 +
                                                   i_node_y*y1 + i_node_t*t1 +
                                                   i_coef_t*8 + i_coef_y*4 +
                                                   i_coef_z*2 + i_coef_x]*
                                            d_aux2*d_dx[i_node_x*2+i_coef_x];
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        /* f_y */
        f_df[2] = 0;
        /* loops to move through the nodes */
        for(int i_node_x=0; i_node_x<2; i_node_x++) {
            for(int i_node_z=0; i_node_z<2; i_node_z++) {
                for(int i_node_y=0; i_node_y<2; i_node_y++) {
                    for(int i_node_t=0; i_node_t<2; i_node_t++) {
                        /* loops to move through the coefficients */
                        for(int i_coef_t=0; i_coef_t<2; i_coef_t++) {
                            for(int i_coef_y=0; i_coef_y<2; i_coef_y++) {
                                d_aux1 = dt[i_node_t*2+i_coef_t]*d_dy[i_node_y*2+i_coef_y];
                                for(int i_coef_z=0; i_coef_z<2; i_coef_z++) {
                                    d_aux2 = d_aux1*dz[i_node_z*2+i_coef_z];
                                    for(int i_coef_x=0; i_coef_x<2; i_coef_x++) {
                                        f_df[2] +=
                                            str->c[n +
                                                   i_node_x*x1 + i_node_z*z1 +
                                                   i_node_y*y1 + i_node_t*t1 +
                                                   i_coef_t*8 + i_coef_y*4 +
                                                   i_coef_z*2 + i_coef_x]*
                                            d_aux2*dx[i_node_x*2+i_coef_x];
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        /* f_z */
        f_df[3] = 0;
        /* loops to move through the nodes */
        for(int i_node_x=0; i_node_x<2; i_node_x++) {
            for(int i_node_z=0; i_node_z<2; i_node_z++) {
                for(int i_node_y=0; i_node_y<2; i_node_y++) {
                    for(int i_node_t=0; i_node_t<2; i_node_t++) {
                        /* loops to move through the coefficients */
                        for(int i_coef_t=0; i_coef_t<2; i_coef_t++) {
                            for(int i_coef_y=0; i_coef_y<2; i_coef_y++) {
                                d_aux1 = dt[i_node_t*2+i_coef_t]*dy[i_node_y*2+i_coef_y];
                                for(int i_coef_z=0; i_coef_z<2; i_coef_z++) {
                                    d_aux2 = d_aux1*d_dz[i_node_z*2+i_coef_z];
                                    for(int i_coef_x=0; i_coef_x<2; i_coef_x++) {
                                        f_df[3] +=
                                            str->c[n +
                                                   i_node_x*x1 + i_node_z*z1 +
                                                   i_node_y*y1 + i_node_t*t1 +
                                                   i_coef_t*8 + i_coef_y*4 +
                                                   i_coef_z*2 + i_coef_x]*
                                            d_aux2*dx[i_node_x*2+i_coef_x];
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        /* f_xx */
        f_df[4] = 0;
        /* loops to move through the nodes */
        for(int i_node_x=0; i_node_x<2; i_node_x++) {
            for(int i_node_z=0; i_node_z<2; i_node_z++) {
                for(int i_node_y=0; i_node_y<2; i_node_y++) {
                    for(int i_node_t=0; i_node_t<2; i_node_t++) {
                        /* loops to move through the coefficients */
                        for(int i_coef_t=0; i_coef_t<2; i_coef_t++) {
                            for(int i_coef_y=0; i_coef_y<2; i_coef_y++) {
                                d_aux1 = dt[i_node_t*2+i_coef_t]*dy[i_node_y*2+i_coef_y];
                                for(int i_coef_z=0; i_coef_z<2; i_coef_z++) {
                                    d_aux2 = d_aux1*dz[i_node_z*2+i_coef_z];
                                    int i_coef_x=1;
                                    f_df[4] +=
                                        str->c[n +
                                               i_node_x*x1 + i_node_z*z1 +
                                               i_node_y*y1 + i_node_t*t1 +
                                               i_coef_t*8 + i_coef_y*4 +
                                               i_coef_z*2 + i_coef_x]*
                                        d_aux2*dd_dx[i_node_x*2+i_coef_x];
                                }
                            }
                        }
                    }
                }
            }
        }

        /* f_yy */
        f_df[5] = 0;
        /* loops to move through the nodes */
        for(int i_node_x=0; i_node_x<2; i_node_x++) {
            for(int i_node_z=0; i_node_z<2; i_node_z++) {
                for(int i_node_y=0; i_node_y<2; i_node_y++) {
                    for(int i_node_t=0; i_node_t<2; i_node_t++) {
                        /* loops to move through the coefficients */
                        for(int i_coef_t=0; i_coef_t<2; i_coef_t++) {
                            int i_coef_y=1;
                            d_aux1 = dt[i_node_t*2+i_coef_t]*dd_dy[i_node_y*2+i_coef_y];
                            for(int i_coef_z=0; i_coef_z<2; i_coef_z++) {
                                d_aux2 = d_aux1*dz[i_node_z*2+i_coef_z];
                                for(int i_coef_x=0; i_coef_x<2; i_coef_x++) {
                                    f_df[5] +=
                                        str->c[n +
                                               i_node_x*x1 + i_node_z*z1 +
                                               i_node_y*y1 + i_node_t*t1 +
                                               i_coef_t*8 + i_coef_y*4 +
                                               i_coef_z*2 + i_coef_x]*
                                        d_aux2*dx[i_node_x*2+i_coef_x];
                                }
                            }
                        }
                    }
                }
            }
        }

        /* f_zz */
        f_df[6] = 0;
        /* loops to move through the nodes */
        for(int i_node_x=0; i_node_x<2; i_node_x++) {
            for(int i_node_z=0; i_node_z<2; i_node_z++) {
                for(int i_node_y=0; i_node_y<2; i_node_y++) {
                    for(int i_node_t=0; i_node_t<2; i_node_t++) {
                        /* loops to move through the coefficients */
                        for(int i_coef_t=0; i_coef_t<2; i_coef_t++) {
                            for(int i_coef_y=0; i_coef_y<2; i_coef_y++) {
                                d_aux1 = dt[i_node_t*2+i_coef_t]*dy[i_node_y*2+i_coef_y];
                                int i_coef_z=1;
                                d_aux2 = d_aux1*dd_dz[i_node_z*2+i_coef_z];
                                for(int i_coef_x=0; i_coef_x<2; i_coef_x++) {
                                    f_df[6] +=
                                        str->c[n +
                                               i_node_x*x1 + i_node_z*z1 +
                                               i_node_y*y1 + i_node_t*t1 +
                                               i_coef_t*8 + i_coef_y*4 +
                                               i_coef_z*2 + i_coef_x]*
                                        d_aux2*dx[i_node_x*2+i_coef_x];
                                }
                            }
                        }
                    }
                }
            }
        }

        /* f_xy */
        f_df[7] = 0;
        /* loops to move through the nodes */
        for(int i_node_x=0; i_node_x<2; i_node_x++) {
            for(int i_node_z=0; i_node_z<2; i_node_z++) {
                for(int i_node_y=0; i_node_y<2; i_node_y++) {
                    for(int i_node_t=0; i_node_t<2; i_node_t++) {
                        /* loops to move through the coefficients */
                        for(int i_coef_t=0; i_coef_t<2; i_coef_t++) {
                            for(int i_coef_y=0; i_coef_y<2; i_coef_y++) {
                                d_aux1 = dt[i_node_t*2+i_coef_t]*d_dy[i_node_y*2+i_coef_y];
                                for(int i_coef_z=0; i_coef_z<2; i_coef_z++) {
                                    d_aux2 = d_aux1*dz[i_node_z*2+i_coef_z];
                                    for(int i_coef_x=0; i_coef_x<2; i_coef_x++) {
                                        f_df[7] +=
                                            str->c[n +
                                                   i_node_x*x1 + i_node_z*z1 +
                                                   i_node_y*y1 + i_node_t*t1 +
                                                   i_coef_t*8 + i_coef_y*4 +
                                                   i_coef_z*2 + i_coef_x]*
                                            d_aux2*d_dx[i_node_x*2+i_coef_x];
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        /* f_xz */
        f_df[8] = 0;
        /* loops to move through the nodes */
        for(int i_node_x=0; i_node_x<2; i_node_x++) {
            for(int i_node_z=0; i_node_z<2; i_node_z++) {
                for(int i_node_y=0; i_node_y<2; i_node_y++) {
                    for(int i_node_t=0; i_node_t<2; i_node_t++) {
                        /* loops to move through the coefficients */
                        for(int i_coef_t=0; i_coef_t<2; i_coef_t++) {
                            for(int i_coef_y=0; i_coef_y<2; i_coef_y++) {
                                d_aux1 = dt[i_node_t*2+i_coef_t]*dy[i_node_y*2+i_coef_y];
                                for(int i_coef_z=0; i_coef_z<2; i_coef_z++) {
                                    d_aux2 = d_aux1*d_dz[i_node_z*2+i_coef_z];
                                    for(int i_coef_x=0; i_coef_x<2; i_coef_x++) {
                                        f_df[8] +=
                                            str->c[n +
                                                   i_node_x*x1 + i_node_z*z1 +
                                                   i_node_y*y1 + i_node_t*t1 +
                                                   i_coef_t*8 + i_coef_y*4 +
                                                   i_coef_z*2 + i_coef_x]*
                                            d_aux2*d_dx[i_node_x*2+i_coef_x];
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        /* f_yz */
        f_df[9] = 0;
        /* loops to move through the nodes */
        for(int i_node_x=0; i_node_x<2; i_node_x++) {
            for(int i_node_z=0; i_node_z<2; i_node_z++) {
                for(int i_node_y=0; i_node_y<2; i_node_y++) {
                    for(int i_node_t=0; i_node_t<2; i_node_t++) {
                        /* loops to move through the coefficients */
                        for(int i_coef_t=0; i_coef_t<2; i_coef_t++) {
                            for(int i_coef_y=0; i_coef_y<2; i_coef_y++) {
                                d_aux1 = dt[i_node_t*2+i_coef_t]*d_dy[i_node_y*2+i_coef_y];
                                for(int i_coef_z=0; i_coef_z<2; i_coef_z++) {
                                    d_aux2 = d_aux1*d_dz[i_node_z*2+i_coef_z];
                                    for(int i_coef_x=0; i_coef_x<2; i_coef_x++) {
                                        f_df[9] +=
                                            str->c[n +
                                                   i_node_x*x1 + i_node_z*z1 +
                                                   i_node_y*y1 + i_node_t*t1 +
                                                   i_coef_t*8 + i_coef_y*4 +
                                                   i_coef_z*2 + i_coef_x]*
                                            d_aux2*dx[i_node_x*2+i_coef_x];
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }


    }

    return err;
}
