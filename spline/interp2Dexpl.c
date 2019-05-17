/**
 * @file interp2Dexpl.c
 * @brief Bicubic spline interpolation in explicit form
 */
#include <stdlib.h>
#include <math.h>
#include "../ascot5.h"
#include "interp.h"
#include "spline.h"

/**
 * @brief Calculate bicubic spline interpolation coefficients for scalar 2D data
 *
 * This function calculates the bicubic spline interpolation coefficients and
 * stores them in a pre-allocated array. Explicit cofficients are calculated.
 *
 * For each data point four coefficients are stored for spline-interpolation.
 *
 * @param c allocated array of length n_y*n_x*16 to store the coefficients
 * @param f 2D data to be interpolated
 * @param n_x number of data points in the x direction
 * @param n_y number of data points in the y direction
 * @param bc_x boundary condition for x axis (0) natural (1) periodic
 * @param bc_y boundary condition for y axis (0) natural (1) periodic
 * @param x_min minimum value of the x axis
 * @param x_max maximum value of the x axis
 * @param y_min minimum value of the y axis
 * @param y_max maximum value of the y axis
 *
 * @return zero if initialization succeeded
 */
int interp2Dexpl_init_coeff(real* c, real* f,
                            int n_x, int n_y, int bc_x, int bc_y,
                            real x_min, real x_max,
                            real y_min, real y_max) {

    /* Allocate helper quantities */
    real* f_x = malloc(n_x*sizeof(real));
    real* f_y = malloc(n_y*sizeof(real));
    real* c_x = malloc((n_x-1*(bc_x==NATURALBC))*NSIZE_EXPL1D*sizeof(real));
    real* c_y = malloc((n_y-1*(bc_y==NATURALBC))*NSIZE_EXPL1D*sizeof(real));
    int i_ct;

    if(f_x == NULL || f_y == NULL || c_x == NULL || c_y == NULL) {
        return 1;
    }

    /* Calculate bicubic spline surface coefficients. For each grid cell
       (i_x, i_y), there are 16 coefficients, one for each variable product
       dx^p_x*dy^p_y in the evaluation formula, where p_x, p_y = 0, 1, 2, 3. */

    /* Cubic spline along x for each y, using f values to get a total of four
       coefficients */
    for(int i_y=0; i_y<n_y; i_y++) {
        for(int i_x=0; i_x<n_x; i_x++) {
            f_x[i_x] = f[i_y*n_x+i_x];
        }
        splineexpl(f_x, n_x, bc_x, c_x);
        for(int i_x=0; i_x<n_x-1; i_x++) {
            for(int i_c=0; i_c<4; i_c++) {
                c[i_y*n_x*16+i_x*16+i_c] = c_x[i_x*4+i_c];
            }
        }
    }

    /* Four cubic splines along y for each x, using the above calculated four
       coefficient values to get a total of 16 coefficients */
    for(int i_x=0; i_x<n_x-1; i_x++) {
        for(int i_s=0; i_s<4; i_s++) {
            for(int i_y=0; i_y<n_y; i_y++) {
                f_y[i_y] = c[i_y*n_x*16+i_x*16+i_s];
            }
            splineexpl(f_y, n_y, bc_y, c_y);
            for(int i_y=0; i_y<n_y-1; i_y++) {
                i_ct = 0;
                for(int i_c=i_s; i_c<16; i_c=i_c+4) {
                    c[i_y*n_x*16+i_x*16+i_c] = c_y[i_y*4+i_ct];
                    i_ct++;
                }
            }
        }
    }

    /* Free allocated memory */
    free(f_x);
    free(f_y);
    free(c_x);
    free(c_y);

    return 0;
}

/**
 * @brief Initialize a bicubic spline
 *
 * @param str pointer to spline to be initialized
 * @param c array where coefficients are stored
 * @param n_x number of data points in the x direction
 * @param n_y number of data points in the y direction
 * @param bc_x boundary condition for x axis
 * @param bc_y boundary condition for y axis
 * @param x_min minimum value of the x axis
 * @param x_max maximum value of the x axis
 * @param y_min minimum value of the y axis
 * @param y_max maximum value of the y axis
 */
void interp2Dexpl_init_spline(interp2D_data* str, real* c,
                              int n_x, int n_y, int bc_x, int bc_y,
                              real x_min, real x_max,
                              real y_min, real y_max) {

    /* Calculate grid intervals. For periodic boundary condition, grid maximum
       value and the last data point are not the same. Take this into account
       in grid intervals. */
    real x_grid = (x_max - x_min) / ( n_x - 1 * (bc_x == NATURALBC) );
    real y_grid = (y_max - y_min) / ( n_y - 1 * (bc_y == NATURALBC) );

    /* Initialize the interp2D_data struct */
    str->n_x    = n_x;
    str->n_y    = n_y;
    str->bc_x   = bc_x;
    str->bc_y   = bc_y;
    str->x_min  = x_min;
    str->x_max  = x_max;
    str->x_grid = x_grid;
    str->y_min  = y_min;
    str->y_max  = y_max;
    str->y_grid = y_grid;
    str->c      = c;
}

/**
 * @brief Evaluate interpolated value of 2D scalar field
 *
 * This function evaluates the interpolated value of a 2D scalar field using
 * bicubic spline interpolation coefficients of the explicit form.
 *
 * @param f variable in which to place the evaluated value
 * @param str data struct for data interpolation
 * @param x x-coordinate
 * @param y y-coordinate
 *
 * @return zero on success and one if (x,y) point is outside the domain.
 */
int interp2Dexpl_eval_f(real* f, interp2D_data* str, real x, real y) {

    /* Make sure periodic coordinates are within [min, max] region. */
    if(str->bc_x == PERIODICBC) {
        x = fmod(x - str->x_min, str->x_max - str->x_min) + str->x_min;
        x = x + (x < str->x_min) * (str->x_max - str->x_min);
    }
    if(str->bc_y == PERIODICBC) {
        y = fmod(y - str->y_min, str->y_max - str->y_min) + str->y_min;
        y = y + (y < str->y_min) * (str->y_max - str->y_min);
    }

    /* Index for x variable. The -1 needed at exactly grid end. */
    int i_x = (x-str->x_min)/str->x_grid - 1*(x==str->x_max);
    /* Normalized x coordinate in current cell */
    real dx = (x-(str->x_min+i_x*str->x_grid))/str->x_grid;
    /* Helper variables */
    real dx2 = dx*dx;
    real dx3 = dx2*dx;

    /* Index for y variable. The -1 needed at exactly grid end. */
    int i_y = (y-str->y_min)/str->y_grid - 1*(y==str->y_max);
    /* Normalized y coordinate in current cell */
    real dy = (y-(str->y_min+i_y*str->y_grid))/str->y_grid;
    /* Helper variables */
    real dy2 = dy*dy;
    real dy3 = dy2*dy;

    int n = i_y*str->n_x*16+i_x*16; /* Index jump to cell */

    int err = 0;

    /* Check that the coordinate is within the domain. */
    if( str->bc_x == NATURALBC && (x < str->x_min || x > str->x_max) ) {
        err = 1;
    }
    if( str->bc_y == NATURALBC && (y < str->y_min || y > str->y_max) ) {
        err = 1;
    }

    if(!err) {
        *f =
                str->c[n+ 0]+dx*str->c[n+ 1]+dx2*str->c[n+ 2]+dx3*str->c[n+ 3]
            +dy*(
                str->c[n+ 4]+dx*str->c[n+ 5]+dx2*str->c[n+ 6]+dx3*str->c[n+ 7])
            +dy2*(
                str->c[n+ 8]+dx*str->c[n+ 9]+dx2*str->c[n+10]+dx3*str->c[n+11])
            +dy3*(
                str->c[n+12]+dx*str->c[n+13]+dx2*str->c[n+14]+dx3*str->c[n+15]);
    }

    return err;
}

/**
 * @brief Evaluate interpolated value and 1st and 2nd derivatives of 2D field
 *
 * This function evaluates the interpolated value of a 2D scalar field and
 * its 1st and 2nd derivatives using bicubic spline interpolation coefficients
 * of the explicit form.
 *
 * The evaluated  values are returned in an array with following elements:
 * - f_df[0] = f
 * - f_df[1] = f_x
 * - f_df[2] = f_y
 * - f_df[3] = f_xx
 * - f_df[4] = f_yy
 * - f_df[5] = f_xy
 *
 * @param f_df array in which to place the evaluated values
 * @param str data struct for data interpolation
 * @param x x-coordinate
 * @param y y-coordinate
 */
int interp2Dexpl_eval_df(real* f_df, interp2D_data* str, real x, real y) {

    /* Make sure periodic coordinates are within [min, max] region. */
    if(str->bc_x == PERIODICBC) {
        x = fmod(x - str->x_min, str->x_max - str->x_min) + str->x_min;
        x = x + (x < str->x_min) * (str->x_max - str->x_min);
    }
    if(str->bc_y == PERIODICBC) {
        y = fmod(y - str->y_min, str->y_max - str->y_min) + str->y_min;
        y = y + (y < str->y_min) * (str->y_max - str->y_min);
    }

    /* Index for x variable. The -1 needed at exactly grid end. */
    int i_x = (x-str->x_min)/str->x_grid - 1*(x==str->x_max);
    /* Normalized x coordinate in current cell */
    real dx = (x-(str->x_min+i_x*str->x_grid))/str->x_grid;
    /* Helper variables */
    real dx2 = dx*dx;
    real dx3 = dx2*dx;
    real xgi = 1.0/str->x_grid;

    /* Index for y variable. The -1 needed at exactly grid end. */
    int i_y = (y-str->y_min)/str->y_grid - 1*(y==str->y_max);
    /* Normalized y coordinate in current cell */
    real dy = (y-(str->y_min+i_y*str->y_grid))/str->y_grid;
    /* Helper variables */
    real dy2 = dy*dy;
    real dy3 = dy2*dy;
    real ygi = 1.0/str->y_grid;

    int n = i_y*str->n_x*16+i_x*16; /* Index jump to cell */

    int err = 0;

    /* Check that the coordinate is within the domain. */
    if( str->bc_x == NATURALBC && (x < str->x_min || x > str->x_max) ) {
        err = 1;
    }
    if( str->bc_y == NATURALBC && (y < str->y_min || y > str->y_max) ) {
        err = 1;
    }

    if(!err) {
        /* f */
        f_df[0] =
                str->c[n+ 0]+dx*str->c[n+ 1]+dx2*str->c[n+ 2]+dx3*str->c[n+ 3]
            +dy*(
                str->c[n+ 4]+dx*str->c[n+ 5]+dx2*str->c[n+ 6]+dx3*str->c[n+ 7])
            +dy2*(
                str->c[n+ 8]+dx*str->c[n+ 9]+dx2*str->c[n+10]+dx3*str->c[n+11])
            +dy3*(
                str->c[n+12]+dx*str->c[n+13]+dx2*str->c[n+14]+dx3*str->c[n+15]);

        /* df/dx */
        f_df[1] =
            xgi*(
                      str->c[n+ 1]+2*dx*str->c[n+ 2]+3*dx2*str->c[n+ 3]
                 +dy*(str->c[n+ 5]+2*dx*str->c[n+ 6]+3*dx2*str->c[n+ 7])
                +dy2*(str->c[n+ 9]+2*dx*str->c[n+10]+3*dx2*str->c[n+11])
                +dy3*(str->c[n+13]+2*dx*str->c[n+14]+3*dx2*str->c[n+15]));

        /* df/dy */
        f_df[2] =
            ygi*(
                         str->c[n+ 4]+dx *str->c[n+ 5]
                    +dx2*str->c[n+ 6]+dx3*str->c[n+ 7]
                +2*dy*(
                         str->c[n+ 8]+dx *str->c[n+ 9]
                    +dx2*str->c[n+10]+dx3*str->c[n+11])
                +3*dy2*(
                        str->c[n+12]+dx *str->c[n+13]
                    +dx2*str->c[n+14]+dx3*str->c[n+15]));

        /* d2f/dx2 */
        f_df[3] =
            xgi*xgi*(
                      2*str->c[n+ 2]+6*dx*str->c[n+ 3]
                 +dy*(2*str->c[n+ 6]+6*dx*str->c[n+ 7])
                +dy2*(2*str->c[n+10]+6*dx*str->c[n+11])
                +dy3*(2*str->c[n+14]+6*dx*str->c[n+15]));

        /* d2f/dy2 */
        f_df[4] =
            ygi*ygi*(
                2*(         str->c[n+ 8]+dx *str->c[n+ 9]
                       +dx2*str->c[n+10]+dx3*str->c[n+11])
                +6*dy*(     str->c[n+12]+dx *str->c[n+13]
                       +dx2*str->c[n+14]+dx3*str->c[n+15]));

        /* d2f/dydx */
        f_df[5] =
            xgi*ygi*(
                        str->c[n+ 5]+2*dx*str->c[n+ 6]+3*dx2*str->c[n+ 7]
                 +2*dy*(str->c[n+ 9]+2*dx*str->c[n+10]+3*dx2*str->c[n+11])
                +3*dy2*(str->c[n+13]+2*dx*str->c[n+14]+3*dx2*str->c[n+15]));
    }

    return err;
}
