/**
 * @file interp3Dexpl.c
 * @brief Tricubic spline interpolation in explicit form
 */
#include <stdlib.h>
#include <math.h>
#include "../ascot5.h"
#include "../consts.h"
#include "interp.h"
#include "spline.h"

/**
 * @brief Calculate tricubic spline interpolation coefficients for 3D data
 *
 * This function calculates the tricubic spline interpolation coefficients for
 * the given data and stores them in an array. Explicit cofficients are
 * calculated.
 *
 * @param c allocated array of length n_z*n_y*n_x*64 to store the coefficients
 * @param f 3D data to be interpolated
 * @param n_x number of data points in the x direction
 * @param n_y number of data points in the y direction
 * @param n_z number of data points in the z direction
 * @param bc_x boundary condition for x axis
 * @param bc_y boundary condition for y axis
 * @param bc_z boundary condition for z axis
 * @param x_min minimum value of the x axis
 * @param x_max maximum value of the x axis
 * @param y_min minimum value of the y axis
 * @param y_max maximum value of the y axis
 * @param z_min minimum value of the z axis
 * @param z_max maximum value of the z axis
 *
 * @return zero if initialization succeeded
 */
int interp3Dexpl_init_coeff(real* c, real* f,
                            int n_x, int n_y, int n_z,
                            int bc_x, int bc_y, int bc_z,
                            real x_min, real x_max,
                            real y_min, real y_max,
                            real z_min, real z_max) {

    /* Allocate helper quantities */
    real* f_x = malloc(n_x*sizeof(real));
    real* f_y = malloc(n_y*sizeof(real));
    real* f_z = malloc(n_z*sizeof(real));
    real* c_x = malloc((n_x-1*(bc_x==NATURALBC))*NSIZE_EXPL1D*sizeof(real));
    real* c_y = malloc((n_y-1*(bc_y==NATURALBC))*NSIZE_EXPL1D*sizeof(real));
    real* c_z = malloc((n_z-1*(bc_z==NATURALBC))*NSIZE_EXPL1D*sizeof(real));
    int i_ct;

    if(f_x == NULL || f_y == NULL || f_z == NULL ||
       c_x == NULL || c_y == NULL || c_z == NULL) {
        return 1;
    }

    /* Calculate tricubic spline volume coefficients. For each grid cell
       (i_x, i_y, i_z), there are 64 coefficients, one for each variable
       product dx^p_x*dy^p_y*dz^p_z in the evaluation formula, where p_x, p_y,
       p_z = 0, 1, 2, 3. Note how we account for normalized grid. */

    /* Bicubic spline surfaces over xy-grid for each z */
    for(int i_z=0; i_z<n_z; i_z++) {

        /* Cubic spline along x for each y, using f values to get a total of
           four coefficients */
        for(int i_y=0; i_y<n_y; i_y++) {
            for(int i_x=0; i_x<n_x; i_x++) {
                f_x[i_x] = f[i_z*n_y*n_x+i_y*n_x+i_x];
            }
            splineexpl(f_x, n_x, bc_x, c_x);
            for(int i_x=0; i_x<n_x-1*(bc_x==NATURALBC); i_x++) {
                for(int i_c=0; i_c<4; i_c++) {
                    c[i_z*n_y*n_x*64+i_y*n_x*64+i_x*64+i_c] = c_x[i_x*4+i_c];
                }
            }
        }

        /* Four cubic splines along y for each x, using the above calulated four
           coefficient values to get a total of 16 coefficients */
        for(int i_x=0; i_x<n_x-1*(bc_x==NATURALBC); i_x++) {
            for(int i_s=0; i_s<4; i_s++) {
                for(int i_y=0; i_y<n_y; i_y++) {
                    f_y[i_y] = c[i_z*n_y*n_x*64+i_y*n_x*64+i_x*64+i_s];
                }
                splineexpl(f_y,n_y,bc_y,c_y);
                for(int i_y=0; i_y<n_y-1*(bc_y==NATURALBC); i_y++) {
                    i_ct = 0;
                    for(int i_c=i_s; i_c<16; i_c=i_c+4) {
                        c[i_z*n_y*n_x*64+i_y*n_x*64+i_x*64+i_c]
                            = c_y[i_y*4+i_ct];
                        i_ct++;
                    }
                }
            }
        }
    }

    /* Cubic splines along z for each xy-pair, using the above calculated 16
       coefficient values to get a total of 64 coefficients */
    for(int i_y=0; i_y<n_y-1*(bc_y==NATURALBC); i_y++) {
        for(int i_x=0; i_x<n_x-1*(bc_x==NATURALBC); i_x++) {
            for(int i_ss=0; i_ss<4; i_ss++) {
                for(int i_s=0; i_s<4; i_s++) {
                    for(int i_z=0; i_z<n_z; i_z++) {
                        f_z[i_z] = c[i_z*n_y*n_x*64+i_y*n_x*64
                                     +i_x*64+(i_ss*4+i_s)];
                    }
                    splineexpl(f_z,n_z,bc_z,c_z);
                    for(int i_z=0; i_z<n_z-1*(bc_z==NATURALBC); i_z++) {
                        i_ct = 0;
                        for(int i_c=4*i_ss+i_s; i_c<64; i_c=i_c+16) {
                            c[i_z*n_y*n_x*64+i_y*n_x*64+i_x*64+i_c]
                                = c_z[i_z*4+i_ct];
                            i_ct++;
                        }
                    }
                }
            }
        }
    }

    /* Free allocated memory */
    free(f_x);
    free(f_y);
    free(f_z);
    free(c_x);
    free(c_y);
    free(c_z);

    return 0;
}

/**
 * @brief Initialize a tricubic spline
 *
 * @param str pointer to spline to be initialized
 * @param c array where coefficients are stored
 * @param n_x number of data points in the x direction
 * @param n_y number of data points in the y direction
 * @param n_z number of data points in the z direction
 * @param bc_x boundary condition for x axis
 * @param bc_y boundary condition for y axis
 * @param bc_y boundary condition for z axis
 * @param x_min minimum value of the x axis
 * @param x_max maximum value of the x axis
 * @param y_min minimum value of the y axis
 * @param y_max maximum value of the y axis
 * @param z_min minimum value of the z axis
 * @param z_max maximum value of the z axis
 */
void interp3Dexpl_init_spline(interp3D_data* str, real* c,
                              int n_x, int n_y, int n_z,
                              int bc_x, int bc_y, int bc_z,
                              real x_min, real x_max,
                              real y_min, real y_max,
                              real z_min, real z_max) {

    /* Calculate grid intervals. For periodic boundary condition, grid maximum
       value and the last data point are not the same. Take this into account
       in grid intervals. */
    real x_grid = (x_max - x_min) / ( n_x - 1 * (bc_x == NATURALBC) );
    real y_grid = (y_max - y_min) / ( n_y - 1 * (bc_y == NATURALBC) );
    real z_grid = (z_max - z_min) / ( n_z - 1 * (bc_z == NATURALBC) );

    /* Initialize the interp3D_data struct */
    str->n_x    = n_x;
    str->n_y    = n_y;
    str->n_z    = n_z;
    str->bc_x   = bc_x;
    str->bc_y   = bc_y;
    str->bc_z   = bc_z;
    str->x_min  = x_min;
    str->x_max  = x_max;
    str->x_grid = x_grid;
    str->y_min  = y_min;
    str->y_max  = y_max;
    str->y_grid = y_grid;
    str->z_min  = z_min;
    str->z_max  = z_max;
    str->z_grid = z_grid;
    str->c      = c;
}

/**
 * @brief Evaluate interpolated value of 3D scalar field
 *
 * This function evaluates the interpolated value of a 3D scalar field using
 * tricubic spline interpolation coefficients of the explicit form.
 *
 * @param f variable in which to place the evaluated value
 * @param str data struct for data interpolation
 * @param x x-coordinate
 * @param y y-coordinate
 * @param z z-coordinate
 *
 * @return zero on success and one if (x,y,z) point is outside the grid.
 */
int interp3Dexpl_eval_f(real* f, interp3D_data* str, real x, real y, real z) {

    /* Make sure periodic coordinates are within [min, max] region. */
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

    /* Index for x variable. The -1 needed at exactly grid end. */
    int i_x = (x-str->x_min)/str->x_grid - 1*(x==str->x_max);
    /* Normalized x coordinate in current cell */
    real dx = (x-(str->x_min+i_x*str->x_grid))/str->x_grid;
    /* Helper varibles */
    real dx2 = dx*dx;
    real dx3 = dx2*dx;

    /* Index for y variable. The -1 needed at exactly grid end. */
    int i_y = (y-str->y_min)/str->y_grid - 1*(y==str->y_max);
    /* Normalized y coordinate in current cell */
    real dy = (y-(str->y_min+i_y*str->y_grid))/str->y_grid;
    /* Helper varibles */
    real dy2 = dy*dy;
    real dy3 = dy2*dy;

    /* Index for z variable. The -1 needed at exactly grid end. */
    int i_z = (z-str->z_min)/str->z_grid - 1*(z==str->z_max);
    /* Normalized z coordinate in current cell */
    real dz = (z-(str->z_min+i_z*str->z_grid))/str->z_grid;
    /* Helper varibles */
    real dz2 = dz*dz;
    real dz3 = dz2*dz;

    /**< Index jump to cell */
    int n = i_z*str->n_y*str->n_x*64+i_y*str->n_x*64+i_x*64;

    int err = 0;

    /* Check that the coordinate is within the domain. */
    if( str->bc_x == NATURALBC && (x < str->x_min || x > str->x_max) ) {
        err = 1;
    }
    if( str->bc_y == NATURALBC && (y < str->y_min || y > str->y_max) ) {
        err = 1;
    }
    if( str->bc_z == NATURALBC && (z < str->z_min || z > str->z_max) ) {
        err = 1;
    }

    if(!err) {

        /* Evaluate spline value */
        *f = (
                    str->c[n+ 0]+dx*str->c[n+ 1]
                    +dx2*str->c[n+ 2]+dx3*str->c[n+ 3]
                +dy*(
                    str->c[n+ 4]+dx*str->c[n+ 5]
                    +dx2*str->c[n+ 6]+dx3*str->c[n+ 7])
                +dy2*(
                    str->c[n+ 8]+dx*str->c[n+ 9]
                    +dx2*str->c[n+10]+dx3*str->c[n+11])
                +dy3*(
                    str->c[n+12]+dx*str->c[n+13]
                    +dx2*str->c[n+14]+dx3*str->c[n+15]))
            +dz*(
                    str->c[n+16]+dx*str->c[n+17]
                    +dx2*str->c[n+18]+dx3*str->c[n+19]
                +dy*(
                    str->c[n+20]+dx*str->c[n+21]
                    +dx2*str->c[n+22]+dx3*str->c[n+23])
                +dy2*(
                    str->c[n+24]+dx*str->c[n+25]
                    +dx2*str->c[n+26]+dx3*str->c[n+27])
                +dy3*(
                    str->c[n+28]+dx*str->c[n+29]
                    +dx2*str->c[n+30]+dx3*str->c[n+31]))
            +dz2*(
                    str->c[n+32]+dx*str->c[n+33]
                    +dx2*str->c[n+34]+dx3*str->c[n+35]
                +dy*(
                    str->c[n+36]+dx*str->c[n+37]
                    +dx2*str->c[n+38]+dx3*str->c[n+39])
                +dy2*(
                    str->c[n+40]+dx*str->c[n+41]
                    +dx2*str->c[n+42]+dx3*str->c[n+43])
                +dy3*(
                    str->c[n+44]+dx*str->c[n+45]
                    +dx2*str->c[n+46]+dx3*str->c[n+47]))
            +dz3*(
                    str->c[n+48]+dx*str->c[n+49]
                    +dx2*str->c[n+50]+dx3*str->c[n+51]
                +dy*(
                    str->c[n+52]+dx*str->c[n+53]
                    +dx2*str->c[n+54]+dx3*str->c[n+55])
                +dy2*(
                    str->c[n+56]+dx*str->c[n+57]
                    +dx2*str->c[n+58]+dx3*str->c[n+59])
                +dy3*(
                    str->c[n+60]+dx*str->c[n+61]
                    +dx2*str->c[n+62]+dx3*str->c[n+63]));

    }

    return err;
}

/**
 * @brief Evaluate interpolated value of 3D field and 1st and 2nd derivatives
 *
 * This function evaluates the interpolated value of a 3D scalar field and
 * its 1st and 2nd derivatives using bicubic spline interpolation coefficients
 * of the explicit form.
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
 *
 * @return zero on success and one if (x,y,z) point is outside the grid.
 */
int interp3Dexpl_eval_df(real* f_df, interp3D_data* str, real x, real y, real z) {

    /* Make sure periodic coordinates are within [min, max] region. */
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

    /* Index for z variable. The -1 needed at exactly grid end. */
    int i_z = (z-str->z_min)/str->z_grid - 1*(z==str->z_max);
    /* Normalized z coordinate in current cell */
    real dz = (z-(str->z_min+i_z*str->z_grid))/str->z_grid;
    /* Helper variables */
    real dz2 = dz*dz;
    real dz3 = dz2*dz;
    real zgi = 1.0/str->z_grid;

    /* Index jump to cell */
    int n = i_z*str->n_y*str->n_x*64+i_y*str->n_x*64+i_x*64;

    int err = 0;

    /* Check that the coordinate is within the domain. */
    if( str->bc_x == NATURALBC && (x < str->x_min || x > str->x_max) ) {
        err = 1;
    }
    if( str->bc_y == NATURALBC && (y < str->y_min || y > str->y_max) ) {
        err = 1;
    }
    if( str->bc_z == NATURALBC && (z < str->z_min || z > str->z_max) ) {
        err = 1;
    }

    if(!err) {

        /* Fetch coefficients explicitly to fetch those that are adjacent
           subsequently and to store in temporary variables coefficients that
           will be used multiple times. This is to decrease computational time,
           by exploiting simultaneous extraction of adjacent memory and by
           avoiding going through the long str->c array repeatedly. */
        real c000 = str->c[n+0];
        real c001 = str->c[n+1];
        real c002 = str->c[n+2];
        real c003 = str->c[n+3];
        real c010 = str->c[n+4];
        real c011 = str->c[n+5];
        real c012 = str->c[n+6];
        real c013 = str->c[n+7];

        real c020 = str->c[n+8];
        real c021 = str->c[n+9];
        real c022 = str->c[n+10];
        real c023 = str->c[n+11];
        real c030 = str->c[n+12];
        real c031 = str->c[n+13];
        real c032 = str->c[n+14];
        real c033 = str->c[n+15];

        real c100 = str->c[n+16];
        real c101 = str->c[n+17];
        real c102 = str->c[n+18];
        real c103 = str->c[n+19];
        real c110 = str->c[n+20];
        real c111 = str->c[n+21];
        real c112 = str->c[n+22];
        real c113 = str->c[n+23];

        real c120 = str->c[n+24];
        real c121 = str->c[n+25];
        real c122 = str->c[n+26];
        real c123 = str->c[n+27];
        real c130 = str->c[n+28];
        real c131 = str->c[n+29];
        real c132 = str->c[n+30];
        real c133 = str->c[n+31];

        real c200 = str->c[n+32];
        real c201 = str->c[n+33];
        real c202 = str->c[n+34];
        real c203 = str->c[n+35];
        real c210 = str->c[n+36];
        real c211 = str->c[n+37];
        real c212 = str->c[n+38];
        real c213 = str->c[n+39];

        real c220 = str->c[n+40];
        real c221 = str->c[n+41];
        real c222 = str->c[n+42];
        real c223 = str->c[n+43];
        real c230 = str->c[n+44];
        real c231 = str->c[n+45];
        real c232 = str->c[n+46];
        real c233 = str->c[n+47];

        real c300 = str->c[n+48];
        real c301 = str->c[n+49];
        real c302 = str->c[n+50];
        real c303 = str->c[n+51];
        real c310 = str->c[n+52];
        real c311 = str->c[n+53];
        real c312 = str->c[n+54];
        real c313 = str->c[n+55];

        real c320 = str->c[n+56];
        real c321 = str->c[n+57];
        real c322 = str->c[n+58];
        real c323 = str->c[n+59];
        real c330 = str->c[n+60];
        real c331 = str->c[n+61];
        real c332 = str->c[n+62];
        real c333 = str->c[n+63];

        /* Evaluate spline value */

        /* f */
        f_df[0] = (
                      c000+dx*c001+dx2*c002+dx3*c003
                 +dy*(c010+dx*c011+dx2*c012+dx3*c013)
                +dy2*(c020+dx*c021+dx2*c022+dx3*c023)
                +dy3*(c030+dx*c031+dx2*c032+dx3*c033))
            +dz*(
                      c100+dx*c101+dx2*c102+dx3*c103
                 +dy*(c110+dx*c111+dx2*c112+dx3*c113)
                +dy2*(c120+dx*c121+dx2*c122+dx3*c123)
                +dy3*(c130+dx*c131+dx2*c132+dx3*c133))
            +dz2*(
                      c200+dx*c201+dx2*c202+dx3*c203
                 +dy*(c210+dx*c211+dx2*c212+dx3*c213)
                +dy2*(c220+dx*c221+dx2*c222+dx3*c223)
                +dy3*(c230+dx*c231+dx2*c232+dx3*c233))
            +dz3*(
                      c300+dx*c301+dx2*c302+dx3*c303
                 +dy*(c310+dx*c311+dx2*c312+dx3*c313)
                +dy2*(c320+dx*c321+dx2*c322+dx3*c323)
                +dy3*(c330+dx*c331+dx2*c332+dx3*c333));

        /* df/dx */
        f_df[1] = xgi*(
            (
                      c001+2*dx*c002+3*dx2*c003
                 +dy*(c011+2*dx*c012+3*dx2*c013)
                +dy2*(c021+2*dx*c022+3*dx2*c023)
                +dy3*(c031+2*dx*c032+3*dx2*c033))
            +dz*(
                      c101+2*dx*c102+3*dx2*c103
                 +dy*(c111+2*dx*c112+3*dx2*c113)
                +dy2*(c121+2*dx*c122+3*dx2*c123)
                +dy3*(c131+2*dx*c132+3*dx2*c133))
            +dz2*(
                      c201+2*dx*c202+3*dx2*c203
                 +dy*(c211+2*dx*c212+3*dx2*c213)
                +dy2*(c221+2*dx*c222+3*dx2*c223)
                +dy3*(c231+2*dx*c232+3*dx2*c233))
            +dz3*(
                      c301+2*dx*c302+3*dx2*c303
                 +dy*(c311+2*dx*c312+3*dx2*c313)
                +dy2*(c321+2*dx*c322+3*dx2*c323)
                +dy3*(c331+2*dx*c332+3*dx2*c333)));

        /* df/dy */
        f_df[2] = ygi*(
            (
                        c010+dx*c011+dx2*c012+dx3*c013
                 +2*dy*(c020+dx*c021+dx2*c022+dx3*c023)
                +3*dy2*(c030+dx*c031+dx2*c032+dx3*c033))
            +dz*(
                        c110+dx*c111+dx2*c112+dx3*c113
                 +2*dy*(c120+dx*c121+dx2*c122+dx3*c123)
                +3*dy2*(c130+dx*c131+dx2*c132+dx3*c133))
            +dz2*(
                        c210+dx*c211+dx2*c212+dx3*c213
                 +2*dy*(c220+dx*c221+dx2*c222+dx3*c223)
                +3*dy2*(c230+dx*c231+dx2*c232+dx3*c233))
            +dz3*(
                        c310+dx*c311+dx2*c312+dx3*c313
                 +2*dy*(c320+dx*c321+dx2*c322+dx3*c323)
                +3*dy2*(c330+dx*c331+dx2*c332+dx3*c333)));

        /* df/dz */
        f_df[3] = zgi*(
            (
                      c100+dx*c101+dx2*c102+dx3*c103
                 +dy*(c110+dx*c111+dx2*c112+dx3*c113)
                +dy2*(c120+dx*c121+dx2*c122+dx3*c123)
                +dy3*(c130+dx*c131+dx2*c132+dx3*c133))
            +2*dz*(
                      c200+dx*c201+dx2*c202+dx3*c203
                 +dy*(c210+dx*c211+dx2*c212+dx3*c213)
                +dy2*(c220+dx*c221+dx2*c222+dx3*c223)
                +dy3*(c230+dx*c231+dx2*c232+dx3*c233))
            +3*dz2*(
                      c300+dx*c301+dx2*c302+dx3*c303
                 +dy*(c310+dx*c311+dx2*c312+dx3*c313)
                +dy2*(c320+dx*c321+dx2*c322+dx3*c323)
                +dy3*(c330+dx*c331+dx2*c332+dx3*c333)));

        /* d2f/dx2 */
        f_df[4] =  xgi*xgi*(
            (
                      2*c002+6*dx*c003
                 +dy*(2*c012+6*dx*c013)
                +dy2*(2*c022+6*dx*c023)
                +dy3*(2*c032+6*dx*c033))
            +dz*(
                      2*c102+6*dx*c103
                 +dy*(2*c112+6*dx*c113)
                +dy2*(2*c122+6*dx*c123)
                +dy3*(2*c132+6*dx*c133))
            +dz2*(
                      2*c202+6*dx*c203
                 +dy*(2*c212+6*dx*c213)
                +dy2*(2*c222+6*dx*c223)
                +dy3*(2*c232+6*dx*c233))
            +dz3*(
                      2*c302+6*dx*c303
                 +dy*(2*c312+6*dx*c313)
                +dy2*(2*c322+6*dx*c323)
                +dy3*(2*c332+6*dx*c333)));

        /* d2f/dy2 */
        f_df[5] = ygi*ygi*(
            (
                    2*(c020+dx*c021+dx2*c022+dx3*c023)
                +6*dy*(c030+dx*c031+dx2*c032+dx3*c033))
            +dz*(
                    2*(c120+dx*c121+dx2*c122+dx3*c123)
                +6*dy*(c130+dx*c131+dx2*c132+dx3*c133))
            +dz2*(
                    2*(c220+dx*c221+dx2*c222+dx3*c223)
                +6*dy*(c230+dx*c231+dx2*c232+dx3*c233))
            +dz3*(
                    2*(c320+dx*c321+dx2*c322+dx3*c323)
                +6*dy*(c330+dx*c331+dx2*c332+dx3*c333)));

        /* d2f/dz2 */
        f_df[6] = zgi*zgi*(
            2*(
                      c200+dx*c201+dx2*c202+dx3*c203
                 +dy*(c210+dx*c211+dx2*c212+dx3*c213)
                +dy2*(c220+dx*c221+dx2*c222+dx3*c223)
                +dy3*(c230+dx*c231+dx2*c232+dx3*c233))
            +6*dz*(
                      c300+dx*c301+dx2*c302+dx3*c303
                 +dy*(c310+dx*c311+dx2*c312+dx3*c313)
                +dy2*(c320+dx*c321+dx2*c322+dx3*c323)
                +dy3*(c330+dx*c331+dx2*c332+dx3*c333)));

        /* d2f/dydx */
        f_df[7] =  xgi*ygi*(
            (
                        c011+2*dx*c012+3*dx2*c013
                 +2*dy*(c021+2*dx*c022+3*dx2*c023)
                +3*dy2*(c031+2*dx*c032+3*dx2*c033))
            +dz*(
                        c111+2*dx*c112+3*dx2*c113
                 +2*dy*(c121+2*dx*c122+3*dx2*c123)
                +3*dy2*(c131+2*dx*c132+3*dx2*c133))
            +dz2*(
                        c211+2*dx*c212+3*dx2*c213
                 +2*dy*(c221+2*dx*c222+3*dx2*c223)
                +3*dy2*(c231+2*dx*c232+3*dx2*c233))
            +dz3*(
                        c311+2*dx*c312+3*dx2*c313
                 +2*dy*(c321+2*dx*c322+3*dx2*c323)
                +3*dy2*(c331+2*dx*c332+3*dx2*c333)));

        /* d2f/dzdx */
        f_df[8] = xgi*zgi*(
            (
                      c101+2*dx*c102+3*dx2*c103
                 +dy*(c111+2*dx*c112+3*dx2*c113)
                +dy2*(c121+2*dx*c122+3*dx2*c123)
                +dy3*(c131+2*dx*c132+3*dx2*c133))
            +2*dz*(
                      c201+2*dx*c202+3*dx2*c203
                 +dy*(c211+2*dx*c212+3*dx2*c213)
                +dy2*(c221+2*dx*c222+3*dx2*c223)
                +dy3*(c231+2*dx*c232+3*dx2*c233))
            +3*dz2*(
                      c301+2*dx*c302+3*dx2*c303
                 +dy*(c311+2*dx*c312+3*dx2*c313)
                +dy2*(c321+2*dx*c322+3*dx2*c323)
                +dy3*(c331+2*dx*c332+3*dx2*c333)));

        /* d2f/dzdy */
        f_df[9] = ygi*zgi*(
            (
                        c110+dx*c111+dx2*c112+dx3*c113
                 +2*dy*(c120+dx*c121+dx2*c122+dx3*c123)
                +3*dy2*(c130+dx*c131+dx2*c132+dx3*c133))
            +2*dz*(
                        c210+dx*c211+dx2*c212+dx3*c213
                 +2*dy*(c220+dx*c221+dx2*c222+dx3*c223)
                +3*dy2*(c230+dx*c231+dx2*c232+dx3*c233))
            +3*dz2*(
                        c310+dx*c311+dx2*c312+dx3*c313
                 +2*dy*(c320+dx*c321+dx2*c322+dx3*c323)
                +3*dy2*(c330+dx*c331+dx2*c332+dx3*c333)));
    }

    return err;
}
