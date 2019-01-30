/**
 * @file interp3Dcomp.c
 * @brief Tricubic spline interpolation in compact form
 */
#include <stdlib.h>
#include <stdio.h> /* Needed for printf debugging purposes */
#include <math.h>
#include "../ascot5.h"
#include "../consts.h"
#include "interp3Dcomp.h"
#include "spline1Dcomp.h"

/**
 * @brief Calculate tricubic spline interpolation coefficients for scalar 3D data
 *
 * This function calculates the tricubic spline interpolation coefficients for
 * the given data and stores them in the data struct. Compact  cofficients are
 * calculated directly.
 *
 * @todo Error checking
 *
 * @param str data struct for data interpolation
 * @param f 3D data to be interpolated
 * @param n_x number of data points in the x direction
 * @param n_y number of data points in the y direction
 * @param n_z number of data points in the z direction
 * @param x_min minimum value of the x axis
 * @param x_max maximum value of the x axis
 * @param y_min minimum value of the y axis
 * @param y_max maximum value of the y axis
 * @param z_min minimum value of the z axis
 * @param z_max maximum value of the z axis
 */
int interp3Dcomp_init(interp3D_data* str, real* f, int n_x, int n_y, int n_z,
                       real x_min, real x_max, real x_grid,
                       real y_min, real y_max, real y_grid,
                       real z_min, real z_max, real z_grid) {

    int err = 0;

    /* Initialize and fill the data struct */
    str->n_x = n_x;
    str->n_y = n_y;
    str->n_z = n_z;
    str->x_min = x_min;
    str->x_max = x_max;
    str->x_grid = x_grid;
    str->y_min = y_min;
    str->y_max = y_max;
    str->y_grid = y_grid;
    str->z_min = z_min;
    str->z_max = z_max;
    str->z_grid = z_grid;
    str->c = malloc(n_y*n_z*n_x*8*sizeof(real));

    /* Declare and allocate the needed variables */
    int i_x;                                    /**< index for x variable */
    int i_y;                                  /**< index for y variable */
    int i_z;                                    /**< index for z variable */
    real* f_x = malloc(n_x*sizeof(real));       /**< Temporary array for data along x */
    real* f_y = malloc(n_y*sizeof(real));   /**< Temp array for data along y */
    real* f_z = malloc(n_z*sizeof(real));       /**< Temp array for data along z */
    real* c_x = malloc(n_x*2*sizeof(real));     /**< Temp array for coefficients along x */
    real* c_y = malloc(n_y*2*sizeof(real)); /**< Temp array for coefs along y */
    real* c_z = malloc(n_z*2*sizeof(real));     /**< Temp array for coefficients along z */

    if(f_x == NULL || f_z == NULL || c_x == NULL || c_z == NULL) {
        err = 1;
    }
    else {

        /* Tricubic spline volume coefficients: For i_x, i_y and i_z the 8 coefficients
           are [f, fxx, fyy, fzz, fxxyy, fxxzz, fyyzz, fxxyyzz].
           Note how we account for normalized grid. */
        /* Bicubic spline surface over xz-grid for each y */
        for(i_y=0; i_y<n_y; i_y++) {
            /* Cubic spline along x for each z to get fxx */
            for(i_z=0; i_z<n_z; i_z++) {
                for(i_x=0; i_x<n_x; i_x++) {
                    f_x[i_x] = f[i_y*n_z*n_x+i_z*n_x+i_x];
                }
                spline1Dcomp(f_x,n_x,0,c_x);
                for(i_x=0; i_x<n_x; i_x++) {
                    str->c[i_y*n_z*n_x*8+i_z*n_x*8+i_x*8] = c_x[i_x*2];
                    str->c[i_y*n_z*n_x*8+i_z*n_x*8+i_x*8+1] = c_x[i_x*2+1]/(x_grid*x_grid);
                }
            }

            /* Two cubic splines along z for each x using f and fxx */
            for(i_x=0; i_x<n_x; i_x++) {
                /* fzz */
                for(i_z=0; i_z<n_z; i_z++) {
                    f_z[i_z] =  f[i_y*n_z*n_x+i_z*n_x+i_x];
                }
                spline1Dcomp(f_z,n_z,0,c_z);
                for(i_z=0; i_z<n_z; i_z++) {
                    str->c[i_y*n_z*n_x*8+i_z*n_x*8+i_x*8+3] = c_z[i_z*2+1]/(z_grid*z_grid);
                }
                /* fxxzz */
                for(i_z=0; i_z<n_z; i_z++) {
                    f_z[i_z] = str->c[i_y*n_z*n_x*8+i_z*n_x*8+i_x*8+1];
                }
                spline1Dcomp(f_z,n_z,0,c_z);
                for(i_z=0; i_z<n_z; i_z++) {
                    str->c[i_y*n_z*n_x*8+i_z*n_x*8+i_x*8+5] = c_z[i_z*2+1]/(z_grid*z_grid);
                }
            }

        }

        /* Cubic spline along y for each xz-pair to find the compact coefficients
           of the tricubic spline volume */
        for(i_z=0; i_z<n_z; i_z++) {
            for(i_x=0; i_x<n_x; i_x++) {
                /* fyy */
                for(i_y=0; i_y<n_y; i_y++) {
                    f_y[i_y] = f[i_y*n_z*n_x+i_z*n_x+i_x];
                }
                spline1Dcomp(f_y,n_y,1,c_y);
                for(i_y=0; i_y<n_y; i_y++) {
                    str->c[i_y*n_z*n_x*8+i_z*n_x*8+i_x*8+2] = c_y[i_y*2+1]/(y_grid*y_grid);
                }
                /* fxxyy */
                for(i_y=0; i_y<n_y; i_y++) {
                    f_y[i_y] = str->c[i_y*n_z*n_x*8+i_z*n_x*8+i_x*8+1];
                }
                spline1Dcomp(f_y,n_y,1,c_y);
                for(i_y=0; i_y<n_y; i_y++) {
                    str->c[i_y*n_z*n_x*8+i_z*n_x*8+i_x*8+4] = c_y[i_y*2+1]/(y_grid*y_grid);
                }
                /* fyyzz */
                for(i_y=0; i_y<n_y; i_y++) {
                    f_y[i_y] = str->c[i_y*n_z*n_x*8+i_z*n_x*8+i_x*8+3];
                }
                spline1Dcomp(f_y,n_y,1,c_y);
                for(i_y=0; i_y<n_y; i_y++) {
                    str->c[i_y*n_z*n_x*8+i_z*n_x*8+i_x*8+6] = c_y[i_y*2+1]/(y_grid*y_grid);
                }
                /* fxxyyzz */
                for(i_y=0; i_y<n_y; i_y++) {
                    f_y[i_y] = str->c[i_y*n_z*n_x*8+i_z*n_x*8+i_x*8+5];
                }
                spline1Dcomp(f_y,n_y,1,c_y);
                for(i_y=0; i_y<n_y; i_y++) {
                    str->c[i_y*n_z*n_x*8+i_z*n_x*8+i_x*8+7] = c_y[i_y*2+1]/(y_grid*y_grid);
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

    return err;
}

/**
 * @brief Evaluate interpolated value of 3D scalar field
 *
 * This function evaluates the interpolated value of a 3D scalar field using
 * tricubic spline interpolation coefficients of the compact form.
 *
 * @todo Seg fault if in last y-sector becaause of compact eval +1-coefs
 * @todo Check discrepency to ascot4 and explicit version
 * @todo Error checking
 *
 * @param f variable in which to place the evaluated value
 * @param str data struct for data interpolation
 * @param x x-coordinate
 * @param y y-coordinate
 * @param z z-coordinate
 */
integer interp3Dcomp_eval_f(real* f, interp3D_data* str, real x, real y, real z) {
    /** Make sure y is in interval [y_min, y_max) */
    real y_range = str->y_max - str->y_min;
    y = fmod(y - str->y_min, y_range) + str->y_min;
    if(y < 0){y = y_range + y;}

    int i_x = (x-str->x_min)/str->x_grid;     /**< index for x variable */
    real dx = (x-(str->x_min+i_x*str->x_grid))/str->x_grid; /**< Normalized x coordinate in
                                                               current cell */
    real dxi = 1.0-dx;
    real dx3 = dx*dx*dx-dx;
    real dxi3 = (1.0-dx)*(1.0-dx)*(1.0-dx)-(1.0-dx);
    real xg2 = str->x_grid*str->x_grid;       /**< Square of cell length in x direction */
    int i_y = (y-str->y_min)/str->y_grid; /**< index for y variable */
    real dy = (y-(str->y_min+i_y*str->y_grid))/str->y_grid; /**< Normalized y
                                                                           coordinate in
                                                                           current cell */
    real dyi = 1.0-dy;
    real dy3 = dy*dy*dy-dy;
    real dyi3 = (1.0-dy)*(1.0-dy)*(1.0-dy)-(1.0-dy);
    real yg2 = str->y_grid*str->y_grid; /**< Square of cell length in y dir */
    int i_z = (z-str->z_min)/str->z_grid;     /**< index for z variable */
    real dz = (z-(str->z_min+i_z*str->z_grid))/str->z_grid; /**< Normalized z coordinate in
                                                               current cell */
    real dzi = 1.0-dz;
    real dz3 = dz*dz*dz-dz;
    real dzi3 = (1.0-dz)*(1.0-dz)*(1.0-dz)-(1.0-dz);
    real zg2 = str->z_grid*str->z_grid;       /**< Square of cell length in z direction */
    int n = i_y*str->n_z*str->n_x*8+i_z*str->n_x*8+i_x*8; /**< Index jump to cell */
    int x1 = 8;                               /**< Index jump one x forward */
    int y1 = str->n_z*str->n_x*8;           /**< Index jump one y forward */
    if(i_y==str->n_y-1) {
        y1 = -(str->n_y-1)*y1;          /**< If last cell, index jump to 1st y */
    }
    int z1 = str->n_x*8;                      /**< Index jump one z forward */

    int err = 0;

    /* Check that the point is not outside the evaluation regime */
    if(x < str->x_min || x > str->x_max
        || z < str->z_min || z > str->z_max) {
        err = 1;
    }
    else {

        *f = (
            dzi*(
                dxi*(dyi*str->c[n+0]+dy*str->c[n+y1+0])+
                dx*(dyi*str->c[n+x1+0]+dy*str->c[n+y1+x1+0]))
            +dz*(
                dxi*(dyi*str->c[n+z1+0]+dy*str->c[n+y1+z1+0])+
                dx*(dyi*str->c[n+x1+z1+0]+dy*str->c[n+y1+z1+x1+0])))
            +xg2/6*(
                dzi*(
                    dxi3*(dyi*str->c[n+1]+dy*str->c[n+y1+1])+
                    dx3*(dyi*str->c[n+x1+1]+dy*str->c[n+y1+x1+1]))
                +dz*(
                    dxi3*(dyi*str->c[n+z1+1]+dy*str->c[n+y1+z1+1])+
                    dx3*(dyi*str->c[n+x1+z1+1]+dy*str->c[n+y1+z1+x1+1])))
            +yg2/6*(
                dzi*(
                    dxi*(dyi3*str->c[n+2]+dy3*str->c[n+y1+2])+
                    dx*(dyi3*str->c[n+x1+2]+dy3*str->c[n+y1+x1+2]))
                +dz*(
                    dxi*(dyi3*str->c[n+z1+2]+dy3*str->c[n+y1+z1+2])+
                    dx*(dyi3*str->c[n+x1+z1+2]+dy3*str->c[n+y1+z1+x1+2])))
            +zg2/6*(
                dzi3*(
                    dxi*(dyi*str->c[n+3]+dy*str->c[n+y1+3])+
                    dx*(dyi*str->c[n+x1+3]+dy*str->c[n+y1+x1+3]))
                +dz3*(
                    dxi*(dyi*str->c[n+z1+3]+dy*str->c[n+y1+z1+3])+
                    dx*(dyi*str->c[n+x1+z1+3]+dy*str->c[n+y1+z1+x1+3])))
            +xg2*yg2/36*(
                dzi*(
                    dxi3*(dyi3*str->c[n+4]+dy3*str->c[n+y1+4])+
                    dx3*(dyi3*str->c[n+x1+4]+dy3*str->c[n+y1+x1+4]))
                +dz*(
                    dxi3*(dyi3*str->c[n+z1+4]+dy3*str->c[n+y1+z1+4])+
                    dx3*(dyi3*str->c[n+x1+z1+4]+dy3*str->c[n+y1+z1+x1+4])))
            +xg2*zg2/36*(
                dzi3*(
                    dxi3*(dyi*str->c[n+5]+dy*str->c[n+y1+5])+
                    dx3*(dyi*str->c[n+x1+5]+dy*str->c[n+y1+x1+5]))
                +dz3*(
                    dxi3*(dyi*str->c[n+z1+5]+dy*str->c[n+y1+z1+5])+
                    dx3*(dyi*str->c[n+x1+z1+5]+dy*str->c[n+y1+z1+x1+5])))
            +yg2*zg2/36*(
                dzi3*(
                    dxi*(dyi3*str->c[n+6]+dy3*str->c[n+y1+6])+
                    dx*(dyi3*str->c[n+x1+6]+dy3*str->c[n+y1+x1+6]))
                +dz3*(
                    dxi*(dyi3*str->c[n+z1+6]+dy3*str->c[n+y1+z1+6])+
                    dx*(dyi3*str->c[n+x1+z1+6]+dy3*str->c[n+y1+z1+x1+6])))
            +xg2*yg2*zg2/216*(
                dzi3*(
                    dxi3*(dyi3*str->c[n+7]+dy3*str->c[n+y1+7])+
                    dx3*(dyi3*str->c[n+x1+7]+dy3*str->c[n+y1+x1+7]))
                +dz3*(
                    dxi3*(dyi3*str->c[n+z1+7]+dy3*str->c[n+y1+z1+7])+
                    dx3*(dyi3*str->c[n+x1+z1+7]+dy3*str->c[n+y1+z1+x1+7])));

    }

    return err;
}

/**
 * @brief Evaluate interpolated value of 3D scalar field and its 1st and 2nd derivatives
 *
 * This function evaluates the interpolated value of a 3D scalar field and
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
 * @todo Seg fault if in last y-sector becaause of compact eval +1-coefs
 * @todo Check discrepency to ascot4 and explicit version
 * @todo Error checking
 *
 * @param B_dB array in which to place the evaluated values
 * @param str data struct for data interpolation
 * @param x x-coordinate
 * @param y y-coordinate
 * @param z z-coordinate
 */
integer interp3Dcomp_eval_df(real* f_df, interp3D_data* str, real x, real y, real z) {
    /** Make sure y is in interval [y_min, y_max) */
    real y_range = str->y_max - str->y_min;
    y = fmod(y - str->y_min, y_range) + str->y_min;
    if(y < 0){y = y_range + y;}

    int i_x = (x-str->x_min)/str->x_grid;       /**< index for x variable */
    real dx = (x-(str->x_min+i_x*str->x_grid))/str->x_grid; /**< Normalized x coordinate in
                                                               current cell */
    real dx3 = dx*dx*dx-dx;
    real dx3dx = 3*dx*dx-1.0;           /**< x-derivative of dx3, not including 1/x_grid */
    real dxi = 1.0-dx;
    real dxi3 = dxi*dxi*dxi-dxi;
    real dxi3dx = -3*dxi*dxi+1.0;      /**< x-derivative of dxi3, not including 1/x_grid */
    real xg = str->x_grid;                      /**< Cell length in x direction */
    real xg2 = xg*xg;
    real xgi = 1.0/xg;
    int i_y = (y-str->y_min)/str->y_grid; /**< index for y variable */
    real dy = (y-(str->y_min+i_y*str->y_grid))/str->y_grid; /**< Normalized y
                                                                           coordinate in
                                                                           current cell */
    real dy3 = dy*dy*dy-dy;
    real dy3dy = 3*dy*dy-1.0; /**< y-derivative of dy3,
                                         not including 1/y_grid */
    real dyi = 1.0-dy;
    real dyi3 = dyi*dyi*dyi-dyi;
    real dyi3dy = -3*dyi*dyi+1.0; /**< y-derivative of dyi3,
                                             not including 1/y_grid */
    real yg = str->y_grid;                  /**< Cell length in y direction */
    real yg2 = yg*yg;
    real ygi = 1.0/yg;
    int i_z = (z-str->z_min)/str->z_grid;       /**< index for z variable */
    real dz = (z-(str->z_min+i_z*str->z_grid))/str->z_grid; /**< Normalized z coordinate in
                                                               current cell */
    real dz3 = dz*dz*dz-dz;
    real dz3dz = 3*dz*dz-1.0;          /**< z-derivative of dz3, not including 1/z_grid */
    real dzi = 1.0-dz;
    real dzi3 = dzi*dzi*dzi-dzi;
    real dzi3dz = -3*dzi*dzi+1.0;      /**< z-derivative of dzi3, not including 1/z_grid */
    real zg = str->z_grid;                      /**< Cell length in z direction */
    real zg2 = zg*zg;
    real zgi = 1.0/zg;
    int n = i_y*str->n_z*str->n_x*8+i_z*str->n_x*8+i_x*8; /**< Index jump to cell */
    int x1 = 8;                                 /**< Index jump one x forward */
    int y1 = str->n_z*str->n_x*8;             /**< Index jump one y forward */
    if(i_y==str->n_y-1) {
        y1 = -(str->n_y-1)*y1;            /**< If last cell, index jump to 1st y */
    }
    int z1 = str->n_x*8;                        /**< Index jump one z forward */


    int err = 0;

    /* Check that the point is not outside the evaluation regime */
    if(x < str->x_min || x > str->x_max
        || z < str->z_min || z > str->z_max) {
        err = 1;
    }
    else {

        real c0000 = str->c[n+0];
        real c0001 = str->c[n+1];
        real c0002 = str->c[n+2];
        real c0003 = str->c[n+3];
        real c0004 = str->c[n+4];
        real c0005 = str->c[n+5];
        real c0006 = str->c[n+6];
        real c0007 = str->c[n+7];

        real c0010 = str->c[n+x1+0];
        real c0011 = str->c[n+x1+1];
        real c0012 = str->c[n+x1+2];
        real c0013 = str->c[n+x1+3];
        real c0014 = str->c[n+x1+4];
        real c0015 = str->c[n+x1+5];
        real c0016 = str->c[n+x1+6];
        real c0017 = str->c[n+x1+7];

        real c0100 = str->c[n+y1+0];
        real c0101 = str->c[n+y1+1];
        real c0102 = str->c[n+y1+2];
        real c0103 = str->c[n+y1+3];
        real c0104 = str->c[n+y1+4];
        real c0105 = str->c[n+y1+5];
        real c0106 = str->c[n+y1+6];
        real c0107 = str->c[n+y1+7];

        real c1000 = str->c[n+z1+0];
        real c1001 = str->c[n+z1+1];
        real c1002 = str->c[n+z1+2];
        real c1003 = str->c[n+z1+3];
        real c1004 = str->c[n+z1+4];
        real c1005 = str->c[n+z1+5];
        real c1006 = str->c[n+z1+6];
        real c1007 = str->c[n+z1+7];

        real c0110 = str->c[n+y1+x1+0];
        real c0111 = str->c[n+y1+x1+1];
        real c0112 = str->c[n+y1+x1+2];
        real c0113 = str->c[n+y1+x1+3];
        real c0114 = str->c[n+y1+x1+4];
        real c0115 = str->c[n+y1+x1+5];
        real c0116 = str->c[n+y1+x1+6];
        real c0117 = str->c[n+y1+x1+7];

        real c1100 = str->c[n+y1+z1+0];
        real c1101 = str->c[n+y1+z1+1];
        real c1102 = str->c[n+y1+z1+2];
        real c1103 = str->c[n+y1+z1+3];
        real c1104 = str->c[n+y1+z1+4];
        real c1105 = str->c[n+y1+z1+5];
        real c1106 = str->c[n+y1+z1+6];
        real c1107 = str->c[n+y1+z1+7];

        real c1010 = str->c[n+x1+z1+0];
        real c1011 = str->c[n+x1+z1+1];
        real c1012 = str->c[n+x1+z1+2];
        real c1013 = str->c[n+x1+z1+3];
        real c1014 = str->c[n+x1+z1+4];
        real c1015 = str->c[n+x1+z1+5];
        real c1016 = str->c[n+x1+z1+6];
        real c1017 = str->c[n+x1+z1+7];

        real c1110 = str->c[n+x1+y1+z1+0];
        real c1111 = str->c[n+x1+y1+z1+1];
        real c1112 = str->c[n+x1+y1+z1+2];
        real c1113 = str->c[n+x1+y1+z1+3];
        real c1114 = str->c[n+x1+y1+z1+4];
        real c1115 = str->c[n+x1+y1+z1+5];
        real c1116 = str->c[n+x1+y1+z1+6];
        real c1117 = str->c[n+x1+y1+z1+7];


        /* f */
        f_df[0] = (
               dzi*(
                   dxi*(dyi*c0000+dy*c0100)+
                   dx*(dyi*c0010+dy*c0110))
               +dz*(
                   dxi*(dyi*c1000+dy*c1100)+
                   dx*(dyi*c1010+dy*c1110)))
        +xg2/6*(
            dzi*(
                dxi3*(dyi*c0001+dy*c0101)+
                dx3*(dyi*c0011+dy*c0111))
            +dz*(
                dxi3*(dyi*c1001+dy*c1101)+
                dx3*(dyi*c1011+dy*c1111)))
        +yg2/6*(
            dzi*(
                dxi*(dyi3*c0002+dy3*c0102)+
                dx*(dyi3*c0012+dy3*c0112))
            +dz*(
                dxi*(dyi3*c1002+dy3*c1102)+
                dx*(dyi3*c1012+dy3*c1112)))
        +zg2/6*(
            dzi3*(
                dxi*(dyi*c0003+dy*c0103)+
                dx*(dyi*c0013+dy*c0113))
            +dz3*(
                dxi*(dyi*c1003+dy*c1103)+
                dx*(dyi*c1013+dy*c1113)))
        +xg2*yg2/36*(
            dzi*(
                dxi3*(dyi3*c0004+dy3*c0104)+
                dx3*(dyi3*c0014+dy3*c0114))
            +dz*(
                dxi3*(dyi3*c1004+dy3*c1104)+
                dx3*(dyi3*c1014+dy3*c1114)))
        +xg2*zg2/36*(
            dzi3*(
                dxi3*(dyi*c0005+dy*c0105)+
                dx3*(dyi*c0015+dy*c0115))
            +dz3*(
                dxi3*(dyi*c1005+dy*c1105)+
                dx3*(dyi*c1015+dy*c1115)))
        +yg2*zg2/36*(
            dzi3*(
                dxi*(dyi3*c0006+dy3*c0106)+
                dx*(dyi3*c0016+dy3*c0116))
            +dz3*(
                dxi*(dyi3*c1006+dy3*c1106)+
                dx*(dyi3*c1016+dy3*c1116)))
        +xg2*yg2*zg2/216*(
            dzi3*(
                dxi3*(dyi3*c0007+dy3*c0107)+
                dx3*(dyi3*c0017+dy3*c0117))
            +dz3*(
                dxi3*(dyi3*c1007+dy3*c1107)+
                dx3*(dyi3*c1017+dy3*c1117)));

    /* df/dx */
    f_df[1] = xgi*(
        dzi*(
            -(dyi*c0000+dy*c0100)
            +(dyi*c0010+dy*c0110))
        +dz*(
            -(dyi*c1000+dy*c1100)
            +(dyi*c1010+dy*c1110)))
        +xg/6*(
            dzi*(
                dxi3dx*(dyi*c0001+dy*c0101)+
                dx3dx*(dyi*c0011+dy*c0111))
            +dz*(
                dxi3dx*(dyi*c1001  +dy*c1101)+
                dx3dx*(dyi*c1011+dy*c1111)))
        +xgi*yg2/6*(
            dzi*(
                -(dyi3*c0002+dy3*c0102)
                +(dyi3*c0012+dy3*c0112))
            +dz*(
                -(dyi3*c1002+dy3*c1102)
                +(dyi3*c1012+dy3*c1112)))
        +xgi*zg2/6*(
            dzi3*(
                -(dyi*c0003+dy*c0103)
                +(dyi*c0013+dy*c0113))
            +dz3*(
                -(dyi*c1003+dy*c1103)
                +(dyi*c1013+dy*c1113)))
        +xg*yg2/36*(
            dzi*(
                dxi3dx*(dyi3*c0004+dy3*c0104)+
                dx3dx*(dyi3*c0014+dy3*c0114))
            +dz*(
                dxi3dx*(dyi3*c1004+dy3*c1104)+
                dx3dx*(dyi3*c1014+dy3*c1114)))
        +xg*zg2/36*(
            dzi3*(
                dxi3dx*(dyi*c0005+dy*c0105)+
                dx3dx*(dyi*c0015+dy*c0115))
            +dz3*(
                dxi3dx*(dyi*c1005+dy*c1105)+
                dx3dx*(dyi*c1015+dy*c1115)))
        +xgi*yg2*zg2/36*(
            dzi3*(
                -(dyi3*c0006+dy3*c0106)
                +(dyi3*c0016+dy3*c0116))
            +dz3*(
                -(dyi3*c1006+dy3*c1106)
                +(dyi3*c1016+dy3*c1116)))
        +xg*yg2*zg2/216*(
            dzi3*(
                dxi3dx*(dyi3*c0007+dy3*c0107)+
                dx3dx*(dyi3*c0017+dy3*c0117))
            +dz3*(
                dxi3dx*(dyi3*c1007+dy3*c1107)+
                dx3dx*(dyi3*c1017+dy3*c1117)));

    /* df/dy */
    f_df[2] = ygi*(
        dzi*(
            dxi*(-c0000+c0100)+
            dx*(-c0010+c0110))
        +dz*(
            dxi*(-c1000+c1100)+
            dx*(-c1010+c1110)))
        +ygi*xg2/6*(
            dzi*(
                dxi3*(-c0001+c0101)+
                dx3*(-c0011+c0111))
            +dz*(
                dxi3*(-c1001+c1101)+
                dx3*(-c1011+c1111)))
        +yg/6*(
            dzi*(
                dxi*(dyi3dy*c0002+dy3dy*c0102)+
                dx*(dyi3dy*c0012+dy3dy*c0112))
            +dz*(
                dxi*(dyi3dy*c1002+dy3dy*c1102)+
                dx*(dyi3dy*c1012+dy3dy*c1112)))
        +ygi*zg2/6*(
            dzi3*(
                dxi*(-c0003+c0103)+
                dx*(-c0013+c0113))
            +dz3*(
                dxi*(-c1003+c1103)+
                dx*(-c1013+c1113)))
        +xg2*yg/36*(
            dzi*(
                dxi3*(dyi3dy*c0004+dy3dy*c0104)+
                dx3*(dyi3dy*c0014+dy3dy*c0114))
            +dz*(
                dxi3*(dyi3dy*c1004+dy3dy*c1104)+
                dx3*(dyi3dy*c1014+dy3dy*c1114)))
        +ygi*xg2*zg2/36*(
            dzi3*(
                dxi3*(-c0005+c0105)+
                dx3*(-c0015+c0115))
            +dz3*(
                dxi3*(-c1005+c1105)+
                dx3*(-c1015+c1115)))
        +yg*zg2/36*(
            dzi3*(
                dxi*(dyi3dy*c0006+dy3dy*c0106)+
                dx*(dyi3dy*c0016+dy3dy*c0116))
            +dz3*(
                dxi*(dyi3dy*c1006+dy3dy*c1106)+
                dx*(dyi3dy*c1016+dy3dy*c1116)))
        +xg2*yg*zg2/216*(
            dzi3*(
                dxi3*(dyi3dy*c0007+dy3dy*c0107)+
                dx3*(dyi3dy*c0017+dy3dy*c0117))
            +dz3*(
                dxi3*(dyi3dy*c1007+dy3dy*c1107)+
                dx3*(dyi3dy*c1017+dy3dy*c1117)));

    /* df/dz */
    f_df[3] = zgi*(
        -(
            dxi*(dyi*c0000+dy*c0100)+
            dx*(dyi*c0010+dy*c0110))
        +(
            dxi*(dyi*c1000+dy*c1100)+
            dx*(dyi*c1010+dy*c1110)))
        +xg2*zgi/6*(
            -(
                dxi3*(dyi*c0001+dy*c0101)+
                dx3*(dyi*c0011+dy*c0111))
            +(
                dxi3*(dyi*c1001+dy*c1101)+
                dx3*(dyi*c1011+dy*c1111)))
        +yg2*zgi/6*(
            -(
                dxi*(dyi3*c0002+dy3*c0102)+
                dx*(dyi3*c0012+dy3*c0112))
            +(
                dxi*(dyi3*c1002+dy3*c1102)+
                dx*(dyi3*c1012+dy3*c1112)))
        +zg/6*(
            dzi3dz*(
                dxi*(dyi*c0003+dy*c0103)+
                dx*(dyi*c0013+dy*c0113))
            +dz3dz*(
                dxi*(dyi*c1003+dy*c1103)+
                dx*(dyi*c1013+dy*c1113)))
        +xg2*yg2*zgi/36*(
            -(
                dxi3*(dyi3*c0004+dy3*c0104)+
                dx3*(dyi3*c0014+dy3*c0114))
            +(
                dxi3*(dyi3*c1004+dy3*c1104)+
                dx3*(dyi3*c1014+dy3*c1114)))
        +xg2*zg/36*(
            dzi3dz*(
                dxi3*(dyi*c0005+dy*c0105)+
                dx3*(dyi*c0015+dy*c0115))
            +dz3dz*(
                dxi3*(dyi*c1005+dy*c1105)+
                dx3*(dyi*c1015+dy*c1115)))
        +yg2*zg/36*(
            dzi3dz*(
                dxi*(dyi3*c0006+dy3*c0106)+
                dx*(dyi3*c0016+dy3*c0116))
            +dz3dz*(
                dxi*(dyi3*c1006+dy3*c1106)+
                dx*(dyi3*c1016+dy3*c1116)))
        +xg2*yg2*zg/216*(
            dzi3dz*(
                dxi3*(dyi3*c0007+dy3*c0107)+
                dx3*(dyi3*c0017+dy3*c0117))
            +dz3dz*(
                dxi3*(dyi3*c1007+dy3*c1107)+
                dx3*(dyi3*c1017+dy3*c1117)));

    /* d2f/dx2 */
    f_df[4] = (
        dzi*(
            dxi*(dyi*c0001+dy*c0101)+
            dx*(dyi*c0011+dy*c0111))
        +dz*(
            dxi*(dyi*c1001+dy*c1101)+
            dx*(dyi*c1011+dy*c1111)))
        +yg2/6*(
            dzi*(
                dxi*(dyi3*c0004+dy3*c0104)+
                dx*(dyi3*c0014+dy3*c0114))
            +dz*(
                dxi*(dyi3*c1004+dy3*c1104)+
                dx*(dyi3*c1014+dy3*c1114)))
        +zg2/6*(
            dzi3*(
                dxi*(dyi*c0005+dy*c0105)+
                dx*(dyi*c0015+dy*c0115))
            +dz3*(
                dxi*(dyi*c1005+dy*c1105)+
                dx*(dyi*c1015+dy*c1115)))
        +yg2*zg2/36*(
            dzi3*(
                dxi*(dyi3*c0007+dy3*c0107)+
                dx*(dyi3*c0017+dy3*c0117))
            +dz3*(
                dxi*(dyi3*c1007+dy3*c1107)+
                dx*(dyi3*c1017+dy3*c1117)));

    /* d2f/dy2 */
    f_df[5] = (
        dzi*(
            dxi*(dyi*c0002+dy*c0102)+
            dx*(dyi*c0012+dy*c0112))
        +dz*(
            dxi*(dyi*c1002+dy*c1102)+
            dx*(dyi*c1012+dy*c1112)))
        +xg2/6*(
            dzi*(
                dxi3*(dyi*c0004+dy*c0104)+
                dx3*(dyi*c0014+dy*c0114))
            +dz*(
                dxi3*(dyi*c1004+dy*c1104)+
                dx3*(dyi*c1014+dy*c1114)))
        +zg2/6*(
            dzi3*(
                dxi*(dyi*c0006+dy*c0106)+
                dx*(dyi*c0016+dy*c0116))
            +dz3*(
                dxi*(dyi*c1006+dy*c1106)+
                dx*(dyi*c1016+dy*c1116)))
        +xg2*zg2/36*(
            dzi3*(
                dxi3*(dyi*c0007+dy*c0107)+
                dx3*(dyi*c0017+dy*c0117))
            +dz3*(
                dxi3*(dyi*c1007+dy*c1107)+
                dx3*(dyi*c1017+dy*c1117)));

    /* d2f/dz2 */
    f_df[6] = (
        dzi*(
            dxi*(dyi*c0003+dy*c0103)+
            dx*(dyi*c0013+dy*c0113))
        +dz*(
            dxi*(dyi*c1003+dy*c1103)+
            dx*(dyi*c1013+dy*c1113)))
        +xg2/6*(
            dzi*(
                dxi3*(dyi*c0005+dy*c0105)+
                dx3*(dyi*c0015+dy*c0115))
            +dz*(
                dxi3*(dyi*c1005+dy*c1105)+
                dx3*(dyi*c1015+dy*c1115)))
        +yg2/6*(
            dzi*(
                dxi*(dyi3*c0006+dy3*c0106)+
                dx*(dyi3*c0016+dy3*c0116))
            +dz*(
                dxi*(dyi3*c1006+dy3*c1106)+
                dx*(dyi3*c1016+dy3*c1116)))
        +xg2*yg2/36*(
            dzi*(
                dxi3*(dyi3*c0007+dy3*c0107)+
                dx3*(dyi3*c0017+dy3*c0117))
            +dz*(
                dxi3*(dyi3*c1007+dy3*c1107)+
                dx3*(dyi3*c1017+dy3*c1117)));

    /* d2f/dxdy */
    f_df[7] = xgi*ygi*(
        dzi*(
            (c0000  -c0100)-
            (c0010-c0110))
        +dz*(
            (c1000  -c1100)-
            (c1010-c1110)))
        +ygi*xg/6*(
            dzi*(
                dxi3dx*(-c0001+c0101)+
                dx3dx*(-c0011+c0111))
            +dz*(
                dxi3dx*(-c1001+c1101)+
                dx3dx*(-c1011+c1111)))
        +xgi*yg/6*(
            dzi*(
                -(dyi3dy*c0002+dy3dy*c0102)
                +(dyi3dy*c0012+dy3dy*c0112))
            +dz*(
                -(dyi3dy*c1002+dy3dy*c1102)
                +(dyi3dy*c1012+dy3dy*c1112)))
        +xgi*ygi*zg2/6*(
            dzi3*(
                (c0003  -c0103)-
                (c0013-c0113))
            +dz3*(
                (c1003  -c1103)-
                (c1013-c1113)))
        +xg*yg/36*(
            dzi*(
                dxi3dx*(dyi3dy*c0004+dy3dy*c0104)+
                dx3dx*(dyi3dy*c0014+dy3dy*c0114))
            +dz*(
                dxi3dx*(dyi3dy*c1004+dy3dy*c1104)+
                dx3dx*(dyi3dy*c1014+dy3dy*c1114)))
        +ygi*xg*zg2/36*(
            dzi3*(
                dxi3dx*(-c0005+c0105)+
                dx3dx*(-c0015+c0115))
            +dz3*(
                dxi3dx*(-c1005+c1105)+
                dx3dx*(-c1015+c1115)))
        +xgi*yg*zg2/36*(
            dzi3*(
                -(dyi3dy*c0006+dy3dy*c0106)
                +(dyi3dy*c0016+dy3dy*c0116))
            +dz3*(
                -(dyi3dy*c1006+dy3dy*c1106)
                +(dyi3dy*c1016+dy3dy*c1116)))
        +xg*yg*zg2/216*(
            dzi3*(
                dxi3dx*(dyi3dy*c0007+dy3dy*c0107)+
                dx3dx*(dyi3dy*c0017+dy3dy*c0117))
            +dz3*(
                dxi3dx*(dyi3dy*c1007+dy3dy*c1107)+
                dx3dx*(dyi3dy*c1017+dy3dy*c1117)));

    /* d2f/dxdz */
    f_df[8] = xgi*zgi*(
        (
            (dyi*c0000+dy*c0100) -
            (dyi*c0010+dy*c0110))
        -(
            (dyi*c1000+dy*c1100) -
            (dyi*c1010+dy*c1110)))
        +xg*zgi/6*(
            -(
                dxi3dx*(dyi*c0001+dy*c0101)+
                dx3dx*(dyi*c0011+dy*c0111))
            +(
                dxi3dx*(dyi*c1001+dy*c1101)+
                dx3dx*(dyi*c1011+dy*c1111)))
        +xgi*yg2*zgi/6*(
            (
                (dyi3*c0002+dy3*c0102) -
                (dyi3*c0012+dy3*c0112))
            -(
                (dyi3*c1002+dy3*c1102) -
                (dyi3*c1012+dy3*c1112)))
        +xgi*zg/6*(
            dzi3dz*(
                -(dyi*c0003+dy*c0103)
                +(dyi*c0013+dy*c0113))
            +dz3dz*(
                -(dyi*c1003+dy*c1103)
                +(dyi*c1013+dy*c1113)))
        +xg*yg2*zgi/36*(
            -(
                dxi3dx*(dyi3*c0004+dy3*c0104)+
                dx3dx*(dyi3*c0014+dy3*c0114))
            +(
                dxi3dx*(dyi3*c1004+dy3*c1104)+
                dx3dx*(dyi3*c1014+dy3*c1114)))
        +xg*zg/36*(
            dzi3dz*(
                dxi3dx*(dyi*c0005+dy*c0105)+
                dx3dx*(dyi*c0015+dy*c0115))
            +dz3dz*(
                dxi3dx*(dyi*c1005+dy*c1105)+
                dx3dx*(dyi*c1015+dy*c1115)))
        +xgi*yg2*zg/36*(
            dzi3dz*(
                -(dyi3*c0006+dy3*c0106)
                +(dyi3*c0016+dy3*c0116))
            +dz3dz*(
                -(dyi3*c1006+dy3*c1106)
                +(dyi3*c1016+dy3*c1116)))
        +xg*yg2*zg/216*(
            dzi3dz*(
                dxi3dx*(dyi3*c0007+dy3*c0107)+
                dx3dx*(dyi3*c0017+dy3*c0117))
            +dz3dz*(
                dxi3dx*(dyi3*c1007+dy3*c1107)+
                dx3dx*(dyi3*c1017+dy3*c1117)));

    /* d2f/dydz */
    f_df[9] = ygi*zgi*(
        (
            dxi*(c0000  -c0100)+
            dx*(c0010-c0110))
        -(
            dxi*(c1000  -c1100)+
            dx*(c1010-c1110)))
        +ygi*xg2*zgi/6*(
            (
                dxi3*(c0001  -c0101)+
                dx3*(c0011-c0111))
            -(
                dxi3*(c1001  -c1101)+
                dx3*(c1011-c1111)))
        +yg*zgi/6*(
            -(
                dxi*(dyi3dy*c0002+dy3dy*c0102)+
                dx*(dyi3dy*c0012+dy3dy*c0112))
            +(
                dxi*(dyi3dy*c1002+dy3dy*c1102)+
                dx*(dyi3dy*c1012+dy3dy*c1112)))
        +ygi*zg/6*(
            dzi3dz*(
                dxi*(-c0003+c0103)+
                dx*(-c0013+c0113))
            +dz3dz*(
                dxi*(-c1003+c1103)+
                dx*(-c1013+c1113)))
        +xg2*yg*zgi/36*(
            -(
                dxi3*(dyi3dy*c0004+dy3dy*c0104)+
                dx3*(dyi3dy*c0014+dy3dy*c0114))
            +(
                dxi3*(dyi3dy*c1004+dy3dy*c1104)+
                dx3*(dyi3dy*c1014+dy3dy*c1114)))
        +ygi*xg2*zg/36*(
            dzi3dz*(
                dxi3*(-c0005+c0105)+
                dx3*(-c0015+c0115))
            +dz3dz*(
                dxi3*(-c1005+c1105)+
                dx3*(-c1015+c1115)))
        +yg*zg/36*(
            dzi3dz*(
                dxi*(dyi3dy*c0006+dy3dy*c0106)+
                dx*(dyi3dy*c0016+dy3dy*c0116))
            +dz3dz*(
                dxi*(dyi3dy*c1006+dy3dy*c1106)+
                dx*(dyi3dy*c1016+dy3dy*c1116)))
        +xg2*yg*zg/216*(
            dzi3dz*(
                dxi3*(dyi3dy*c0007+dy3dy*c0107)+
                dx3*(dyi3dy*c0017+dy3dy*c0117))
            +dz3dz*(
                dxi3*(dyi3dy*c1007+dy3dy*c1107)+
                dx3*(dyi3dy*c1017+dy3dy*c1117)));

    }

    return err;
}

/**
 * @brief Free allocated memory in interpolation data struct
 *
 * This function frees the memory allocated for interpolation coefficients
 * in the interpolation data struct
 *
 * @todo Error checking
 *
 * @param str data struct for data interpolation
 */
void interp3Dcomp_free(interp3D_data* str) {
    free(str->c);
}
