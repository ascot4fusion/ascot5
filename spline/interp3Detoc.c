/**
 * @file interp3Detoc.c
 * @brief Tricubic spline interpolation, i.e. cubic spline interpolation of 3D scalar data
 */
#include <stdlib.h>
#include <stdio.h> /* Needed for printf debugging purposes */
#include <math.h>
#include "../ascot5.h"
#include "../consts.h"
#include "interp3Detoc.h"
#include "spline1D.h"

/**
 * @brief Calculate tricubic spline interpolation coefficients for scalar 3D data
 *
 * This function calculates the tricubic spline interpolation coefficients for
 * the given data and stores them in the data struct. The explicit cofficients
 * are first calculated and then compact coefficients are derived from these.
 *
 * @todo Directly calculate compact coefficients
 * @todo Error checking
 *
 * @param str data struct for data interpolation
 * @param f 3D data to be interpolated
 * @param n_r number of data points in the r direction
 * @param n_phi number of data points in the phi direction
 * @param n_z number of data points in the z direction
 * @param r_min minimum value of the r axis
 * @param r_max maximum value of the r axis
 * @param phi_min minimum value of the phi axis
 * @param phi_max maximum value of the phi axis
 * @param z_min minimum value of the z axis
 * @param z_max maximum value of the z axis
 */
void interp3Detoc_init(interp3D_data* str, real* f, int n_r, int n_phi, int n_z,
                       real r_min, real r_max, real r_grid,
                       real phi_min, real phi_max, real phi_grid,
                       real z_min, real z_max, real z_grid) {

    /* Initialize and fill the data struct */
    str->n_r = n_r;
    str->n_phi = n_phi;
    str->n_z = n_z;
    str->r_min = r_min;
    str->r_max = r_max;
    str->r_grid = r_grid;
    str->phi_min = phi_min;
    str->phi_max = phi_max;
    str->phi_grid = phi_grid;
    str->z_min = z_min;
    str->z_max = z_max;
    str->z_grid = z_grid;
    str->c = malloc(n_phi*n_z*n_r*64*sizeof(real));

    /* Declare and allocate the needed variables */
    int i_r;                                  /**< index for r variable */
    int i_phi;                                /**< index for phi variable */
    int i_z;                                  /**< index for z variable */
    int i_c;                                  /**< index for coefficient for data struct */
    real* f_r = malloc(n_r*sizeof(real));       /**< Temporary array for data along r */
    real* f_phi = malloc(n_phi*sizeof(real));   /**< Temp array for data along phi */
    real* f_z = malloc(n_z*sizeof(real));       /**< Temp array for data along z */
    real* c_r = malloc((n_r-1)*4*sizeof(real)); /**< Temp array for coefficients along r */
    real* c_phi = malloc(n_phi*4*sizeof(real)); /**< Temp array for coefs along phi */
    real* c_z = malloc((n_z-1)*4*sizeof(real)); /**< Temp array for coefficients along z */
    int i_s;                                 /**< index for spline regarding degree in r */
    int i_ss;                                /**< index for spline regarding degree in z */
    int i_ct;                                   /**< index for coefficient array */

    /* Bicubic spline surface over rz-grid for each phi */
    for(i_phi=0; i_phi<n_phi; i_phi++) {
        /* Cubic spline along r for each z */
        for(i_z=0; i_z<n_z; i_z++) {
            for(i_r=0; i_r<n_r; i_r++) {
                f_r[i_r] = f[i_phi*n_z*n_r+i_z*n_r+i_r];
            }
            spline1D(f_r,n_r,0,c_r);
            for(i_r=0; i_r<n_r-1; i_r++) {
                for(i_c=0; i_c<4; i_c++) {
                    i_ct = i_c;
                    str->c[i_phi*n_z*n_r*64+i_z*n_r*64+i_r*64+i_c] = c_r[i_r*4+i_ct];
                }
            }
        }

        /* Four cubic splines along z for each r using four different data sets */
        for(i_r=0; i_r<n_r-1; i_r++) {
            /* s0, s1, s2, s3 */
            for(i_s=0; i_s<4; i_s++) {
                for(i_z=0; i_z<n_z; i_z++) {
                    f_z[i_z] = str->c[i_phi*n_z*n_r*64+i_z*n_r*64+i_r*64+i_s];
                }
                spline1D(f_z,n_z,0,c_z);
                for(i_z=0; i_z<n_z-1; i_z++) {
                    i_ct = 0;
                    for(i_c=i_s; i_c<16; i_c=i_c+4) {
                        str->c[i_phi*n_z*n_r*64+i_z*n_r*64+i_r*64+i_c] = c_z[i_z*4+i_ct];
                        i_ct++;
                    }
                }
            }
        }
    }

    /* Cubic spline along phi for each rz-pair to find the coefficients
       of the tricubic spline volume */
    for(i_z=0; i_z<n_z-1; i_z++) {
        for(i_r=0; i_r<n_r-1; i_r++) {
            for(i_ss=0; i_ss<4; i_ss++) {
                for(i_s=0; i_s<4; i_s++) {
                    for(i_phi=0; i_phi<n_phi; i_phi++) {
                        f_phi[i_phi] = str->c[i_phi*n_z*n_r*64+i_z*n_r*64
                                            +i_r*64+(i_ss*4+i_s)];
                    }
                    spline1D(f_phi,n_phi,1,c_phi);
                    for(i_phi=0; i_phi<n_phi; i_phi++) {
                        i_ct = 0;
                        for(i_c=4*i_ss+i_s; i_c<64; i_c=i_c+16) {
                            str->c[i_phi*n_z*n_r*64+i_z*n_r*64+i_r*64+i_c]
                                = c_phi[i_phi*4+i_ct];
                            i_ct++;
                        }
                    }
                }
            }
        }
    }

    /* Free allocated memory */
    free(f_r);
    free(f_phi);
    free(f_z);
    free(c_r);
    free(c_phi);
    free(c_z);

    /* Transform from explicit to compact */
    real* cc = malloc(n_phi*n_z*n_r*8*sizeof(real)); /**< Temporary coefficient array */
    for(i_phi=0; i_phi<n_z; i_phi++) {
        for(i_z=0; i_z<n_z-1; i_z++) {
            for(i_r=0; i_r<n_r-1; i_r++) {
                cc[i_phi*n_z*n_r*8+i_z*n_r*8+i_r*8+0] =
                    str->c[i_phi*n_z*n_r*64+i_z*n_r*64+i_r*64+0];
                cc[i_phi*n_z*n_r*8+i_z*n_r*8+i_r*8+1] = (1.0/(r_grid*r_grid))*
                    2*str->c[i_phi*n_z*n_r*64+i_z*n_r*64+i_r*64+2];
                cc[i_phi*n_z*n_r*8+i_z*n_r*8+i_r*8+2] = (1.0/(phi_grid*phi_grid))*
                    2*str->c[i_phi*n_z*n_r*64+i_z*n_r*64+i_r*64+32];
                cc[i_phi*n_z*n_r*8+i_z*n_r*8+i_r*8+3] = (1.0/(z_grid*z_grid))*
                    2*str->c[i_phi*n_z*n_r*64+i_z*n_r*64+i_r*64+8];
                cc[i_phi*n_z*n_r*8+i_z*n_r*8+i_r*8+4] = (1.0/(r_grid*r_grid*phi_grid*phi_grid))*
                    4*str->c[i_phi*n_z*n_r*64+i_z*n_r*64+i_r*64+34];
                cc[i_phi*n_z*n_r*8+i_z*n_r*8+i_r*8+5] = (1.0/(r_grid*r_grid*z_grid*z_grid))*
                    4*str->c[i_phi*n_z*n_r*64+i_z*n_r*64+i_r*64+10];
                cc[i_phi*n_z*n_r*8+i_z*n_r*8+i_r*8+6] = (1.0/(r_grid*r_grid*phi_grid*phi_grid*z_grid*z_grid))*
                    4*str->c[i_phi*n_z*n_r*64+i_z*n_r*64+i_r*64+40];
                cc[i_phi*n_z*n_r*8+i_z*n_r*8+i_r*8+7] = (1.0/(r_grid*r_grid*phi_grid*phi_grid*
                    z_grid*z_grid))*
                    8*str->c[i_phi*n_z*n_r*64+i_z*n_r*64+i_r*64+9+42];

            }
            cc[i_phi*n_z*n_r*8+i_z*n_r*8+(n_r-1)*8+0] = f[i_phi*n_z*n_r+i_z*n_r+n_r-1]; // From coefs?
            cc[i_phi*n_z*n_r*8+i_z*n_r*8+(n_r-1)*8+1] = 0; // Are all derivs 0?
            cc[i_phi*n_z*n_r*8+i_z*n_r*8+(n_r-1)*8+2] = 0;
            cc[i_phi*n_z*n_r*8+i_z*n_r*8+(n_r-1)*8+3] = 0;
            cc[i_phi*n_z*n_r*8+i_z*n_r*8+(n_r-1)*8+4] = 0;
            cc[i_phi*n_z*n_r*8+i_z*n_r*8+(n_r-1)*8+5] = 0;
            cc[i_phi*n_z*n_r*8+i_z*n_r*8+(n_r-1)*8+6] = 0;
            cc[i_phi*n_z*n_r*8+i_z*n_r*8+(n_r-1)*8+7] = 0;
        }
    }
    for(i_r=0; i_r<n_r; i_r++) {
        cc[i_phi*n_z*n_r*8+(n_z-1)*n_r*8+i_r*8+0] = f[i_phi*n_z*n_r+(n_z-1)*n_r+i_r]; // From coefs?
        cc[i_phi*n_z*n_r*8+(n_z-1)*n_r*8+i_r*8+0] = 0; // Are all derivs 0?
        cc[i_phi*n_z*n_r*8+(n_z-1)*n_r*8+i_r*8+0] = 0;
        cc[i_phi*n_z*n_r*8+(n_z-1)*n_r*8+i_r*8+0] = 0;
        cc[i_phi*n_z*n_r*8+(n_z-1)*n_r*8+i_r*8+0] = 0;
        cc[i_phi*n_z*n_r*8+(n_z-1)*n_r*8+i_r*8+0] = 0;
        cc[i_phi*n_z*n_r*8+(n_z-1)*n_r*8+i_r*8+0] = 0;
        cc[i_phi*n_z*n_r*8+(n_z-1)*n_r*8+i_r*8+0] = 0;
    }
    free(str->c);
    str->c = cc;
}

/**
 * @brief Evaluate interpolated value of 3D scalar field
 *
 * This function evaluates the interpolated value of a 3D scalar field using
 * tricubic spline interpolation coefficients of the compact form.
 *
 * @todo Check discrepency to ascot4 and explicit version
 * @todo Error checking
 *
 * @param B variable in which to place the evaluated value
 * @param str data struct for data interpolation
 * @param r r-coordinate
 * @param phi phi-coordinate
 * @param z z-coordinate
 */
void interp3Detoc_eval_B(real* B, interp3D_data* str, real r, real phi, real z) {
    /** Make sure phi is in interval [phi_min, phi_max + phi_grid) */
    phi_range = (str->phi_max + str->phi_grid - str->phi_min);
    phi = fmod(phi - str->phi_min, phi_range) + str->phi_min;
    if(phi < 0){phi = phi_range + phi;}

    int i_r = (r-str->r_min)/str->r_grid;     /**< index for r variable */
    real dr = (r-(str->r_min+i_r*str->r_grid))/str->r_grid; /**< Normalized r coordinate in
                                                               current cell */
    real dri = 1.0-dr;
    real dr3 = dr*(dr*dr-1);
    real dri3 = (1.0-dr)*(1.0-dr)*(1.0-dr)-(1.0-dr);
    real rg2 = str->r_grid*str->r_grid;       /**< Square of cell length in r direction */
    int i_phi = (phi-str->phi_min)/str->phi_grid; /**< index for phi variable */
    real dphi = (phi-(str->phi_min+i_phi*str->phi_grid))/str->phi_grid; /**< Normalized phi
                                                                           coordinate in
                                                                           current cell */
    real dphii = 1.0-dphi;
    real dphi3 = dphi*dphi*dphi-dphi;
    real dphii3 = (1.0-dphi)*(1.0-dphi)*(1.0-dphi)-(1.0-dphi);
    real phig2 = str->phi_grid*str->phi_grid; /**< Square of cell length in phi dir */
    int i_z = (z-str->z_min)/str->z_grid;     /**< index for z variable */
    real dz = (z-(str->z_min+i_z*str->z_grid))/str->z_grid; /**< Normalized z coordinate in
                                                               current cell */
    real dzi = 1.0-dz;
    real dz3 = dz*dz*dz-dz;
    real dzi3 = (1.0-dz)*(1.0-dz)*(1.0-dz)-(1.0-dz);
    real zg2 = str->z_grid*str->z_grid;       /**< Square of cell length in z direction */
    int n = i_phi*str->n_z*str->n_r*8+i_z*str->n_r*8+i_r*8; /**< Index jump to cell */
    int r1 = 8;                               /**< Index jump one r forward */
    int phi1 = str->n_z*str->n_r*8;           /**< Index jump one phi forward */
    if(i_phi==str->n_phi-1) {
        phi1 = -(str->n_phi-1)*phi1;          /**< If last cell, index jump to 1st phi */
    }
    int z1 = str->n_r*8;                      /**< Index jump one z forward */

    *B = (
          dzi*(
               dri*(dphii*str->c[n+0]+dphi*str->c[n+phi1+0])+
               dr*(dphii*str->c[n+r1+0]+dphi*str->c[n+phi1+r1+0]))
          +dz*(
               dri*(dphii*str->c[n+z1+0]+dphi*str->c[n+phi1+z1+0])+
               dr*(dphii*str->c[n+r1+z1+0]+dphi*str->c[n+phi1+z1+r1+0])))
          +rg2/6*(
                  dzi*(
                       dri3*(dphii*str->c[n+1]+dphi*str->c[n+phi1+1])+
                       dr3*(dphii*str->c[n+r1+1]+dphi*str->c[n+phi1+r1+1]))
                  +dz*(
                       dri3*(dphii*str->c[n+z1+1]+dphi*str->c[n+phi1+z1+1])+
                       dr3*(dphii*str->c[n+r1+z1+1]+dphi*str->c[n+phi1+z1+r1+1])))
          +phig2/6*(
                    dzi*(
                         dri*(dphii3*str->c[n+2]+dphi3*str->c[n+phi1+2])+
                         dr*(dphii3*str->c[n+r1+2]+dphi3*str->c[n+phi1+r1+2]))
                    +dz*(
                         dri*(dphii3*str->c[n+z1+2]+dphi3*str->c[n+phi1+z1+2])+
                         dr*(dphii3*str->c[n+r1+z1+2]+dphi3*str->c[n+phi1+z1+r1+2])))
          +zg2/6*(
                  dzi3*(
                        dri*(dphii*str->c[n+3]+dphi*str->c[n+phi1+3])+
                        dr*(dphii*str->c[n+r1+3]+dphi*str->c[n+phi1+r1+3]))
                  +dz3*(
                        dri*(dphii*str->c[n+z1+3]+dphi*str->c[n+phi1+z1+3])+
                        dr*(dphii*str->c[n+r1+z1+3]+dphi*str->c[n+phi1+z1+r1+3])))
          +rg2*phig2/36*(
                         dzi*(
                              dri3*(dphii3*str->c[n+4]+dphi3*str->c[n+phi1+4])+
                              dr3*(dphii3*str->c[n+r1+4]+dphi3*str->c[n+phi1+r1+4]))
                         +dz*(
                              dri3*(dphii3*str->c[n+z1+4]+dphi3*str->c[n+phi1+z1+4])+
                              dr3*(dphii3*str->c[n+r1+z1+4]+dphi3*str->c[n+phi1+z1+r1+4])))
          +rg2*zg2/36*(
                       dzi3*(
                             dri3*(dphii*str->c[n+5]+dphi*str->c[n+phi1+5])+
                             dr3*(dphii*str->c[n+r1+5]+dphi*str->c[n+phi1+r1+5]))
                       +dz3*(
                             dri3*(dphii*str->c[n+z1+5]+dphi*str->c[n+phi1+z1+5])+
                             dr3*(dphii*str->c[n+r1+z1+5]+dphi*str->c[n+phi1+z1+r1+5])))
          +phig2*zg2/36*(
                         dzi3*(
                               dri*(dphii3*str->c[n+6]+dphi3*str->c[n+phi1+6])+
                               dr*(dphii3*str->c[n+r1+6]+dphi3*str->c[n+phi1+r1+6]))
                         +dz3*(
                               dri*(dphii3*str->c[n+z1+6]+dphi3*str->c[n+phi1+z1+6])+
                               dr*(dphii3*str->c[n+r1+z1+6]+dphi3*str->c[n+phi1+z1+r1+6])))
          +rg2*phig2*zg2/216*(
                              dzi3*(
                                    dri3*(dphii3*str->c[n+7]+dphi3*str->c[n+phi1+7])+
                                    dr3*(dphii3*str->c[n+r1+7]+dphi3*str->c[n+phi1+r1+7]))
                              +dz3*(
                                    dri3*(dphii3*str->c[n+z1+7]+dphi3*str->c[n+phi1+z1+7])+
                                    dr3*(dphii3*str->c[n+r1+z1+7]+dphi3*str->c[n+phi1+z1+r1+7])));
}

/**
 * @brief Evaluate interpolated value of 3D scalar field and its 1st and 2nd derivatives
 *
 * This function evaluates the interpolated value of a 3D scalar field and
 * its 1st and 2nd derivatives using bicubic spline interpolation coefficients
 * of the compact form.
 *
 * @todo Check discrepency to ascot4 and explicit version
 * @todo Error checking
 *
 * @param B_dB array in which to place the evaluated values
 * @param str data struct for data interpolation
 * @param r r-coordinate
 * @param phi phi-coordinate
 * @param z z-coordinate
 */
void interp3Detoc_eval_dB(real* B_dB, interp3D_data* str, real r, real phi, real z) {
    /** Make sure phi is in interval [phi_min, phi_max + phi_grid) */
    phi_range = (str->phi_max + str->phi_grid - str->phi_min);
    phi = fmod(phi - str->phi_min, phi_range) + str->phi_min;
    if(phi < 0){phi = phi_range + phi;}

    int i_r = (r-str->r_min)/str->r_grid;       /**< index for r variable */
    real dr = (r-(str->r_min+i_r*str->r_grid))/str->r_grid; /**< Normalized r coordinate in
                                                               current cell */
    real dr3 = dr*dr*dr-dr;
    real dr3dr = 3*dr*dr-1.0;           /**< r-derivative of dr3, not including 1/r_grid */
    real dri = 1.0-dr;
    real dri3 = dri*dri*dri-dri;
    real dri3dr = -3*dri*dri+1.0;      /**< r-derivative of dri3, not including 1/r_grid */
    real rg = str->r_grid;                      /**< Cell length in r direction */
    real rg2 = rg*rg;
    real rgi = 1.0/rg;
    int i_phi = (phi-str->phi_min)/str->phi_grid; /**< index for phi variable */
    real dphi = (phi-(str->phi_min+i_phi*str->phi_grid))/str->phi_grid; /**< Normalized phi
                                                                           coordinate in
                                                                           current cell */
    real dphi3 = dphi*dphi*dphi-dphi;
    real dphi3dphi = 3*dphi*dphi-1.0; /**< phi-derivative of dphi3,
                                         not including 1/phi_grid */
    real dphii = 1.0-dphi;
    real dphii3 = dphii*dphii*dphii-dphii;
    real dphii3dphi = -3*dphii*dphii+1.0; /**< phi-derivative of dphii3,
                                             not including 1/r_grid */
    real phig = str->phi_grid;                  /**< Cell length in phi direction */
    real phig2 = phig*phig;
    real phigi = 1.0/phig;
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
    int n = i_phi*str->n_z*str->n_r*8+i_z*str->n_r*8+i_r*8; /**< Index jump to cell */
    int r1 = 8;                                 /**< Index jump one r forward */
    int phi1 = str->n_z*str->n_r*8;             /**< Index jump one phi forward */
    if(i_phi==str->n_phi-1) {
        phi1 = -(str->n_phi-1)*phi1;            /**< If last cell, index jump to 1st phi */
    }
    int z1 = str->n_r*8;                        /**< Index jump one z forward */

    /* f */
    B_dB[0] = (
               dzi*(
                    dri*(dphii*str->c[n+0]+dphi*str->c[n+phi1+0])+
                    dr*(dphii*str->c[n+r1+0]+dphi*str->c[n+phi1+r1+0]))
               +dz*(
                    dri*(dphii*str->c[n+z1+0]+dphi*str->c[n+phi1+z1+0])+
                    dr*(dphii*str->c[n+r1+z1+0]+dphi*str->c[n+phi1+z1+r1+0])))
        +rg2/6*(
                dzi*(
                     dri3*(dphii*str->c[n+1]+dphi*str->c[n+phi1+1])+
                     dr3*(dphii*str->c[n+r1+1]+dphi*str->c[n+phi1+r1+1]))
                +dz*(
                     dri3*(dphii*str->c[n+z1+1]+dphi*str->c[n+phi1+z1+1])+
                     dr3*(dphii*str->c[n+r1+z1+1]+dphi*str->c[n+phi1+z1+r1+1])))
        +phig2/6*(
                  dzi*(
                       dri*(dphii3*str->c[n+2]+dphi3*str->c[n+phi1+2])+
                       dr*(dphii3*str->c[n+r1+2]+dphi3*str->c[n+phi1+r1+2]))
                  +dz*(
                       dri*(dphii3*str->c[n+z1+2]+dphi3*str->c[n+phi1+z1+2])+
                       dr*(dphii3*str->c[n+r1+z1+2]+dphi3*str->c[n+phi1+z1+r1+2])))
        +zg2/6*(
                dzi3*(
                      dri*(dphii*str->c[n+3]+dphi*str->c[n+phi1+3])+
                      dr*(dphii*str->c[n+r1+3]+dphi*str->c[n+phi1+r1+3]))
                +dz3*(
                      dri*(dphii*str->c[n+z1+3]+dphi*str->c[n+phi1+z1+3])+
                      dr*(dphii*str->c[n+r1+z1+3]+dphi*str->c[n+phi1+z1+r1+3])))
        +rg2*phig2/36*(
                       dzi*(
                            dri3*(dphii3*str->c[n+4]+dphi3*str->c[n+phi1+4])+
                            dr3*(dphii3*str->c[n+r1+4]+dphi3*str->c[n+phi1+r1+4]))
                       +dz*(
                            dri3*(dphii3*str->c[n+z1+4]+dphi3*str->c[n+phi1+z1+4])+
                            dr3*(dphii3*str->c[n+r1+z1+4]+dphi3*str->c[n+phi1+z1+r1+4])))
        +rg2*zg2/36*(
                     dzi3*(
                           dri3*(dphii*str->c[n+5]+dphi*str->c[n+phi1+5])+
                           dr3*(dphii*str->c[n+r1+5]+dphi*str->c[n+phi1+r1+5]))
                     +dz3*(
                           dri3*(dphii*str->c[n+z1+5]+dphi*str->c[n+phi1+z1+5])+
                           dr3*(dphii*str->c[n+r1+z1+5]+dphi*str->c[n+phi1+z1+r1+5])))
        +phig2*zg2/36*(
                       dzi3*(
                             dri*(dphii3*str->c[n+6]+dphi3*str->c[n+phi1+6])+
                             dr*(dphii3*str->c[n+r1+6]+dphi3*str->c[n+phi1+r1+6]))
                       +dz3*(
                             dri*(dphii3*str->c[n+z1+6]+dphi3*str->c[n+phi1+z1+6])+
                             dr*(dphii3*str->c[n+r1+z1+6]+dphi3*str->c[n+phi1+z1+r1+6])))
        +rg2*phig2*zg2/216*(
                            dzi3*(
                                  dri3*(dphii3*str->c[n+7]+dphi3*str->c[n+phi1+7])+
                                  dr3*(dphii3*str->c[n+r1+7]+dphi3*str->c[n+phi1+r1+7]))
                            +dz3*(
                                  dri3*(dphii3*str->c[n+z1+7]+dphi3*str->c[n+phi1+z1+7])+
                                  dr3*(dphii3*str->c[n+r1+z1+7]+dphi3*str->c[n+phi1+z1+r1+7])));

    /* df/dr */
    B_dB[1] = rgi*(
                   dzi*(
                        -(dphii*str->c[n+0]+dphi*str->c[n+phi1+0])
                        +(dphii*str->c[n+r1+0]+dphi*str->c[n+phi1+r1+0]))
                   +dz*(
                        -(dphii*str->c[n+z1+0]+dphi*str->c[n+phi1+z1+0])
                        +(dphii*str->c[n+r1+z1+0]+dphi*str->c[n+phi1+z1+r1+0])))
        +rg/6*(
               dzi*(
                    dri3dr*(dphii*str->c[n+1]+dphi*str->c[n+phi1+1])+
                    dr3dr*(dphii*str->c[n+r1+1]+dphi*str->c[n+phi1+r1+1]))
               +dz*(
                    dri3dr*(dphii*str->c[n+z1+1]  +dphi*str->c[n+phi1+z1+1])+
                    dr3dr*(dphii*str->c[n+r1+z1+1]+dphi*str->c[n+phi1+z1+r1+1])))
        +rgi*phig2/6*(
                      dzi*(
                           -(dphii3*str->c[n+2]+dphi3*str->c[n+phi1+2])
                           +(dphii3*str->c[n+r1+2]+dphi3*str->c[n+phi1+r1+2]))
                      +dz*(
                           -(dphii3*str->c[n+z1+2]+dphi3*str->c[n+phi1+z1+2])
                           +(dphii3*str->c[n+r1+z1+2]+dphi3*str->c[n+phi1+z1+r1+2])))
        +rgi*zg2/6*(
                    dzi3*(
                          -(dphii*str->c[n+3]+dphi*str->c[n+phi1+3])
                          +(dphii*str->c[n+r1+3]+dphi*str->c[n+phi1+r1+3]))
                    +dz3*(
                          -(dphii*str->c[n+z1+3]+dphi*str->c[n+phi1+z1+3])
                          +(dphii*str->c[n+r1+z1+3]+dphi*str->c[n+phi1+z1+r1+3])))
        +rg*phig2/36*(
                      dzi*(
                           dri3dr*(dphii3*str->c[n+4]+dphi3*str->c[n+phi1+4])+
                           dr3dr*(dphii3*str->c[n+r1+4]+dphi3*str->c[n+phi1+r1+4]))
                      +dz*(
                           dri3dr*(dphii3*str->c[n+z1+4]+dphi3*str->c[n+phi1+z1+4])+
                           dr3dr*(dphii3*str->c[n+r1+z1+4]+dphi3*str->c[n+phi1+z1+r1+4])))
        +rg*zg2/36*(
                    dzi3*(
                          dri3dr*(dphii*str->c[n+5]+dphi*str->c[n+phi1+5])+
                          dr3dr*(dphii*str->c[n+r1+5]+dphi*str->c[n+phi1+r1+5]))
                    +dz3*(
                          dri3dr*(dphii*str->c[n+z1+5]+dphi*str->c[n+phi1+z1+5])+
                          dr3dr*(dphii*str->c[n+r1+z1+5]+dphi*str->c[n+phi1+z1+r1+5])))
        +rgi*phig2*zg2/36*(
                           dzi3*(
                                 -(dphii3*str->c[n+6]+dphi3*str->c[n+phi1+6])
                                 +(dphii3*str->c[n+r1+6]+dphi3*str->c[n+phi1+r1+6]))
                           +dz3*(
                                 -(dphii3*str->c[n+z1+6]+dphi3*str->c[n+phi1+z1+6])
                                 +(dphii3*str->c[n+r1+z1+6]+dphi3*str->c[n+phi1+z1+r1+6])))
        +rg*phig2*zg2/216*(
                           dzi3*(
                                 dri3dr*(dphii3*str->c[n+7]+dphi3*str->c[n+phi1+7])+
                                 dr3dr*(dphii3*str->c[n+r1+7]+dphi3*str->c[n+phi1+r1+7]))
                           +dz3*(
                                 dri3dr*(dphii3*str->c[n+z1+7]+dphi3*str->c[n+phi1+z1+7])+
                                 dr3dr*(dphii3*str->c[n+r1+z1+7]+dphi3*str->c[n+phi1+z1+r1+7])));

    /* df/dphi */
    B_dB[2] = phigi*(
                     dzi*(
                          dri*(-str->c[n+0]+str->c[n+phi1+0])+
                          dr*(-str->c[n+r1+0]+str->c[n+phi1+r1+0]))
                     +dz*(
                          dri*(-str->c[n+z1+0]+str->c[n+phi1+z1+0])+
                          dr*(-str->c[n+r1+z1+0]+str->c[n+phi1+z1+r1+0])))
        +phigi*rg2/6*(
                      dzi*(
                           dri3*(-str->c[n+1]+str->c[n+phi1+1])+
                           dr3*(-str->c[n+r1+1]+str->c[n+phi1+r1+1]))
                      +dz*(
                           dri3*(-str->c[n+z1+1]+str->c[n+phi1+z1+1])+
                           dr3*(-str->c[n+r1+z1+1]+str->c[n+phi1+z1+r1+1])))
        +phig/6*(
                 dzi*(
                      dri*(dphii3dphi*str->c[n+2]+dphi3dphi*str->c[n+phi1+2])+
                      dr*(dphii3dphi*str->c[n+r1+2]+dphi3dphi*str->c[n+phi1+r1+2]))
                 +dz*(
                      dri*(dphii3dphi*str->c[n+z1+2]+dphi3dphi*str->c[n+phi1+z1+2])+
                      dr*(dphii3dphi*str->c[n+r1+z1+2]+dphi3dphi*str->c[n+phi1+z1+r1+2])))
        +phigi*zg2/6*(
                      dzi3*(
                            dri*(-str->c[n+3]+str->c[n+phi1+3])+
                            dr*(-str->c[n+r1+3]+str->c[n+phi1+r1+3]))
                      +dz3*(
                            dri*(-str->c[n+z1+3]+str->c[n+phi1+z1+3])+
                            dr*(-str->c[n+r1+z1+3]+str->c[n+phi1+z1+r1+3])))
        +rg2*phig/36*(
                      dzi*(
                           dri3*(dphii3dphi*str->c[n+4]+dphi3dphi*str->c[n+phi1+4])+
                           dr3*(dphii3dphi*str->c[n+r1+4]+dphi3dphi*str->c[n+phi1+r1+4]))
                      +dz*(
                           dri3*(dphii3dphi*str->c[n+z1+4]+dphi3dphi*str->c[n+phi1+z1+4])+
                           dr3*(dphii3dphi*str->c[n+r1+z1+4]+dphi3dphi*str->c[n+phi1+z1+r1+4])))
        +phigi*rg2*zg2/36*(
                           dzi3*(
                                 dri3*(-str->c[n+5]+str->c[n+phi1+5])+
                                 dr3*(-str->c[n+r1+5]+str->c[n+phi1+r1+5]))
                           +dz3*(
                                 dri3*(-str->c[n+z1+5]+str->c[n+phi1+z1+5])+
                                 dr3*(-str->c[n+r1+z1+5]+str->c[n+phi1+z1+r1+5])))
        +phig*zg2/36*(
                      dzi3*(
                            dri*(dphii3dphi*str->c[n+6]+dphi3dphi*str->c[n+phi1+6])+
                            dr*(dphii3dphi*str->c[n+r1+6]+dphi3dphi*str->c[n+phi1+r1+6]))
                      +dz3*(
                            dri*(dphii3dphi*str->c[n+z1+6]+dphi3dphi*str->c[n+phi1+z1+6])+
                            dr*(dphii3dphi*str->c[n+r1+z1+6]+dphi3dphi*str->c[n+phi1+z1+r1+6])))
        +rg2*phig*zg2/216*(
                           dzi3*(
                                 dri3*(dphii3dphi*str->c[n+7]+dphi3dphi*str->c[n+phi1+7])+
                                 dr3*(dphii3dphi*str->c[n+r1+7]+dphi3dphi*str->c[n+phi1+r1+7]))
                           +dz3*(
                                 dri3*(dphii3dphi*str->c[n+z1+7]+dphi3dphi*str->c[n+phi1+z1+7])+
                                 dr3*(dphii3dphi*str->c[n+r1+z1+7]+dphi3dphi*str->c[n+phi1+z1+r1+7])));

    /* df/dz */
    B_dB[3] = zgi*(
                   -(
                     dri*(dphii*str->c[n+0]+dphi*str->c[n+phi1+0])+
                     dr*(dphii*str->c[n+r1+0]+dphi*str->c[n+phi1+r1+0]))
                   +(
                     dri*(dphii*str->c[n+z1+0]+dphi*str->c[n+phi1+z1+0])+
                     dr*(dphii*str->c[n+r1+z1+0]+dphi*str->c[n+phi1+z1+r1+0])))
        +rg2*zgi/6*(
                    -(
                      dri3*(dphii*str->c[n+1]+dphi*str->c[n+phi1+1])+
                      dr3*(dphii*str->c[n+r1+1]+dphi*str->c[n+phi1+r1+1]))
                    +(
                      dri3*(dphii*str->c[n+z1+1]+dphi*str->c[n+phi1+z1+1])+
                      dr3*(dphii*str->c[n+r1+z1+1]+dphi*str->c[n+phi1+z1+r1+1])))
        +phig2*zgi/6*(
                      -(
                        dri*(dphii3*str->c[n+2]+dphi3*str->c[n+phi1+2])+
                        dr*(dphii3*str->c[n+r1+2]+dphi3*str->c[n+phi1+r1+2]))
                      +(
                        dri*(dphii3*str->c[n+z1+2]+dphi3*str->c[n+phi1+z1+2])+
                        dr*(dphii3*str->c[n+r1+z1+2]+dphi3*str->c[n+phi1+z1+r1+2])))
        +zg/6*(
               dzi3dz*(
                       dri*(dphii*str->c[n+3]+dphi*str->c[n+phi1+3])+
                       dr*(dphii*str->c[n+r1+3]+dphi*str->c[n+phi1+r1+3]))
               +dz3dz*(
                       dri*(dphii*str->c[n+z1+3]+dphi*str->c[n+phi1+z1+3])+
                       dr*(dphii*str->c[n+r1+z1+3]+dphi*str->c[n+phi1+z1+r1+3])))
        +rg2*phig2*zgi/36*(
                           -(
                             dri3*(dphii3*str->c[n+4]+dphi3*str->c[n+phi1+4])+
                             dr3*(dphii3*str->c[n+r1+4]+dphi3*str->c[n+phi1+r1+4]))
                           +(
                             dri3*(dphii3*str->c[n+z1+4]+dphi3*str->c[n+phi1+z1+4])+
                             dr3*(dphii3*str->c[n+r1+z1+4]+dphi3*str->c[n+phi1+z1+r1+4])))
        +rg2*zg/36*(
                    dzi3dz*(
                            dri3*(dphii*str->c[n+5]+dphi*str->c[n+phi1+5])+
                            dr3*(dphii*str->c[n+r1+5]+dphi*str->c[n+phi1+r1+5]))
                    +dz3dz*(
                            dri3*(dphii*str->c[n+z1+5]+dphi*str->c[n+phi1+z1+5])+
                            dr3*(dphii*str->c[n+r1+z1+5]+dphi*str->c[n+phi1+z1+r1+5])))
        +phig2*zg/36*(
                      dzi3dz*(
                              dri*(dphii3*str->c[n+6]+dphi3*str->c[n+phi1+6])+
                              dr*(dphii3*str->c[n+r1+6]+dphi3*str->c[n+phi1+r1+6]))
                      +dz3dz*(
                              dri*(dphii3*str->c[n+z1+6]+dphi3*str->c[n+phi1+z1+6])+
                              dr*(dphii3*str->c[n+r1+z1+6]+dphi3*str->c[n+phi1+z1+r1+6])))
        +rg2*phig2*zg/216*(
                           dzi3dz*(
                                   dri3*(dphii3*str->c[n+7]+dphi3*str->c[n+phi1+7])+
                                   dr3*(dphii3*str->c[n+r1+7]+dphi3*str->c[n+phi1+r1+7]))
                           +dz3dz*(
                                   dri3*(dphii3*str->c[n+z1+7]+dphi3*str->c[n+phi1+z1+7])+
                                   dr3*(dphii3*str->c[n+r1+z1+7]+dphi3*str->c[n+phi1+z1+r1+7])));

    /* d2f/dr2 */
    B_dB[4] = (
               dzi*(
                    dri*(dphii*str->c[n+1]+dphi*str->c[n+phi1+1])+
                    dr*(dphii*str->c[n+r1+1]+dphi*str->c[n+phi1+r1+1]))
               +dz*(
                    dri*(dphii*str->c[n+z1+1]+dphi*str->c[n+phi1+z1+1])+
                    dr*(dphii*str->c[n+r1+z1+1]+dphi*str->c[n+phi1+z1+r1+1])))
        +phig2/6*(
                  dzi*(
                       dri*(dphii3*str->c[n+4]+dphi3*str->c[n+phi1+4])+
                       dr*(dphii3*str->c[n+r1+4]+dphi3*str->c[n+phi1+r1+4]))
                  +dz*(
                       dri*(dphii3*str->c[n+z1+4]+dphi3*str->c[n+phi1+z1+4])+
                       dr*(dphii3*str->c[n+r1+z1+4]+dphi3*str->c[n+phi1+z1+r1+4])))
        +zg2/6*(
                dzi3*(
                      dri*(dphii*str->c[n+5]+dphi*str->c[n+phi1+5])+
                      dr*(dphii*str->c[n+r1+5]+dphi*str->c[n+phi1+r1+5]))
                +dz3*(
                      dri*(dphii*str->c[n+z1+5]+dphi*str->c[n+phi1+z1+5])+
                      dr*(dphii*str->c[n+r1+z1+5]+dphi*str->c[n+phi1+z1+r1+5])))
        +phig2*zg2/36*(
                       dzi3*(
                             dri*(dphii3*str->c[n+7]+dphi3*str->c[n+phi1+7])+
                             dr*(dphii3*str->c[n+r1+7]+dphi3*str->c[n+phi1+r1+7]))
                       +dz3*(
                             dri*(dphii3*str->c[n+z1+7]+dphi3*str->c[n+phi1+z1+7])+
                             dr*(dphii3*str->c[n+r1+z1+7]+dphi3*str->c[n+phi1+z1+r1+7])));

    /* d2f/dphi2 */
    B_dB[5] = (
               dzi*(
                    dri*(dphii*str->c[n+2]+dphi*str->c[n+phi1+2])+
                    dr*(dphii*str->c[n+r1+2]+dphi*str->c[n+phi1+r1+2]))
               +dz*(
                    dri*(dphii*str->c[n+z1+2]+dphi*str->c[n+phi1+z1+2])+
                    dr*(dphii*str->c[n+r1+z1+2]+dphi*str->c[n+phi1+z1+r1+2])))
        +rg2/6*(
                dzi*(
                     dri3*(dphii*str->c[n+4]+dphi*str->c[n+phi1+4])+
                     dr3*(dphii*str->c[n+r1+4]+dphi*str->c[n+phi1+r1+4]))
                +dz*(
                     dri3*(dphii*str->c[n+z1+4]+dphi*str->c[n+phi1+z1+4])+
                     dr3*(dphii*str->c[n+r1+z1+4]+dphi*str->c[n+phi1+z1+r1+4])))
        +zg2/6*(
                dzi3*(
                      dri*(dphii*str->c[n+6]+dphi*str->c[n+phi1+6])+
                      dr*(dphii*str->c[n+r1+6]+dphi*str->c[n+phi1+r1+6]))
                +dz3*(
                      dri*(dphii*str->c[n+z1+6]+dphi*str->c[n+phi1+z1+6])+
                      dr*(dphii*str->c[n+r1+z1+6]+dphi*str->c[n+phi1+z1+r1+6])))
        +rg2*zg2/36*(
                     dzi3*(
                           dri3*(dphii*str->c[n+7]+dphi*str->c[n+phi1+7])+
                           dr3*(dphii*str->c[n+r1+7]+dphi*str->c[n+phi1+r1+7]))
                     +dz3*(
                           dri3*(dphii*str->c[n+z1+7]+dphi*str->c[n+phi1+z1+7])+
                           dr3*(dphii*str->c[n+r1+z1+7]+dphi*str->c[n+phi1+z1+r1+7])));

    /* d2f/dz2 */
    B_dB[6] = (
               dzi*(
                    dri*(dphii*str->c[n+3]+dphi*str->c[n+phi1+3])+
                    dr*(dphii*str->c[n+r1+3]+dphi*str->c[n+phi1+r1+3]))
               +dz*(
                    dri*(dphii*str->c[n+z1+3]+dphi*str->c[n+phi1+z1+3])+
                    dr*(dphii*str->c[n+r1+z1+3]+dphi*str->c[n+phi1+z1+r1+3])))
        +rg2/6*(
                dzi*(
                     dri3*(dphii*str->c[n+5]+dphi*str->c[n+phi1+5])+
                     dr3*(dphii*str->c[n+r1+5]+dphi*str->c[n+phi1+r1+5]))
                +dz*(
                     dri3*(dphii*str->c[n+z1+5]+dphi*str->c[n+phi1+z1+5])+
                     dr3*(dphii*str->c[n+r1+z1+5]+dphi*str->c[n+phi1+z1+r1+5])))
        +phig2/6*(
                  dzi*(
                       dri*(dphii3*str->c[n+6]+dphi3*str->c[n+phi1+6])+
                       dr*(dphii3*str->c[n+r1+6]+dphi3*str->c[n+phi1+r1+6]))
                  +dz*(
                       dri*(dphii3*str->c[n+z1+6]+dphi3*str->c[n+phi1+z1+6])+
                       dr*(dphii3*str->c[n+r1+z1+6]+dphi3*str->c[n+phi1+z1+r1+6])))
        +rg2*phig2/36*(
                       dzi*(
                            dri3*(dphii3*str->c[n+7]+dphi3*str->c[n+phi1+7])+
                            dr3*(dphii3*str->c[n+r1+7]+dphi3*str->c[n+phi1+r1+7]))
                       +dz*(
                            dri3*(dphii3*str->c[n+z1+7]+dphi3*str->c[n+phi1+z1+7])+
                            dr3*(dphii3*str->c[n+r1+z1+7]+dphi3*str->c[n+phi1+z1+r1+7])));

    /* d2f/drdphi */
    B_dB[7] = rgi*phigi*(
                         dzi*(
                              (str->c[n+0]  -str->c[n+phi1+0])-
                              (str->c[n+r1+0]-str->c[n+phi1+r1+0]))
                         +dz*(
                              (str->c[n+z1+0]  -str->c[n+phi1+z1+0])-
                              (str->c[n+r1+z1+0]-str->c[n+phi1+z1+r1+0])))
        +phigi*rg/6*(
                     dzi*(
                          dri3dr*(-str->c[n+1]+str->c[n+phi1+1])+
                          dr3dr*(-str->c[n+r1+1]+str->c[n+phi1+r1+1]))
                     +dz*(
                          dri3dr*(-str->c[n+z1+1]+str->c[n+phi1+z1+1])+
                          dr3dr*(-str->c[n+r1+z1+1]+str->c[n+phi1+z1+r1+1])))
        +rgi*phig/6*(
                     dzi*(
                          -(dphii3dphi*str->c[n+2]+dphi3dphi*str->c[n+phi1+2])
                          +(dphii3dphi*str->c[n+r1+2]+dphi3dphi*str->c[n+phi1+r1+2]))
                     +dz*(
                          -(dphii3dphi*str->c[n+z1+2]+dphi3dphi*str->c[n+phi1+z1+2])
                          +(dphii3dphi*str->c[n+r1+z1+2]+dphi3dphi*str->c[n+phi1+z1+r1+2])))
        +rgi*phigi*zg2/6*(
                          dzi3*(
                                (str->c[n+3]  -str->c[n+phi1+3])-
                                (str->c[n+r1+3]-str->c[n+phi1+r1+3]))
                          +dz3*(
                                (str->c[n+z1+3]  -str->c[n+phi1+z1+3])-
                                (str->c[n+r1+z1+3]-str->c[n+phi1+z1+r1+3])))
        +rg*phig/36*(
                     dzi*(
                          dri3dr*(dphii3dphi*str->c[n+4]+dphi3dphi*str->c[n+phi1+4])+
                          dr3dr*(dphii3dphi*str->c[n+r1+4]+dphi3dphi*str->c[n+phi1+r1+4]))
                     +dz*(
                          dri3dr*(dphii3dphi*str->c[n+z1+4]+dphi3dphi*str->c[n+phi1+z1+4])+
                          dr3dr*(dphii3dphi*str->c[n+r1+z1+4]+dphi3dphi*str->c[n+phi1+z1+r1+4])))
        +phigi*rg*zg2/36*(
                          dzi3*(
                                dri3dr*(-str->c[n+5]+str->c[n+phi1+5])+
                                dr3dr*(-str->c[n+r1+5]+str->c[n+phi1+r1+5]))
                          +dz3*(
                                dri3dr*(-str->c[n+z1+5]+str->c[n+phi1+z1+5])+
                                dr3dr*(-str->c[n+r1+z1+5]+str->c[n+phi1+z1+r1+5])))
        +rgi*phig*zg2/36*(
                          dzi3*(
                                -(dphii3dphi*str->c[n+6]+dphi3dphi*str->c[n+phi1+6])
                                +(dphii3dphi*str->c[n+r1+6]+dphi3dphi*str->c[n+phi1+r1+6]))
                          +dz3*(
                                -(dphii3dphi*str->c[n+z1+6]+dphi3dphi*str->c[n+phi1+z1+6])
                                +(dphii3dphi*str->c[n+r1+z1+6]+dphi3dphi*str->c[n+phi1+z1+r1+6])))
        +rg*phig*zg2/216*(
                          dzi3*(
                                dri3dr*(dphii3dphi*str->c[n+7]+dphi3dphi*str->c[n+phi1+7])+
                                dr3dr*(dphii3dphi*str->c[n+r1+7]+dphi3dphi*str->c[n+phi1+r1+7]))
                          +dz3*(
                                dri3dr*(dphii3dphi*str->c[n+z1+7]+dphi3dphi*str->c[n+phi1+z1+7])+
                                dr3dr*(dphii3dphi*str->c[n+r1+z1+7]+dphi3dphi*str->c[n+phi1+z1+r1+7])));

    /* d2f/drdz */
    B_dB[8] = rgi*zgi*(
                       (
                        (dphii*str->c[n+0]+dphi*str->c[n+phi1+0]) -
                        (dphii*str->c[n+r1+0]+dphi*str->c[n+phi1+r1+0]))
                       -(
                         (dphii*str->c[n+z1+0]+dphi*str->c[n+phi1+z1+0]) -
                         (dphii*str->c[n+r1+z1+0]+dphi*str->c[n+phi1+z1+r1+0])))
        +rg*zgi/6*(
                   -(
                     dri3dr*(dphii*str->c[n+1]+dphi*str->c[n+phi1+1])+
                     dr3dr*(dphii*str->c[n+r1+1]+dphi*str->c[n+phi1+r1+1]))
                   +(
                     dri3dr*(dphii*str->c[n+z1+1]+dphi*str->c[n+phi1+z1+1])+
                     dr3dr*(dphii*str->c[n+r1+z1+1]+dphi*str->c[n+phi1+z1+r1+1])))
        +rgi*phig2*zgi/6*(
                          (
                           (dphii3*str->c[n+2]+dphi3*str->c[n+phi1+2]) -
                           (dphii3*str->c[n+r1+2]+dphi3*str->c[n+phi1+r1+2]))
                          -(
                            (dphii3*str->c[n+z1+2]+dphi3*str->c[n+phi1+z1+2]) -
                            (dphii3*str->c[n+r1+z1+2]+dphi3*str->c[n+phi1+z1+r1+2])))
        +rgi*zg/6*(
                   dzi3dz*(
                           -(dphii*str->c[n+3]+dphi*str->c[n+phi1+3])
                           +(dphii*str->c[n+r1+3]+dphi*str->c[n+phi1+r1+3]))
                   +dz3dz*(
                           -(dphii*str->c[n+z1+3]+dphi*str->c[n+phi1+z1+3])
                           +(dphii*str->c[n+r1+z1+3]+dphi*str->c[n+phi1+z1+r1+3])))
        +rg*phig2*zgi/36*(
                          -(
                            dri3dr*(dphii3*str->c[n+4]+dphi3*str->c[n+phi1+4])+
                            dr3dr*(dphii3*str->c[n+r1+4]+dphi3*str->c[n+phi1+r1+4]))
                          +(
                            dri3dr*(dphii3*str->c[n+z1+4]+dphi3*str->c[n+phi1+z1+4])+
                            dr3dr*(dphii3*str->c[n+r1+z1+4]+dphi3*str->c[n+phi1+z1+r1+4])))
        +rg*zg/36*(
                   dzi3dz*(
                           dri3dr*(dphii*str->c[n+5]+dphi*str->c[n+phi1+5])+
                           dr3dr*(dphii*str->c[n+r1+5]+dphi*str->c[n+phi1+r1+5]))
                   +dz3dz*(
                           dri3dr*(dphii*str->c[n+z1+5]+dphi*str->c[n+phi1+z1+5])+
                           dr3dr*(dphii*str->c[n+r1+z1+5]+dphi*str->c[n+phi1+z1+r1+5])))
        +rgi*phig2*zg/36*(
                          dzi3dz*(
                                  -(dphii3*str->c[n+6]+dphi3*str->c[n+phi1+6])
                                  +(dphii3*str->c[n+r1+6]+dphi3*str->c[n+phi1+r1+6]))
                          +dz3dz*(
                                  -(dphii3*str->c[n+z1+6]+dphi3*str->c[n+phi1+z1+6])
                                  +(dphii3*str->c[n+r1+z1+6]+dphi3*str->c[n+phi1+z1+r1+6])))
        +rg*phig2*zg/216*(
                          dzi3dz*(
                                  dri3dr*(dphii3*str->c[n+7]+dphi3*str->c[n+phi1+7])+
                                  dr3dr*(dphii3*str->c[n+r1+7]+dphi3*str->c[n+phi1+r1+7]))
                          +dz3dz*(
                                  dri3dr*(dphii3*str->c[n+z1+7]+dphi3*str->c[n+phi1+z1+7])+
                                  dr3dr*(dphii3*str->c[n+r1+z1+7]+dphi3*str->c[n+phi1+z1+r1+7])));

    /* d2f/dphidz */
    B_dB[9] = phigi*zgi*(
                         (
                          dri*(str->c[n+0]  -str->c[n+phi1+0])+
                          dr*(str->c[n+r1+0]-str->c[n+phi1+r1+0]))
                         -(
                           dri*(str->c[n+z1+0]  -str->c[n+phi1+z1+0])+
                           dr*(str->c[n+r1+z1+0]-str->c[n+phi1+z1+r1+0])))
        +phigi*rg2*zgi/6*(
                          (
                           dri3*(str->c[n+1]  -str->c[n+phi1+1])+
                           dr3*(str->c[n+r1+1]-str->c[n+phi1+r1+1]))
                          -(
                            dri3*(str->c[n+z1+1]  -str->c[n+phi1+z1+1])+
                            dr3*(str->c[n+r1+z1+1]-str->c[n+phi1+z1+r1+1])))
        +phig*zgi/6*(
                     -(
                       dri*(dphii3dphi*str->c[n+2]+dphi3dphi*str->c[n+phi1+2])+
                       dr*(dphii3dphi*str->c[n+r1+2]+dphi3dphi*str->c[n+phi1+r1+2]))
                     +(
                       dri*(dphii3dphi*str->c[n+z1+2]+dphi3dphi*str->c[n+phi1+z1+2])+
                       dr*(dphii3dphi*str->c[n+r1+z1+2]+dphi3dphi*str->c[n+phi1+z1+r1+2])))
        +phigi*zg/6*(
                     dzi3dz*(
                             dri*(-str->c[n+3]+str->c[n+phi1+3])+
                             dr*(-str->c[n+r1+3]+str->c[n+phi1+r1+3]))
                     +dz3dz*(
                             dri*(-str->c[n+z1+3]+str->c[n+phi1+z1+3])+
                             dr*(-str->c[n+r1+z1+3]+str->c[n+phi1+z1+r1+3])))
        +rg2*phig*zgi/36*(
                          -(
                            dri3*(dphii3dphi*str->c[n+4]+dphi3dphi*str->c[n+phi1+4])+
                            dr3*(dphii3dphi*str->c[n+r1+4]+dphi3dphi*str->c[n+phi1+r1+4]))
                          +(
                            dri3*(dphii3dphi*str->c[n+z1+4]+dphi3dphi*str->c[n+phi1+z1+4])+
                            dr3*(dphii3dphi*str->c[n+r1+z1+4]+dphi3dphi*str->c[n+phi1+z1+r1+4])))
        +phigi*rg2*zg/36*(
                          dzi3dz*(
                                  dri3*(-str->c[n+5]+str->c[n+phi1+5])+
                                  dr3*(-str->c[n+r1+5]+str->c[n+phi1+r1+5]))
                          +dz3dz*(
                                  dri3*(-str->c[n+z1+5]+str->c[n+phi1+z1+5])+
                                  dr3*(-str->c[n+r1+z1+5]+str->c[n+phi1+z1+r1+5])))
        +phig*zg/36*(
                     dzi3dz*(
                             dri*(dphii3dphi*str->c[n+6]+dphi3dphi*str->c[n+phi1+6])+
                             dr*(dphii3dphi*str->c[n+r1+6]+dphi3dphi*str->c[n+phi1+r1+6]))
                     +dz3dz*(
                             dri*(dphii3dphi*str->c[n+z1+6]+dphi3dphi*str->c[n+phi1+z1+6])+
                             dr*(dphii3dphi*str->c[n+r1+z1+6]+dphi3dphi*str->c[n+phi1+z1+r1+6])))
        +rg2*phig*zg/216*(
                          dzi3dz*(
                                  dri3*(dphii3dphi*str->c[n+7]+dphi3dphi*str->c[n+phi1+7])+
                                  dr3*(dphii3dphi*str->c[n+r1+7]+dphi3dphi*str->c[n+phi1+r1+7]))
                          +dz3dz*(
                                  dri3*(dphii3dphi*str->c[n+z1+7]+dphi3dphi*str->c[n+phi1+z1+7])+
                                  dr3*(dphii3dphi*str->c[n+r1+z1+7]+dphi3dphi*str->c[n+phi1+z1+r1+7])));
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
void interp3Detoc_free(interp3D_data* str) {
    free(str->c);
}
