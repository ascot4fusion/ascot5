/**
 * @file interp2D.c
 * @file Cubic spline interpolation of 2D magnetic field component, tricubic
 */
#include <stdlib.h>
#include <stdio.h>
#include "../ascot5.h"
#include "interp2D.h"
#include "spline1D.h"

/* This function calculates the interpolation coefficients for a
   bicubic spline interpolation of a 2D magnetic field component.*/
void interp2D_init(interp2D_data* str, real* f, int n_r, int n_z,
		   real r_min, real r_max,
		   real z_min, real z_max) {
    
    /* We initialize the struct */
    str->n_r = n_r;
    str->n_z = n_z;
    str->r_min = r_min;
    str->r_max = r_max;
    str->r_grid = (r_max-r_min)/(n_r-1);
    str->z_min = z_min;
    str->z_max = z_max;
    str->z_grid = (z_max-z_min)/(n_z-1);
    str->c = malloc(n_z*n_r*16*sizeof(real));
    printf("Here n_r = %d, n_z = %d\n", n_r, n_z);

    /* We declare, initialize and allocate the needed variables */
    int i_r;
    int i_z;
    int i_c;
    real* f_r = malloc(n_r*sizeof(real));
    real* f_z = malloc(n_z*sizeof(real));
    real* c_r = malloc(n_r*4*sizeof(real)); // How many slots needed?
    real* c_z = malloc(n_z*4*sizeof(real)); // How many slots needed?
    int i_s;
    int i_ct;

    /* Bicubic spline surface over z and r */
    /* Cubic spline along r for each z */
    for(i_z=0; i_z<n_z; i_z++) {
	for(i_r=0; i_r<n_r; i_r++) {
	    f_r[i_r] = f[i_z*n_r+i_r];
	}
	spline1D(f_r,n_r,0,c_r);
	for(i_r=0; i_r<n_r-1; i_r++) {
	    for(i_c=0; i_c<4; i_c++) { // Number of indices?
		i_ct = i_c;
		str->c[i_z*n_r*16+i_r*16+i_c] = c_r[i_r*4+i_ct];
	    }
	}
    }
    
    /* Four cubic splines along z for each r using four different data sets */
    for(i_r=0; i_r<n_r-1; i_r++) {
	/* s0, s1, s2, s3 */
	for(i_s=0; i_s<4; i_s++) {
	    for(i_z=0; i_z<n_z; i_z++) {
		f_z[i_z] = str->c[i_z*n_r*16+i_r*16+i_s];
	    }
	    spline1D(f_z,n_z,0,c_z);
	    for(i_z=0; i_z<n_z-1; i_z++) {
		i_ct = 0;
		for(i_c=i_s; i_c<16; i_c=i_c+4) { // Number of indices?
		    str->c[i_z*n_r*16+i_r*16+i_c] = c_z[i_z*4+i_ct];
		    i_ct++;
		}
	    }
	}
    }
    
    /* We free allocated memory */
    free(f_r);
    free(f_z);
    free(c_r);
    free(c_z);

    /* We transform from explicit to compact */
    real* cc = malloc(n_z*n_r*4*sizeof(real));
    for(i_z=0; i_z<n_z-1; i_z++) {
    	for(i_r=0; i_r<n_r-1; i_r++) {
    	    cc[i_z*n_r*4+i_r*4  ] =   str->c[i_z*n_r*16+i_r*16   ];
    	    cc[i_z*n_r*4+i_r*4+1] = 2*str->c[i_z*n_r*16+i_r*16+ 2];
    	    cc[i_z*n_r*4+i_r*4+2] = 2*str->c[i_z*n_r*16+i_r*16+ 8];
    	    cc[i_z*n_r*4+i_r*4+3] = 4*str->c[i_z*n_r*16+i_r*16+10];
    	}
    	cc[i_z*n_r*4+(n_r-1)*4  ] = f[i_z*n_r+n_r-1];
    	cc[i_z*n_r*4+(n_r-1)*4+1] = 0;
    	cc[i_z*n_r*4+(n_r-1)*4+2] = 0; // Is this correct?
    	cc[i_z*n_r*4+(n_r-1)*4+3] = 0;
    }
    for(i_r=0; i_r<n_r; i_r++) {
    	cc[(n_z-1)*n_r*4+i_r*4  ] = f[(n_z-1)*n_r+i_r];
    	cc[(n_z-1)*n_r*4+i_r*4+1] = 0; // Is this correct?
    	cc[(n_z-1)*n_r*4+i_r*4+2] = 0;
    	cc[(n_z-1)*n_r*4+i_r*4+3] = 0;
    }
    free(str->c);
    str->c = cc;
}



void interp2D_eval_B(real* B, interp2D_data* str, real r, real z) {
    int i_r = (r-str->r_min)/str->r_grid;
    real dr = (r-(str->r_min+i_r*str->r_grid))/str->r_grid;
    real dr3 = dr*dr*dr-dr;
    real dri = 1.0-dr;
    real dri3 = dri*dri*dri-dri;
    real rg2 = str->r_grid*str->r_grid;
    int i_z = (z-str->z_min)/str->z_grid;
    real dz = (z-(str->z_min+i_z*str->z_grid))/str->z_grid;
    real dz3 = dz*dz*dz-dz;
    real dzi = 1.0-dz;
    real dzi3 = dzi*dzi*dzi-dzi;
    real zg2 = str->z_grid*str->z_grid;
    int n = i_z*str->n_r*4+i_r*4;
    int z1 = str->n_r*4;

    *B = (
	  dri*( dzi*str->c[n]+ dz*str->c[n+str->n_r*4] )
	  +dr*( dzi*str->c[n+4]+ dz*str->c[n+str->n_r*4+4] ) )
    	+rg2/6*(
		dri3*( dzi*str->c[n+1]+ dz*str->c[n+str->n_r*4+1] )
		+dr3*( dzi*str->c[n+4+1]+ dz*str->c[n+str->n_r*4+4+1] ) )
    	+zg2/6*(
		dri*( dzi3*str->c[n+2]+dz3*str->c[n+str->n_r*4+2] )
		+dr*( dzi3*str->c[n+4+2]+dz3*str->c[n+str->n_r*4+4+2] ) )
    	+rg2*zg2/36*(
		     dri3*( dzi3*str->c[n+3]+dz3*str->c[n+str->n_r*4+3] )
		     +dr3*( dzi3*str->c[n+4+3]+dz3*str->c[n+str->n_r*4+4+3] ) );
}



void interp2D_eval_dB(real* B_dB, interp2D_data* str, real r, real z) {
    int i_r = (r-str->r_min)/str->r_grid;
    real dr = (r-(str->r_min+i_r*str->r_grid))/str->r_grid;
    real dr3 = dr*dr*dr-dr;
    real dr3dr = 3*dr*dr-1;
    real dri = 1.0-dr;
    real dri3 = dri*dri*dri-dri;
    real dri3dr = -3*(1-dr)*(1-dr)+1;
    real rg = str->r_grid;
    real rg2 = rg*rg;
    real rgi = 1.0/rg;
    int i_z = (z-str->z_min)/str->z_grid;
    real dz = (z-(str->z_min+i_z*str->z_grid))/str->z_grid;
    real dz3 = dz*dz*dz-dz;
    real dz3dz = 3*dz*dz-1;
    real dzi = 1.0-dz;
    real dzi3 = dzi*dzi*dzi-dzi;
    real dzi3dz = -3*(1-dz)*(1-dz)+1;
    real zg = str->z_grid;
    real zg2 = zg*zg;
    real zgi = 1.0/zg;
    int n = i_z*str->n_r*4+i_r*4;
    int z1 = str->n_r*4;

    /* f */
    B_dB[0] = (
	      dri*( dzi*str->c[n]+dz*str->c[n+z1] )+
	      dr*( dzi*str->c[n+4]+dz*str->c[n+z1+4] ) )
	+rg2/6*(
		dri3*( dzi*str->c[n+1]+ dz*str->c[n+z1+1] )+
		dr3*( dzi*str->c[n+4+1]+ dz*str->c[n+z1+4+1] ) )
	+zg2/6*(
		dri*( dzi3*str->c[n+2]+dz3*str->c[n+z1+2] )+
		dr*( dzi3*str->c[n+4+2]+dz3*str->c[n+z1+4+2] ) )
	+rg2*zg2/36*(
		     dri3*( dzi3*str->c[n+3]+dz3*str->c[n+z1+3] )+
		     dr3*( dzi3*str->c[n+4+3]+dz3*str->c[n+z1+4+3] ) );

    /* df/dr */
    B_dB[1] = rgi*(
		  -(dzi*str->c[n]  +dz*str->c[n+z1])
		  +(dzi*str->c[n+4]+dz*str->c[n+z1+4]))
	+rg/6*(
	       dri3dr*(dzi*str->c[n+1]  +dz*str->c[n+z1+1])+
	       dr3dr*(dzi*str->c[n+4+1]+dz*str->c[n+z1+4+1]))
	+rgi*zg2/6*(
		    -(dzi3*str->c[n+2]  +dz3*str->c[n+z1+2])
		    +(dzi3*dzi3*str->c[n+4+2]+dz3*str->c[n+z1+4+2]))
	+rg*zg2/36*(
		    dri3dr*(dzi3*str->c[n+3]  +dz3*str->c[n+z1+3])+
		    dr3dr*(dzi3*str->c[n+4+3]+dz3*str->c[n+z1+4+3]));

    /* dr/dz */
    B_dB[2] = zgi*(
		  dri*(-str->c[n]  +str->c[n+z1])+
		  dr*(-str->c[n+4]+str->c[n+z1+4]))
	+rg2*zgi/6*(
                    dri3*(-str->c[n+1]  +str->c[n+z1+1])+
                    dr3*(-str->c[n+4+1]+str->c[n+z1+4+1]))
	+zg/6*(
	       dri*(dzi3dz*str->c[n+2]  +dz3dz*str->c[n+z1+2])+
	       dr*(dzi3dz*dzi3*str->c[n+4+2]+dz3dz*str->c[n+z1+4+2]))
	+rg2*zg/36*(
                    dri3*(dzi3dz*str->c[n+3]  +dz3dz*str->c[n+z1+3])+
                    dr3*(dzi3dz*str->c[n+4+3]+dz3dz*str->c[n+z1+4+3]));

    /* d2f/dr2 */
    B_dB[3] = (
	      dri*(dzi*str->c[n+1]  +dz*str->c[n+z1+1])+
	      dr*(dzi*str->c[n+4+1]+dz*str->c[n+z1+4+1]))
	+zg2/6*(
		dri*(dzi3*str->c[n+3]  +dz3*str->c[n+z1+3])+
		dr*(dzi3*str->c[n+4+3]+dz3*str->c[n+z1+4+3]));

    /* d2f/dz2 */
    B_dB[4] = (
	      dri*(dzi*str->c[n+2]  +dz*str->c[n+z1+2])+
	      dr*(dzi*dzi3*str->c[n+4+2]+dz*str->c[n+z1+4+2]))
	+rg2/6*(
		dri3*(dzi*str->c[n+3]  +dz*str->c[n+z1+3])+
		dr3*(dzi*str->c[n+4+3]+dz*str->c[n+z1+4+3]));

    /* d2f/dzdr */
    B_dB[5] = rgi*zgi*(
		      str->c[n]  -str->c[n+z1]
		      -str->c[n+4]+str->c[n+z1+4])
	+rg/6*zgi*(
		   dri3dr*(-str->c[n+1]  +str->c[n+z1+1])+
		   dr3dr*(-str->c[n+4+1]+str->c[n+z1+4+1]))
	+rgi/6*zg*(
		   -(dzi3dz*str->c[n+2]  +dz3dz*str->c[n+z1+2])
		   +(dzi3dz*dzi3*str->c[n+4+2]+dz3dz*str->c[n+z1+4+2]))
	+rg*zg/36*(
		   dri3dr*(dzi3dz*str->c[n+3]  +dz3dz*str->c[n+z1+3])+
		   dr3dr*(dzi3dz*str->c[n+4+3]+dz3dz*str->c[n+z1+4+3]));
}



void interp2D_free(interp2D_data* str) {
    free(str->c);
}
