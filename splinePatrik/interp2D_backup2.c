/**
 * @file interp2D.c
 * @file Cubic spline interpolation of 2D magnetic field component, tricubic
 */
#include <stdlib.h>
#include <stdio.h>
#include "ascot5.h"
#include "interp2D.h"
#include "spline1D.h"

/* This function calculates the interpolation coefficients for a
   bicubic spline interpolation of a 2D magnetic field component.*/
void interp2D_init(interp2D_data* str, real* f, int n_r, int n_z,
		   real r_min, real r_max,
		   real z_min, real z_max) {
    
    printf("The first data value is %le\n", f[0]);
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
    printf("n_r = %d, r_min = %le, r_max = %le, r_grid = %le\nn_z = %d, z_min = %le, z_max %le, z_grid = %le\n",
	   str->n_r, str->r_min, str->r_max, str->r_grid,
	   str->n_z, str->z_min, str->z_max, str->z_grid);

    /* We declare, initialize and allocate the needed variables */
    int i_r;
    int i_z;
    int i_c;
    real* f_t = malloc(n_r*sizeof(real)); // Maybe 2 custom length arrays instead? (vv)
    real* c_t = malloc(n_r*4*sizeof(real));
    int i_s;
    int i_ct;

    /* Bicubic spline surface over z and r */
    /* Cubic spline along r for each z */
    for(i_z=0; i_z<n_z; i_z++) {
	f_t = realloc(f_t, n_r*sizeof(real));
	for(i_r=0; i_r<n_r; i_r++) {
	    f_t[i_r] = f[i_z*n_r+i_r];
	}
	c_t = realloc(c_t, n_r*4*sizeof(real)); // How many slots needed?
	spline1D(f_t,n_r,0,c_t);
	for(i_r=0; i_r<n_r-1; i_r++) {
	    for(i_c=0; i_c<4; i_c++) { // Number of indices?
		i_ct = i_c;
		str->c[i_z*n_r*16+i_r*16+i_c] = c_t[i_r*4+i_ct];
	    }
	}
    }
    
    /* Four cubic splines along z for each r using four different data sets */
    for(i_r=0; i_r<n_r-1; i_r++) {
	/* s0, s1, s2, s3 */
	for(i_s=0; i_s<4; i_s++) {
	    f_t = realloc(f_t, n_z*sizeof(real));
	    for(i_z=0; i_z<n_z; i_z++) {
		f_t[i_z] = str->c[i_z*n_r*16+i_r*16+i_s];
	    }
	    c_t = realloc(c_t, n_z*4*sizeof(real)); // How many slots needed?
	    spline1D(f_t,n_z,0,c_t);
	    for(i_z=0; i_z<n_z-1; i_z++) {
		i_ct = 0;
		for(i_c=i_s; i_c<16; i_c=i_c+4) { // Number of indices?
		    str->c[i_z*n_r*16+i_r*16+i_c] = c_t[i_z*4+i_ct];
		    i_ct++;
		}
	    }
	}
    }
    
    /* We free allocated memory */
    free(f_t);
    free(c_t);

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



real interp2D_eval(interp2D_data* str, real r, real z) {
    /* Explicit */
    /* int i_r = (r-str->r_min)/str->r_grid; */
    /* real dr = (r-(str->r_min+i_r*str->r_grid))/str->r_grid; */
    /* real dr2 = dr*dr; */
    /* real dr3 = dr2*dr; */
    /* int i_z = (z-str->z_min)/str->z_grid; */
    /* real dz = (z-(str->z_min+i_z*str->z_grid))/str->z_grid; */
    /* real dz2 = dz*dz; */
    /* real dz3 = dz2*dz; */
    /* int n = i_z*str->n_r*16+i_r*16; */
    /* printf("i_r = %d, dr = %le\ni_z = %d, dz = %le\n", */
    /* 	   i_r, dr, i_z, dz); */

    /* real val = */
    /* 	              str->c[n+ 0]+dr*str->c[n+ 1]+dr2*str->c[n+ 2]+dr3*str->c[n+ 3] */
    /* 	         +dz*(str->c[n+ 4]+dr*str->c[n+ 5]+dr2*str->c[n+ 6]+dr3*str->c[n+ 7]) */
    /* 	        +dz2*(str->c[n+ 8]+dr*str->c[n+ 9]+dr2*str->c[n+10]+dr3*str->c[n+11]) */
    /* 	        +dz3*(str->c[n+12]+dr*str->c[n+13]+dr2*str->c[n+14]+dr3*str->c[n+15]); */

    /* return val; */

    /* Compact */
    int i_r = (r-str->r_min)/str->r_grid;
    real dr = (r-(str->r_min+i_r*str->r_grid))/str->r_grid;
    real r_r = dr*dr*dr-dr;
    real s_r = (1-dr)*(1-dr)*(1-dr)-(1-dr);
    int i_z = (z-str->z_min)/str->z_grid;
    real dz = (z-(str->z_min+i_z*str->z_grid))/str->z_grid;
    real r_z = dz*dz*dz-dz;
    real s_z = (1-dz)*(1-dz)*(1-dz)-(1-dz);
    int n = i_z*str->n_r*4+i_r*4;
    printf("i_r = %d, dr = %le\ni_z = %d, dz = %le\n",
    	   i_r, dr, i_z, dz);

    real val =   (1-dr)*( (1-dz)*str->c[n    ]+ dz*str->c[n+str->n_r*4    ] )
    	            +dr*( (1-dz)*str->c[n+4  ]+ dz*str->c[n+str->n_r*4+4  ] )
    	+str->r_grid*str->r_grid/6
    	      *(    s_r*( (1-dz)*str->c[n  +1]+ dz*str->c[n+str->n_r*4  +1] )
    		   +r_r*( (1-dz)*str->c[n+4+1]+ dz*str->c[n+str->n_r*4+4+1] ) )
    	+str->z_grid*str->z_grid/6
    	      *( (1-dr)*(    s_z*str->c[n  +2]+r_z*str->c[n+str->n_r*4  +2] )
    		    +dr*(    s_z*str->c[n+4+2]+r_z*str->c[n+str->n_r*4+4+2] ) )
    	+str->r_grid*str->r_grid*str->z_grid*str->z_grid/36
    	      *(    s_r*(    s_z*str->c[n  +3]+r_z*str->c[n+str->n_r*4  +3] )
    		   +r_r*(    s_z*str->c[n+4+3]+r_z*str->c[n+str->n_r*4+4+3] ) );
		  
    return val;
}



void interp2D_free(interp2D_data* str) {
    free(str->c);
}
