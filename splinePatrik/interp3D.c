/**
 * @file interp3D.c
 * @file Cubic spline interpolation of 3D magnetic field, tricubic
 */
#include <stdlib.h>
#include <stdio.h>
#include "ascot5.h"
#include "interp3D.h"
#include "spline1D.h"

/* This function calculates the interpolation coefficients for a
   tricubic spline interpolation of a 3D magnetic field.*/
void interp3D_init(interp3D_data* str, real* f, int n_r, int n_phi, int n_z,
		   real r_min, real r_max,
		   real phi_min, real phi_max,
		   real z_min, real z_max) {
    /* We initialize the struct */
    str->n_r = n_r;
    str->n_phi = n_phi;
    str->n_z = n_z;
    str->r_min = r_min;
    str->r_max = r_max;
    str->r_grid = (r_max-r_min)/(n_r-1);
    str->phi_min = phi_min;
    str->phi_max = phi_max;
    str->phi_grid = (phi_max-phi_min)/n_phi;
    str->z_min = z_min;
    str->z_max = z_max;
    str->z_grid = (z_max-z_min)/(n_z-1);
    str->c = malloc(n_phi*n_z*n_r*64*sizeof(real));
    printf("n_r = %d, r_min = %le, r_max = %le, r_grid = %le\nn_phi = %d, phi_min = %le, phi_max %le, phi_grid = %le\nn_z = %d, z_min = %le, z_max %le, z_grid = %le\n",
	   str->n_r, str->r_min, str->r_max, str->r_grid,
	   str->n_phi, str->phi_min, str->phi_max, str->phi_grid,
	   str->n_z, str->z_min, str->z_max, str->z_grid);

    /* We declare, initialize and allocate the needed variables */
    int i_r;
    int i_phi;
    int i_z;
    int i_c;
    real* f_t = malloc(n_r*sizeof(real)); // Maybe 3 custom length arrays instead? (vv)
    real* c_t = malloc(n_r*4*sizeof(real));
    int i_s;
    int i_ss;
    int i_ct;

    /* Bicubic spline surface over z and r for each phi */
    for(i_phi=0; i_phi<n_phi; i_phi++) {
	/* Cubic spline along r for each z */
	for(i_z=0; i_z<n_z; i_z++) {
	    f_t = realloc(f_t, n_r*sizeof(real));
	    for(i_r=0; i_r<n_r; i_r++) {
		/* f_t[i_r] = f[i_phi*n_z*n_r+i_z*n_r+i_r]; */
		f_t[i_r] = f[i_phi*168*92+i_z*92+i_r]; // When only part for tests
	    }
	    c_t = realloc(c_t, n_r*4*sizeof(real)); // How many slots needed?
	    spline1D(f_t,n_r,0,c_t);
	    for(i_r=0; i_r<n_r-1; i_r++) {
		for(i_c=0; i_c<4; i_c++) { // Number of indices?
		    i_ct = i_c;
		    str->c[i_phi*n_z*n_r*64+i_z*n_r*64+i_r*64+i_c] = c_t[i_r*4+i_ct];
		}
	    }
	}
	
	/* Four cubic splines along z for each r using four different data sets */
	for(i_r=0; i_r<n_r-1; i_r++) {
	    /* s0, s1, s2, s3 */
	    for(i_s=0; i_s<4; i_s++) {
		f_t = realloc(f_t, n_z*sizeof(real));
		for(i_z=0; i_z<n_z; i_z++) {
		    f_t[i_z] = str->c[i_phi*n_z*n_r*64+i_z*n_r*64+i_r*64+i_s];
		}
		c_t = realloc(c_t, n_z*4*sizeof(real)); // How many slots needed?
		spline1D(f_t,n_z,0,c_t);
		for(i_z=0; i_z<n_z-1; i_z++) {
		    i_ct = 0;
		    for(i_c=i_s; i_c<16; i_c=i_c+4) { // Number of indies?
			str->c[i_phi*n_z*n_r*64+i_z*n_r*64+i_r*64+i_c] = c_t[i_z*4+i_ct];
			i_ct++;
		    }
		}
	    }
	}
    }

    /* Cubic spline along phi for each pair of z and r to find the coefficients
       of the tricubic spline volume */
    for(i_z=0; i_z<n_z-1; i_z++) {
	for(i_r=0; i_r<n_r-1; i_r++) {
	    for(i_ss=0; i_ss<4; i_ss++) {
		for(i_s=0; i_s<4; i_s++) {
		    f_t = realloc(f_t, n_phi*sizeof(real));
		    for(i_phi=0; i_phi<n_phi; i_phi++) {
			f_t[i_phi] = str->c[i_phi*n_z*n_r*64+i_z*n_r*64
					    +i_r*64+(i_ss*4+i_s)];
		    }
		    c_t = realloc(c_t, n_phi*4*sizeof(real));
		    spline1D(f_t,n_phi,1,c_t);
		    for(i_phi=0; i_phi<n_phi; i_phi++) {
			i_ct = 0;
			for(i_c=4*i_ss+i_s; i_c<64; i_c=i_c+16) {
			    str->c[i_phi*n_z*n_r*64+i_z*n_r*64+i_r*64+i_c] 
				= c_t[i_phi*4+i_ct];
			    i_ct++;
			}
		    }
		}
	    }
	}
    }

    /* We free allocated memory */
    free(f_t);
    free(c_t);

    /* We transform from explicit to compact */
    real* cc = malloc(n_phi*n_z*n_r*8*sizeof(real));
    for(i_phi=0; i_phi<n_z; i_phi++) {
	for(i_z=0; i_z<n_z-1; i_z++) {
	    for(i_r=0; i_r<n_r-1; i_r++) {
		cc[i_phi*n_z*n_r*8+i_z*n_r*8+i_r*8+0] =
		    str->c[i_phi*n_z*n_r*64+i_z*n_r*64+i_r*64+0];
		cc[i_phi*n_z*n_r*8+i_z*n_r*8+i_r*8+1] =
		    2*str->c[i_phi*n_z*n_r*64+i_z*n_r*64+i_r*64+2];
		cc[i_phi*n_z*n_r*8+i_z*n_r*8+i_r*8+2] =
		    2*str->c[i_phi*n_z*n_r*64+i_z*n_r*64+i_r*64+32]; // 1/r2?
		cc[i_phi*n_z*n_r*8+i_z*n_r*8+i_r*8+3] =
		    2*str->c[i_phi*n_z*n_r*64+i_z*n_r*64+i_r*64+8];
		cc[i_phi*n_z*n_r*8+i_z*n_r*8+i_r*8+4] =
		    4*str->c[i_phi*n_z*n_r*64+i_z*n_r*64+i_r*64+34]; // 1/r2?
		cc[i_phi*n_z*n_r*8+i_z*n_r*8+i_r*8+5] =
		    4*str->c[i_phi*n_z*n_r*64+i_z*n_r*64+i_r*64+10];
		cc[i_phi*n_z*n_r*8+i_z*n_r*8+i_r*8+6] =
		    4*str->c[i_phi*n_z*n_r*64+i_z*n_r*64+i_r*64+40]; // 1/r2?
		cc[i_phi*n_z*n_r*8+i_z*n_r*8+i_r*8+7] =
		    8*str->c[i_phi*n_z*n_r*64+i_z*n_r*64+i_r*64+9+42]; // 1/r2?
		
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



real interp3D_eval(interp3D_data* str, real r, real phi, real z) {
    /* Excplicit */
    /* int i_r = (r-str->r_min)/str->r_grid; */
    /* real dr = (r-(str->r_min+i_r*str->r_grid))/str->r_grid; */
    /* real dr2 = dr*dr; */
    /* real dr3 = dr2*dr; */
    /* int i_phi = (phi-str->phi_min)/str->phi_grid; */ //Degrees and radians mixed!
    /* real dphi = (phi-(str->phi_min+i_phi*str->phi_grid))/str->phi_grid; */
    /* real dphi2 = dphi*dphi; */
    /* real dphi3 = dphi2*dphi; */
    /* int i_z = (z-str->z_min)/str->z_grid; */
    /* real dz = (z-(str->z_min+i_z*str->z_grid))/str->z_grid; */
    /* real dz2 = dz*dz; */
    /* real dz3 = dz2*dz; */
    /* int n = i_phi*str->n_z*str->n_r*64+i_z*str->n_r*64+i_r*64; */
    /* printf("i_r = %d, dr = %le\ni_phi = %d, dphi = %le\ni_z = %d, dz = %le\n", */
    /* 	   i_r, dr, i_phi, dphi, i_z, dz); */

    /* real val = */
    /* 	              str->c[n+ 0]+dr*str->c[n+ 1]+dr2*str->c[n+ 2]+dr3*str->c[n+ 3] */
    /* 	         +dz*(str->c[n+ 4]+dr*str->c[n+ 5]+dr2*str->c[n+ 6]+dr3*str->c[n+ 7]) */
    /* 	        +dz2*(str->c[n+ 8]+dr*str->c[n+ 9]+dr2*str->c[n+10]+dr3*str->c[n+11]) */
    /* 	        +dz3*(str->c[n+12]+dr*str->c[n+13]+dr2*str->c[n+14]+dr3*str->c[n+15]) */
    /* 	 +dphi*(      str->c[n+16]+dr*str->c[n+17]+dr2*str->c[n+18]+dr3*str->c[n+19] */
    /* 	         +dz*(str->c[n+20]+dr*str->c[n+21]+dr2*str->c[n+22]+dr3*str->c[n+23]) */
    /* 	        +dz2*(str->c[n+24]+dr*str->c[n+25]+dr2*str->c[n+26]+dr3*str->c[n+27]) */
    /* 	        +dz3*(str->c[n+28]+dr*str->c[n+29]+dr2*str->c[n+30]+dr3*str->c[n+31])) */
    /* 	+dphi2*(      str->c[n+32]+dr*str->c[n+33]+dr2*str->c[n+34]+dr3*str->c[n+35] */
    /* 	         +dz*(str->c[n+36]+dr*str->c[n+37]+dr2*str->c[n+38]+dr3*str->c[n+39]) */
    /* 	        +dz2*(str->c[n+40]+dr*str->c[n+41]+dr2*str->c[n+42]+dr3*str->c[n+43]) */
    /* 	        +dz3*(str->c[n+44]+dr*str->c[n+45]+dr2*str->c[n+46]+dr3*str->c[n+47])) */
    /* 	+dphi3*(      str->c[n+48]+dr*str->c[n+49]+dr2*str->c[n+50]+dr3*str->c[n+51] */
    /* 	         +dz*(str->c[n+52]+dr*str->c[n+53]+dr2*str->c[n+54]+dr3*str->c[n+55]) */
    /* 	        +dz2*(str->c[n+56]+dr*str->c[n+57]+dr2*str->c[n+58]+dr3*str->c[n+59]) */
    /* 	        +dz3*(str->c[n+60]+dr*str->c[n+61]+dr2*str->c[n+62]+dr3*str->c[n+63])); */

    /* return val; */


    /* Compact */
    int i_r = (r-str->r_min)/str->r_grid;
    real dr = (r-(str->r_min+i_r*str->r_grid))/str->r_grid;
    real dri = 1.0-dr;
    real dr3 = dr*dr*dr-dr;
    real dri3 = (1.0-dr)*(1.0-dr)*(1.0-dr)-(1.0-dr);
    real rg2 = str->r_grid*str->r_grid;
    phi = phi*6.283185307179586/360; // Degrees to radians
    int i_phi = (phi-str->phi_min)/str->phi_grid;
    real dphi = (phi-(str->phi_min+i_phi*str->phi_grid))/str->phi_grid;
    real dphii = 1.0-dphi;
    real dphi3 = dphi*dphi*dphi-dphi;
    real dphii3 = (1.0-dphi)*(1.0-dphi)*(1.0-dphi)-(1.0-dphi);
    real phig2 = str->phi_grid*str->phi_grid;
    int i_z = (z-str->z_min)/str->z_grid;
    real dz = (z-(str->z_min+i_z*str->z_grid))/str->z_grid;
    real dzi = 1.0-dz;
    real dz3 = dz*dz*dz-dz;
    real dzi3 = (1.0-dz)*(1.0-dz)*(1.0-dz)-(1.0-dz);
    real zg2 = str->z_grid*str->z_grid;
    int n = i_phi*str->n_z*str->n_r*8+i_z*str->n_r*8+i_r*8;
    int r1 = 8;
    int phi1 = str->n_z*str->n_r*8;
    int z1 = str->n_r*8;
    printf("i_r = %d, dr = %le\ni_phi = %d, dphi = %le\ni_z = %d, dz = %le\n",
	   i_r, dr, i_phi, dphi, i_z, dz);
	   
    real val =
	        dzi*(
		     dri*(dphii*str->c[n+0]+dphi*str->c[n+phi1+0])+
		     dr*(dphii*str->c[n+r1+0]+dphi*str->c[n+phi1+r1+0]))
	        +dz*(
		     dri*(dphii*str->c[n+z1+0]+dphi*str->c[n+phi1+z1+0])+
		     dr*(dphii*str->c[n+r1+z1+0]+dphi*str->c[n+phi1+z1+r1+0]))
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

    return val;
}



void interp3D_free(interp3D_data* str) {
    free(str->c);
}
