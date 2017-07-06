#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <hdf5.h>
#include <hdf5_hl.h>
#include "../wall.h"
#include "../wall_2d.h"
#include "../wall_3d.h"
#include "hdf5_wall.h"

void hdf5_wall_init_offload(hid_t f, wall_offload_data* offload_data, real** offload_array) {    
    herr_t err;

    err = H5LTfind_dataset(f, "/wall/");
    char type[32];
    err = H5LTget_attribute_string(f, "/wall/", "type", type);
    if(err < 0) {
        return;
    }
    if(strncmp(type,"2D",2) == 0) {
        offload_data->type = 1;
        hdf5_wall_init_offload_2D(f, &(offload_data->w2d), offload_array);
        offload_data->offload_array_length = offload_data->w2d.offload_array_length;
	
	#if VERBOSE > 0
	    printf("Loaded 2D wall (w2d)\n");
	    printf("with parameters:\n");
	#endif

    }
    else if (strcmp(type, "3D") == 0) {
        offload_data->type = 3;
        hdf5_wall_init_offload_3D(f, &(offload_data->w3d), offload_array);
        offload_data->offload_array_length = offload_data->w3d.offload_array_length;
    }
}

void hdf5_wall_init_offload_2D(hid_t f, wall_2d_offload_data* offload_data, real** offload_array) {    
    herr_t err;

    /* Read number of wall elements */
    err = H5LTget_attribute_int(f, "/wall/2D", "n", &(offload_data->n));
    offload_data->offload_array_length = 2 * offload_data->n;
    *offload_array = (real*) malloc(2 * offload_data->n * sizeof(real));

    /* Pointers to beginning of different data series to make code more
     * readable */
    real* r = &(*offload_array)[0];
    real* z = &(*offload_array)[offload_data->n];
        
    err = H5LTread_dataset_double(f,"wall/2D/r", r);
    err = H5LTread_dataset_double(f,"wall/2D/z", z);

}

/**
 * @brief Load 3D wall data and prepare parameters
 *
 * Reads a 3D wall from hdf5 file (not the same format as ASCOT4!), stores
 * parameters in offload struct and allocates and fills the offload array.
 *
 */
void hdf5_wall_init_offload_3D(hid_t f, wall_3d_offload_data* offload_data, real** offload_array) {    
    herr_t err;

    /* Read number of wall elements */
    err = H5LTget_attribute_int(f, "/wall/3D", "n", &(offload_data->n));

    *offload_array = (real*) malloc(9 * offload_data->n * sizeof(real));
    offload_data->offload_array_length = 9 * offload_data->n;

    H5LTget_attribute_double(f, "/wall/3D", "min_x", &(offload_data->xmin));
    H5LTget_attribute_double(f, "/wall/3D", "max_x", &(offload_data->xmax));
    H5LTget_attribute_double(f, "/wall/3D", "min_y", &(offload_data->ymin));
    H5LTget_attribute_double(f, "/wall/3D", "max_y", &(offload_data->ymax));
    H5LTget_attribute_double(f, "/wall/3D", "min_z", &(offload_data->zmin));
    H5LTget_attribute_double(f, "/wall/3D", "max_z", &(offload_data->zmax));

    /* Add a little bit of padding so we don't need to worry about triangles
       clipping the edges */
    offload_data->xmin -= 0.1;
    offload_data->xmax += 0.1;
    offload_data->ymin -= 0.1;
    offload_data->ymax += 0.1;
    offload_data->zmin -= 0.1;
    offload_data->zmax += 0.1;

    /* Allocate space for x1x2x3, y1y2y3, z1z2z3 for each element */
    *offload_array = (real*) malloc(9 * offload_data->n * sizeof(real));

    real x1x2x3[3*offload_data->n];
    real y1y2y3[3*offload_data->n];
    real z1z2z3[3*offload_data->n];
    
    err = H5LTread_dataset_double(f,"wall/3D/x1x2x3", x1x2x3);
    err = H5LTread_dataset_double(f,"wall/3D/y1y2y3", y1y2y3);
    err = H5LTread_dataset_double(f,"wall/3D/z1z2z3", z1z2z3);

    /* The data in the offload array is to be in the format
     *  [x1 y1 z1 x2 y2 z2 x3 y3 z3; ... ]
     */
    int i, j;
    for(i = 0; i < offload_data->n; i++) {
        for(j = 0; j < 3; j++) {
            (*offload_array)[(i+j)*3] = x1x2x3[3*i+j];
            (*offload_array)[(i+j)*3+1] = y1y2y3[3*i+j];
            (*offload_array)[(i+j)*3+2] = z1z2z3[3*i+j];
        }
    }
    
    /* Depth of the octree in which the triangles are sorted */
    offload_data->depth = 7;
    offload_data->ngrid = 1;
    for(i = 0; i < offload_data->depth - 1; i++) {
        offload_data->ngrid *= 2;
    }
    offload_data->xgrid = (offload_data->xmax - offload_data->xmin)
                          / offload_data->ngrid;
    offload_data->ygrid = (offload_data->ymax - offload_data->ymin)
                          / offload_data->ngrid;
    offload_data->zgrid = (offload_data->zmax - offload_data->zmin)
                          / offload_data->ngrid;  
}
