#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <hdf5.h>
#include <hdf5_hl.h>
#include "../wall.h"
#include "../wall/wall_2d.h"
#include "../wall/wall_3d.h"
#include "hdf5_wall.h"
#include "hdf5_helpers.h"

int hdf5_wall_init_offload(hid_t f, wall_offload_data* offload_data, real** offload_array) {    
    herr_t err;

    #if VERBOSE > 0
        printf("Reading wall input from the HDF5 file...\n");
    #endif

    err = hdf5_find_group(f, "/wall/");
    if(err < 0) {
        return -1;
    }

    char active[11];
    err = H5LTget_attribute_string(f, "/wall/", "active", active);
    if(err < 0) {
        return -1;
    }
    active[10] = '\0';
    
    #if VERBOSE > 0
        printf("Active qid is %s\n", active);
    #endif
    
    /* Go through all different input types and see which one the active qid corresponds to.
     * Then read this input. */
    char path[256];

    hdf5_generate_qid_path("/wall/wall_2D-XXXXXXXXXX", active, path);
    if(hdf5_find_group(f, path) == 0) {
	
        offload_data->type = 1;
        hdf5_wall_init_offload_2D(f, &(offload_data->w2d), offload_array, active);
        offload_data->offload_array_length = offload_data->w2d.offload_array_length;
	
	#if VERBOSE > 0
	    printf("Loaded 2D wall (wall_2D)\n");
	    printf("with parameters:\n");
	    printf("- number of wall elements = %d\n",
		   offload_data->w2d.n);
	#endif
	return 1;
    }

    hdf5_generate_qid_path("/wall/wall_3D-XXXXXXXXXX", active, path);
    if (hdf5_find_group(f, path) == 0) {
        offload_data->type = 3;
        hdf5_wall_init_offload_3D(f, &(offload_data->w3d), offload_array, active);
        offload_data->offload_array_length = offload_data->w3d.offload_array_length;
	#if VERBOSE > 0
	    printf("\nLoaded 3D wall (wall_3D)\n");
	    printf("with parameters:\n");
	    printf("- number of wall elements = %d\n",
		   offload_data->w3d.n);
	    printf("- xmin = %le and xmax = %le\n",
		   offload_data->w3d.xmin,offload_data->w3d.xmax);
	    printf("- ymin = %le and ymax = %le\n",
		   offload_data->w3d.ymin,offload_data->w3d.ymax);
	    printf("- zmin = %le and zmax = %le\n",
		   offload_data->w3d.zmin,offload_data->w3d.zmax);

        #endif
       return 1;
    }
    
    return -1;
}

void hdf5_wall_init_offload_2D(hid_t f, wall_2d_offload_data* offload_data, real** offload_array, char* qid) {    
    herr_t err;
    char path[256];

    /* Read number of wall elements */
    err = H5LTread_dataset_int(f, hdf5_generate_qid_path("/wall/wall_2D-XXXXXXXXXX/n", qid, path), &(offload_data->n));
    offload_data->offload_array_length = 2 * offload_data->n;
    *offload_array = (real*) malloc(2 * offload_data->n * sizeof(real));

    /* Pointers to beginning of different data series to make code more
     * readable */
    real* r = &(*offload_array)[0];
    real* z = &(*offload_array)[offload_data->n];
        
    err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/wall/wall_2D-XXXXXXXXXX/r", qid, path), r);
    err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/wall/wall_2D-XXXXXXXXXX/z", qid, path), z);

}

/**
 * @brief Load 3D wall data and prepare parameters
 *
 * Reads a 3D wall from hdf5 file (not the same format as ASCOT4!), stores
 * parameters in offload struct and allocates and fills the offload array.
 *
 */
void hdf5_wall_init_offload_3D(hid_t f, wall_3d_offload_data* offload_data, real** offload_array, char* qid) {    
    herr_t err;
    char path[256];

    /* Read number of wall elements */
    err = H5LTread_dataset_int(f, hdf5_generate_qid_path("/wall/wall_3D-XXXXXXXXXX/n", qid, path), &(offload_data->n));

    err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/wall/wall_3D-XXXXXXXXXX/min_x", qid, path), &(offload_data->xmin));
    err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/wall/wall_3D-XXXXXXXXXX/max_x", qid, path), &(offload_data->xmax));
    err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/wall/wall_3D-XXXXXXXXXX/min_y", qid, path), &(offload_data->ymin));
    err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/wall/wall_3D-XXXXXXXXXX/max_y", qid, path), &(offload_data->ymax));
    err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/wall/wall_3D-XXXXXXXXXX/min_z", qid, path), &(offload_data->zmin));
    err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/wall/wall_3D-XXXXXXXXXX/max_z", qid, path), &(offload_data->zmax));

    /* Add a little bit of padding so we don't need to worry about triangles
       clipping the edges */
    offload_data->xmin -= 0.1;
    offload_data->xmax += 0.1;
    offload_data->ymin -= 0.1;
    offload_data->ymax += 0.1;
    offload_data->zmin -= 0.1;
    offload_data->zmax += 0.1;

    /* Allocate space for x1x2x3, y1y2y3, z1z2z3 for each element */

    offload_data->offload_array_length = 9 * offload_data->n;
    *offload_array = (real*) malloc(9 * offload_data->n * sizeof(real));

    
    real* x1x2x3 = (real*)malloc(3 * offload_data->n * sizeof(real));
    real* y1y2y3 = (real*)malloc(3 * offload_data->n * sizeof(real));
    real* z1z2z3 = (real*)malloc(3 * offload_data->n * sizeof(real));
    
    err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/wall/wall_3D-XXXXXXXXXX/x1x2x3", qid, path), x1x2x3);
    err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/wall/wall_3D-XXXXXXXXXX/y1y2y3", qid, path), y1y2y3);
    err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/wall/wall_3D-XXXXXXXXXX/z1z2z3", qid, path), z1z2z3);

    /* The data in the offload array is to be in the format
     *  [x1 y1 z1 x2 y2 z2 x3 y3 z3; ... ]
     */
    int i, j;
    for(i = 0; i < offload_data->n; i++) {
        for(j = 0; j < 3; j++) {
            (*offload_array)[i*9 + j*3 + 0] = x1x2x3[3*i+j];
            (*offload_array)[i*9 + j*3 + 1] = y1y2y3[3*i+j];
            (*offload_array)[i*9 + j*3 + 2] = z1z2z3[3*i+j];
        }
    }
    
    free(x1x2x3);
    free(y1y2y3);
    free(z1z2z3);
    
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
