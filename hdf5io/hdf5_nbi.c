#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <hdf5.h>
#include <hdf5_hl.h>
#include "../ascot5.h"
#include "../nbi.h"
#include "hdf5_helpers.h"

int hdf5_nbi_read(hid_t f, int* n_inj_out, nbi_injector** inj_out) {
    herr_t err;

    #if VERBOSE > 0
        printf("Reading NBI input from the HDF5 file...\n");
    #endif

    err = hdf5_find_group(f, "/nbi/");
    if(err < 0) {
        return -1;
    }

    char active[11];
    err = H5LTget_attribute_string(f, "/nbi/", "active", active);
    if(err < 0) {
        return -1;
    }
    active[10] = '\0';

    #if VERBOSE > 0
        printf("Active qid is %s\n", active);
    #endif

    /* Go through all different input types and see which one the active qid
     *corresponds to. Then read this input. */
    char path[256];

    hdf5_generate_qid_path("/nbi/nbi_XXXXXXXXXX", active, path);
    if(hdf5_find_group(f, path) == 0) {
        int n_inj;

        err = H5LTread_dataset_int(f, hdf5_generate_qid_path("/nbi/nbi_XXXXXXXXXX/ninj", active, path), &n_inj);
        if(err) {printf("Error while reading HDF5 data at %s line %d\n", __FILE__, __LINE__); return -1;}

        printf("Reading %d injectors...\n", n_inj);

        nbi_injector* inj = malloc(n_inj * sizeof(nbi_injector));

        for(int i = 0; i < n_inj; i++) {
            sprintf(path, "/nbi/nbi_XXXXXXXXXX/inj%d/id", i+1);
            err = H5LTread_dataset_int(f, hdf5_generate_qid_path(path, active, path), &(inj[i].id));

            sprintf(path, "/nbi/nbi_XXXXXXXXXX/inj%d/nbeamlet", i+1);
            err = H5LTread_dataset_int(f, hdf5_generate_qid_path(path, active, path), &(inj[i].n_beamlet));

            inj[i].beamlet_x = malloc(inj[i].n_beamlet * sizeof(real));
            inj[i].beamlet_y = malloc(inj[i].n_beamlet * sizeof(real));
            inj[i].beamlet_z = malloc(inj[i].n_beamlet * sizeof(real));
            inj[i].beamlet_dx = malloc(inj[i].n_beamlet * sizeof(real));
            inj[i].beamlet_dy = malloc(inj[i].n_beamlet * sizeof(real));
            inj[i].beamlet_dz = malloc(inj[i].n_beamlet * sizeof(real));
            
            sprintf(path, "/nbi/nbi_XXXXXXXXXX/inj%d/beamletx", i+1);
            err = H5LTread_dataset_double(f, hdf5_generate_qid_path(path, active, path), inj[i].beamlet_x);
            sprintf(path, "/nbi/nbi_XXXXXXXXXX/inj%d/beamlety", i+1);
            err = H5LTread_dataset_double(f, hdf5_generate_qid_path(path, active, path), inj[i].beamlet_y);
            sprintf(path, "/nbi/nbi_XXXXXXXXXX/inj%d/beamletz", i+1);
            err = H5LTread_dataset_double(f, hdf5_generate_qid_path(path, active, path), inj[i].beamlet_z);
            sprintf(path, "/nbi/nbi_XXXXXXXXXX/inj%d/beamletdx", i+1);
            err = H5LTread_dataset_double(f, hdf5_generate_qid_path(path, active, path), inj[i].beamlet_dx);
            sprintf(path, "/nbi/nbi_XXXXXXXXXX/inj%d/beamletdy", i+1);
            err = H5LTread_dataset_double(f, hdf5_generate_qid_path(path, active, path), inj[i].beamlet_dy);
            sprintf(path, "/nbi/nbi_XXXXXXXXXX/inj%d/beamletdz", i+1);
            err = H5LTread_dataset_double(f, hdf5_generate_qid_path(path, active, path), inj[i].beamlet_dz);

            sprintf(path, "/nbi/nbi_XXXXXXXXXX/inj%d/power", i+1);
            err = H5LTread_dataset_double(f, hdf5_generate_qid_path(path, active, path), &(inj[i].power));
            sprintf(path, "/nbi/nbi_XXXXXXXXXX/inj%d/energy", i+1);
            err = H5LTread_dataset_double(f, hdf5_generate_qid_path(path, active, path), &(inj[i].energy));
            sprintf(path, "/nbi/nbi_XXXXXXXXXX/inj%d/efrac", i+1);
            err = H5LTread_dataset_double(f, hdf5_generate_qid_path(path, active, path), inj[i].efrac);
            sprintf(path, "/nbi/nbi_XXXXXXXXXX/inj%d/divergence", i+1);
            err = H5LTread_dataset_double(f, hdf5_generate_qid_path(path, active, path), inj[i].divergence);
            sprintf(path, "/nbi/nbi_XXXXXXXXXX/inj%d/anum", i+1);
            err = H5LTread_dataset_int(f, hdf5_generate_qid_path(path, active, path), &(inj[i].anum));
            sprintf(path, "/nbi/nbi_XXXXXXXXXX/inj%d/znum", i+1);
            err = H5LTread_dataset_int(f, hdf5_generate_qid_path(path, active, path), &(inj[i].znum));
        }

        *n_inj_out = n_inj;
        *inj_out = inj;
        return 1;
    }

    return -1;
}
