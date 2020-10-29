/**
 * @file hdf5_helpers.c
 * @brief HDF5 wrappers and shortcuts
 */
#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <hdf5.h>
#include <hdf5_hl.h>
#include "../ascot5.h"
#include "../print.h"

/**
 * @brief Initialize hdf5, right now just disables automatic error messages.
 */
void hdf5_init(void) {
    H5Eset_auto1(NULL, NULL);
}

/**
 * @brief Create an hdf5 file, fail if file exists. A negative value is
 * returned on failure.
 */
hid_t hdf5_create(const char* filename) {
    hid_t file;
    file = H5Fcreate(filename, H5F_ACC_EXCL, H5P_DEFAULT, H5P_DEFAULT);
    return file;
}

/**
 * @brief Open a hdf5 file for reading and writing. A negative value is
 *returned on failure.
 */
hid_t hdf5_open(const char *filename) {
    hid_t file;
    file = H5Fopen(filename, H5F_ACC_RDWR, H5P_DEFAULT);
    return file;
}

/**
 * @brief Open a hdf5 file for reading and writing. A negative value is
 *returned on failure.
 */
hid_t hdf5_open_ro(const char *filename) {
    hid_t file;
    file = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
    return file;
}

/**
 * @brief Close access to given hdf5 file identifier. A negative value is
 *returned on failure.
 */
herr_t hdf5_close(hid_t file_id) {
    herr_t err;
    err = H5Fclose(file_id);
    return err;
}

/**
 * @brief Create a group (with parent groups if necessary). Returns a handle to
 * the group. Negative on failure.
 */
hid_t hdf5_create_group(hid_t loc, const char* path) {
    const char* start;
    if(path[0] == '/') {
        start = path + 1;
    } else {
        start = path;
    }

    const char* end = strstr(start, "/");

    if(end == NULL) {
        hid_t g = H5Gcreate2(loc, start, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        return g;
    } else {
        char* group = (char*) malloc((strlen(path) + 1) * sizeof(char));
        strncpy(group, start, end-start);
        group[end-start] = '\0';
        hid_t g = H5Gcreate2(loc, group, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        free(group);
        return hdf5_create_group(g, end);
    }
}

/**
 * @brief Checks if given group exists within given hdf5 file. Negative value is
 * returned if the group doesn't exist.
 */
herr_t hdf5_find_group(hid_t loc, const char* path) {
    return H5Gget_objinfo (loc, path, 0, NULL);
}

/**
 * @brief Generate a valid path from a given template and qid.
 *
 * Data in ASCOT5 HDF5 files is stored in groups, where each group is assigned
 * a unique identifier. The paths are then in format such as "bfield/B_2D-0123456789".
 * This function turns a template e.g. "bfield/B_2D-XXXXXXXXXX" to a valid path.
 */
char* hdf5_generate_qid_path(const char* original, char* qid, char* path) {
    strcpy(path, original);
    char* ptr = strstr(path,"XXXXXXXXXX");

    for(int i = 0; i < 10; i++) {
        ptr[i] = qid[i];
    }

    return path;
}

/**
 * @brief Generate a valid path from a given template and qid.
 *
 * Data in ASCOT5 HDF5 files is stored in groups, where each group is assigned
 * a unique identifier. The paths are then in format such as "bfield/B_2D-0123456789".
 * This function turns a template e.g. "bfield/B_2D-XXXXXXXXXX" to a valid path.
 */
char* hdf5_gen_path(const char* original, char* qid, char* path) {
    strcpy(path, original);
    char* ptr = strstr(path,"XXXXXXXXXX");

    for(int i = 0; i < 10; i++) {
        ptr[i] = qid[i];
    }

    return path;
}

/**
 * @brief Read double-valued data from ASCOT5 HDF5 file.
 *
 * Datasets in ASCOT5 are located in the HDF5 file at paths that look like this:
 * "/options/opt-XXXXXXXXXX/dataset", where X's is the QID value that is used
 * differentiate between inputs of same type. This function constructs the path
 * for a given dummy path, dataset, and qid, reads the dataset, and writes error
 * message if the dataset could not be read.
 *
 * This function reads data that has type double.
 *
 * The file is opened and closed outside of this function.
 *
 * @param var "dummy" (otherwise valid but with X's) path to variable
 * @param ptr pointer where data will be stored
 * @param file HDF5 file
 * @param qid QID value
 * @param errfile use macro __FILE__ here to indicate the file this function
 *        was called
 * @param errline use macro __LINE__ here to indicate the line this function
 *        was called
 *
 * @return zero if success
 */
int hdf5_read_double(const char* var, real* ptr, hid_t file, char* qid,
                     const char* errfile, int errline) {
    char temp[256];
    if( H5LTread_dataset_double(file, hdf5_gen_path(var, qid, temp), ptr) < 0 ){
        print_err("Error: could not read HDF5 dataset %s FILE %s LINE %d\n",
                  temp, errfile, errline);
        return 1;
    }
    return 0;
}

/**
 * @brief Read int-valued data from ASCOT5 HDF5 file.
 *
 * Datasets in ASCOT5 are located in the HDF5 file at paths that look like this:
 * "/options/opt-XXXXXXXXXX/dataset", where X's is the QID value that is used
 * differentiate between inputs of same type. This function constructs the path
 * for a given dummy path, dataset, and qid, reads the dataset, and writes error
 * message if the dataset could not be read.
 *
 * This function reads data that has type double.
 *
 * The file is opened and closed outside of this function.
 *
 * @param var "dummy" (otherwise valid but with X's) path to variable
 * @param ptr pointer where data will be stored
 * @param file HDF5 file
 * @param qid QID value
 * @param errfile use macro __FILE__ here to indicate the file this function
 *        was called
 * @param errline use macro __LINE__ here to indicate the line this function
 *        was called
 *
 * @return zero if success
 */
int hdf5_read_int(const char* var, int* ptr, hid_t file, char* qid,
                  const char* errfile, int errline) {
    char temp[256];
    if( H5LTread_dataset_int(file, hdf5_gen_path(var, qid, temp), ptr) < 0 ){
        print_err("Error: could not read HDF5 dataset %s FILE %s LINE %d\n",
                  temp, errfile, errline);
        return 1;
    }
    return 0;
}

/**
 * @brief Read long-valued data from ASCOT5 HDF5 file.
 *
 * Datasets in ASCOT5 are located in the HDF5 file at paths that look like this:
 * "/options/opt-XXXXXXXXXX/dataset", where X's is the QID value that is used
 * differentiate between inputs of same type. This function constructs the path
 * for a given dummy path, dataset, and qid, reads the dataset, and writes error
 * message if the dataset could not be read.
 *
 * This function reads data that has type double.
 *
 * The file is opened and closed outside of this function.
 *
 * @param var "dummy" (otherwise valid but with X's) path to variable
 * @param ptr pointer where data will be stored
 * @param file HDF5 file
 * @param qid QID value
 * @param errfile use macro __FILE__ here to indicate the file this function
 *        was called
 * @param errline use macro __LINE__ here to indicate the line this function
 *        was called
 *
 * @return zero if success
 */
int hdf5_read_long(const char* var, long* ptr, hid_t file, char* qid,
                   const char* errfile, int errline) {
    char temp[256];
    if( H5LTread_dataset_long(file, hdf5_gen_path(var, qid, temp), ptr) < 0 ){
        print_err("Error: could not read HDF5 dataset %s FILE %s LINE %d\n",
                  temp, errfile, errline);
        return 1;
    }
    return 0;
}

/**
 * @brief Write string attribute with null-padding.
 *
 * There is a H5LTset_attribute_string function but it writes strings as null-terminated. However, string
 * attributes in ASCOT5 HDF5 file are assumed to be null-padded.
 */
herr_t hdf5_write_string_attribute(hid_t loc, const char* path, const char* attrname,  const char* string) {
    herr_t err;

    hid_t grp = H5Gopen(loc, path, H5P_DEFAULT);
    hid_t aid = H5Screate(H5S_SCALAR);
    hid_t atype = H5Tcopy(H5T_C_S1);
    H5Tset_size(atype, strlen(string));
    H5Tset_strpad(atype,H5T_STR_NULLPAD);

    hid_t attr = H5Aopen(grp, attrname, H5P_DEFAULT);
    if(attr > 0) {
        H5Adelete(grp, attrname);
    }

    attr = H5Acreate2(grp, attrname, atype, aid, H5P_DEFAULT, H5P_DEFAULT);

    err = H5Awrite(attr, atype, string);
    if(err) {return err;}

    err = H5Sclose(aid);
    if(err) {return err;}

    err = H5Tclose(atype);
    if(err) {return err;}

    err = H5Aclose(attr);
    if(err) {return err;}

    err = H5Gclose(grp);
    if(err) {return err;}

    return 0;
}

/**
 * @brief Create and write to an extendible dataset for double data.
 */
herr_t hdf5_write_extendible_dataset_double(hid_t group,
                                            const char* datasetname,
                                            int length, double* data) {
    /* Create the data space with unlimited dimensions. */
    hsize_t dim[1]    = {length};
    hsize_t maxdim[1] = {H5S_UNLIMITED};
    hid_t dataspace   = H5Screate_simple(1, dim, maxdim);

    /* Modify dataset creation properties, i.e. enable chunking  */
    hsize_t chunk_dim[1] = {(int)ceil(length/2.0)};
    hid_t prop   = H5Pcreate(H5P_DATASET_CREATE);
    H5Pset_chunk (prop, 1, chunk_dim);

    /* Create a new dataset within the file using chunk creation properties. */
    hid_t dataset = H5Dcreate2(group, datasetname, H5T_IEEE_F64LE, dataspace,
                               H5P_DEFAULT, prop, H5P_DEFAULT);

    /* Write data to dataset */
    if(H5Dwrite(dataset, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data)) {
        return -1;
    }
    return 0;
}


/**
 * @brief Create and write to an extendible dataset for long data.
 */
herr_t hdf5_write_extendible_dataset_long(hid_t group,
                                          const char* datasetname,
                                          int length, long* data) {
    /* Create the data space with unlimited dimensions. */
    hsize_t dim[1]    = {length};
    hsize_t maxdim[1] = {H5S_UNLIMITED};
    hid_t dataspace   = H5Screate_simple(1, dim, maxdim);

    /* Modify dataset creation properties, i.e. enable chunking  */
    hsize_t chunk_dim[1] = {(int)ceil(length/2.0)};
    hid_t prop   = H5Pcreate(H5P_DATASET_CREATE);
    H5Pset_chunk (prop, 1, chunk_dim);

    /* Create a new dataset within the file using chunk creation properties. */
    hid_t dataset = H5Dcreate2(group, datasetname, H5T_STD_I64LE, dataspace,
                               H5P_DEFAULT, prop, H5P_DEFAULT);

    /* Write data to dataset */
    if(H5Dwrite(dataset, H5T_STD_I64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data)) {
        return -1;
    }
    return 0;
}

/**
 * @brief Create and write to an extendible dataset int data.
 */
herr_t hdf5_write_extendible_dataset_int(hid_t group,
                                         const char* datasetname,
                                         int length, int* data) {
    /* Create the data space with unlimited dimensions. */
    hsize_t dim[1]    = {length};
    hsize_t maxdim[1] = {H5S_UNLIMITED};
    hid_t dataspace   = H5Screate_simple(1, dim, maxdim);

    /* Modify dataset creation properties, i.e. enable chunking  */
    hsize_t chunk_dim[1] = {(int)ceil(length/2.0)};
    hid_t prop   = H5Pcreate(H5P_DATASET_CREATE);
    H5Pset_chunk (prop, 1, chunk_dim);

    /* Create a new dataset within the file using chunk creation properties. */
    hid_t dataset = H5Dcreate2(group, datasetname, H5T_STD_I32LE, dataspace,
                               H5P_DEFAULT, prop, H5P_DEFAULT);

    /* Write data to dataset */
    if(H5Dwrite(dataset, H5T_STD_I32LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data)) {
        return -1;
    }
    return 0;
}
