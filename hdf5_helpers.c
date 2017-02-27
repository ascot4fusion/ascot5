/**
 * @file hdf5_helpers.c
 * @brief HDF5 wrappers and shortcuts
 */
#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <stdlib.h>
#include <hdf5.h>

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
