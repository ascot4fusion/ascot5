/**
 * @file hdf5_histogram.h
 * @brief Header file for hdf5_histogram.c
 */
#ifndef HDF5_HISTOGRAM
#define HDF5_HISTOGRAM

#include <hdf5.h>

int hdf5_histogram_write_uniform_double(hid_t f, const char *path,
                                        int abscissaDim, int ordinateDim,
                                        int *abscissaNslots,
                                        double *abscissaMin,
                                        double *abscissaMax,
                                        char **abscissaUnits,
                                        char **abscissaNames,
                                        char **ordinateUnits,
                                        char **ordinateNames,
                                        double *ordinate);
#endif
