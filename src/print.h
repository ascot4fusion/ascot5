/**
 * @file print.h
 * @brief Macros for printing console output
 */
#ifndef PRINT_H
#define PRINT_H

#include <stdio.h>

enum {
    VERBOSE_DEBUG   = 2,
    VERBOSE_NORMAL  = 1,
    VERBOSE_MINIMAL = 0,
    VERBOSE_IO      = 1
};

extern const char VERBOSE_LEVEL;

/**
 * @brief Print to standard output
 */
#define print_out(v,...) { if(VERBOSE_LEVEL >= (v)) printf(__VA_ARGS__); }

/**
 * @brief Print to standard output only for rank 0
 */
#define print_out0(v,rank,...) { \
    if(VERBOSE_LEVEL >= (v) &&(rank)==0) printf(__VA_ARGS__); }

/**
 * @brief Print to standard error
 */
#define print_err(...) fprintf(stderr,__VA_ARGS__)

#endif
