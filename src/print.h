/**
 * @file print.h
 * @brief Macros for printing console output
 */
#ifndef PRINT_H
#define PRINT_H

#include <stdio.h>

/**
 * @brief Versbosity levels and definitions
 *
 * DEBUG, NORMAL, and MINIMAL are levels that can be specified and the rest are
 * used to specify to which levels operations belongs to.
 */
enum VERBOSE_RANK {
    VERBOSE_DEBUG   = 2, /**< For debugging purposes                   */
    VERBOSE_NORMAL  = 1, /**< For normal output                        */
    VERBOSE_MINIMAL = 0, /**< Only print version and simulation status */
    VERBOSE_IO      = 1  /**< Flag for IO operations                   */
};

/**
 * @brief Verbose level
 */
extern const char VERBOSE_LEVEL;

/**
 * @brief Print to standard output
 */
#define print_out(v,...) { if(VERBOSE_LEVEL >= (v)) printf(__VA_ARGS__); }

/**
 * @brief Print to standard output only for root process
 */
#define print_out0(v,rank,root,...) { \
    if(VERBOSE_LEVEL >= (v) &&(rank)==(root)) printf(__VA_ARGS__); }

/**
 * @brief Print to standard error
 */
#define print_err(...) fprintf(stderr,__VA_ARGS__)

#endif
