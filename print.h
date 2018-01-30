/**
 * @file print.h
 * @brief Macros for printing console output
 */
#ifndef PRINT_H
#define PRINT_H

#include <stdio.h>

/**
 * @brief Print to standard output
 */
#define print_out(...) printf(__VA_ARGS__)

/**
 * @brief Print to standard output only for rank 0
 */
#define print_out0(rank,...) { if((rank)==0) printf(__VA_ARGS__); }

/**
 * @brief Print to standard error
 */
#define print_err(...) fprintf(stderr,__VA_ARGS__)

/**
 * @brief Print to standard error only for rank 0
 */
#define print_err0(...) { if((rank)==0) fprintf(stderr,__VA_ARGS__); }

#endif
