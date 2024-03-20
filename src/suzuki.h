/**
 * @file suzuki.h
 * @brief Header file for suzuki.c
 */
#ifndef SUZUKI_H
#define SUZUKI_H

#include "ascot5.h"
#include "error.h"

a5err suzuki_sigmav(real* sigmav, real EperAmu, real vnorm, real ne, real te,
                    integer nion, real* ni, const int* Anum, const int* Znum);

#endif
