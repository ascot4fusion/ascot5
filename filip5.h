/**
 * @file filip5.h
 * @brief Header file for filip5.h
 */
#ifndef FILIP5_H
#define FILIP5_H

#include "ascot5.h"

typedef struct {
    real r;
    real phi;
    real z;
    real ini_r;
    real ini_phi;
    real ini_z;
    real dir;
    integer hit_wall;
    real length;
} fieldline;

#endif
