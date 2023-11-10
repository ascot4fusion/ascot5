/**
 * @file boschhale.h
 * @brief Header file for boschdale.c
 */
#ifndef BOSCHHALE_H
#define BOSCHHALE_H

#include "ascot5.h"

/**
 * @brief Available reactions
 */
typedef enum Reaction {
    DT_He4n   = 1,
    DHe3_He4p = 2,
    DD_Tp     = 3,
    DD_He3n   = 4,
} Reaction;

void boschhale_reaction(
    Reaction reaction, real* m1, real* q1, real* m2, real* q2,
    real* mprod1, real* qprod1, real* mprod2, real* qprod2, real* Q);
real boschhale_sigma(Reaction reaction, real E);
real boschhale_sigmav(Reaction reaction, real Ti);

#endif
