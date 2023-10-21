/**
 * @file boschhale.h
 * @brief Header file for boschdale.c
 */
#ifndef BOSCHHALE_H
#define BOSCHHALE_H

#include "ascot5.h"

void boschhale_reaction(int reaction, real* m1, real* q1, real* m2, real* q2,
                        real* mprod1, real* qprod1, real* mprod2, real* qprod2,
                        real* Q);
real boschhale_sigma(int reaction, real E);
real boschhale_sigmav(int reaction, real Ti);

#endif
