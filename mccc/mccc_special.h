/**
 * @file mccc_special.h
 * @brief Header file for mccc_special.c
*/
#ifndef MCCC_SPECIAL_H
#define MCCC_SPECIAL_H
#include "../ascot5.h"

void mccc_special_G(real x, real* G, int exact);

void mccc_special_GdG(real x, real* GdG, int exact);

void mccc_special_fo(real x, real* fdf, int exact);

void mccc_special_mu(real u, real th, real* mu, int exact);

void mccc_special_mudmu(real u, real th, real* mudmu, int exact);
#endif
