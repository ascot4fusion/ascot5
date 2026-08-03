/**
 * @file boschhale.h
 * @brief Header file for boschdale.c
 */
#ifndef BOSCHHALE_H
#define BOSCHHALE_H

#include "ascot5.h"
#include "plasma.h"
#include "linint/linint.h"

/**
 * @brief Available reactions
 */
typedef enum Reaction {
    DT_He4n   = 1,
    DHe3_He4p = 2,
    DD_Tp     = 3,
    DD_He3n   = 4,
} Reaction;

/**
 * @brief pre-calculated <sigma*v> for fusion
 */
typedef struct {
    linint2D_data* DT_He4n;    /**< D+T ==> He4+n                             */
    linint2D_data* DHe3_He4p;  /**< D+He3 ==> He4+p                           */
    linint2D_data* DD_Tp;      /**< D+D ==> T+p                               */
    linint2D_data* DD_He3n;    /**< D+D ==> He3+n                             */
} tabulated_sigmav_data;

void boschhale_reaction(
    Reaction reaction, real* m1, real* q1, real* m2, real* q2,
    real* mprod1, real* qprod1, real* mprod2, real* qprod2, real* Q);
real boschhale_sigma(Reaction reaction, real E);
real boschhale_sigmav(Reaction reaction, real Ti);
real boschhale_sigmav_beam_bulk(Reaction reaction, real vt, real vf, int N);
a5err boschhale_sigmav_beam_bulk_tabulated(real* sigmav, Reaction reaction, real vt, real vf, tabulated_sigmav_data* data);
int tabulated_sigmav_init(tabulated_sigmav_data* data, plasma_data* plasma_data);
void tabulated_sigmav_free(tabulated_sigmav_data* data);

#endif
