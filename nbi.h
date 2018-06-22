#ifndef NBI_H
#define NBI_H

#include "ascot5.h"

typedef struct {
    int id;
    int n_beamlet;
    real* beamlet_x;
    real* beamlet_y;
    real* beamlet_z;
    real* beamlet_dirx;
    real* beamlet_diry;
    real* beamlet_dirz;
    real power;
    real energy;
    real efrac[3];
    real divergence[3];
    int anum;
    int znum;
} nbi_injector;

#endif
