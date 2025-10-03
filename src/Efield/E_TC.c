#include <stdlib.h>
#include "ascot5.h"
#include "math.h"
#include "error.h"
#include "B_field.h"
#include "E_TC.h"


int EfieldCartesian_init(EfieldCartesian* efield, real exyz[3]) {
    efield->exyz[0] = exyz[0];
    efield->exyz[1] = exyz[1];
    efield->exyz[2] = exyz[2];
    return 0;
}



void EfieldCartesian_free(EfieldCartesian* efield) {
    // No resources were allocated
    (void)efield;
}


void EfieldCartesian_offload(EfieldCartesian* efield) {
    (void)efield;
    GPU_MAP_TO_DEVICE( data->exyz[0:3] )
}


a5err EfieldCartesian_eval_e(
    real e[3], real r, real phi, real z, EfieldCartesian* efield,
    B_field_data* bfield
) {
    math_vec_xyz2rpz(efield->exyz, e, phi);

    return 0;
}
