#include <stdio.h>
#include <stdlib.h>
#include "ascot5.h"
#include "particle.h"
#include "diag.h"


void diag_init(diag_data* data, diag_offload_data* offload_data){
    data->diag_orb_collect = 1;
    data->orbits = diag_orb_init();
}

void diag_update_gc(diag_data* d, particle_simd_gc* p_f, particle_simd_gc* p_i){
    if(d->diag_orb_collect){
	diag_orb_updategc(p_f, p_i, d->orbits);
    }
    if(d->diag_debug_collect){
	
    }
    if(d->diag_dist4D_collect){
	
    }
}

void diag_update_fo(diag_data* d, particle_simd_fo* p_f, particle_simd_fo* p_i){
    if(d->diag_orb_collect){
	diag_orb_updatefo(p_f, p_i, d->orbits);
    }
    if(d->diag_debug_collect){
	
    }
    if(d->diag_dist4D_collect){
	
    }
}

void diag_write(diag_data* d){
    if(d->diag_orb_collect){
	diag_orb_write(d->orbits);
    }
    if(d->diag_debug_collect){
	
    }
    if(d->diag_dist4D_collect){
	
    }
}

void diag_clean(diag_data* d){
    if(d->diag_orb_collect){
	diag_orb_clean(d->orbits);
    }
    if(d->diag_debug_collect){
	
    }
    if(d->diag_dist4D_collect){
	
    }
}

