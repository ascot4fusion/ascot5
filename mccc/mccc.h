#ifndef MCCC_H
#define MCCC_H

void mccc_init();

void mccc_update_fo();

void mccc_update_gc();

void mccc_step_fo_fixed();

void mccc_step_gc_fixed();

void mccc_step_gc_adaptive();

void mccc_printerror(int err);

#endif
