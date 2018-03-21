#include <stdlib.h>

/**
 * @brief Returns the interpolated y-value.
 * Saturates to y0 or y1 if x outside interval [x0, x1].
 */
int linint1D_eval(real* val, real x0, real y0, real x1, real y1, real x){
  real t;
  int err = 0;

  t =  (x-x0);
  t /= (x1-x0);

  *val = y0 + t*(y1-y0);
  return err;
}
