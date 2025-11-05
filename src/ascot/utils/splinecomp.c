/**
 * Cubic spline interpolation coefficients of a 1D data set (see interp.h).
 */
#include "defines.h"
#include "utils/interp.h"
#include <stdlib.h>

int solve_compact_cubic_spline(size_t n, real c[2 * n], real f[n], int bc)
{
    int err = 0;
    real *rhs = (real *)xmalloc(&err, n * sizeof(real));
    if (err)
        return 1;
    real *second_derivative = (real *)xmalloc(&err, n * sizeof(real));
    if (err)
    {
        free(rhs);
        return 1;
    }
    real *p = (real *)xmalloc(&err, (n-1) * sizeof(real));
    if (err)
    {
        free(rhs);
        free(second_derivative);
        return 1;
    }

    if (bc == NATURALBC)
    {
        /* Initialize RHS of equation */
        for (size_t i = 1; i < n - 1; i++)
            rhs[i] = 6 * (f[i + 1] - 2 * f[i] + f[i - 1]);

        /* Enforce boundary condition */
        rhs[0] = 0.0;
        rhs[n - 1] = 0.0;

        /* Forward sweep */
        p[0] = 0.0;
        for (size_t i = 1; i < n - 1; i++)
        {
            p[i] = 1 / (4 - p[i - 1]);
            rhs[i] = (rhs[i] - rhs[i - 1]) / (4 - p[i - 1]);
        }
        rhs[n - 1] = 0.0;

        /* Back substitution */
        second_derivative[n - 1] = rhs[n - 1];
        for (size_t i = n - 1; i > 0; i--)
        {
            second_derivative[i - 1] =
                rhs[i - 1] - p[i - 1] * second_derivative[i];
        }
    }
    else if (bc == PERIODICBC)
    {
        /** PERIODIC (Function has same value and derivatives at both ends) **/

        /* Initialize some additional helper variables */

        /* Value that starts from lower left corner and moves right */
        real l = 1.0;
        /* Last diagonal value        */
        real dlast = 4.0;
        /* Matrix right column values */
        real *r = malloc((n - 2) * sizeof(real));
        /* Last subdiagonal value     */
        real blast;

        /* Initialize RHS of equation */
        rhs[0] = 6 * (f[1] - 2 * f[0] + f[n - 1]);
        for (size_t i = 1; i < n - 1; i++)
        {
            rhs[i] = 6 * (f[i + 1] - 2 * f[i] + f[i - 1]);
        }
        rhs[n - 1] = 6 * (f[0] - 2 * f[n - 1] + f[n - 2]);

        /* Simplified Gauss elimination is used (own algorithm) */

        /* Forward sweep */
        p[0] = 1.0 / 4;
        r[0] = 1.0 / 4;
        rhs[0] = rhs[0] / 4;
        for (size_t i = 1; i < n - 2; i++)
        {
            dlast = dlast - l * r[i - 1];
            rhs[n - 1] = rhs[n - 1] - l * rhs[i - 1];
            l = -l * p[i - 1];
            p[i] = 1 / (4 - p[i - 1]);
            r[i] = -r[i - 1] / (4 - p[i - 1]);
            rhs[i] = (rhs[i] - rhs[i - 1]) / (4 - p[i - 1]);
        }
        blast = 1.0 - l * p[n - 3];
        dlast = dlast - l * r[n - 3];
        rhs[n - 1] = rhs[n - 1] - l * rhs[n - 3];

        p[n - 2] = (1 - r[n - 3]) / (4 - p[n - 3]);
        rhs[n - 2] = (rhs[n - 2] - rhs[n - 3]) / (4 - p[n - 3]);
        rhs[n - 1] =
            (rhs[n - 1] - blast * rhs[n - 2]) / (dlast - blast * p[n - 2]);

        /* Back substitution */
        second_derivative[n - 1] = rhs[n - 1];
        second_derivative[n - 2] =
            rhs[n - 2] - p[n - 2] * second_derivative[n - 1];
        for (size_t i = n - 2; i > 0; i--)
            second_derivative[i - 1] = rhs[i - 1] -
                                       p[i - 1] * second_derivative[i] -
                                       r[i - 1] * second_derivative[n - 1];

        free(r);
    }

    for (size_t i = 0; i < n; i++)
    {
        c[2 * i] = f[i];
        c[2 * i + 1] = second_derivative[i];
    }

    free(rhs);
    free(second_derivative);
    free(p);
    return 0;
}
