#include "ascot.h"
#include "data/bfield.h"
#include "defines.h"
#include "parallel.h"


size_t ascot_sizeof_real(void) {
    return sizeof(real);
}

void ascot_map_rhotheta_to_rz(
    Bfield *bfield, size_t npnt, size_t maxiter, real tol, real t, real rho[npnt],
    real theta[npnt], real phi[npnt], real rz[2][npnt])
{
    (void)t;
    OMP_PARALLEL_CPU_ONLY
    for (size_t j = 0; j < npnt; j++)
    {
        real axisrz[2], rho_drho[4];
        if (Bfield_eval_axis_rz(axisrz, bfield, phi[j]))
        {
            continue;
        }
        if (Bfield_eval_rho_drho(
                rho_drho, axisrz[0], phi[j], axisrz[1], bfield))
        {
            continue;
        }
        if (rho_drho[0] > rho[j])
        {
            rz[0][j] = axisrz[0];
            rz[1][j] = axisrz[1];
            continue;
        }

        real a = 0.0, b = 5.0;
        real costh = cos(theta[j]);
        real sinth = sin(theta[j]);
        for (size_t i = 0; i < maxiter; i++)
        {
            real c = 0.5 * (a + b);
            real rj = axisrz[0] + c * costh;
            real zj = axisrz[1] + c * sinth;
            if (rj < 0)
            {
                b = c;
                continue;
            }
            if (Bfield_eval_rho_drho(rho_drho, rj, phi[j], zj, bfield))
            {
                b = c;
                continue;
            }
            if (fabs(rho[j] - rho_drho[0]) < tol)
            {
                rz[0][j] = rj;
                rz[1][j] = zj;
                break;
            }
            if (rho[j] < rho_drho[0])
            {
                b = c;
            }
            else
            {
                a = c;
            }
        }
    }
}


void ascot_find_psi_on_axis_2d(
    Bfield *bfield, size_t maxiter, size_t ascent, real step, real tol,
    real psi[1], real rz[2])
{

    if (ascent)
    {
        step = -1 * step;
    }

    real phi = 0.0, time = 0.0;
    real psidpsi[4], nextrz[2];
    Bfield_eval_psi_dpsi(psidpsi, rz[0], phi, rz[1], time, bfield);

    size_t iter = 0;
    while (1)
    {
        if (Bfield_eval_psi_dpsi(psidpsi, rz[0], phi, rz[1], time, bfield))
        {
            break;
        }
        nextrz[0] = rz[0] - step * psidpsi[1];
        nextrz[1] = rz[1] - step * psidpsi[3];

        // Check convergence
        if (sqrt(
                (nextrz[0] - rz[0]) * (nextrz[0] - rz[0]) +
                (nextrz[1] - rz[1]) * (nextrz[1] - rz[1])) < tol)
        {
            psi[0] = psidpsi[0];
            rz[0] = nextrz[0];
            rz[1] = nextrz[1];

            // Add a bit of padding
            Bfield_eval_psi_dpsi(psidpsi, rz[0], phi, rz[1], time, bfield);
            psi[0] = psi[0] + (tol * psidpsi[1] + tol * psidpsi[3]);
            break;
        }

        rz[0] = nextrz[0];
        rz[1] = nextrz[1];
        iter++;

        if (iter == maxiter)
        {
            break;
        }
    }
}

void ascot_find_psi_on_axis_3d(
    Bfield *bfield, size_t maxiter, int ascent, real phimin, real phimax,
    real step, real tol, real psi[1], real rzphi[3])
{

    if (ascent)
    {
        step = -1 * step;
    }

    double time = 0.0;
    real psidpsi[4], nextrzphi[3];
    Bfield_eval_psi_dpsi(psidpsi, rzphi[0], rzphi[2], rzphi[1], time, bfield);

    size_t iter = 0;
    while (1)
    {
        if (Bfield_eval_psi_dpsi(
                psidpsi, rzphi[0], rzphi[2], rzphi[1], time, bfield))
        {
            break;
        }
        nextrzphi[0] = rzphi[0] - step * psidpsi[1]; // R
        nextrzphi[1] = rzphi[1] - step * psidpsi[3]; // z
        nextrzphi[2] =
            rzphi[2] - step / rzphi[0] * psidpsi[2]; // phi

        /* Check that phi remained inside the sector. If not, use the value on
           the sector boundary. */
        if (nextrzphi[2] > phimax)
        {
            nextrzphi[2] = phimax;
        }
        if (nextrzphi[2] < phimin)
        {
            nextrzphi[2] = phimin;
        }

        /* Check convergence (phi difference must be multiplied by R to get
        the arc length which has dimensions of L) */
        if (sqrt(
                (nextrzphi[0] - rzphi[0]) * (nextrzphi[0] - rzphi[0]) +
                (nextrzphi[1] - rzphi[1]) * (nextrzphi[1] - rzphi[1]) +
                rzphi[0] * (nextrzphi[2] - rzphi[2]) * rzphi[0] *
                    (nextrzphi[2] - rzphi[2])) < tol)
        {
            psi[0] = psidpsi[0];
            rzphi[0] = nextrzphi[0];
            rzphi[1] = nextrzphi[1];
            rzphi[2] = nextrzphi[2];

            // Add a bit of padding
            Bfield_eval_psi_dpsi(
                psidpsi, rzphi[0], rzphi[2], rzphi[1], time, bfield);
            psi[0] = psi[0] +
                     (tol * (psidpsi[1] + psidpsi[2] / rzphi[0] + psidpsi[3]));
            break;
        }

        rzphi[0] = nextrzphi[0];
        rzphi[1] = nextrzphi[1];
        rzphi[2] = nextrzphi[2];
        iter++;

        if (iter == maxiter)
        {
            break;
        }
    }
}
