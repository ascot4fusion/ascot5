#include "defines.h"
#include "bfield.h"

/**
 * @brief Map (rho, theta, phi) to (R,z) coordinates.
 *
 * This function implements the Newton method. If the function fails on
 * a given position, the corresponding (R,z) values in the output arrays are
 * not altered.
 *
 * @param bfield magnetic field data
 * @param Neval number of query points.
 * @param rho the square root of the normalized poloidal flux values.
 * @param theta poloidal angles [rad].
 * @param phi toroidal angles [rad].
 * @param t time coordinate (same for all) [s].
 * @param maxiter maximum number of iterations in Newton algorithm.
 * @param tol algorithm is stopped when |rho - rho(r,z)| < tol
 * @param r output array for R coordinates [m].
 * @param z output array for z coordinates [m].
 */
void ascot_map_rhotheta_to_rz(
    Bfield *bfield, int Neval, int maxiter, real tol, real t, real *rho,
    real *theta, real *phi, real *r, real *z)
{
    (void)t;
#pragma omp parallel for
    for (int j = 0; j < Neval; j++)
    {
        real axisrz[2];
        real rhodrho[4];
        if (B_field_get_axis_rz(axisrz, bfield, phi[j]))
        {
            continue;
        }
        if (B_field_eval_rho_drho(
                rhodrho, axisrz[0], phi[j], axisrz[1], bfield))
        {
            continue;
        }
        if (rhodrho[0] > rho[j])
        {
            /* Due to padding, rho might not be exactly zero on the axis so we
             * return the axis position for small values of queried rho */
            r[j] = axisrz[0];
            z[j] = axisrz[1];
            continue;
        }

        real a = 0.0, b = 5.0;
        real costh = cos(theta[j]);
        real sinth = sin(theta[j]);
        for (int i = 0; i < maxiter; i++)
        {
            real c = 0.5 * (a + b);
            real rj = axisrz[0] + c * costh;
            real zj = axisrz[1] + c * sinth;
            if (rj < 0)
            {
                b = c;
                continue;
            }
            if (B_field_eval_rho_drho(rhodrho, rj, phi[j], zj, bfield))
            {
                b = c;
                continue;
            }
            if (fabs(rho[j] - rhodrho[0]) < tol)
            {
                r[j] = rj;
                z[j] = zj;
                break;
            }
            if (rho[j] < rhodrho[0])
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

/**
 * @brief Find psi on axis using the gradient descent method
 *
 * Note that the psi value is not returned in case this algorithm fails.
 *
 * @param bfield magnetic field data
 * @param psi value of psi on axis if this function did not fail
 * @param rz initial (R,z) position where also the result is stored
 * @param step the step size
 * @param tol the current position is accepted if the distance (in meters)
 * between this and the previous point is below this value
 * @param maxiter maximum number of iterations before failure
 * @param ascent if true the algorithm instead ascends to find psi0 (> psi1)
 */
void ascot_find_psi_on_axis_2d(
    Bfield *bfield, int maxiter, int ascent, real step, real tol,
    real psi[1], real rz[2])
{

    if (ascent)
    {
        step = -1 * step;
    }

    real phi = 0.0, time = 0.0;
    real psidpsi[4], nextrz[2];
    B_field_eval_psi_dpsi(psidpsi, rz[0], phi, rz[1], time, bfield);

    int iter = 0;
    while (1)
    {
        if (B_field_eval_psi_dpsi(psidpsi, rz[0], phi, rz[1], time, bfield))
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
            B_field_eval_psi_dpsi(psidpsi, rz[0], phi, rz[1], time, bfield);
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
    Bfield *bfield, int maxiter, int ascent, real phimin, real phimax,
    real step, real tol, real psi[1], real rzphi[3])
{

    if (ascent)
    {
        step = -1 * step;
    }

    real time = 0.0;
    real psidpsi[4], nextrzphi[3];
    B_field_eval_psi_dpsi(psidpsi, rzphi[0], rzphi[2], rzphi[1], time, bfield);

    int iter = 0;
    while (1)
    {
        if (B_field_eval_psi_dpsi(
                psidpsi, rzphi[0], rzphi[2], rzphi[1], time, bfield))
        {
            break;
        }
        nextrzphi[0] = rzphi[0] - step * psidpsi[1]; // R
        nextrzphi[1] = rzphi[1] - step * psidpsi[3]; // z
        nextrzphi[2] =
            rzphi[2] - step / rzphi[0] * psidpsi[2]; /* phi. phidpsi[2]
                                                      is dimensionless,
                                                      must divide by R
                                                      because in
                                                      cylindrical
                                                      co-ordinates   */

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
            B_field_eval_psi_dpsi(
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
