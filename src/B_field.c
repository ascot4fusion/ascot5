/**
 * @file B_field.c
 * @brief Magnetic field interface
 *
 * This is an interface through which magnetic field data is initialized and
 * accessed. Reading e.g. from disk is done elsewhere.
 *
 * To add a new magnetic field instance, make sure these functions are
 * implemented and called from this interface, and that B_field.h contains enum
 * type for the new instance.
 *
 * The interface checks which instance given data corresponds to from
 * B_field_offload_data.type and B_field_data.type from the structure that is
 * given as an argument, and calls the relevant function for that instance.
 */
#include <stdio.h>
#include "ascot5.h"
#include "error.h"
#include "print.h"
#include "B_field.h"
#include "Bfield/B_GS.h"
#include "Bfield/B_2DS.h"
#include "Bfield/B_3DS.h"
#include "Bfield/B_STS.h"
#include "Bfield/B_TC.h"

/**
 * @brief Free allocated resources
 *
 * @param data pointer to data struct
 */
void B_field_free(B_field_data* data) {
    switch(data->type) {
        case B_field_type_GS:
            B_GS_free(&(data->BGS));
            break;

        case B_field_type_2DS:
            B_2DS_free(&(data->B2DS));
            break;

        case B_field_type_3DS:
            B_3DS_free(&(data->B3DS));
            break;

        case B_field_type_STS:
            B_STS_free(&(data->BSTS));
            break;

        case B_field_type_TC:
            B_TC_free(&(data->BTC));
            break;
    }
}

/**
 * @brief Evaluate poloidal flux psi
 *
 * This function evaluates the poloidal flux psi at the given coordinates. The
 * psi is exactly as it appears in Grad-Shafranov equation.
 *
 * This is a SIMD function.
 *
 * @param psi pointer where psi [V*s*m^-1] value will be stored
 * @param r R coordinate [m]
 * @param phi phi coordinate [rad]
 * @param z z coordinate [m]
 * @param t time coordinate [s]
 * @param Bdata pointer to magnetic field data struct
 *
 * @return Non-zero a5err value if evaluation failed, zero otherwise
 */
a5err B_field_eval_psi(real* psi, real r, real phi, real z, real t,
                      B_field_data* Bdata) {
    a5err err = 0;

    switch(Bdata->type) {
        case B_field_type_GS:
            err = B_GS_eval_psi(psi, r, phi, z, &(Bdata->BGS));
            break;

        case B_field_type_2DS:
            err = B_2DS_eval_psi(psi, r, phi, z, &(Bdata->B2DS));
            break;

        case B_field_type_3DS:
            err = B_3DS_eval_psi(psi, r, phi, z, &(Bdata->B3DS));
            break;

        case B_field_type_STS:
            err = B_STS_eval_psi(psi, r, phi, z, &(Bdata->BSTS));
            break;

        case B_field_type_TC:
            err = B_TC_eval_psi(psi, r, phi, z, &(Bdata->BTC));
            break;

        default:
            /* Unregonized input. Produce error. */
            err = error_raise( ERR_UNKNOWN_INPUT, __LINE__, EF_B_FIELD );
            break;
    }

    if(err) {
        /* In case of error, return some reasonable value to avoid further
           complications */
        psi[0] = 1;
    }

    return err;
}

/**
 * @brief Evaluate poloidal flux psi and its derivatives
 *
 * This function evaluates the poloidal flux psi and its derivatives at the
 * given coordinates. The psi is exactly as it appears in Grad-Shafranov
 * equation.
 *
 * The values are stored in the given array as:
 * - psi_dpsi[0] = psi
 * - psi_dpsi[1] = dpsi/dr
 * - psi_dpsi[2] = dpsi/dphi
 * - psi_dpsi[3] = dpsi/dz
 *
 * This is a SIMD function.
 *
 * @param psi_dpsi pointer for storing psi [V*s*m^-1] and its derivatives
 * @param r R coordinate [m]
 * @param phi phi coordinate [rad]
 * @param z z coordinate [m]
 * @param t time coordinate [s]
 * @param Bdata pointer to magnetic field data struct
 *
 * @return Non-zero a5err value if evaluation failed, zero otherwise
 */
a5err B_field_eval_psi_dpsi(real psi_dpsi[4], real r, real phi, real z, real t,
                            B_field_data* Bdata) {
    a5err err = 0;

    switch(Bdata->type) {
        case B_field_type_GS:
            err = B_GS_eval_psi_dpsi(psi_dpsi, r, phi, z, &(Bdata->BGS));
            break;

        case B_field_type_2DS:
            err = B_2DS_eval_psi_dpsi(psi_dpsi, r, phi, z, &(Bdata->B2DS));
            break;

        case B_field_type_3DS:
            err = B_3DS_eval_psi_dpsi(psi_dpsi, r, phi, z, &(Bdata->B3DS));
            break;

        case B_field_type_STS:
            err = B_STS_eval_psi_dpsi(psi_dpsi, r, phi, z, &(Bdata->BSTS));
            break;

        case B_field_type_TC:
            err = B_TC_eval_psi_dpsi(psi_dpsi, r, phi, z, &(Bdata->BTC));
            break;

        default:
            /* Unregonized input. Produce error. */
            err = error_raise( ERR_UNKNOWN_INPUT, __LINE__, EF_B_FIELD );
            break;
    }

    if(err) {
        /* In case of error, return some reasonable values to avoid further
           complications */
        psi_dpsi[0] = 1;
        for(int k=1; k<4; k++) {psi_dpsi[k] = 0;}
    }

    return err;
}

/**
 * @brief Evaluate normalized poloidal flux rho and its psi derivative
 *
 * This function evaluates the normalized poloidal flux rho at the given
 * coordinates. The rho is evaluated from psi as:
 *
 * \f{equation*}{
 * \rho = \sqrt{ \frac{\psi - \psi_0}{\psi_1 - \psi_0} },
 * \f}
 *
 * where \f$\psi_0\f$ is psi at magnetic axis and \f$\psi_1\f$ is psi at
 * separatrix.
 *
 * This is a SIMD function.
 *
 * @param rho pointer where rho value will be stored
 * @param psi poloidal flux from which rho is evaluated
 * @param Bdata pointer to magnetic field data struct
 *
 * @return Non-zero a5err value if evaluation failed, zero otherwise
 */
a5err B_field_eval_rho(real rho[2], real psi, B_field_data* Bdata) {
    a5err err = 0;

    real psi0 = 0.0, psi1 = 1.0;
    switch(Bdata->type) {
        case B_field_type_GS:
            psi0 = Bdata->BGS.psi0;
            psi1 = Bdata->BGS.psi1;
            break;

        case B_field_type_2DS:
            psi0 = Bdata->B2DS.psi0;
            psi1 = Bdata->B2DS.psi1;
            break;

        case B_field_type_3DS:
            psi0 = Bdata->B3DS.psi0;
            psi1 = Bdata->B3DS.psi1;
            break;

        case B_field_type_STS:
            psi0 = Bdata->BSTS.psi0;
            psi1 = Bdata->BSTS.psi1;
            break;

        case B_field_type_TC:
            psi0 = Bdata->BTC.psival - Bdata->BTC.rhoval * Bdata->BTC.rhoval;
            psi1 = Bdata->BTC.psival - Bdata->BTC.rhoval * Bdata->BTC.rhoval
                   + 1.0;
            break;

        default:
            /* Unregonized input. Produce error. */
            err = error_raise( ERR_UNKNOWN_INPUT, __LINE__, EF_B_FIELD );
            break;
    }

    real delta = (psi1 - psi0);
    if( (psi - psi0) / delta < 0 ) {
         err = error_raise( ERR_INPUT_UNPHYSICAL, __LINE__, EF_B_FIELD );
    } else {
        rho[0] = sqrt( (psi - psi0) / delta );
        rho[1] = 1.0 / (2*delta*rho[0]);
    }

    if(err) {
        /* In case of error, return some reasonable value to avoid further
           complications */
        rho[0] = 1;
        rho[1] = 0;
    }

    return err;
}

/**
 * @brief Evaluate normalized poloidal flux rho and its derivatives
 *
 * This function evaluates the normalized poloidal flux rho at the given
 * coordinates. The rho is evaluated from psi as:
 *
 * \f{equation*}{
 * \rho = \sqrt{ \frac{\psi - \psi_0}{\psi_1 - \psi_0} },
 * \f}
 *
 * where \f$\psi_0\f$ is psi at magnetic axis and \f$\psi_1\f$ is psi at
 * separatrix.
 *
 * The values are stored in the given array as:
 * - rho_drho[0] = rho
 * - rho_drho[1] = drho/dr
 * - rho_drho[2] = drho/dphi
 * - rho_drho[3] = drho/dz
 *
 * This is a SIMD function.
 *
 * @param rho_drho pointer where rho and its derivatives will be stored
 * @param r R coordinate [m]
 * @param phi phi coordinate [rad]
 * @param z z coordinate [m]
 * @param Bdata pointer to magnetic field data struct
 *
 * @return Non-zero a5err value if evaluation failed, zero otherwise
 */
a5err B_field_eval_rho_drho(real rho_drho[4], real r, real phi, real z,
                            B_field_data* Bdata) {
    a5err err = 0;

    switch(Bdata->type) {
        case B_field_type_GS:
            err = B_GS_eval_rho_drho(rho_drho, r, phi, z, &(Bdata->BGS));
            break;

        case B_field_type_2DS:
            err = B_2DS_eval_rho_drho(rho_drho, r, phi, z, &(Bdata->B2DS));
            break;

        case B_field_type_3DS:
            err = B_3DS_eval_rho_drho(rho_drho, r, phi, z, &(Bdata->B3DS));
            break;

        case B_field_type_STS:
            err = B_STS_eval_rho_drho(rho_drho, r, phi, z, &(Bdata->BSTS));
            break;

        case B_field_type_TC:
            err = B_TC_eval_rho_drho(rho_drho, r, phi, z, &(Bdata->BTC));
            break;

        default:
            /* Unregonized input. Produce error. */
            err = error_raise( ERR_UNKNOWN_INPUT, __LINE__, EF_B_FIELD );
            break;
    }

    if(err) {
        /* In case of error, return some reasonable values to avoid further
           complications */
        rho_drho[0] = 1;
        for(int k=1; k<4; k++) {rho_drho[k] = 0;}
    }

    return err;
}

/**
 * @brief Evaluate magnetic field
 *
 * This function evaluates the magnetic field at the given coordinates.
 *
 * The values are stored in the given array as:
 * - B[0] = BR
 * - B[1] = Bphi
 * - B[2] = Bz
 *
 * This is a SIMD function.
 *
 * @param B pointer to array where magnetic field values are stored
 * @param r R coordinate [m]
 * @param phi phi coordinate [deg]
 * @param z z coordinate [m]
 * @param t time coordinate [s]
 * @param Bdata pointer to magnetic field data struct
 *
 * @return Non-zero a5err value if evaluation failed, zero otherwise
 */
a5err B_field_eval_B(real B[3], real r, real phi, real z, real t,
                     B_field_data* Bdata) {
    a5err err = 0;

    switch(Bdata->type) {
        case B_field_type_GS:
            err = B_GS_eval_B(B, r, phi, z, &(Bdata->BGS));
            break;

        case B_field_type_2DS:
            err = B_2DS_eval_B(B, r, phi, z, &(Bdata->B2DS));
            break;

        case B_field_type_3DS:
            err = B_3DS_eval_B(B, r, phi, z, &(Bdata->B3DS));
            break;

        case B_field_type_STS:
            err = B_STS_eval_B(B, r, phi, z, &(Bdata->BSTS));
            break;

        case B_field_type_TC:
            err = B_TC_eval_B(B, r, phi, z, &(Bdata->BTC));
            break;

        default:
            /* Unregonized input. Produce error. */
            err = error_raise( ERR_UNKNOWN_INPUT, __LINE__, EF_B_FIELD );
            break;
    }

    if(err) {
        /* In case of error, return some reasonable values to avoid further
           complications */
        B[0] = 1;
        for(int k=1; k<3; k++) {B[k] = 0;}
    }

    return err;
}

/**
 * @brief Evaluate magnetic field and its derivatives
 *
 * This function evaluates the magnetic field and its derivatives at the given
 * coordinates.
 *
 * The values are stored in the given array as:
 * - B[0]  = BR
 * - B[1]  = dBR/dR
 * - B[2]  = dBR/dphi
 * - B[3]  = dBR/dz
 * - B[4]  = Bphi
 * - B[5]  = dBphi/dR
 * - B[6]  = dBphi/dphi
 * - B[7]  = dBphi/dz
 * - B[8]  = Bz
 * - B[9]  = dBz/dR
 * - B[10] = dBz/dphi
 * - B[11] = dBz/dz
 * - B[12] = dBR/dt
 * - B[13] = dBphi/dt
 * - B[14] = dBz/dt
 *
 * This is a SIMD function.
 *
 * @param B_dB pointer to array where the field and its derivatives are stored
 * @param r R coordinate [m]
 * @param phi phi coordinate [deg]
 * @param z z coordinate [m]
 * @param t time coordinate [s]
 * @param Bdata pointer to magnetic field data struct
 *
 * @return Non-zero a5err value if evaluation failed, zero otherwise
 */
a5err B_field_eval_B_dB(real B_dB[15], real r, real phi, real z, real t,
                        B_field_data* Bdata) {
    a5err err = 0;

    switch(Bdata->type) {
        case B_field_type_GS:
            err = B_GS_eval_B_dB(B_dB, r, phi, z, &(Bdata->BGS));
            break;

        case B_field_type_2DS:
            err = B_2DS_eval_B_dB(B_dB, r, phi, z, &(Bdata->B2DS));
            break;

        case B_field_type_3DS:
            err = B_3DS_eval_B_dB(B_dB, r, phi, z, &(Bdata->B3DS));
            break;

        case B_field_type_STS:
            err = B_STS_eval_B_dB(B_dB, r, phi, z, &(Bdata->BSTS));
            break;

        case B_field_type_TC:
            err = B_TC_eval_B_dB(B_dB, r, phi, z, &(Bdata->BTC));
            break;

        default:
            /* Unregonized input. Produce error. */
            err = error_raise( ERR_UNKNOWN_INPUT, __LINE__, EF_B_FIELD );
            break;
    }

    if(err) {
        /* In case of error, return some reasonable values to avoid further
           complications */
        B_dB[0] = 1;
        for(int k=1; k<12; k++) {B_dB[k] = 0;}
    }

    return err;
}

/**
 * @brief Return magnetic axis Rz-coordinates
 *
 * Returns magnetic axis Rz-coordinates at given toroidal angle.
 *
 * @param rz pointer where axis R and z [m] values will be stored
 * @param Bdata pointer to magnetic field data struct
 * @param phi phi coordinate [deg]
 *
 * @return Magnetic axis R-coordinate [m]
 */
a5err B_field_get_axis_rz(real rz[2], B_field_data* Bdata, real phi) {
    a5err err = 0;

    switch(Bdata->type) {
        case B_field_type_GS:
            err = B_GS_get_axis_rz(rz, &(Bdata->BGS));
            break;

        case B_field_type_2DS:
            err = B_2DS_get_axis_rz(rz, &(Bdata->B2DS));
            break;

        case B_field_type_3DS:
            err = B_3DS_get_axis_rz(rz, &(Bdata->B3DS));
            break;

        case B_field_type_STS:
            err = B_STS_get_axis_rz(rz, &(Bdata->BSTS), phi);
            break;

        case B_field_type_TC:
            err = B_TC_get_axis_rz(rz, &(Bdata->BTC));
            break;

        default:
            /* Unregonized input. Produce error. */
            err = error_raise( ERR_UNKNOWN_INPUT, __LINE__, EF_B_FIELD );
            break;
    }

    if(err) {
      /* In case of error, return some reasonable values to avoid further
         complications */
      rz[0] = 1.0;
      rz[1] = 0.0;
    }

    return err;
}
