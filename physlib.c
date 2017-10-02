/**
 * @brief Physics library
 *
 * Contains elementary quantities and coordinate transformations.
 */
#include <math.h>
#include "ascot5.h"
#include "physlib.h"
#include "consts.h"
#include "math.h"

/**
 * @brief Lorentz factor
 *
 * gamma = sqrt(1 / (1 - v^2/c^2) )
 *
 * @param vnorm particle velocity norm
 */
real physlib_relfactorv_fo(real vnorm) {
    return sqrt( 1.0 / (1 - (vnorm * vnorm) / CONST_C2) );
}

/**
 * @brief Lorentz factor
 *
 * gamma = sqrt(1 + (p/mc)^2 )
 *
 * @param mass  particle mass
 * @param pnorm particle momentum norm
 */
real physlib_relfactorp_fo(real mass, real pnorm) {
    return sqrt( 1 + pow( pnorm/mass ,2) / CONST_C2 );
}

/**
 * @brief Lorentz factor
 *
 * gamma = sqrt( (1 + 2*mu*B/mc^2) / (1 - vpar^2/c^2) )
 *
 * @param mass  guiding center mass
 * @param mu    guiding center magnetic moment
 * @param vpar  guiding center parallel velocity
 * @param Bnorm magnetic field norm
 */
real physlib_relfactorv_gc(real mass, real mu, real vpar, real Bnorm) {
    return sqrt( ( 1 + 2*mu*Bnorm/(mass*CONST_C2) ) / (1 - vpar*vpar/CONST_C2) );
}

/**
 * @brief Lorentz factor
 *
 * gamma = sqrt(1 + 2*mu*B/mc^2 + (ppar/mc)^2 )
 *
 * @param mass  guiding center mass
 * @param mu    guiding center magnetic moment
 * @param ppar  guiding center parallel momentum
 * @param Bnorm magnetic field norm
 */
real physlib_relfactorp_gc(real mass, real mu, real ppar, real Bnorm) {
    return sqrt( 1 + 2*mu*Bnorm/(mass*CONST_C2) + pow(ppar/mass,2)/CONST_C2 );
}

/**
 * @brief Kinetic energy
 *
 * Ekin = (gamma - 1) * mc^2
 *
 * @param mass  particle mass
 * @param vnorm particle velocity norm
 */
real physlib_Ekin_fo(real mass, real vnorm) {
    return (physlib_relfactorv_fo(vnorm) - 1) * mass*CONST_C2;
}

/**
 * @brief Kinetic energy
 *
 * Ekin = (gamma - 1) * mc^2
 *
 * @param mass  guiding center mass
 * @param mu    guiding center magnetic moment
 * @param vpar  guiding center parallel velocity
 * @param Bnorm magnetic field norm
 */
real physlib_Ekin_gc(real mass, real mu, real vpar, real Bnorm) {
    return (physlib_relfactorv_gc(mass, mu, vpar, Bnorm) - 1) * mass*CONST_C2;
}

/**
 * @brief Transform guiding center (v, xi) coordinates to (mu, vpar)
 *
 * @param mass  guiding center mass
 * @param Bnorm magnetic field norm
 * @param v     guiding center total velocity
 * @param xi    guiding center pitch
 * @param mu    pointer to returned magnetic moment
 * @param vpar  pointer to returned parallel velocity
 */
void physlib_gc_vxi2muvpar(real mass, real Bnorm, real v, real xi, real* mu, real* vpar) {
    vpar[0] = xi*v;
    real v2 = v * v;
    real gamma2 = ( 1.0 / (1 - v2 / CONST_C2) );
    mu[0] = mass * gamma2 * v2 * (1 - xi*xi) / (2*Bnorm);
}

/**
 * @brief Transform guiding center (v, xi) coordinates to (mu, vpar)
 *
 * @param mass  guiding center mass
 * @param Bnorm magnetic field norm
 * @param v     guiding center total velocity
 * @param xi    guiding center pitch
 * @param mu    pointer to returned magnetic moment
 * @param vpar  pointer to returned parallel velocity
 */
void physlib_gc_muvpar2vxi(real mass, real Bnorm, real mu, real vpar, real* v, real* xi) {
    real gamma2 = ( 1 + 2*mu*Bnorm/(mass*CONST_C2) ) / (1 - vpar*vpar/CONST_C2);
    real vperp2 = 2*mu*Bnorm/(gamma2*mass);
    
    v[0]  = sqrt(vperp2 + vpar*vpar);
    xi[0] = vpar/v[0];
}

/**
 * @brief First order guiding center transformation from fo to gc
 *
 * @param mass   particle mass
 * @param charge particle charge
 * @param B_dB   magnetic field and Jacobian at particle location
 * @param Rprt, phiprt, zprt particle location
 * @param pR, pphi, pz particle momentum components
 * @param R, phi, z pointers to returned guiding center location
 * @param mu, ppar, theta pointers to returned guiding center magnetic moment,
 *        parallel momentum and gyroangle
 */
void physlib_fo2gc(real mass, real charge, real* B_dB,
		   real Rprt, real phiprt, real zprt, real pR, real pphi, real pz,
		   real* R, real* phi, real* z, real* mu, real* ppar, real* theta) {
    
    /* Temporary re-used variables */
    real TEMP_S1;
    real TEMP_V1[3];

    /* Helper variables */
    real prtrpz[3];
    prtrpz[0] = Rprt;
    prtrpz[1] = phiprt;
    prtrpz[2] = zprt;

    real prpz[3];
    prpz[0] = pR;
    prpz[1] = pphi;
    prpz[2] = pz;
    real ptot  = math_norm(prpz);
    

    /* Transform into the cartesian system */
    real prtxyz[3];
    math_rpz2xyz(prtrpz, prtxyz);

    real p_unit[3];
    math_vec_rpz2xyz(prpz, p_unit, phiprt);
    p_unit[0] /= ptot;
    p_unit[1] /= ptot;
    p_unit[2] /= ptot;

    real B_dBc[12];
    math_jac_rpz2xyz(B_dB, B_dBc, Rprt, phiprt);
    real jacB[9] = {B_dBc[1], B_dBc[2], B_dBc[3],
		    B_dBc[5], B_dBc[6], B_dBc[7],
		    B_dBc[9], B_dBc[10], B_dBc[11]};
    real B[3] = {B_dBc[0], B_dBc[4], B_dBc[8]};
    real B_norm = math_norm(B);
    real B_unit[3];
    B_unit[0] = B[0]/B_norm;
    B_unit[1] = B[1]/B_norm;
    B_unit[2] = B[2]/B_norm;

    real pitch = math_dot(p_unit,B_unit);

    /* Evaluate some required magnetic field quantities */
    real gradB[3];
    math_matmul(jacB,B_unit,3,3,1,gradB);

    real curlB[3] = {jacB[7]-jacB[5], jacB[2]-jacB[6], jacB[3]-jacB[1]};

    real tau_B = math_dot(B_unit, curlB)/B_norm;
    real nablabhat[9];
    nablabhat[0]=(jacB[0] - gradB[0]*B_unit[0])/B_norm;
    nablabhat[1]=(jacB[1] - gradB[0]*B_unit[1])/B_norm;
    nablabhat[2]=(jacB[2] - gradB[0]*B_unit[2])/B_norm;
    nablabhat[3]=(jacB[3] - gradB[1]*B_unit[0])/B_norm;
    nablabhat[4]=(jacB[4] - gradB[1]*B_unit[1])/B_norm;
    nablabhat[5]=(jacB[5] - gradB[1]*B_unit[2])/B_norm;
    nablabhat[6]=(jacB[6] - gradB[2]*B_unit[0])/B_norm;
    nablabhat[7]=(jacB[7] - gradB[2]*B_unit[1])/B_norm;
    nablabhat[8]=(jacB[8] - gradB[2]*B_unit[2])/B_norm; 

    real kappa[3];
    math_matmul(nablabhat,B_unit,3,3,1,kappa);
    
    /* Zeroth order momentum terms */
    real p_para0 = pitch*ptot;
    real mu_0 = ( 1 - pitch*pitch )*(ptot*ptot)/(2*mass*B_norm);

    /* Make the spatial transformation */
    real rho[3];
    real pperp[3];
    TEMP_V1[0] = p_unit[0]*(1 - B_unit[0]);
    TEMP_V1[1] = p_unit[1]*(1 - B_unit[1]);
    TEMP_V1[2] = p_unit[2]*(1 - B_unit[2]);
    math_unit(TEMP_V1,pperp);

    math_cross(B_unit, pperp, rho);
    math_prod(rho, sqrt(2*mass*mu_0/B_norm)/charge );

    real rho_unit[3];
    math_unit(rho,rho_unit);

    real xyz[3];
    xyz[0] = prtxyz[0] - rho[0];
    xyz[1] = prtxyz[1] - rho[1];
    xyz[2] = prtxyz[2] - rho[2];

    /* First order momentum terms */
    real perphat[3];
    math_cross(rho_unit, B_unit, perphat);

    real a1ddotgradb = -0.5*(2*(rho_unit[0]*perphat[0]*nablabhat[0]+
				rho_unit[1]*perphat[1]*nablabhat[4]+
				rho_unit[2]*perphat[2]*nablabhat[8])
			     +(rho_unit[0]*perphat[1]+rho_unit[1]*perphat[0])*
			     (nablabhat[1]+nablabhat[3])
			     +(rho_unit[0]*perphat[2]+rho_unit[2]*perphat[0])*
			     (nablabhat[2]+nablabhat[6])
			     +(rho_unit[1]*perphat[2]+rho_unit[2]*perphat[1])*
			     (nablabhat[5]+nablabhat[7]));

    real p_para1 = -p_para0*math_dot(rho,kappa)+((mass*mu_0)/charge)*(tau_B+a1ddotgradb);
    TEMP_S1 = pow(p_para0,2)/(mass*B_norm);
    TEMP_V1[0] = TEMP_S1 * kappa[0] + mu_0*gradB[0]/B_norm;
    TEMP_V1[1] = TEMP_S1 * kappa[1] + mu_0*gradB[1]/B_norm;
    TEMP_V1[2] = TEMP_S1 * kappa[2] + mu_0*gradB[2]/B_norm;
    real mu_1 = math_dot(rho,TEMP_V1)-((mu_0*p_para0)/(charge*B_norm))*(tau_B+a1ddotgradb);

    /* Make the momentum transformation */
    ppar[0] = p_para0 + p_para1;
    mu[0]   = mu_0 + mu_1;

    /* Back to cylindrical system */
    real rpz[3];
    math_xyz2rpz(xyz,rpz);
    R[0]   = rpz[0];
    phi[0] = rpz[1];
    z[0]   = rpz[2];

    /* Calculate gyroangle (this is zeroth order) */
    real a1[3];
    real z_unit[3];
    z_unit[0] = 0.0;
    z_unit[1] = 0.0;
    z_unit[2] = 1.0;
    math_cross(B_unit,z_unit,a1);
    math_unit(a1,a1);

    real a2[3];
    math_cross(a1,B_unit,a2);
    math_unit(a2,a2);

    theta[0] = atan2(math_dot(rho_unit,a2),math_dot(rho_unit,a1));
    theta[0] = fmod(theta[0] + CONST_2PI, CONST_2PI);
}

/**
 * @brief First order guiding center transformation from gc to fo
 *
 * @param mass   guiding center mass
 * @param charge guiding center charge
 * @param B_dB   magnetic field and Jacobian at guiding center location
 * @param R, phi, z guiding center location
 * @param mu, ppar, theta guiding center magnetic moment,
 *        parallel momentum and gyroangle
 * @param Rprt, phiprt, zprt pointers to returned particle location
 * @param pR, pphi, pz pointers to returned particle momentum components
 */
void physlib_gc2fo(real mass, real charge, real* B_dB,
		   real R, real phi, real z, real mu, real ppar, real theta,
		   real* Rprt, real* phiprt, real* zprt, real* pR, real* pphi, real* pz) {
     /* Temporary re-used variables */
    real TEMP_S1;
    real TEMP_V1[3];

    /* Helper variables */
    real B_dBc[12];
    math_jac_rpz2xyz(B_dB,B_dBc,R,phi);
    real xyz[3];
    real rpz[3] = {R, phi, z};
    math_rpz2xyz(rpz,xyz);

    real B_vec[3]  = {B_dBc[0], B_dBc[4], B_dBc[8]};
    real B_norm = math_norm(B_vec);
    real B_unit[3];
    math_unit(B_vec, B_unit);

    /* Magnetic field jacobian, gradient and curl in cartesian coordinates */
    real jacB[9] = {B_dBc[1], B_dBc[2],  B_dBc[3],
		    B_dBc[5], B_dBc[6],  B_dBc[7],
		    B_dBc[9], B_dBc[10], B_dBc[11]};

    real gradB[3];
    gradB[0] = B_unit[0]*jacB[0] + B_unit[0]*jacB[3] + B_unit[0]*jacB[6];
    gradB[1] = B_unit[1]*jacB[1] + B_unit[1]*jacB[4] + B_unit[1]*jacB[7];
    gradB[2] = B_unit[2]*jacB[2] + B_unit[2]*jacB[5] + B_unit[2]*jacB[8];
    math_matmul(jacB,B_unit,3,3,1,gradB);

    real curlB[3] = {jacB[7]-jacB[5], jacB[2]-jacB[6], jacB[3]-jacB[1]};

    real tau_B = math_dot(B_unit, curlB)/B_norm;
    real nablabhat[9];
    nablabhat[0]=(jacB[0] - gradB[0]*B_unit[0])/B_norm;
    nablabhat[1]=(jacB[1] - gradB[0]*B_unit[1])/B_norm;
    nablabhat[2]=(jacB[2] - gradB[0]*B_unit[2])/B_norm;
    nablabhat[3]=(jacB[3] - gradB[1]*B_unit[0])/B_norm;
    nablabhat[4]=(jacB[4] - gradB[1]*B_unit[1])/B_norm;
    nablabhat[5]=(jacB[5] - gradB[1]*B_unit[2])/B_norm;
    nablabhat[6]=(jacB[6] - gradB[2]*B_unit[0])/B_norm;
    nablabhat[7]=(jacB[7] - gradB[2]*B_unit[1])/B_norm;
    nablabhat[8]=(jacB[8] - gradB[2]*B_unit[2])/B_norm; 

    real kappa[3];
    math_matmul(nablabhat,B_unit,3,3,1,kappa);

    /* Zeroth order momentum terms */
    real p_para0 = ppar;
    real mu_0 = mu;

    /* Calculate gyroradius and its unit vector. Choose basis vectors a1,
       which is perpendicular to B and z (since we always have Bphi),
       and a2 which is perpendicular to a1 and B */
    real a1[3];
    real z_unit[3];
    z_unit[0] = 0.0;
    z_unit[1] = 0.0;
    z_unit[2] = 1.0;
    math_cross(B_unit,z_unit,a1);
    math_unit(a1,a1);

    real a2[3];
    math_cross(a1,B_unit,a2);
    math_unit(a2,a2);

    real rho_unit[3];
    rho_unit[0] = cos(theta)*a1[0]+sin(theta)*a2[0]; 
    rho_unit[1] = cos(theta)*a1[1]+sin(theta)*a2[1];
    rho_unit[2] = cos(theta)*a1[2]+sin(theta)*a2[2];

    real rho[3];
    rho[0] = rho_unit[0];
    rho[1] = rho_unit[1];
    rho[2] = rho_unit[2];
    math_prod(rho,sqrt(2*mass*mu/B_norm)/charge);

    /* First order momentum terms */
    real perphat[3];
    math_cross(rho_unit,B_unit,perphat); 

/*    real a1ddotgradb = -0.5*(2*(rho_unit[0]*perphat[0]*nablabhat[0]+
				rho_unit[1]*perphat[1]*nablabhat[4]+
				rho_unit[2]*perphat[2]*nablabhat[8])
			     +(rho_unit[0]*perphat[1]+rho_unit[1]*perphat[0])*
			      (nablabhat[1]+nablabhat[3])
			     +(rho_unit[0]*perphat[2]+rho_unit[2]*perphat[0])*
			      (nablabhat[2]+nablabhat[6])
			     +(rho_unit[1]*perphat[2]+rho_unit[2]*perphat[1])*
			      (nablabhat[5]+nablabhat[7]));
*/
    real A1[9] = {rho_unit[0]*perphat[0], rho_unit[1]*perphat[0], rho_unit[2]*perphat[0],
		  rho_unit[0]*perphat[1], rho_unit[1]*perphat[1], rho_unit[2]*perphat[1],
		  rho_unit[0]*perphat[2], rho_unit[1]*perphat[2], rho_unit[2]*perphat[2]};
    real a1ddotgradb = 0;
    for(int i=0; i<9; i++) {
	a1ddotgradb += A1[i]*nablabhat[i];
    }
    a1ddotgradb *= -0.5;

    real p_para1 = -p_para0*math_dot(rho,kappa)+((mass*mu_0)/charge)*(tau_B+a1ddotgradb);
    TEMP_S1 = pow(p_para0,2)/(mass*B_norm);
    TEMP_V1[0] = TEMP_S1 * kappa[0] + mu_0*gradB[0]/B_norm;
    TEMP_V1[1] = TEMP_S1 * kappa[1] + mu_0*gradB[1]/B_norm;
    TEMP_V1[2] = TEMP_S1 * kappa[2] + mu_0*gradB[2]/B_norm;
    real mu_1 = math_dot(rho,TEMP_V1)-((mu_0*p_para0)/(charge*B_norm))*(tau_B+a1ddotgradb);

    /* Calculate the parallel momentum vector */
    real p_para[3];
    p_para[0] = B_unit[0];
    p_para[1] = B_unit[1];
    p_para[2] = B_unit[2];
    math_prod(p_para,p_para0-p_para1);

    /* Calculate perpendicular momentum vector. */
    real p_perp[3];
    p_perp[0] = perphat[0];
    p_perp[1] = perphat[1];
    p_perp[2] = perphat[2];
    math_prod(p_perp,sqrt(2*mass*(mu_0-mu_1)*B_norm));

    /* Make the spatial transformation */
    real xyzprt[3];
    xyzprt[0] = xyz[0]+rho[0];
    xyzprt[2] = xyz[1]+rho[1];
    xyzprt[1] = xyz[2]+rho[2];

    real rpzprt[3];
    math_xyz2rpz(xyzprt,rpzprt);
    Rprt[0]   = rpzprt[0];
    phiprt[0] = rpzprt[1];
    zprt[0]   = rpzprt[2];

    /* Calculate the momentum vector */
    real pxyz[3];
    pxyz[0]   = p_para[0]+p_perp[0];
    pxyz[1]   = p_para[1]+p_perp[1];
    pxyz[2]   = p_para[2]+p_perp[2];
    
    real prpz[3];
    math_vec_xyz2rpz(pxyz,prpz,phiprt[0]);
    pR[0]   = prpz[0];
    pphi[0] = prpz[1];
    pz[0]   = prpz[2];
}
