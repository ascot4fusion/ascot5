/** @author Konsta Sarkimaki konsta.sarkimaki@aalto.fi
*
* A module for testing mccc package
*/
#define _XOPEN_SOURCE
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../../math.h"
#include "../../ascot5.h"
#include "../../consts.h"
#include "mccc_wiener.h"
#include "mccc_special.h"
#include "mccc_coefs.h"
#include "mccc_push.h"
#include "mccc.h"

void mccc_initdist(int nbin, real* dist){
    for(int i=0; i<nbin; i=i+1){
	dist[i] = 0;
    }
}

void mccc_updatedist(real val, real dt, real xmin, real xmax, int nbin, real* dist){
    real dx = (xmax-xmin)/nbin;
    int bin = floor((val-xmin)/dx);
    if(bin >= nbin){bin=nbin-1;};
    dist[bin] = dist[bin] + dt;
}

void mccc_writedist(FILE* fn, real xmin, real xmax, int nbin, real* distA, real* distB, real* distC){
    real dx = (xmax-xmin)/nbin;
    for(int i=0; i<nbin; i=i+1){
	fprintf(fn,"%g %g %g %g\n",xmin+dx*i,distA[i],distB[i],distC[i]);
    }
}

void mccc_writeonedist(FILE* fn, real xmin, real xmax, int nbin, real* dist){
    real dx = (xmax-xmin)/nbin;
    for(int i=0; i<nbin; i=i+1){
	fprintf(fn,"%g %g\n",xmin+dx*i,dist[i]);
    }
}

void mccc_push_test(){

    real me = 9.10938356e-31;
    real mp = 1.6726219e-27;
    real q = 1.60217662e-19;

    real ma, qa, va;

    int nspec = 2;
    real clogab[2];
    real Fb[2],Dparab[2],Dperpb[2],Kb[2],nub[2],DXb[2];
    real F, Dpara, Dperp, dt, rnd[3], vin[3], vout[3];
    real K,nu,rndgc[5], DX, xiin, xiiout, vgcin, vgcout, Xin[3], Xout[3];

    int mode;
    int err,i,j,k;

    real mb[2] = {me, mp};
    real qb[2] = {-q, q};
    real nb[2] = {1e20, 1e20};
    real Tb[2] = {10e3*q, 10e3*q};
    real B[3] = {1, 0, 0};

    real tprt,tmax,Eprt,Emin,cutoff;

    real Ethmin=0,Ethmax=1e5;
    int Ethnbin = 1000;
    real Esdmin=0,Esdmax=1.5e6;
    int Esdnbin = 100;
    real ximin=-1,ximax=1;
    int xibin = 50;
    real rmin=0,rmax=0.1;
    int rbin = 50;
    real Dmin=0,Dmax=0.005;
    int Dbin = 50;
    real* distA;
    real* distB;
    real* distC;
    real* distxi;
    real* distr;
    real* distD;

    real Ekin, pitch, r, D;

 

    srand48(1);

    printf("\n");
    printf("This test program tests functionality of mccc push operators.\n");
    printf("Tests are done for fo-fixed, gc-fixed, and gc-adaptive modes.\n");
    printf("\n");
    printf("First we check that a thermal test electron in a electron-proton.\n");
    printf("plasma obeys Maxwellian distribution.\n");
    printf("\n");
    printf("Calculating...\n");

    distA = malloc(Ethnbin*sizeof(real));
    distB = malloc(Ethnbin*sizeof(real));
    distC = malloc(Ethnbin*sizeof(real));
    distxi = malloc(xibin*sizeof(real));
    distr = malloc(rbin*sizeof(real));
    distD = malloc(Dbin*sizeof(real));

    mccc_initdist(Ethnbin,distA);
    mccc_initdist(Ethnbin,distB);
    mccc_initdist(Ethnbin,distC);
    mccc_initdist(xibin,distxi);
    mccc_initdist(rbin,distr);
    mccc_initdist(Dbin,distD);

    for(mode=1; mode<2; mode=mode+1){

	// Options
	tmax = 1e-1;
	dt = 1e-7;

	// Init a particle
	ma = me;
	qa = -q;
	tprt = 0;
	vin[0] = 0;
	vin[1] = 0;
	vin[2] = sqrt(2*Tb[0]/ma);

	vgcin = math_norm(vin);
	xiin = 0;// 1.0-2*drand48();//pitch = [-1,1]
	Xin[0] = 0;
	Xin[1] = 0;
	Xin[2] = 0;

	cutoff = 0.1*sqrt(2*Tb[0]/ma);

	tprt = 0;
	while(tprt < tmax){
	    
	    switch(mode) {
	    case 0 :
		// Update coefficients
		va = math_norm(vin);
		mccc_coefs_clog(ma,qa,va,mb,qb,nb,Tb,clogab,nspec);
		mccc_coefs_fo(ma,qa,va,mb,qb,nb,Tb,clogab,nspec, 
			      Fb,Dparab,Dperpb,Kb,nub);
		F=0;
		Dpara=0;
		Dperp=0;
		for(i=0; i<nspec; i=i+1){
		    F=F+Fb[i];
		    Dpara=Dpara+Dparab[i];
		    Dperp=Dperp+Dperpb[i];
		}
		rnd[0] = 2*(drand48() > 0.5)-1.0;
		rnd[1] = 2*(drand48() > 0.5)-1.0;
		rnd[2] = 2*(drand48() > 0.5)-1.0;

		// Take a time step
		mccc_push_foEM(F,Dpara,Dperp,dt,rnd,vin,vout,&err);
		mccc_printerror(err);
	   

		// Update distributions
		Ekin = 0.5*ma*va*va/q;
		mccc_updatedist(Ekin,dt,Ethmin,Ethmax,Ethnbin,distA);

		// Update particle
		vin[0] = vout[0];
		vin[1] = vout[1];
		vin[2] = vout[2];
		break;

	    case 1 :
		// Update coefficients
		mccc_coefs_clog(ma,qa,vgcin,mb,qb,nb,Tb,clogab,nspec);
		mccc_coefs_gcfixed(ma,qa,vgcin,xiin,mb,qb,nb,Tb,math_norm(B),clogab,nspec, 
			Dparab,DXb,Kb,nub);
		K = 0;
		nu = 0;
		Dpara = 0;
		DX = 0;
		for(i=0; i<nspec; i=i+1){
		    K=K+Kb[i];
		    nu=nu+nub[i];
		    Dpara=Dpara+Dparab[i];
		    DX=DX+DXb[i];
		}
		rndgc[0] = 1-2*drand48();//2*(drand48() > 0.5)-1.0;
		rndgc[1] = 1-2*drand48();//2*(drand48() > 0.5)-1.0;
		rndgc[2] = 1-2*drand48();//2*(drand48() > 0.5)-1.0;
		rndgc[3] = 2*(drand48() > 0.5)-1.0;
		rndgc[4] = 2*(drand48() > 0.5)-1.0;
		Ekin = 0.5*ma*vgcout*vgcout/q;

		// Take a time step
		mccc_push_gcEM(K,nu,Dpara,DX,B,dt,rndgc, 
			       vgcin,&vgcout,xiin,&xiiout,Xin,Xout,cutoff,&err);
		mccc_printerror(err);

		// Update distributions
		Ekin = 0.5*ma*vgcin*vgcin/q;
		r = sqrt(pow(Xin[1],2)+pow(Xin[2],2))/dt;
		D = (pow(Xin[1]-Xout[1],2)+pow(Xin[2]-Xout[2],2))/(2*dt);
		mccc_updatedist(Ekin,dt,Ethmin,Ethmax,Ethnbin,distB);
		mccc_updatedist(xiin,dt,ximin,ximax,xibin,distxi);
		mccc_updatedist(r,dt,rmin,rmax,rbin,distr);
		mccc_updatedist(D,dt,Dmin,Dmax,Dbin,distD);
		
		Xin[0] = Xout[0];
		Xin[1] = Xout[1];
		Xin[2] = Xout[2];
		//vgcin = vgcout;
		//xiin = xiiout;
		break;
	    case 2 :
		// Update coefficients
		mccc_coefs_clog(ma,qa,vgcin,mb,qb,nb,Tb,clogab,nspec);
		mccc_coefs_gcfixed(ma,qa,vgcin,xiin,mb,qb,nb,Tb,math_norm(B),clogab,nspec, 
				   Dparab,DXb,Kb,nub);
		K = 0;
		nu = 0;
		Dpara = 0;
		DX = 0;
		for(i=0; i<nspec; i=i+1){
		    K=K+Kb[i];
		    nu=nu+nub[i];
		    Dpara=Dpara+Dparab[i];
		    DX=DX+DXb[i];
		}
		rndgc[0] = 2*(drand48() > 0.5)-1.0;
		rndgc[1] = 2*(drand48() > 0.5)-1.0;
		rndgc[2] = 2*(drand48() > 0.5)-1.0;
		rndgc[3] = 2*(drand48() > 0.5)-1.0;
		rndgc[4] = 2*(drand48() > 0.5)-1.0;
		Ekin = 0.5*ma*vgcout*vgcout/q;

		// Take a time step
		mccc_push_gcEM(K,nu,Dpara,DX,B,dt,rndgc, 
			       vgcin,&vgcout,xiin,&xiiout,Xin,Xout,cutoff,&err);
		mccc_printerror(err);

		// Update distributions
		Ekin = 0.5*ma*vgcin*vgcin/q;
		mccc_updatedist(Ekin,dt,Ethmin,Ethmax,Ethnbin,distB);

		Xin[0] = Xout[0];
		Xin[1] = Xout[1];
		Xin[2] = Xout[2];
		vgcin = vgcout;
		xiin = xiiout;
		break;
	    }
	
        
	    if(err != 0){exit(0);}

	
	    tprt = tprt+dt;
	}
    }
    
    // Write distributions
    FILE* fn_th = fopen("mccc_push_thermal.test","w");
    mccc_writedist(fn_th,Ethmin,Ethmax,Ethnbin,distA,distB,distC);
    fclose(fn_th);
    FILE* fn_th_xi = fopen("mccc_thermal_pitch.test","w");
    mccc_writeonedist(fn_th_xi,ximin,ximax,xibin,distxi);
    fclose(fn_th_xi);
    FILE* fn_th_r = fopen("mccc_thermal_r.test", "w");
    mccc_writeonedist(fn_th_r,rmin,rmax,rbin,distr);
    fclose(fn_th_r);
    FILE* fn_th_D = fopen("mccc_thermal_diffusion.test","w");
    mccc_writeonedist(fn_th_D,Dmin,Dmax,Dbin,distD);
    fclose(fn_th_D);
    free(distA);
    free(distB);
    free(distC);
    free(distxi);
    free(distr);
    free(distD);

    printf("\n");
    printf("Done. Results are written in mccc_push_thermal.test\n");
    printf("\n");
    printf("Next we check that fast test protons in a electron-proton.\n");
    printf("plasma obeys slowing-down distribution.\n");
    printf("\n");
    printf("Calculating...\n");


    distA = malloc(Esdnbin*sizeof(real));
    distB = malloc(Esdnbin*sizeof(real));
    distC = malloc(Esdnbin*sizeof(real));
    distxi = malloc(xibin*sizeof(real));
    distr = malloc(rbin*sizeof(real));
    distD = malloc(Dbin*sizeof(real));

    mccc_initdist(Esdnbin,distA);
    mccc_initdist(Esdnbin,distB);
    mccc_initdist(Esdnbin,distC);
    mccc_initdist(xibin,distxi);
    mccc_initdist(rbin,distr);
    mccc_initdist(Dbin,distD);

    
    // Options
    int Nprt = 100;
    dt = 1e-6;
    for(mode=1; mode<2; mode=mode+1){
	// Init a particle
	ma = mp;
	qa = q;
	cutoff = 0;
    
	Emin = 1*Tb[0];
	for(j = 0; j < Nprt; j = j+1){
	    // Init a new particle
	    Eprt = 1e6*q;
	    vin[0] = 0;
	    vin[1] = 0;
	    vin[2] = sqrt(2*Eprt/ma);
	
	    vgcin = math_norm(vin);
	    xiin =  1.0-2*drand48(); 
	    Xin[0] = 0;
	    Xin[1] = 0;
	    Xin[2] = 0;
	    
	    while(Eprt > Emin){
		switch(mode){
		case 0 :
		    // Update coefficients
		    va = math_norm(vin);
		    mccc_coefs_clog(ma,qa,va,mb,qb,nb,Tb,clogab,nspec);
		    mccc_coefs_fo(ma,qa,va,mb,qb,nb,Tb,clogab,nspec, 
				  Fb,Dparab,Dperpb,Kb,nub);
		    F=0;
		    Dpara=0;
		    Dperp=0;
		    for(i=0; i<nspec; i=i+1){
			F=F+Fb[i];
			Dpara=Dpara+0*Dparab[i];
			Dperp=Dperp+0*Dperpb[i];
		    }
		    rnd[0] = 2*(drand48() > 0.5)-1.0;
		    rnd[1] = 2*(drand48() > 0.5)-1.0;
		    rnd[2] = 2*(drand48() > 0.5)-1.0;

		    // Take a time step
		    mccc_push_foEM(F,Dpara,Dperp,dt,rnd,vin,vout,&err);
		    mccc_printerror(err);
	    
		    if(err != 0){exit(0);}

		    // Update distributions
		    Eprt = 0.5*ma*va*va;
		    mccc_updatedist(Eprt/q,dt,Esdmin,Esdmax,Esdnbin,distA);

		    // Update particle
		    vin[0] = vout[0];
		    vin[1] = vout[1];
		    vin[2] = vout[2];

		    break;

		case 1 :
		    mccc_coefs_clog(ma,qa,vgcin,mb,qb,nb,Tb,clogab,nspec);
		    mccc_coefs_gcfixed(ma,qa,vgcin,xiin,mb,qb,nb,Tb,math_norm(B),clogab,nspec, 
				       Dparab,DXb,Kb,nub);
		    K = 0;
		    nu = 0;
		    Dpara = 0;
		    DX = 0;
		    for(i=0; i<nspec; i=i+1){
			K=K+Kb[i];
			nu=nu+nub[i];
			Dpara=Dpara+Dparab[i];
			DX=DX+DXb[i];
		    }
		    rndgc[0] = 2*(drand48() > 0.5)-1.0;
		    rndgc[1] = 2*(drand48() > 0.5)-1.0;
		    rndgc[2] = 2*(drand48() > 0.5)-1.0;
		    rndgc[3] = 2*(drand48() > 0.5)-1.0;
		    rndgc[4] = 2*(drand48() > 0.5)-1.0;

		    // Take a time step
		    mccc_push_gcEM(K,nu,Dpara,DX,B,dt,rndgc, 
				   vgcin,&vgcout,xiin,&xiiout,Xin,Xout,cutoff,&err);
		    mccc_printerror(err);

		    // Update distributions
		    Eprt = 0.5*ma*vgcin*vgcin/q;
		    r = sqrt(pow(Xin[1],2)+pow(Xin[2],2))/(2*dt);
		    D = (pow(Xin[1]-Xout[1],2)+pow(Xin[2]-Xout[2],2))/dt;
		    mccc_updatedist(Eprt,dt,Esdmin,Esdmax,Esdnbin,distB);
		    mccc_updatedist(xiin,dt,ximin,ximax,xibin,distxi);
		    mccc_updatedist(r,dt,rmin,rmax,rbin,distr);
		    mccc_updatedist(D,dt,Dmin,Dmax,Dbin,distD);

		    Xin[0] = Xout[0];
		    Xin[1] = Xout[1];
		    Xin[2] = Xout[2];
		    vgcin = vgcout;
		    xiin = xiiout;
		    Eprt = 0.5*ma*vgcin*vgcin;
		    break;
		}
	
        
		if(err != 0){exit(0);}
	
	    }
	}
    }
    // Write distributions
    FILE* fn_sd = fopen("mccc_push_fast.test","w");
    mccc_writedist(fn_sd,Esdmin,Esdmax,Esdnbin,distA,distB,distC);
    fclose(fn_sd);
    FILE* fn_sd_xi = fopen("mccc_fast_pitch.test","w");
    mccc_writeonedist(fn_sd_xi,ximin,ximax,xibin,distxi);
    fclose(fn_sd_xi);
    FILE* fn_sd_r = fopen("mccc_fast_r.test","w");
    mccc_writeonedist(fn_sd_r,rmin,rmax,rbin,distr);
    fclose(fn_sd_r);
    FILE* fn_sd_D = fopen("mccc_fast_diffusion.test","w");
    mccc_writeonedist(fn_sd_D,Dmin,Dmax,Dbin,distD);
    fclose(fn_sd_D);
    free(distA);
    free(distB);
    free(distC);
    free(distxi);
    free(distr);
    free(distD);


    printf("\n");
    printf("Done. Results are written in mccc_push_fast.test\n");
    printf("\n");
}

void mccc_coefs_test(){

    real umin = 0;
    real umax = 10;
    int usteps = 1000;
    real du, u;

    real thmin = 0.01;
    real thmax = 0.1;
    int thsteps = 10;
    real dth, th;

    real density = 1e20;
    real me = 9.10938356e-31;
    real mp = 1.6726219e-27;
    real q = 1.60217662e-19;

    int nspec = 2;
    real clogab[2];
    real Tb[2];
    real F[2],Dpara[2],Dperp[2],K[2],nu[2];
    real x, va, vth;

    // Define plasma
    real nb[2] = {density, density};
    real mb[2] = {me, mp};
    real qb[2] = {-q, q};

    printf("\n");
    printf("This test program tests functionality of mccc special functions.\n");
    printf("The results are printed on mccc_coefs.test.\n");
    printf("Format: clog F Dpar Dperp Q mu.\n");
    printf("b={ep} x u x th x a={ep}\n");
    printf("u = gamma*v/c %d momentum steps\n",usteps);
    printf("th = T/mc^2 %d temperature steps\n",thsteps);
    printf("e = electron, p = proton, a = test particle, b=plasma \n");
    printf("\n");

    FILE* fn = fopen("mccc_coefs.test","w");

    printf("Computing...\n");
    for(int particle = 0; particle < 2; particle = particle + 1){

	// Define particle
	real ma;
	real qa;
	if(particle){
	    // electron
	    ma = me;
	    qa = -q;
	}
	else{
	    // proton
	    ma = mp;
	    qa = q;
	}

	// Iterate over u and th
	th = thmin;
	dth = (thmax-thmin)/thsteps;
	for(int thi = 0; thi < thsteps; thi = thi + 1){

	    u = umin;
	    du = (umax-umin)/usteps;
	    for(int ui = 0; ui < usteps; ui = ui + 1){
		
		va = CONST_C*u/sqrt(1+u*u);
		vth = sqrt(2*th*CONST_C2);
		Tb[0] = th*CONST_C2*mb[0];
		Tb[1] = th*CONST_C2*mb[1];
		x = va/vth;
		// Evaluate
		mccc_coefs_clog(ma,qa,va,mb,qb,nb,Tb,clogab,nspec);
		mccc_coefs_fo(ma,qa,va,mb,qb,nb,Tb,clogab,nspec, 
			      F,Dpara,Dperp,K,nu);
		// Write
		fprintf(fn,"%g %g %g %g %g %g\n",clogab[0],F[0],Dpara[0],Dperp[0],K[0],nu[0]);
		fprintf(fn,"%g %g %g %g %g %g\n",clogab[1],F[1],Dpara[1],Dperp[1],K[1],nu[1]);
		// Proceed
		u = u + du;
	    }
	    th = th + dth;
	}
	
    }


    fclose(fn);
    printf("Done.\n");
    printf("\n");

}

void mccc_special_test(){

    real umin = 0;
    real umax = 10;
    int usteps = 1000;
    real du, u;

    real thmin = 0.01;
    real thmax = 0.1;
    int thsteps = 10;
    real dth, th;

    real density = 1e20;
    real me = 9.10938356e-31;
    real mp = 1.6726219e-27;
    real q = 1.60217662e-19;

    real GdG[4];
    real mudmu[6];
    real x, va, vth;

    // Define plasma
    real nb[2] = {density, density};
    real mb[2] = {me, mp};
    real qb[2] = {-q, q};

    printf("\n");
    printf("This test program tests functionality of mccc special functions.\n");
    printf("The results are printed on mccc_special.test.\n");
    printf("Format: G G1 dG dG1 mu0 mu1 mu2 dmu0 dmu1 dmu2.\n");
    printf("u x th x a={ep} x 2\n");
    printf("u = gamma*v/c %d momentum steps\n",usteps);
    printf("th = T/mc^2 %d temperature steps\n",thsteps);
    printf("e = electron, p = proton, a = test particle, b=plasma \n");
    printf("2: exact and interpolated coefficients\n");
    printf("\n");

    FILE* fn = fopen("mccc_special.test","w");

    printf("Computing...\n");
    for(int exact = 1; exact > -1; exact = exact - 1){
	for(int particle = 0; particle < 2; particle = particle + 1){

	    // Define particle
	    real ma;
	    real qa;
	    if(particle){
		// electron
		ma = me;
		qa = -q;
	    }
	    else{
		// proton
		ma = mp;
		qa = q;
	    }

	    // Iterate over u and th
	    th = thmin;
	    dth = (thmax-thmin)/thsteps;
	    for(int thi = 0; thi < thsteps; thi = thi + 1){

		u = umin;
		du = (umax-umin)/usteps;
		for(int ui = 0; ui < usteps; ui = ui + 1){
		
		    va = CONST_C*u/sqrt(1+u*u);
		    vth = sqrt(2*th*CONST_C2);
		    x = va/vth;
		    // Evaluate
		    mccc_special_GdG(x, GdG, exact);
		    mccc_special_mudmu(u, th, mudmu, exact);
		    // Write
		    fprintf(fn,"%g %g %g %g %g %g %g %g %g %g\n",GdG[0],GdG[1],GdG[2],GdG[3],mudmu[0],mudmu[1],mudmu[2],mudmu[3],mudmu[4],mudmu[5]);
		    // Proceed
		    u = u + du;
		}
		th = th + dth;
	    }
	}
    }


    fclose(fn);
    printf("Done.\n");
    printf("\n");

}

/** Tests functionality of wiener process handling.
 *
 */
void mccc_wiener_test(){

    mccc_wienarr *wienarr, *temparr;
    int Ndim = 5, Nslot = 5, windex;
    int j,k;
    int Niter = 10000, iter;
    int err;
    real W[30000];
    real initime = 0.0;
    real meanCom[5], meanThr[5], varCom[5], varThr[5];

    srand48(1);

    printf("\n");
    printf("This test program tests functionality of wiener module.\n");
    printf("The results are printed directly on the screen.\n");
    printf("\n");

    printf("\n");
    printf("Initializing a 5 slot array of 3D Wiener processes...\n");
    mccc_wiener_initialize(wienarr,initime);
    mccc_wiener_initialize(temparr,initime);
    printf("Done. The wiener, time, and index values are:\n");
    printf("\n");
    printf("%g %g %g %g %g\n",wienarr->wiener[0],wienarr->wiener[1*Ndim],wienarr->wiener[2*Ndim],wienarr->wiener[3*Ndim],wienarr->wiener[4*Ndim]);
    printf("%g %g %g %g %g\n",wienarr->time[0],wienarr->time[1],wienarr->time[2],wienarr->time[3],wienarr->time[4]);
    printf("%d %d %d %d %d\n",wienarr->nextslot[0],wienarr->nextslot[1],wienarr->nextslot[2],wienarr->nextslot[3],wienarr->nextslot[4]);
    printf("\n");
    printf("Expected: [0 -999 -999 -999 -999]\n");
    printf("\n");

    printf("\n");
    printf("Generating a new Wiener process for t = 1.0.\n");
    mccc_wiener_generate(wienarr, 1.0, &windex, &err);
    mccc_printerror(err);
    printf("Done. The wiener, time, and index values are:\n");
    printf("\n");
    printf("%g %g %g %g %g\n",wienarr->wiener[0],wienarr->wiener[1*Ndim],wienarr->wiener[2*Ndim],wienarr->wiener[3*Ndim],wienarr->wiener[4*Ndim]);
    printf("%g %g %g %g %g\n",wienarr->time[0],wienarr->time[1],wienarr->time[2],wienarr->time[3],wienarr->time[4]);
    printf("%d %d %d %d %d\n",wienarr->nextslot[0],wienarr->nextslot[1],wienarr->nextslot[2],wienarr->nextslot[3],wienarr->nextslot[4]);
    printf("\n");	
    printf("Expected: [1 1 -999 -999 -999]");	
    printf("\n");

    printf("\n");
    printf("Generating a new Wiener process for t = 0.5.\n");
    mccc_wiener_generate(wienarr, 0.5, &windex, &err);
    mccc_printerror(err);
    printf("Done. The wiener, time, and index values are:\n");
    printf("\n");
    printf("%g %g %g %g %g\n",wienarr->wiener[0],wienarr->wiener[1*Ndim],wienarr->wiener[2*Ndim],wienarr->wiener[3*Ndim],wienarr->wiener[4*Ndim]);
    printf("%g %g %g %g %g\n",wienarr->time[0],wienarr->time[1],wienarr->time[2],wienarr->time[3],wienarr->time[4]);
    printf("%d %d %d %d %d\n",wienarr->nextslot[0],wienarr->nextslot[1],wienarr->nextslot[2],wienarr->nextslot[3],wienarr->nextslot[4]);
    printf("\n");	
    printf("Expected: [2 1 1 -999 -999]");	
    printf("\n");

    printf("\n");
    printf("Generating a new Wiener process for t = 0.75.\n");
    mccc_wiener_generate(wienarr, 0.75, &windex, &err);
    mccc_printerror(err);
    printf("Done. The wiener, time, and index values are:\n");
    printf("\n");
    printf("%g %g %g %g %g\n",wienarr->wiener[0],wienarr->wiener[1*Ndim],wienarr->wiener[2*Ndim],wienarr->wiener[3*Ndim],wienarr->wiener[4*Ndim]);
    printf("%g %g %g %g %g\n",wienarr->time[0],wienarr->time[1],wienarr->time[2],wienarr->time[3],wienarr->time[4]);
    printf("%d %d %d %d %d\n",wienarr->nextslot[0],wienarr->nextslot[1],wienarr->nextslot[2],wienarr->nextslot[3],wienarr->nextslot[4]);
    printf("\n");	
    printf("Expected: [2 1 3 1 -999]");	
    printf("\n");

    printf("\n");
    printf("Generating a new Wiener process for t = 2.0.\n");
    mccc_wiener_generate(wienarr, 2.0, &windex, &err);
    mccc_printerror(err);
    printf("Done. The wiener, time, and index values are:\n");
    printf("\n");
    printf("%g %g %g %g %g\n",wienarr->wiener[0],wienarr->wiener[1*Ndim],wienarr->wiener[2*Ndim],wienarr->wiener[3*Ndim],wienarr->wiener[4*Ndim]);
    printf("%g %g %g %g %g\n",wienarr->time[0],wienarr->time[1],wienarr->time[2],wienarr->time[3],wienarr->time[4]);
    printf("%d %d %d %d %d\n",wienarr->nextslot[0],wienarr->nextslot[1],wienarr->nextslot[2],wienarr->nextslot[3],wienarr->nextslot[4]);
    printf("\n");	
    printf("Expected: [2 4 3 1 4]");	
    printf("\n");

    printf("\n");
    printf("Next we try to generate value at t = 1.5.\n");
    mccc_wiener_generate(wienarr, 1.5, &windex, &err);
    mccc_printerror(err);
    printf("This should have resulted in error.\n");
    printf("\n");

    printf("\n");
    printf("Assuming simulation time is 1.0, clean the array.\n");
    mccc_wiener_clean(wienarr,1.0, &err);
    mccc_printerror(err);
    printf("Done. The wiener, time, and index values are:\n");
    printf("\n");
    printf("%g %g %g %g %g\n",wienarr->wiener[2],wienarr->wiener[1*Ndim],wienarr->wiener[2*Ndim],wienarr->wiener[3*Ndim],wienarr->wiener[4*Ndim]);
    printf("%g %g %g %g %g\n",wienarr->time[0],wienarr->time[1],wienarr->time[2],wienarr->time[3],wienarr->time[4]);
    printf("%d %d %d %d %d\n",wienarr->nextslot[0],wienarr->nextslot[1],wienarr->nextslot[2],wienarr->nextslot[3],wienarr->nextslot[4]);
    printf("\n");
    printf("Expected: [4 -999 -999 -999 4]\n");
    printf("\n");

    printf("\n");
    printf("Finally, check that Wiener processes are indeed generated correctly.\n");
    printf("With simulation time 1.0, mean and var for W(3.0) calculated N = 1000\n");
    printf("iterations is (expected values are in parenthesis)\n");
    printf("\n");

    mccc_wiener_clean(wienarr,2.0, &err);
    mccc_printerror(err);
    for(j = 0; j < Nslot; j = j+1){
	temparr->time[j] = wienarr->time[j];
	temparr->nextslot[j] = wienarr->nextslot[j];
	for(k = 0; k < Ndim; k = k+1){
	    temparr->wiener[j*Ndim + k] = wienarr->wiener[j*Ndim + k];
	}
    }
    
    for (iter = 0; iter < Niter; iter = iter+1){
	
	mccc_wiener_generate(temparr, 3.0, &windex, &err);
        mccc_printerror(err);
	W[iter*Ndim+0] = temparr->wiener[windex*Ndim+0];
	W[iter*Ndim+1] = temparr->wiener[windex*Ndim+1];
	W[iter*Ndim+2] = temparr->wiener[windex*Ndim+2];

	for(j = 0; j < Nslot; j = j+1){
	    temparr->time[j] = wienarr->time[j];
	    temparr->nextslot[j] = wienarr->nextslot[j];
	}
    }

    //mccc_wiener_generate(temparr, 2.0, lastindex,err);
    //mccc_wiener_error(err);
    meanCom[0] = 0;
    meanCom[1] = 0;
    meanCom[2] = 0;
    for (iter = 0; iter < Niter; iter = iter+1){
	meanCom[0] = meanCom[0] + W[iter*Ndim+0]/Niter;
	meanCom[1] = meanCom[1] + W[iter*Ndim+1]/Niter;
	meanCom[2] = meanCom[2] + W[iter*Ndim+2]/Niter;
    }
    meanThr[0] = wienarr->wiener[0];
    meanThr[1] = wienarr->wiener[1];
    meanThr[2] = wienarr->wiener[2];

    varCom[0] = 0;
    varCom[1] = 0;
    varCom[2] = 0;
    for (iter = 0; iter < Niter; iter = iter+1){
	varCom[0] = varCom[0] + pow(W[iter*Ndim+0]-meanCom[0],2)/(Niter-1);
	varCom[1] = varCom[1] + pow(W[iter*Ndim+1]-meanCom[1],2)/(Niter-1);
	varCom[2] = varCom[2] + pow(W[iter*Ndim+2]-meanCom[2],2)/(Niter-1);
    }
    varThr[0] = (3.0-2.0);
    varThr[1] = (3.0-2.0);
    varThr[2] = (3.0-2.0);

    printf("\n");
    printf("- When W(t>3.0) does not exist:\n");
    printf("%g, %g",meanCom[0],varCom[0]);
    printf("  (%g, %g)\n",meanThr[0],varThr[0]);
    printf("%g, %g",meanCom[1],varCom[1]);
    printf("  (%g, %g)\n",meanThr[1],varThr[1]);
    printf("%g, %g",meanCom[2],varCom[2]);
    printf("  (%g, %g)\n",meanThr[2],varThr[2]);
    printf("\n");

    mccc_wiener_generate(wienarr, 3.1, &windex, &err);
    mccc_printerror(err);
    for(j = 0; j < Nslot; j = j+1){
	temparr->time[j] = wienarr->time[j];
	temparr->nextslot[j] = wienarr->nextslot[j];
	for(k = 0; k < wienarr->Ndim; k = k+1){
	    temparr->wiener[j*Ndim + k] = wienarr->wiener[j*Ndim + k];
	}
    }
    
    for (iter = 0; iter < Niter; iter = iter+1){
	
	mccc_wiener_generate(temparr, 3.0, &windex, &err);
        mccc_printerror(err);
	W[iter*Ndim+0] = temparr->wiener[windex*Ndim+0];
	W[iter*Ndim+1] = temparr->wiener[windex*Ndim+1];
	W[iter*Ndim+2] = temparr->wiener[windex*Ndim+2];

	for(j = 0; j < Nslot; j = j+1){
	    temparr->time[j] = wienarr->time[j];
	    temparr->nextslot[j] = wienarr->nextslot[j];
	}
    }

    // Analytical values for Brownian bridge:
    // mean = W(tm) + ( W(ip)-W(im) )*(t-tm)/(tp-tm)
    // variance = (t-tm)*(tp-t)/(tp-tm)
    meanCom[0] = 0;
    meanCom[1] = 0;
    meanCom[2] = 0;
    for (iter = 0; iter < Niter; iter = iter+1){
	meanCom[0] = meanCom[0] + W[iter*Ndim+0]/Niter;
	meanCom[1] = meanCom[1] + W[iter*Ndim+1]/Niter;
	meanCom[2] = meanCom[2] + W[iter*Ndim+2]/Niter;
    }
    mccc_wiener_generate(temparr, 3.1, &windex, &err);
    meanThr[0] = wienarr->wiener[0] + (wienarr->wiener[windex*Ndim+0] - wienarr->wiener[0])*(3.0-2.0)/(3.1-2.0);
    meanThr[1] = wienarr->wiener[1] + (wienarr->wiener[windex*Ndim+1] - wienarr->wiener[1])*(3.0-2.0)/(3.1-2.0);
    meanThr[2] = wienarr->wiener[2] + (wienarr->wiener[windex*Ndim+2] - wienarr->wiener[2])*(3.0-2.0)/(3.1-2.0);

    varCom[0] = 0;
    varCom[1] = 0;
    varCom[2] = 0;
    for (iter = 0; iter < Niter; iter = iter+1){
	varCom[0] = varCom[0] + pow(W[iter*Ndim+0]-meanCom[0],2)/(Niter-1);
	varCom[1] = varCom[1] + pow(W[iter*Ndim+1]-meanCom[1],2)/(Niter-1);
	varCom[2] = varCom[2] + pow(W[iter*Ndim+2]-meanCom[2],2)/(Niter-1);
    }
    varThr[0] = (3.0-2.0)*(3.1-3.0)/(3.1-2.0);
    varThr[1] = (3.0-2.0)*(3.1-3.0)/(3.1-2.0);
    varThr[2] = (3.0-2.0)*(3.1-3.0)/(3.1-2.0);

    printf("\n");
    printf("- When W(3.1) exists:\n");
    printf("%g, %g",meanCom[0],varCom[0]);
    printf("  (%g, %g)\n",meanThr[0],varThr[0]);
    printf("%g, %g",meanCom[1],varCom[1]);
    printf("  (%g, %g)\n",meanThr[1],varThr[1]);
    printf("%g, %g",meanCom[2],varCom[2]);
    printf("  (%g, %g)\n",meanThr[2],varThr[2]);
    printf("\n");

    printf("Wiener test complete.\n");
    printf("\n");

}

int main(int argc, char** argv) {
    
    mccc_coefs_test();
    mccc_push_test();
    //mccc_special_test();
    //mccc_wiener_test();
    return 0;
}
