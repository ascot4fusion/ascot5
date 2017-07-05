/**
 * @file ascot4_interface.c
 * @brief Functions for exchanging data with ASCOT4
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "ascot5.h"
#include "ascot4_interface.h"
#include "math.h"
#include "consts.h"
#include "phys_orbit.h"
#include "B_field.h"
#include "hdf5io/hdf5_helpers.h"
#include "hdf5io/hdf5_histogram.h"
#include "hdf5io/hdf5_particlestate.h"
#include "distributions.h"
#include "particle.h"

/**
 * @brief Write magnetic field component
 *
 * This function writes the magnetic field components in a hardcoded grid to
 * BR.out, Bphi.out and Bz.out for further processing with prepare_magn_bkg.m.
 *
 * @param Bdata pointer to magnetic field data struct
 */
void ascot4_write_B(B_field_data* Bdata) {
    FILE* fBr = fopen("BR.out", "w");
    FILE* fBphi = fopen("Bphi.out", "w");
    FILE* fBz = fopen("Bz.out", "w");

    real B[3];

    real r, z;

    int n_r = 100;
    real r_min = 3;
    real r_max = 9;
    real r_inc = (r_max-r_min)/(n_r-1);
    
    int n_z = 150;
    real z_min = 5;
    real z_max = -5;
    real z_inc = (z_max-z_min)/(n_z-1);
    
    int i, j;
    for(i = n_z-1; i >= 0; i--) {
        for(j = 0; j < n_r; j++) {
            r = r_min + j*r_inc;
            z = z_min + i*z_inc;
            B_field_eval_B(B, r, 0.0, z, Bdata);
            fprintf(fBr, "%lf ", B[0]);
            fprintf(fBphi, "%lf ", B[1]);
            fprintf(fBz, "%lf ", B[2]);
        }
        fprintf(fBr, "\n");
        fprintf(fBphi, "\n");
        fprintf(fBz, "\n");
    }

    fclose(fBr);
    fclose(fBphi);
    fclose(fBz);
}

/**
 * @brief Output particle data
 *
 * This function outputs the particles in the specified structure in
 * input.particles format for ASCOT4.
 *
 * @param p pointer to the particle group array
 * @param n number of particle groups
 * @param filename name of the output file
 * @param comment single line comment added to the output file
 */
void ascot4_write_particles(particle* p, int n, char* filename, 
                            char* comment) {
    int i;
    FILE* f = fopen(filename, "w");

    fprintf(f, " PARTICLE INITIAL DATA FOR ASCOT\n");
    fprintf(f, " 4 VERSION =====================\n");
    fprintf(f, "\n");
    fprintf(f, " 1  # Number of comment lines, max length 256 char, (defined in prtRead_lineLength)\n");
    fprintf(f, "%s\n", comment);
    fprintf(f, "\n");
    fprintf(f, "   %d # Number of particles (-1 means unknown number)\n", n*NSIMD);
    fprintf(f, "\n");
    fprintf(f, "14 # Number of different fields for each particle [10 first letters are significant]\n");
    fprintf(f, "Anum      - mass number of particle        (integer)\n");
    fprintf(f, "mass      - mass of the particle           (amu)\n");
    fprintf(f, "Znum      - charge number of particle      (integer)\n");
    fprintf(f, "charge    - charge of the particle         (elemental charge)\n");
    fprintf(f, "energy    - kinetic energy of particle     (eV)\n");
    fprintf(f, "phiprt    - toroidal angle of particle     (deg)\n");
    fprintf(f, "Rprt      - R of particle                  (m)\n");
    fprintf(f, "zprt      - z of particle                  (m)\n");
    fprintf(f, "vphi      - toroidal velocity of particle  (m/s)\n");
    fprintf(f, "vR        - radial velocity of particle    (m/s)\n");
    fprintf(f, "vz        - vertical velocity of particle  (m/s)\n");
    fprintf(f, "weight    - weight factor of particle      (particle/second)\n");
    fprintf(f, "id        - unique identifier of particle  (integer)\n");
    fprintf(f, "origin    - origin of the particle         ()\n");
    fprintf(f, "\n");

    for(i = 0; i < n; i++) {
        int Anum = (int) lround(p[i].mass / CONST_U);
        double mass = p[i].mass / CONST_U;

        int Znum = (int) lround(p[i].charge / CONST_E);
        double charge = p[i].charge / CONST_E;
        
	double gamma = phys_gammaprtv(p[i].v_r*p[i].v_r + p[i].v_phi*p[i].v_phi + p[i].v_z*p[i].v_z);
        double energy = (gamma - 1) * p[i].mass *CONST_C2  / CONST_E;
        
        double phiprt = p[i].phi;
        double Rprt = p[i].r;
        double zprt = p[i].z;

        double vphi = p[i].v_phi;
        double vR = p[i].v_r;
        double vz = p[i].v_z;

        double weight = 1.0;
        int id = p[i].id;
        int origin = 0;

        fprintf(f, "%d %le %d %le %le %le %le %le %le %le %le %le %d %d\n",
                Anum, mass, Znum, charge, energy, phiprt, Rprt, zprt,
                vphi, vR, vz, weight, id, origin);
    }
    fprintf(f, "#EOF\n");
}

/**
 * @brief Read particle data
 *
 * This function reads the particles from a file in the ASCOT4 input.particles
 * format into a particle group array, allocating memory in the process.
 * The function reads particles into complete NSIMD structs. Any remanining
 * particles that would partially fill a struct are discarded.
 *
 * @param p pointer to pointer to the particle group array
 * @param n number of particle groups
 * @param filename name of the input file
 */
void ascot4_read_particles(input_particle** p, int *n, char* filename) {
    int i, j;

    FILE* f = fopen(filename, "r");

    /* Skip first 3 lines */
    while(fgetc(f) != '\n');
    while(fgetc(f) != '\n');
    while(fgetc(f) != '\n');

    /* Read number of comment lines */
    int num_comment_lines;
    fscanf(f, "%d", &num_comment_lines);

    /* Skip to the end of line + num_comment_lines + empty line */
    while(fgetc(f) != '\n');
    for(i = 0; i < num_comment_lines; i++) {
        while(fgetc(f) != '\n');
    }
    while(fgetc(f) != '\n');

    /* Read number of particles */
    fscanf(f, "%d", n);

    *p = (input_particle*) malloc(*n * sizeof(input_particle));

    while(fgetc(f) != '\n');
    while(fgetc(f) != '\n');

    /* Read number of fields for each particle */
    int num_fields;
    fscanf(f, "%d", &num_fields);
    while(fgetc(f) != '\n');

    /* Find indexes for the fields we are interested in */
    int i_Rprt = -1;
    int i_phiprt = -1;
    int i_zprt = -1;
    int i_vR = -1;
    int i_vphi = -1;
    int i_vz = -1;
    int i_R = -1;
    int i_phi = -1;
    int i_z = -1;
    int i_energy = -1;
    int i_pitch = -1;
    int i_mass = -1;
    int i_charge = -1;
    int i_weight = -1;
    int i_id = -1;
    for(i = 0; i < num_fields; i++) {
        #define FIELD_NAME_LENGTH 10
        char field_name[FIELD_NAME_LENGTH+1];
        fgets(field_name, FIELD_NAME_LENGTH+1, f);
        if(!strcmp(field_name, "Rprt      "))
            i_Rprt = i;
        else if(!strcmp(field_name, "phiprt    "))
            i_phiprt = i;
        else if(!strcmp(field_name, "zprt      "))
            i_zprt = i;
        else if(!strcmp(field_name, "vR        "))
            i_vR = i;
        else if(!strcmp(field_name, "vphi      "))
            i_vphi = i;
        else if(!strcmp(field_name, "vz        "))
            i_vz = i;
        else if(!strcmp(field_name, "R         "))
            i_R = i;
        else if(!strcmp(field_name, "phi       "))
            i_phi = i;
        else if(!strcmp(field_name, "z         "))
            i_z = i;
        else if(!strcmp(field_name, "energy    "))
            i_energy = i;
        else if(!strcmp(field_name, "pitch     "))
            i_pitch = i;
        else if(!strcmp(field_name, "mass      "))
            i_mass = i;
        else if(!strcmp(field_name, "charge    "))
            i_charge = i;
        else if(!strcmp(field_name, "weight    "))
            i_weight = i;
        else if(!strcmp(field_name, "id        "))
            i_id = i;
        while(fgetc(f) != '\n');
    }
    while(fgetc(f) != '\n');
    
    /* Read all the fields in an array and fill the particle struct with
     * some of the fields */
    real* fields = (real*) malloc(num_fields * sizeof(real));
    if (i_Rprt == -1) {
        for(i = 0; i < *n; i++) {
            for(j = 0; j < num_fields; j++) {
                fscanf(f, "%lf", &fields[j]);
            }
            (*p)[i].p_gc.r = fields[i_R];
            (*p)[i].p_gc.phi = fields[i_phi] * math_pi / 180;
            (*p)[i].p_gc.z = fields[i_z];
            (*p)[i].p_gc.energy = fields[i_energy];
            (*p)[i].p_gc.pitch = fields[i_pitch];
            (*p)[i].p_gc.mass = fields[i_mass] * CONST_U;
            (*p)[i].p_gc.charge = fields[i_charge] * CONST_E;
            (*p)[i].p_gc.weight = fields[i_weight];
            (*p)[i].p_gc.id = (integer) fields[i_id];
            (*p)[i].p_gc.running = 1;
            (*p)[i].type = input_particle_type_gc;
        }
    }
    else {
        for(i = 0; i < *n; i++) {
            for(j = 0; j < num_fields; j++) {
                fscanf(f, "%lf", &fields[j]);
            }
            (*p)[i].p.r = fields[i_Rprt];
            (*p)[i].p.phi = fields[i_phiprt] * math_pi / 180;
            (*p)[i].p.z = fields[i_zprt];
            (*p)[i].p.v_r = fields[i_vR];
            (*p)[i].p.v_phi = fields[i_vphi];
            (*p)[i].p.v_z = fields[i_vz];
            (*p)[i].p.mass = fields[i_mass] * CONST_U;
            (*p)[i].p.charge = fields[i_charge] * CONST_E;
            (*p)[i].p.weight = fields[i_weight];
            (*p)[i].p.id = (integer) fields[i_id];
            (*p)[i].p.running = 1;
            (*p)[i].type = input_particle_type_p;
        }
    }
    free(fields);
    
    fclose(f);
}

void ascot4_write_dist_rzvv(dist_rzvv_offload_data* dist, real* hist,
        char* filename) {
    int abscissa_dim = 6;
    int ordinate_length = 1;

    /* transpose the histogram data for ascot4 rzVDist ordinate */
    double* ordinate = (double*) malloc(dist->n_r * dist->n_z * dist->n_vpara
                                        * dist->n_vperp * sizeof(double));
    int i, j, k, l;
    for(i = 0; i < dist->n_r; i++) {
        for(j = 0; j < dist->n_z; j++) {
            for(k = 0; k < dist->n_vpara; k++) {
                for(l = 0; l < dist->n_vperp; l++) {
                    ordinate[  l * (dist->n_vpara * dist->n_z * dist->n_r)
                             + k * (dist->n_z * dist->n_r)
                             + j * (dist->n_r)
                             + i] =
                    hist[  i * (dist->n_z * dist->n_vpara * dist->n_vperp)
                         + j * (dist->n_vpara * dist->n_vperp)
                         + k * (dist->n_vperp)
                         + l];
                }
            }
        }
    }

    int abscissa_n_slots[6];
    abscissa_n_slots[0] = dist->n_r;
    abscissa_n_slots[1] = dist->n_z;
    abscissa_n_slots[2] = dist->n_vpara;
    abscissa_n_slots[3] = dist->n_vperp;
    abscissa_n_slots[5] = 1;
    abscissa_n_slots[4] = 1;

    double abscissa_min[6];
    abscissa_min[0] = dist->min_r;
    abscissa_min[1] = dist->min_z;
    abscissa_min[2] = dist->min_vpara;
    abscissa_min[3] = dist->min_vperp;
    abscissa_min[4] = 0;
    abscissa_min[5] = 0.5;

    double abscissa_max[6];
    abscissa_max[0] = dist->max_r;
    abscissa_max[1] = dist->max_z;
    abscissa_max[2] = dist->max_vpara;
    abscissa_max[3] = dist->max_vperp;
    abscissa_max[4] = 100;
    abscissa_max[5] = 1.5;

    char* abscissa_names[] = { "R", "z", "vpa", "vpe", "time", "species" };
    char* abscissa_units[] = { "m", "m", "m/s", "m/s", "s", "" };
    char* ordinate_names[] = { "density" };
    char* ordinate_units[] = { "s^2/m^5" };

    int retval;
    retval =  hdf5_histogram_write_uniform_double(
		      filename, 
		      "rzVDist",
		      abscissa_dim,
		      ordinate_length,
		      abscissa_n_slots,
		      abscissa_min,
		      abscissa_max,
		      abscissa_units,
		      abscissa_names,
		      ordinate_units,
		      ordinate_names,
              ordinate
			         );
}

void ascot4_write_inistate(int n, input_particle* p, char* filename) {
    hid_t f = hdf5_create(filename);
    hdf5_particlestate_write(f, "inistate", n, p);
}

void ascot4_write_endstate(int n, input_particle* p, char* filename) {
    hid_t f = hdf5_open(filename);
    hdf5_particlestate_write(f, "endstate", n, p);
}

