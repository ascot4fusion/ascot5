
///////////////////////////////////////////////////////////////////
//
//      For a brief summary/tutorial,
//      please see: https://www.youtube.com/watch?v=xvFZjo5PgG0
//
/////////////////////////////////////////////////////////////////////


// Here you have the declarations of the functions and subroutines of RFOF. 
// Also, the types needed when calling these functions are declared here.

// NOTE: the names of the fortran functions and subroutines are now different. 
//     1. All upper case letters have been converted to lower case letters.
//     2. A prefix of the type "__modulename_MOD_" has been added. 
//
//     Example:
//         In RFOF_kick.f90 we have "tau_RF_phase_integral" but we write now in C
//         "__rfof_kick_MOD_tau_rf_phase_integral" instead of 
//         "tau_RF_phase_integral"
//
//     note (inside a note): please note that also the name of the module 
//                           in the prefix is written in lower case 
//                           (__rfof_kick_MOD and NOT __RFOF_kick_MOD)

// ROPLEEMAT
//
// RFOF_types puuttuu (?)
// resonance_memoryssa pitää ottaa muistin allokointi huomioon


#include <complex.h>
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <stdbool.h>
#include <string.h>
//#include <mpi.h>        
//#include <H5Part.h>
//#include <Dierckx.h>



///////////////////////////////////////////////////////////////////////             structs
/////////////////////////////////////////////////////////////////////

// From rfof/itm/euitm_schemas.f90
typedef struct {
    char** parameters;
    char** default_param;
    char** schema;
} type_param;

struct type_whatref {
    char** user;  // /whatref/user - Name of the user if private data, public if public ITM database.
    char** machine;  // /whatref/machine - Name of the device
    int shot;  // /whatref/shot - Shot number
    int run;  // /whatref/run - Run number
    int occurrence;  // /whatref/occurrence - Occurrence number of the CPO in the reference entry
};

// type_putinfo struct
struct type_putinfo {
    char** putmethod;  // /putinfo/putmethod - Storage method for this data
    char** putaccess;  // /putinfo/putaccess - Instructions to access the data using this method
    char** putlocation;  // /putinfo/putlocation - Name of this data under this method
    char** rights;  // /putinfo/rights - Access rights to this data
};

struct type_datainfo {
    char** dataprovider;  // /datainfo/dataprovider - Name of the actual data provider (the person who filled the data)
    char** putdate;  // /datainfo/putdate - Date at which the data has been put in the DB
    char** source;  // /datainfo/source - Exact reference of the data source (e.g. original reference in the native machine data base)
    char** comment;  // /datainfo/comment - Any additional comment
    int cocos;  // /datainfo/cocos - COordinates COnventionS followed by this CPO
    int id;  // /datainfo/id - CPO id for checking its provenance in the workflow
    int isref;  // /datainfo/isref - 1 if the data can be found in the present data base entry; 2 if the data can be found in a parent data base
    struct type_whatref whatref;  // /datainfo/whatref - 
    struct type_putinfo putinfo;  // /datainfo/putinfo - 
};

struct type_identifier {
    char** id;  // /id - Short string identifier
    int flag;  // /flag - Integer identifier
    char** description;  // /description - Verbose description of identifier
};

struct type_enum_instance {
    struct type_identifier type;  // /type - Identify the type of the object or process.
    char** name;  // /name - The name of the object or process. Here the object should be an instans of the type specified in the
    int index;  // /index - Index the separating objects or processes with the same name.
};

struct type_composition {
    double* amn;  // /composition/amn - Atomic mass number (lumped ions are allowed); Vector (nion)
    double* zn;  // /composition/zn - Nuclear charge (lumped ions are allowed); Vector (nion)
    double* zion;  // /composition/zion - Ion charge (of the dominant ionisation state; lumped ions are allowed); Vector (nion)
    int* imp_flag;  // /composition/imp_flag - Multiple charge state calculation flag : 0-Only one charge state is considered; 1-Multiple charge states are considered
    char** label;  // /composition/label - Label for the ions - note the charge state is not included; String Vector (nion)
};

struct type_compositions_type {
    struct type_nuclei** nuclei;  // /nuclei(i) - Array of nuclei considered.
    struct type_ions** ions;  // /ions(i) - Array of main plasma ions.
    struct type_impurities** impurities;  // /impurities(i) - Array of impurities.
    struct type_composition_neutralscomp** neutralscomp;  // /neutralscomp(i) - Array of neutrals.
    struct type_edgespecies** edgespecies;  // /edgespecies(i) - Array of edge species.
    struct type_identifier signature;  // /signature - Identifier for species choices. The goal of this is to uniquely capture the species blocks so that i
};

typedef struct {
    int* range;
    int* ind;
} type_complexgrid_indexlist;
/*
typedef struct {
    int* cls;
    type_complexgrid_indexlist** indset;
    int** ind;
} type_complexgrid_objectlist;
*/
/*
typedef struct {
    char* id;
    type_complexgrid_objectlist** list;
} type_complexgrid_subgrid;
*/
/*
typedef struct {
    double complex scalar;
    double complex* vector;
    double complex (*matrix)[3][3];
} type_complexgrid_scalar_cplx;
*/

typedef struct {
    int griduid;
    int subgrid;
    char** label;
    struct type_complexgrid_scalar_cplx** comp;
    int* align;
    char** alignid;
    int basis;
} type_complexgrid_vector;

typedef struct {
    int griduid;
    int subgrid;
    complex double* scalar;
    complex double (*vector)[2];
    complex double (*matrix)[3][3];
} type_complexgrid_scalar_cplx;

/*
typedef struct {
    type_complexgrid_scalar_cplx* measure;
    type_complexgrid_scalar_cplx** g11;
    type_complexgrid_scalar_cplx** g12;
    type_complexgrid_scalar_cplx** g13;
    type_complexgrid_scalar_cplx** g22;
    type_complexgrid_scalar_cplx** g23;
    type_complexgrid_scalar_cplx** g33;
    type_complexgrid_scalar_cplx** jacobian;
} type_complexgrid_metric;
*/
typedef struct {
    int geotype;
    char** geotypeid;
    int* coordtype;
    struct type_complexgrid_scalar_cplx** geo_matrix;
    struct type_complexgrid_scalar_cplx* measure;
} type_complexgrid_geo_global;

typedef struct {
    double r;
    double z;
} type_rz0D;

typedef struct {
    double r0;
    double b0;
} type_b0r0;

typedef struct {
    char** name;
    char** type;
    int* f_assumption;
    int code_type;
    double frequency;
    int* ntor;
    double power_tot;
    double* p_frac_ntor;
    double pow_e;
    double* pow_i;
    double complex (*pow_z)[3];
    double pow_fe;
    double* pow_fi;
    double complex (*pow_fz)[3];
    double* pow_ntor_e;
    double complex (*pow_ntor_i)[3];
    double complex (*pow_ntor_z)[3][3];
    double* pow_ntor_fe;
    double complex (*pow_ntor_fi)[3];
    double complex (*pow_ntor_fz)[3][3];
    double cur_tor;
    double* cur_tor_ntor;
    struct type_rz0D mag_axis;
    struct type_b0r0 toroid_field;
} type_waves_global_param;

typedef struct {
    double* rho_tor;
    double* rho_tor_norm;
    double* psi;
    double* volume;
    double* area;
} type_waves_grid_1d;
/*
typedef struct {
    char** name;
    char** type;
    int* f_assumption;
    int code_type;
    double frequency;
    int* ntor;
    double power_tot;
    double* p_frac_ntor;
    double pow_e;
    double* pow_i;
    double complex (*pow_z)[3];
    double pow_fe;
    double* pow_fi;
    double complex (*pow_fz)[3];
    double* pow_ntor_e;
    double complex (*pow_ntor_i)[3];
    double complex (*pow_ntor_z)[3][3];
    double* pow_ntor_fe;
    double complex (*pow_ntor_fi)[3];
    double complex (*pow_ntor_fz)[3][3];
    double cur_tor;
    double* cur_tor_ntor;
    type_rz0D mag_axis;
    type_b0r0 toroid_field;
} type_waves_global_param;
*/
typedef struct {
    double* powd_tot;
    double* powd_e;
    double (*powd_i)[2];
    double (*powd_z)[3][4];
    double* powd_fe;
    double (*powd_fi)[3][4];
    double (*powd_fz)[3][4][5];
    double (*powd_ntor)[2];
    double (*powd_ntor_e)[2];
    double (*powd_ntor_i)[3][4];
    double (*powd_ntor_z)[3][4][5];
    double (*powd_ntor_fe)[2];
    double (*powd_ntor_fi)[3][4];
    double (*powd_ntor_fz)[3][4][5];
    double* curd_tor;
    double (*curd_torntor)[2];
    double* pow_tot;
    double* pow_e;
    double (*pow_i)[2];
    double (*pow_z)[3][4];
    double* pow_fe;
    double (*pow_fi)[3][4];
    double (*pow_fz)[3][4][5];
    double (*pow_ntor)[2];
    double (*pow_ntor_e)[2];
    double (*pow_ntor_i)[3][4];
    double (*pow_ntor_z)[3][4][5];
    double (*pow_ntor_fe)[2];
    double (*pow_ntor_fi)[3][4];
    double (*pow_ntor_fz)[3][4][5];
    double* curd_par;
    double (*curd_parntor)[2];
    double* cur_tor;
    double (*cur_tor_ntor)[2];
    double (*e_plus_ave)[2];
    double (*e_minus_ave)[2];
    double (*e_para_ave)[2];
    double (*k_perp_ave)[2];
} type_waves_profiles_1d;

typedef struct {
    double (*powd_tot)[2];
    double (*powd_e)[2];
    double (*powd_i)[3][4];
    double (*powd_z)[4][5][6];
    double (*powd_fe)[2];
    double (*powd_fi)[3][4];
    double (*powd_fz)[4][5][6];
    double (*powd_ntor)[3][4];
    double (*powd_ntor_e)[3][4];
    double (*powd_ntor_i)[4][5];
    double (*powd_ntor_z)[4][5][6];
    double (*powd_ntor_fe)[3][4];
    double (*powd_ntor_fi)[4][5];
    double (*powd_ntor_fz)[4][5][6];
    double (*powd_iharm)[5][6][7][8];
} type_waves_profiles_2d;

typedef struct {
    int angl_type;
    double (*th2th_pol)[2];
} type_theta_info;

typedef struct {
    double *r;
    double *z;
    double *phi;
    double *psi;
    double *theta;
} type_waves_rtposition;

typedef struct {
    double *kr;
    double *kz;
    double *kphi;
    double *npar;
    double *nperp;
    double *ntor;
    int var_ntor;
} type_waves_rtwavevector;


typedef struct {
    double *epol_p_re;
    double *epol_p_im;
    double *epol_m_re;
    double *epol_m_im;
    double *epol_par_re;
    double *epol_par_im;
} type_polarization;


typedef struct {
    double *phi_perp;
    double *phi_par;
    double *power_e;
    double **power_i;
} type_powerflow;

typedef struct {
    int **boundary;
    int ***neighbour;
    double ****geo;
    double **measure;
} type_objects;

typedef struct {
    int *geotype;
    char **geotypeid;
    int **coordtype;
    struct type_objects **objects;
    int *xpoints;
} type_complexgrid_space;

typedef struct {
    char **id;
    struct type_complexgrid_objectlist **list;
} type_complexgrid_subgrid;

/*
typedef struct {
    int *range;
    int *ind;
} type_complexgrid_indexlist;
*/
typedef struct {
    int *cls;
    struct type_complexgrid_indexlist **indset;
    int **ind;
} type_complexgrid_objectlist;

typedef struct {
    int griduid;
    int subgrid;
    double *scalar;
    double **vector;
    double ***matrix;
} type_complexgrid_scalar;
/*
typedef struct {
    int geotype;
    char **geotypeid;
    int *coordtype;
    struct type_complexgrid_scalar **geo_matrix;
    struct type_complexgrid_scalar **measure;
} type_complexgrid_geo_global;
*/
typedef struct {
    struct type_complexgrid_scalar **measure;
    struct type_complexgrid_scalar **g11;
    struct type_complexgrid_scalar **g12;
    struct type_complexgrid_scalar **g13;
    struct type_complexgrid_scalar **g22;
    struct type_complexgrid_scalar **g23;
    struct type_complexgrid_scalar **g33;
    struct type_complexgrid_scalar **jacobian;
} type_complexgrid_metric;
/*
typedef struct {
    int griduid;
    char **label;
    struct type_complexgrid_scalar **comp;
    int *align;
    char **alignid;
    int basis;
} type_complexgrid_vector;
*/
typedef struct {
    int uid;
    char **id;
    struct type_complexgrid_space **spaces;
    struct type_complexgrid_subgrid **subgrids;
    struct type_complexgrid_metric metric;
    struct type_complexgrid_geo_global **geo;
    struct type_complexgrid_vector **bases;
} type_complexgrid;

typedef struct {
    int* mpol;
    double*** e_plus;
    double*** e_plus_ph;
    double*** e_minus;
    double*** e_minus_ph;
    double*** e_norm;
    double*** e_norm_ph;
    double*** e_binorm;
    double*** e_binorm_ph;
    double*** e_para;
    double*** e_para_ph;
    double*** b_norm;
    double*** b_norm_ph;
    double*** b_binorm;
    double*** b_binorm_ph;
    double*** b_para;
    double*** b_para_ph;
    double*** k_perp;
} type_pol_decomp;

typedef struct {
    double*** e_plus;
    double*** e_plus_ph;
    double*** e_minus;
    double*** e_minus_ph;
    int*** e_norm;
    double*** enorm_ph;
    double*** e_binorm;
    double*** e_binorm_ph;
    double*** e_para;
    double*** e_para_ph;
    double*** b_norm;
    double*** b_norm_ph;
    double*** b_binorm;
    double*** b_binorm_ph;
    double*** b_para;
    double*** b_para_ph;
    double*** k_perp;
} type_local;

typedef struct {
    struct type_complexgrid_scalar_cplx e_plus;
    struct type_complexgrid_scalar_cplx e_minus;
    struct type_complexgrid_scalar_cplx e_para;
    struct type_complexgrid_scalar_cplx e_norm;
    struct type_complexgrid_scalar_cplx e_binorm;
    struct type_complexgrid_scalar_cplx b_norm;
    struct type_complexgrid_scalar_cplx b_binorm;
    struct type_complexgrid_scalar_cplx b_para;
    struct type_complexgrid_scalar_cplx k_perp;
} type_e_components;



typedef struct {
    struct type_complexgrid grid;
    struct type_e_components e_components;
    struct type_pol_decomp pol_decomp;
    struct type_local local;
} type_fullwave;

typedef struct {
    int npoints;
    double power;
    double *dnpar;
    double *length;
    struct type_waves_rtposition position;
    struct type_waves_rtwavevector wavevector;
    struct type_polarization polarization;
    struct type_powerflow powerflow;
} type_beamtracing;

struct type_coherentwave {
    struct type_enum_instance wave_id;  // /waves/coherentwave(i)/wave_id - List of identifiers for the coherent-wave, in terms of the type and name of the antenna driving the 
    struct type_composition composition;  // /waves/coherentwave(i)/composition - 
    struct type_compositions_type compositions;  // /waves/coherentwave(i)/compositions - Contains detailed information on the plasma composition (main ions, impurities, neutrals, edge speci
    struct type_waves_global_param global_param;  // /waves/coherentwave(i)/global_param - Global wave deposition parameters
    struct type_waves_grid_1d grid_1d;  // /waves/coherentwave(i)/grid_1d - Grid points for 1D profiles.
    struct type_waves_grid_2d grid_2d;  // /waves/coherentwave(i)/grid_2d - Grid points for 2D profiles and for full wave solutions.
    struct type_waves_profiles_1d profiles_1d;  // /waves/coherentwave(i)/profiles_1d - 1D radial profiles
    struct type_waves_profiles_2d profiles_2d;  // /waves/coherentwave(i)/profiles_2d - 2D profiles in poloidal cross-section
    struct type_beamtracing* beamtracing;  // /waves/coherentwave(i)/beamtracing(i) - Beam-tracing or ray-tracing solver. Vector(nbeams). Time-dependent
    struct type_fullwave fullwave;  // /waves/coherentwave(i)/fullwave - Solution by full wave code
    struct type_codeparam codeparam;  // /waves/coherentwave(i)/codeparam - Code parameters of physics code, i.e. codes calculating a wave field.
};

struct type_codeparam {
    char** codename;  // /codeparam/codename - Name of the code
    char** codeversion;  // /codeparam/codeversion - Version of the code (as in the ITM repository)
    char** parameters;  // /codeparam/parameters - List of the code specific parameters, string expected to be in XML format.
    char** output_diag;  // /codeparam/output_diag - List of the code specific diagnostic/output, string expected to be in XML format.
    int output_flag;  // /codeparam/output_flag - Output flag : 0 means the run is successful, other values meaning some difficulty has been encountered
};

struct type_waves {
    struct type_datainfo datainfo;  // /waves/datainfo - 
    struct type_coherentwave* coherentwave;  // /waves/coherentwave(i) - Wave description for each frequency. Time-dependent. Structure array(nfreq)
    struct type_codeparam codeparam;  // /waves/codeparam - Code parameters of datajoiners, i.e. codes that merge the wave field of two or more physics codes.
    double time;  // /waves/time - Time [s]; Time-dependent; Scalar
};

//From RFOF_input.f90
typedef struct {
    // 1. rfof_core_param
    // 1.1 assumptions
    _Bool assume_static_resonance_position_during_RF_kick;
    _Bool use_drift_velocity_in_doppler_shift;
    _Bool use_parallel_velocity_in_doppler_shift;
    _Bool assume_zero_larmor_radius_in_KPERPxRHO;
    _Bool assume_kpar_is_nphi_over_R;
    _Bool assume_zero_order_FLR_for_Pphi;
    double width_of_rf_resonance_layer;

    // 1.2 bounding_box
    double plasma_boundingbox_Rmin;
    double plasma_boundingbox_Rmax;
    double plasma_boundingbox_Zmin;
    double plasma_boundingbox_Zmax;

    // 1.3 resonance_memory
    int n_store_times_in_resonance_memory;

    // 1.4 IO_control
    double start_time_event_output;
    _Bool output__2D_RZ_out;
    int NRedges_2DgridRZ;
    int NZedges_2DgridRZ;
    _Bool output__Orbit;
    int MAX_number_of_points_stored_in_the_Orbit;
    _Bool output__rf_kicks;
    int MAX_number_of_points_stored_in_rf_kick;
    _Bool output__resonace_predictions;
    int MAX_number_of_points_stored_in_resonance_memory;
    _Bool output__efield_normalization;
    int MAX_number_of_points_stored_in_the_efield_normalization;

    // 1.5 Quasilinear model parameters
    double MAX_relative_energy_kick;

    // 2. rfof_plasma_param
    int n_species;
    double* amn;
    double* zn;
    double* zion;
    double* densities;
    double ion_temperature;

    // 3. rfof_wave_param
    int select_wave_from;

    // 3.1 rfof_wave_param/wave_modifiers
    _Bool append_kperp;
    _Bool map_efield_from_pol_decomp_to_local;

    // 3.2 rfof_wave_param/parametric_wave_field
    int nfreq;
    int nnphi;
    double* RFpower;
    double* EfieldNormalisation;
    double* ratioEPlusOverEMinus;
    double* freq;
    int* nphi;
    double* kperp;
    double* verticalCentre;
    double* verticalWidth;
    int nGridHorisontal;
    int nGridVertical;
    char dielectric_response_model[125];
    char filename_lion_fields[125];

    // 3.3 ascii_itm_wave
    char filename_ascii_itm_wave[125];

    // 4. rfof_wrapper_param
    // 4.1 rfof_wrapper_param/time_stepping
    int NtimeSteps;
    double dt;
    int nStoreOutTimes;

    // 4.2 rfof_wrapper_param/magnetic_field
    double R0;
    double aminor;
    double B0;
    double q;

    // 4.3 rfof_wrapper_param/input_marker
    int species_index;
    double weight;
    double R;
    double z;
    double phi;
    double charge;
    double mass;
    double E;
    double xi;

} RFOF_input_param;



//From RFOF_resonance_memory.f90
struct resonance_memory {
    int Number_points_in_memory;
    double* time;
    double* NACC;
    double* r;
    double* phi;
    double* z;
    double* omega_res;
    double* omega_c;
    int* nharm;
    bool previous_resonance_exists;
    double time_of_last_resonance_crossing;
    int sign_omega_dot_at_last_resonance_crossing;
};



//This should be initialised somewhere?
//From RFOF_constants.f90
typedef struct {
    double PI;
    double C;    // speed of light [m/s]
    double ME;   // electron mass [kg]
    double MP;   // proton mass [kg]
    double MN;   // neutron mass [kg]
    double MD;   // deuteron mass [kg]
    double MT;   // triton mass [kg]
    double MA;   // alpha mass [kg]
    double AMU;  // atomic mass unit [kg]
    double EV;   // electron volt (eV)
    double QE;   // elementary charge [coulomb]
    double MU0;  // vacuum permeability
    double EPS0;
    double AVOGR;
    double KBOLT; // Boltzmann's constant
} RFOF_CONST_TYPE;



//From dierckx_wrapper.f90 (located in rfof/dierckx/ and not rfof/src/ like the other files)
typedef struct {
    double min_x;
    double max_x;
    double min_y;
    double max_y;
    double fp;
    double smoothing_factor;
    int iopt;
    int nx;
    int ny;
    int kx;
    int ky;
    int nxest;
    int nyest;
    int* iwrk;
    double* wrk;
    double* x;
    double* y;
    double** z;
    double* tx;
    double* ty;
    double* c;
} spline_fcn_2d_rect;



//From RFOF_parameters.f90
typedef struct {
    double R0;
    double aminor;
    double B0;
    int nfreq;
    int nnphi;
    double RFpower;
    double EfieldNormalisation;
    double freq;
    int nphi;
    double kperp;
    char filename_lion_fields[125];  // Assuming 124 characters plus null terminator
} rf_wave_wrapper_parameters;



//From RFOF_waves.f90
typedef struct {
    double normalisation_factor;
    double power;
    double frequency;
    int n_phi;
    struct spline_fcn_2d_rect Eplus;
    struct spline_fcn_2d_rect Eplus_phase;
    struct spline_fcn_2d_rect Eminus;
    struct spline_fcn_2d_rect Eminus_phase;
    struct spline_fcn_2d_rect k_par;
    struct spline_fcn_2d_rect k_perp;
    struct spline_fcn_2d_rect k_rho;
    struct spline_fcn_2d_rect k_dia;
    double bounding_box_x1_min;
    double bounding_box_x1_max;
    double bounding_box_x2_min;
    double bounding_box_x2_max;
} rf_wave_mode_global;


typedef struct {
    struct rf_wave_mode_global *mode;
} rf_wave_coherent_global;

typedef struct {
    int index_present_wave;
    int index_present_mode;
    int max_number_modes;
    int coord_name_x1;
    int coord_name_x2;
    struct rf_wave_coherent_global *wave;
} rf_wave_global;

typedef struct {
    int jfreq;
    int jnphi;
    complex double Erho;
    complex double Edia;
    complex double Epar;
    complex double Eplus;
    complex double Eminus;
    double krho;
    double kdia;
    double kpar;
    double kperp;
    double omega;
    double nphi;
} rf_wave_local;

typedef struct {
    int jfreq;
    int jnphi;
    double krho;
    double kdia;
    double kpar;
    double kperp;
    double omega;
    double nphi;
} rf_wave_numbers_local;



//From RFOF_magnetic_field
typedef struct {
    double R;
    double phi;
    double z;
    double psi;
    double theta;
    double rho_tor;
    double Bmod;
    double F;
    double psi_Estatic;
    double dBmod_dpsi;
    double dF_dpsi;
    double dBmod_dtheta;
    double dF_dtheta;
} magnetic_field_local;



// From RFOF_markers.f90
typedef struct {
    int *Id;                  // Number to identify particle
    double *weight;           // Number of real particles represented by the marker
    double *charge;           // Charge (SI units; everything is in SI)
    double *mass;             // Mass (SI units; everything is in SI)
    double *R;                // Major radius
    double *phi;              // Toroidal angle
    double *z;                // Vertical position
    double *psi;              // Poloidal flux function (following the ITM and ITER conventions; COCOS 13/11; not divided by 2π)
    double *rho_tor;          // The square root of the toroidal flux function normalised according to the ITM convention
    double *theta;            // Poloidal angle; generic angle that could be geometrical, or straight field line, etc.
    double *energy;           // Total energy, i.e. sum of kinetic and potential energy
    double *energy_kinetic;   // Kinetic energy
    double *velocity;         // Magnitude of the velocity
    double *magneticMoment;   // Magnetic moment
    double *Pphi;             // The canonical toroidal angular momentum Pφ = Zeψ/2π + mv‖F/B
    double *vpar;             // Parallel velocity
    double *vperp;            // Perpendicular velocity
    double *omega_gyro;       // The angular gyro (Larmor) frequency
    double *tauBounce;        // Estimate for the time of a poloidal orbit (bounce time for trapped / transit time for passing particles)
    double *vDrift;           // The magnitude of the drift velocity
    double *vDriftRho;        // The component of the drift velocity w.r.t. the radial direction
    double *vDriftDia;        // The component of the drift velocity w.r.t. the diamagnetic direction (diamagnetic=poloidal)
    double *d_vpar_d_rho;     // Partial derivative of the parallel velocity w.r.t. the radial direction
    double *d_vpar_d_dia;     // Partial derivative of the parallel velocity w.r.t. the diamagnetic direction (diamagnetic=poloidal)
    double *d_vperp_d_rho;    // Partial derivative of the perpendicular velocity w.r.t. the radial direction
    double *d_vperp_d_dia;    // Partial derivative of the perpendicular velocity w.r.t. the diamagnetic direction (diamagnetic=poloidal)
    double *d_vDriftRho_d_rho;// Partial derivative of the radial component of the drift velocity w.r.t. the radial direction
    double *d_vDriftRho_d_dia;// Partial derivative of the radial component of the drift velocity w.r.t. the diamagnetic direction (diamagnetic=poloidal)
    double *d_vDriftDia_d_rho;// Partial derivative of the diamagnetic component of the drift velocity w.r.t. the radial direction (diamagnetic=poloidal)
    double *d_vDriftDia_d_dia;// Factor by which the orbit integration has been accelerated (number of of real bounce times per calculated bounce in the orbit following.
    double *time_acceleration; // Time-acceleration factor
    int isOrbitTimeAccelerated; // When using time-acceleration one can choose to supply an orbit time in two ways
} particle;


//This is probably not even needed actually.
typedef struct {
    int Id;
    double weight;
    double charge;
    double mass;
    double R;
    double phi;
    double z;
    double psi;
    double rho_tor;
    double theta;
    double energy;
    double energy_kinetic;
    double velocity;
    double magneticMoment;
    double Pphi;
    double vpar;
    double vperp;
    double omega_gyro;
    double tauBounce;
    double vDrift;
    double vDriftRho;
    double vDriftDia;
    double d_vpar_d_rho;
    double d_vpar_d_dia;
    double d_vperp_d_rho;
    double d_vperp_d_dia;
    double d_vDriftRho_d_rho;
    double d_vDriftRho_d_dia;
    double d_vDriftDia_d_rho;
    double d_vDriftDia_d_dia;
    double time_acceleration;
    int isOrbitTimeAccelerated; // Using int to represent logical in C (0 for false, 1 for true)
} particle_static;



//From RFOF_resonance_memory
typedef struct {
    int Number_points_in_memory;
    double *time;
    double *NACC;
    double *r;
    double *phi;
    double *z;
    double *omega_res;
    double *omega_c;
    int *nharm;
    int previous_resonance_exists;  // Using int to represent logical in C (0 for false, 1 for true)
    double time_of_last_resonance_crossing;
    int sign_omega_dot_at_last_resonance_crossing;
} resonance_memory;



//In RFOF_resonance_condition.f90 there are no types so there is nothing to be declared here.



//In RFOF_random_numbers.f90 there are no types so there is nothing to be declared here.



//From RFOF_diagnostics.f90
typedef struct {
    int diagnostic_active;  // Using int to represent logical in C (0 for false, 1 for true)
    long long int *kick_counter;
    double time_start_sum_diagnostics;
    double time_end_sum_diagnostics;
    double energy_to_markers;
    double toroidal_momentum_to_markers;
    double *energy_drift_from_mode;
    double *energy_to_markers_from_mode;
    double *energy_square_to_markers_from_mode;
    double *toroidal_momentum_to_markers_from_mode;
    double *sum_weight_at_resonance_with_mode;
    double *R_2DgridRZ;
    double *Z_2DgridRZ;
    double *energy_2DgridRZ;
    double *momentum_2DgridRZ;
} RFOF_cumlative_diagnostics;




//From RFOF_kick.f90 itself we get finally this one last type (struct)
typedef struct {
    double d_mu_dI;
    double d_vpar_dI;
    double d_psi_dI;
    double d_theta_dI;
} coeff_for_guiding_centre_kick;






/////////////////////////////////////////////////////////////////////
// FUNCTIONS AND SUBROUTINES
////////////////////////////////////////////////////////////////////

// These are given roughly in the order in which they are used, one module at a time.



// From RFOF_input_cpo.f90. 
//BOTH OF THESE WERE SUPPOSED TO BE FOUNF IN THE "RFOF_input module", but that only contained RFOF_parse_xml_input (with the 
// last "param" missing)!!

//  Read RFOF code-parameters from file and store in the type_param structure
void __rfof_input_cpo_MOD_fill_rfof_type_param(type_param* param, const char* fn_xml, const char* fn_xml_default, const char* fn_xsd, int* error_flag, int* mpirank);

// This function is used after the fill_rfof_type_param function. The resulting rfof_param is 
//then later used when calling different RFOF-constructors.
void __rfof_input_cpo_MOD_rfof_parse_xml_input_param(void* xml_node, RFOF_input_param* rfof_param, int* return_status);




// From RFOF_markers

//  "The logical, optional, intent(in) parameter isOrbitTimeAccelerated in the 
// Fortran subroutine is not directly translatable to C. In C, you can use a 
// pointer to an integer to simulate a similar behavior. If the pointer is not 
// NULL, you can consider it as true, otherwise, it can be considered as false."


void __rfof_markers_MOD_set_marker_pointers(struct particle* marker,
    int* Id, double* weight, double* R, double* phi, double* z,
    double* psi, double* rho_tor, double* theta,
    double* charge, double* mass, double* energy, double* energy_kinetic,
    double* velocity, double* magneticMoment, double* Pphi, double* vpar,
    double* vperp, double* omega_gyro, double* tauBounce,
    double* vDrift, double* vDriftRho, double* vDriftDia,
    double* d_vpar_d_rho, double* d_vpar_d_dia, double* d_vperp_d_rho, double* d_vperp_d_dia,
    double* d_vDriftRho_d_rho, double* d_vDriftRho_d_dia, double* d_vDriftDia_d_rho, double* d_vDriftDia_d_dia,
    double* time_acceleration, int isOrbitTimeAccelerated);



//From RFOF_resonance_memory.f90

/* 
 * Constructor of a matrix of the type resonance_memory.
 *
 * Allocates memory and initializes counters and state variables.
 *
 * @param mem      Resonance memory information to be initialized.
 * @param nWaves1  Size of matrix in the first dimension.
 * @param nWaves2  Size of matrix in the second dimension.
 * @param n_store_times_in_resonance_memory Number of time points to be stored in the resonance memory.
 */
void __rfof_resonance_memory_MOD_constructor_rf_resonance_memory_matrix(struct resonance_memory** mem, int nWaves1, int nWaves2, int n_store_times_in_resonance_memory);



void __rfof_resonance_memory_MOD_reset_rf_resonance_memory_matrix(struct resonance_memory** mem, int nWaves1, int nWaves2);



// From RFOF_waves_cpo.f90

void __rfof_waves_cpo_MOD_rf_wave_global_constructor(struct type_waves* itm_waves_in, struct rf_wave_global* RFglobal, struct RFOF_input_param rfof_param, int* use_compositions_ions, int* ierror, char** error_message, int mpirank);



// From RFOF_diagnostics.f90
void __rfof_diagnostics_MOD_contructor_RFOF_cumlative_diagnostics(struct RFOF_cumlative_diagnostics* diagno, int nfreq, int nnphi, int* error_flag, char** error_message);



// From RFOF_parameters.f90

//This should be called from the orbit code to set the global RFOF parameters. This routine will set the (R,Z) bounding box inside which the particles are expected to be located.
void __rfof_parameters_MOD_initialise_RFOF_parameters(double* Rmin_bound, double* Rmax_bound, double* Zmin_bound, double* Zmax_bound, char* filename_in, struct RFOF_input_param* rfof_param, int* error_flag);



// From RFOF_main.f90 

/*
The subroutine which kicks the markers in the ass thus changing their velocity. This routine is called after every time step of the orbit solver. This routine has two parts:
    1. It checks the resonance condition (with respect to all RF wave field wavemodes).
    2.
        2.1. If (close enough to) resonance, gives a kick.
        2.2. If not in resonance, check if particle has crossed a resonance without receiving a kick. 
            2.2.1 If yes, blow into the whizzle and release the hounds. (Hounds are yet to be implemented, so as a placeholder, an error flag is raised.) The routine returns an estimate when the crossing happened. You should then refine the time step and redo this step. In case the step is successful, then the \c RFOF_master will return the time at which the next resonance should occure. You may then use this time to decide the length of the folloing time step.

Note: If also the marker position is to be changed inside RFOF, then local_magnetic_field_interf... must be written to provide the magnetic field data to RFOF.
*/
void __rfof_main_MOD_rfof_master(double time, double time_previous_step, 
                 struct RFglobal* RFglobal, struct particle* marker, 
                 struct resonance_memory** resonanceMemory, 
                 struct RFOF_cumlative_diagnostics* diagno, int MPI_nod_Id, 
                 int* interaction_failed_particle_overshot_resonance, 
                 double* time_at_resonance, double* normalised_diffusion_rate, 
                 int* do_kick_marker, int* force_kick_marker, int* verbosity);



// From RFOF_Efield_update.f90

//To make sure that the RF power absorbed in RFOF the wave field strength has to be renormalised during the simulation. The is done using this routine. The call frequency to this routine is not trivial to determine. In principle it should be given by the time scale on which the absorption changes.
void __rfof_efield_update_MOD_update_efield_normalisation(double dt, struct rf_wave_global* RFglobal, int MPI_node_Id, int* ierr, double RFpower[][2], struct RFOF_cumlative_diagnostics* diagno);



// At the end of the program all allocated memory should be deallocated using 
// the routine \c RFOF_destructor in the module \c RFOF_main.
void __rfof_main_MOD_rfof_destructor(rf_wave_global* wave, resonance_memory** mem, RFOF_cumlative_diagnostics* diagno);



// From RFOF_kick.f90 

void __rfof_kick_MOD_quasilinear_rf_kick_steinbrecher_integrator(struct particle* marker, struct resonance_memory* mem, 
        struct magnetic_field_local* Blocal, struct rf_wave_local* RFlocal, struct RFOF_cumlative_diagnostics* diagno, 
        double diffusion_normalised, int MPI_node_Id);


//This is not in the librfof.so and not needed?? Originally declared inside the quasilinear_RF_steinbrecher_integrator (see above)
void __rfof_kick_MOD_MC_kick_steinbrecher_integrator(double dt, double* x, void (*calc_diffusion_coeff)(double, RFOF_state*, double*, int*),
                                    RFOF_state* state, double relative_tolerance_MonteCarlo, double* out_x_evaluations,
                                    double* out_diffusion_evaluations, int* errorFlag, int MPI_node_Id);

typedef void (*calc_diffusion_coeff_func)(double, RFOF_state*, double*, int*);                      

void __rfof_kick_MOD_wrap_rf_diffusion_and_move_marker(double Iperp, struct RFOF_state* state, double* diffusion, int* errorFlag);

double __rfof_kick_MOD_quasilinear_rf_diffusion_coeff(struct particle marker, struct rf_wave_local RFlocal, struct resonance_memory mem);

double __rfof_kick_MOD_max_rfkick_in_single_pass(struct particle marker, struct rf_wave_local wave, struct resonance_memory mem);

double __rfof_kick_MOD_tau_RF_phase_integral_acc(struct resonance_memory mem, struct particle marker);

double __rfof_kick_MOD_tau_rf_phase_integral(double time[], double omega_res[], int Nelements);

struct coeff_for_guiding_centre_kick __rfof_kick_MOD_coeff_rf_characteristic(struct particle marker, struct rf_wave_local RFlocal, struct magnetic_field_local Blocal, int nharm);

void __rfof_kick_MOD_wrap_move_marker_on_characteristic(double Iperp, struct RFOF_state* state, int* errorFlag);

void __rfof_kick_MOD_move_marker_on_characteristic(particle* marker, double dIperp, magnetic_field_local Blocal, coeff_for_guiding_centre_kick coeff);

double __rfof_kick_MOD_energy_drift(double Iperp[2], double diffusion[2], particle marker, rf_wave_local RFlocal, magnetic_field_local Blocal, int nharm);

double __rfof_kick_MOD_get_Bmod(particle marker, rf_wave_local RFlocal, magnetic_field_local Blocal, double rel_err);
