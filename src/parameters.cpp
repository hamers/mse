#include "parameters.h"
#include <math.h>

// Constants //
double CONST_G = 4.0*M_PI*M_PI; 
double CONST_G_P2 = CONST_G*CONST_G;
double CONST_G_P3 = CONST_G_P2*CONST_G;
double CONST_C_LIGHT = 63239.72638679138;
double CONST_C_LIGHT_P2 = CONST_C_LIGHT*CONST_C_LIGHT;
double CONST_C_LIGHT_P3 = CONST_C_LIGHT_P2*CONST_C_LIGHT;
double CONST_C_LIGHT_P4 = CONST_C_LIGHT_P3*CONST_C_LIGHT;
double CONST_C_LIGHT_P5 = CONST_C_LIGHT_P4*CONST_C_LIGHT;
double CONST_C_LIGHT_P6 = CONST_C_LIGHT_P5*CONST_C_LIGHT;
double CONST_MSUN = 1.0;
double CONST_R_SUN = 0.004649130343817401;
double CONST_L_SUN = 0.0002710404109745588;
double CONST_KM_PER_S = 0.210862;
double CONST_PER_PC3 = 1.14059e-16;
double CONST_MJUP = 0.000954248;

// Default parameters //
double relative_tolerance = 1.0e-10;
double absolute_tolerance_eccentricity_vectors = 1.0e-8;
double absolute_tolerance_spin_vectors = 1.0e-4;
double absolute_tolerance_angular_momentum_vectors = 1.0e-2;
bool include_quadrupole_order_terms = true;
bool include_octupole_order_binary_pair_terms = true;
bool include_octupole_order_binary_triplet_terms = true;
bool include_hexadecupole_order_binary_pair_terms = true;
bool include_dotriacontupole_order_binary_pair_terms = true;
bool include_double_averaging_corrections = false;
int random_seed = 0;
int verbose_flag = 0;
bool stop_after_root_found = false;

double secular_integration_exclusion_safety_factor = 1.0e-3;
double ODE_min_dt = 1.0;

double epsilon = 1.0e-12; /* used for tiny numbers close to machine precision */

double effective_radius_multiplication_factor_for_collisions_stars = 3.0;
double effective_radius_multiplication_factor_for_collisions_compact_objects = 1.0e3;

bool include_flybys = true;
bool flybys_correct_for_gravitational_focussing = true;
int flybys_velocity_distribution = 0;
int flybys_mass_distribution = 0;
int flybys_N_enc = 0;
int flybys_N_not_impulsive = 0;
int flybys_reference_binary = -1;
double flybys_t_next_encounter = 0.0;
double flybys_mass_distribution_lower_value = 0.1;
double flybys_mass_distribution_upper_value = 100.0;
double flybys_encounter_sphere_radius = 1.0e5;
double flybys_stellar_density = 0.1*CONST_PER_PC3; /* density at infinity */
double flybys_stellar_relative_velocity_dispersion = 30.0*CONST_KM_PER_S;
double flybys_W_max = 0.0;
double flybys_total_encounter_rate_at_R_enc = 0.0;
double flybys_stellar_density_at_R_enc = 0.0; /* density at R_enc */
double flybys_internal_mass = 0.0;
double flybys_internal_semimajor_axis = 0.0;

double MSTAR_gbs_tolerance_default = 1.0e-12;
double MSTAR_gbs_tolerance_kick = 1.0e-8;
double MSTAR_stopping_condition_tolerance = 1.0e-10;
double MSTAR_output_time_tolerance = 1.0e-4;

double nbody_analysis_fractional_semimajor_axis_change_parameter = 0.01;
double nbody_analysis_fractional_integration_time = 0.05;
double nbody_analysis_minimum_integration_time = 1.0e1;
double nbody_analysis_maximum_integration_time = 1.0e5;

double nbody_dynamical_instability_direct_integration_time_multiplier = 1.5;
double nbody_semisecular_direct_integration_time_multiplier = 1.0e2;
double nbody_supernovae_direct_integration_time_multiplier = 1.5;    
double nbody_other_direct_integration_time_multiplier = 1.5;
double nbody_maximum_separation_for_inclusion = 1.0e5;

int binary_evolution_CE_energy_flag = 0;
int binary_evolution_CE_spin_flag = 1;
int binary_evolution_CE_mass_loss_Nsteps = 100;
double binary_evolution_mass_transfer_timestep_parameter = 0.05;

double chandrasekhar_mass = 1.44;
double eddington_accretion_factor = 10.0;
double nova_accretion_factor = 1.0e-3;
double alpha_wind_accretion = 1.5;
double beta_wind_accretion = 0.125;

double triple_mass_transfer_primary_star_accretion_efficiency_no_disk = 0.1;
double triple_mass_transfer_secondary_star_accretion_efficiency_no_disk = 0.1;
double triple_mass_transfer_primary_star_accretion_efficiency_disk = 0.9;
double triple_mass_transfer_secondary_star_accretion_efficiency_disk = 0.9;
double triple_mass_transfer_inner_binary_alpha_times_lambda = 5.0;

double kroupa_alpha1 = -1.3;
double kroupa_alpha2 = -2.2;
double kroupa_alpha3 = -2.7;
double kroupa_m1 = 0.1;
double kroupa_m2 = 0.5;
double kroupa_m3 = 1.0;
double kroupa_m4 = 100.0;

double kroupa_alpha1_plus_1 = 1.0 + kroupa_alpha1;
double kroupa_alpha2_plus_1 = 1.0 + kroupa_alpha2;
double kroupa_alpha3_plus_1 = 1.0 + kroupa_alpha3;
double kroupa_alpha1_plus_1_pm1 = 1.0/kroupa_alpha1_plus_1;
double kroupa_alpha2_plus_1_pm1 = 1.0/kroupa_alpha2_plus_1;
double kroupa_alpha3_plus_1_pm1 = 1.0/kroupa_alpha3_plus_1;

double kroupa_m1_pow_alpha1_plus_one = pow(kroupa_m1,kroupa_alpha1_plus_1);
double kroupa_m2_pow_alpha1_plus_one = pow(kroupa_m2,kroupa_alpha1_plus_1);
double kroupa_m2_pow_alpha2_plus_one = pow(kroupa_m2,kroupa_alpha2_plus_1);
double kroupa_m3_pow_alpha2_plus_one = pow(kroupa_m3,kroupa_alpha2_plus_1);
double kroupa_m3_pow_alpha3_plus_one = pow(kroupa_m3,kroupa_alpha3_plus_1);
double kroupa_m4_pow_alpha3_plus_one = pow(kroupa_m4,kroupa_alpha3_plus_1);

double kroupa_C1 = 1.0/( (1.0/kroupa_alpha1_plus_1)*(kroupa_m2_pow_alpha1_plus_one - kroupa_m1_pow_alpha1_plus_one) + pow(kroupa_m2,kroupa_alpha1-kroupa_alpha2)*(1.0/kroupa_alpha2_plus_1)*(kroupa_m3_pow_alpha2_plus_one - kroupa_m2_pow_alpha2_plus_one) + pow(kroupa_m2,kroupa_alpha1-kroupa_alpha2)*pow(kroupa_m3,kroupa_alpha2-kroupa_alpha3)*(1.0/kroupa_alpha3_plus_1)*(kroupa_m4_pow_alpha3_plus_one - kroupa_m3_pow_alpha3_plus_one) );
double kroupa_C2 = kroupa_C1 * pow(kroupa_m2,kroupa_alpha1-kroupa_alpha2);
double kroupa_C3 = kroupa_C2 * pow(kroupa_m3,kroupa_alpha2-kroupa_alpha3);

double kroupa_x1 = (kroupa_C1/kroupa_alpha1_plus_1)*( kroupa_m2_pow_alpha1_plus_one - kroupa_m1_pow_alpha1_plus_one);
double kroupa_x2 = kroupa_x1 + (kroupa_C2/kroupa_alpha2_plus_1)*( kroupa_m3_pow_alpha2_plus_one - kroupa_m2_pow_alpha2_plus_one);
double kroupa_x3 = kroupa_x2 + (kroupa_C3/kroupa_alpha3_plus_1)*( kroupa_m4_pow_alpha3_plus_one - kroupa_m3_pow_alpha3_plus_one);



/* Used in MSTAR only */
double SPEEDOFLIGHT = CONST_C_LIGHT;
double GCONST = CONST_G;

double MSTAR_maximum_separation_for_inclusion = 1.0e5;
bool MSTAR_verbose = true;

bool MSTAR_include_PN_acc_10 = true;
bool MSTAR_include_PN_acc_20 = true;
bool MSTAR_include_PN_acc_25 = true;
bool MSTAR_include_PN_acc_30 = true;
bool MSTAR_include_PN_acc_35 = true;

bool MSTAR_include_PN_acc_SO = true;
bool MSTAR_include_PN_acc_SS = true;
bool MSTAR_include_PN_acc_Q = true;

bool MSTAR_include_PN_spin_SO = true;
bool MSTAR_include_PN_spin_SS = true;
bool MSTAR_include_PN_spin_Q = true;

