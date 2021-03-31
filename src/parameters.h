#ifndef __PARAMETERS
#define __PARAMETERS
#include <stdbool.h>
#include <time.h>

// Constants //
extern double CONST_G;
extern double CONST_G_P2;
extern double CONST_G_P3;
extern double CONST_C_LIGHT;
extern double CONST_C_LIGHT_P2;
extern double CONST_C_LIGHT_P3;
extern double CONST_C_LIGHT_P4;
extern double CONST_C_LIGHT_P5;
extern double CONST_C_LIGHT_P6;
extern double CONST_MSUN;
extern double CONST_R_SUN;
extern double CONST_L_SUN;
extern double CONST_KM_PER_S;
extern double CONST_PER_PC3;
extern double CONST_MJUP;

// Parameters //
extern double relative_tolerance;
extern double absolute_tolerance_eccentricity_vectors;
extern double absolute_tolerance_spin_vectors;
extern double absolute_tolerance_angular_momentum_vectors;
extern bool include_quadrupole_order_terms;
extern bool include_octupole_order_binary_pair_terms;
extern bool include_octupole_order_binary_triplet_terms;
extern bool include_hexadecupole_order_binary_pair_terms;
extern bool include_dotriacontupole_order_binary_pair_terms;
extern bool include_double_averaging_corrections;
extern double epsilon;
extern int random_seed;
extern int verbose_flag;
extern bool stop_after_root_found;
extern bool check_numbers_internally;
extern int error_code;
extern time_t wall_time_start;
extern double wall_time_max_s;

extern double secular_integration_exclusion_safety_factor;
extern double ODE_min_dt;

extern double effective_radius_multiplication_factor_for_collisions_stars;
extern double effective_radius_multiplication_factor_for_collisions_compact_objects;

extern bool include_flybys;
extern bool flybys_correct_for_gravitational_focussing;
extern int flybys_velocity_distribution;
extern int flybys_mass_distribution;
extern int flybys_N_enc;
extern int flybys_N_not_impulsive;
extern int flybys_reference_binary;
extern double flybys_t_next_encounter;
extern double flybys_mass_distribution_lower_value;
extern double flybys_mass_distribution_upper_value;
extern double flybys_encounter_sphere_radius;
extern double flybys_stellar_density;
extern double flybys_stellar_relative_velocity_dispersion;
extern double flybys_W_max;
extern double flybys_total_encounter_rate_at_R_enc;
extern double flybys_stellar_density_at_R_enc;
extern double flybys_internal_mass;
extern double flybys_internal_semimajor_axis;

extern int binary_evolution_CE_energy_flag;
extern int binary_evolution_CE_spin_flag;
extern int binary_evolution_CE_mass_loss_Nsteps;
extern double binary_evolution_CE_recombination_fraction;
extern double binary_evolution_mass_transfer_timestep_parameter;
extern bool binary_evolution_use_eCAML_model;

extern double chandrasekhar_mass;
extern double eddington_accretion_factor;
extern double nova_accretion_factor;
extern double alpha_wind_accretion;
extern double beta_wind_accretion;

extern double MSTAR_gbs_tolerance_default;
extern double MSTAR_gbs_tolerance_kick;
extern double MSTAR_stopping_condition_tolerance;
extern double MSTAR_output_time_tolerance;

extern double triple_mass_transfer_primary_star_accretion_efficiency_no_disk;
extern double triple_mass_transfer_secondary_star_accretion_efficiency_no_disk;
extern double triple_mass_transfer_primary_star_accretion_efficiency_disk;
extern double triple_mass_transfer_secondary_star_accretion_efficiency_disk;
extern double triple_mass_transfer_inner_binary_alpha_times_lambda;

extern double nbody_analysis_fractional_semimajor_axis_change_parameter;
extern double nbody_analysis_fractional_integration_time;
extern double nbody_analysis_minimum_integration_time;
extern double nbody_analysis_maximum_integration_time;

extern double nbody_dynamical_instability_direct_integration_time_multiplier;
extern double nbody_semisecular_direct_integration_time_multiplier;
extern double nbody_supernovae_direct_integration_time_multiplier;
extern double nbody_other_direct_integration_time_multiplier;
extern double nbody_maximum_separation_for_inclusion;

extern double kroupa_alpha1;
extern double kroupa_alpha2;
extern double kroupa_alpha3;
extern double kroupa_m1;
extern double kroupa_m2;
extern double kroupa_m3;
extern double kroupa_m4;

extern double kroupa_alpha1_plus_1;
extern double kroupa_alpha2_plus_1;
extern double kroupa_alpha3_plus_1;
extern double kroupa_alpha1_plus_1_pm1;
extern double kroupa_alpha2_plus_1_pm1;
extern double kroupa_alpha3_plus_1_pm1;

extern double kroupa_m1_pow_alpha1_plus_one;
extern double kroupa_m2_pow_alpha1_plus_one;
extern double kroupa_m2_pow_alpha2_plus_one;
extern double kroupa_m3_pow_alpha2_plus_one;
extern double kroupa_m3_pow_alpha3_plus_one;
extern double kroupa_m4_pow_alpha3_plus_one;

extern double kroupa_C1;
extern double kroupa_C2;
extern double kroupa_C3;

extern double kroupa_x1;
extern double kroupa_x2;
extern double kroupa_x3;


/* Used in MSTAR only */
extern double SPEEDOFLIGHT;
extern double GCONST;

extern double MSTAR_maximum_separation_for_inclusion;
extern bool MSTAR_verbose;

extern bool MSTAR_include_PN_acc_10;
extern bool MSTAR_include_PN_acc_20;
extern bool MSTAR_include_PN_acc_25;
extern bool MSTAR_include_PN_acc_30;
extern bool MSTAR_include_PN_acc_35;

extern bool MSTAR_include_PN_acc_SO;
extern bool MSTAR_include_PN_acc_SS;
extern bool MSTAR_include_PN_acc_Q;

extern bool MSTAR_include_PN_spin_SO;
extern bool MSTAR_include_PN_spin_SS;
extern bool MSTAR_include_PN_spin_Q;


#endif
