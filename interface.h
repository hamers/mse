#include "src/types.h"
extern "C"
{

/*******************
/* basic interface *
 ******************/

int add_particle(int *index, bool is_binary, bool is_external);
int delete_particle(int index);

int set_children(int index, int child1, int child2);
int get_children(int index, int *child1, int *child2);

int get_number_of_particles();
int get_internal_index_in_particlesMap(int absolute_index);
bool get_is_binary(int index);
int get_is_bound(int index, bool *is_bound);

int set_mass(int index, double mass);
int get_mass(int index, double *mass);

int set_mass_transfer_terms(int index, bool include_mass_transfer_terms);
int get_mass_transfer_terms(int index, bool *include_mass_transfer_terms);

int get_mass_dot(int index, double *mass_dot);

int set_radius(int index, double radius, double radius_dot);
int get_radius(int index, double *radius, double *radius_dot);

int get_level(int index, int *level);

int set_stellar_evolution_properties(int index, int stellar_type, int object_type, double sse_initial_mass, double metallicity, double sse_time_step, double epoch, double age, 
    double convective_envelope_mass, double convective_envelope_radius, double core_mass, double core_radius, double luminosity, double apsidal_motion_constant, double gyration_radius, double tides_viscous_time_scale, int tides_viscous_time_scale_prescription);
int get_stellar_evolution_properties(int index, int *stellar_type, int *object_type, double *sse_initial_mass, double *metallicity, double *sse_time_step, double *epoch, double *age, 
    double *convective_envelope_mass, double *convective_envelope_radius, double *core_mass, double *core_radius, double *luminosity, double *apsidal_motion_constant, double *gyration_radius, double *tides_viscous_time_scale, double *roche_lobe_radius_pericenter);  
int set_kick_properties(int index, int kick_distribution, bool include_WD_kicks, double kick_distribution_sigma_km_s_NS, double kick_distribution_sigma_km_s_BH, double kick_distribution_sigma_km_s_WD, double kick_distribution_2_m_NS, double kick_distribution_4_m_NS, double kick_distribution_4_m_ej, \
    double kick_distribution_5_v_km_s_NS, double kick_distribution_5_v_km_s_BH, double kick_distribution_5_sigma, double kick_distribution_sigma_km_s_NS_ECSN);
int get_kick_properties(int index, int *kick_distribution, bool *include_WD_kicks, double *kick_distribution_sigma_km_s_NS, double *kick_distribution_sigma_km_s_BH,  double *kick_distribution_sigma_km_s_WD, double *kick_distribution_2_m_NS, double *kick_distribution_4_m_NS, double *kick_distribution_4_m_ej, \
    double *kick_distribution_5_v_km_s_NS, double *kick_distribution_5_v_km_s_BH, double *kick_distribution_5_sigma, double *kick_distribution_sigma_km_s_NS_ECSN);

int set_binary_evolution_properties(int index, double dynamical_mass_transfer_low_mass_donor_timescale, double dynamical_mass_transfer_WD_donor_timescale, double compact_object_disruption_mass_loss_timescale, \
    double common_envelope_alpha, double common_envelope_lambda, double common_envelope_timescale, double triple_common_envelope_alpha);
int get_binary_evolution_properties(int index, double *dynamical_mass_transfer_low_mass_donor_timescale, double *dynamical_mass_transfer_WD_donor_timescale, double *compact_object_disruption_mass_loss_timescale, \
    double *common_envelope_alpha, double *common_envelope_lambda, double *common_envelope_timescale, double *triple_common_envelope_alpha);

int set_true_anomaly(int index, double value);
int get_true_anomaly(int index, double *value);

int set_sample_orbital_phases_randomly(int index, int value);
int get_sample_orbital_phases_randomly(int index, int *value);


/*******************************
 * instantaneous perturbations *
 * ****************************/

int set_instantaneous_perturbation_properties(int index, double delta_mass, double delta_X, double delta_Y, double delta_Z, double delta_VX, double delta_VY, double delta_VZ);


/************
 * external *
 * *********/
 
int set_external_particle_properties(int index, double t_ref, double e, double r_p, double INCL, double AP, double LAN);


/****************
/* spin vectors *
 ****************/
int set_spin_vector(int index, double spin_vec_x, double spin_vec_y, double spin_vec_z);
int get_spin_vector(int index, double *spin_vec_x, double *spin_vec_y, double *spin_vec_z);

int set_spin_vector_dot(int index, double spin_vec_x_dot, double spin_vec_y_dot, double spin_vec_z_dot);
int get_spin_vector_dot(int index, double *spin_vec_x_dot, double *spin_vec_y_dot, double *spin_vec_z_dot);

/****************************
/* orbital vectors/elements *
 ****************************/

int set_orbital_vectors(int index, double e_vec_x, double e_vec_y, double e_vec_z, \
    double h_vec_x, double h_vec_y, double h_vec_z);
int get_orbital_vectors(int index, double *e_vec_x, double *e_vec_y, double *e_vec_z, \
    double *h_vec_x, double *h_vec_y, double *h_vec_z);
    
int set_orbital_elements(int index, double semimajor_axis, double eccentricity, double true_anomaly, \
    double inclination, double argument_of_pericenter, double longitude_of_ascending_node, bool sample_orbital_phase_randomly);
int get_orbital_elements(int index, double *semimajor_axis, double *eccentricity, double *true_anomaly, \
    double *inclination, double *argument_of_pericenter, double *longitude_of_ascending_node);


int get_relative_position_and_velocity(int index, double *x, double *y, double *z, double *vx, double *vy, double *vz);
int get_absolute_position_and_velocity(int index, double *X, double *Y, double *Z, double *VX, double *VY, double *VZ);

    
/************
/* PN terms *
 ************/
int set_PN_terms(int index, bool include_pairwise_1PN_terms, bool include_pairwise_25PN_terms, bool include_spin_orbit_1PN_terms, bool exclude_1PN_precession_in_case_of_isolated_binary);
int get_PN_terms(int index, bool *include_pairwise_1PN_terms, bool *include_pairwise_25PN_terms, bool *include_spin_orbit_1PN_terms, bool *exclude_1PN_precession_in_case_of_isolated_binary);


/*********
/* tides *
 *********/
int set_tides_terms(int index, bool include_tidal_friction_terms, int tides_method, bool include_tidal_bulges_precession_terms, bool include_rotation_precession_terms, double minimum_eccentricity_for_tidal_precession, bool exclude_rotation_and_bulges_precession_in_case_of_isolated_binary);
int get_tides_terms(int index, bool *include_tidal_friction_terms, int *tides_method, bool *include_tidal_bulges_precession_terms, bool *include_rotation_precession_terms, double *minimum_eccentricity_for_tidal_precession, bool *exclude_rotation_and_bulges_precession_in_case_of_isolated_binary);

/*******
 * VRR *
 * *****/
 
int set_VRR_properties(int index, int VRR_model, int VRR_include_mass_precession, double VRR_mass_precession_rate, 
    double VRR_Omega_vec_x, double VRR_Omega_vec_y, double VRR_Omega_vec_z, 
    double VRR_eta_20_init, double VRR_eta_a_22_init, double VRR_eta_b_22_init, double VRR_eta_a_21_init, double VRR_eta_b_21_init,
    double VRR_eta_20_final, double VRR_eta_a_22_final, double VRR_eta_b_22_final, double VRR_eta_a_21_final, double VRR_eta_b_21_final,
	double VRR_initial_time, double VRR_final_time);

/****************
/* root finding *
 ****************/

int set_root_finding_terms(int index, bool check_for_secular_breakdown, bool check_for_dynamical_instability, int dynamical_instability_criterion, int dynamical_instability_central_particle, double dynamical_instability_K_parameter,
    bool check_for_physical_collision_or_orbit_crossing, bool check_for_minimum_periapse_distance, double check_for_minimum_periapse_distance_value, bool check_for_RLOF_at_pericentre, bool check_for_RLOF_at_pericentre_use_sepinsky_fit, bool check_for_GW_condition);
int get_root_finding_terms(int index, bool *check_for_secular_breakdown, bool *check_for_dynamical_instability, int *dynamical_instability_criterion, int *dynamical_instability_central_particle, double *dynamical_instability_K_parameter,
    bool *check_for_physical_collision_or_orbit_crossing, bool *check_for_minimum_periapse_distance, double *check_for_minimum_periapse_distance_value, bool *check_for_RLOF_at_pericentre, bool *check_for_RLOF_at_pericentre_use_sepinsky_fit, bool *check_for_GW_condition);
int set_root_finding_state(int index, bool secular_breakdown_has_occurred, bool dynamical_instability_has_occurred, bool physical_collision_or_orbit_crossing_has_occurred, bool minimum_periapse_distance_has_occurred, bool RLOF_at_pericentre_has_occurred, bool GW_condition_has_occurred);
int get_root_finding_state(int index, bool *secular_breakdown_has_occurred, bool *dynamical_instability_has_occurred, bool *physical_collision_or_orbit_crossing_has_occurred, bool *minimum_periapse_distance_has_occurred, bool *RLOF_at_pericentre_has_occurred, bool *GW_condition_has_occurred);


/***********************
/* interface functions *
 ***********************/
int evolve_interface(double start_time, double end_time, double *output_time, double *hamiltonian, int *state, int *CVODE_flag, int *CVODE_error_code, int *integration_flag);
int determine_binary_parents_levels_and_masses_interface();
int apply_external_perturbation_assuming_integrated_orbits_interface();
int apply_external_perturbation_assuming_integrated_orbits_single_perturber_interface(double M_per, double e_per, double Q_per, double e_vec_unit_per_x, double e_vec_unit_per_y, double e_vec_unit_per_z, double h_vec_unit_per_x, double h_vec_unit_per_y, double h_vec_unit_per_z);
int apply_user_specified_instantaneous_perturbation_interface();
int set_positions_and_velocities_interface();
int reset_interface();
int initialize_code_interface();

int get_inclination_relative_to_parent_interface(int index, double *inclination_relative_to_parent);
int get_de_dt(int index, double *de_dt);


/************************
/* parameters *
 ************************/

int set_constants(double CONST_G_, double CONST_C_, double CONST_MSUN_, double CONST_R_SUN_, double CONST_L_SUN_, double CONST_KM_PER_S_, double CONST_PER_PC3_, double CONST_MJUP_);
int set_parameters(double relative_tolerance_, double absolute_tolerance_eccentricity_vectors_, double absolute_tolerance_spin_vectors_, double absolute_tolerance_angular_momentum_vectors_,
    bool include_quadrupole_order_terms_, bool include_octupole_order_binary_pair_terms_, bool include_octupole_order_binary_triplet_terms_,
    bool include_hexadecupole_order_binary_pair_terms_, bool include_dotriacontupole_order_binary_pair_terms_,  bool include_double_averaging_corrections_,
    bool include_flybys_, int flybys_reference_binary_, bool flybys_correct_for_gravitational_focussing_, int flybys_velocity_distribution_, int flybys_mass_distribution_, bool flybys_include_secular_encounters_,
    double flybys_mass_distribution_lower_value_, double flybys_mass_distribution_upper_value_, double flybys_encounter_sphere_radius_, 
    double flybys_stellar_density_, double flybys_stellar_relative_velocity_dispersion_,
    int binary_evolution_CE_energy_flag_, int binary_evolution_CE_spin_flag_, double binary_evolution_mass_transfer_timestep_parameter_, double binary_evolution_CE_recombination_fraction_, bool binary_evolution_use_eCAML_model_, \
    double MSTAR_gbs_tolerance_default_, double MSTAR_gbs_tolerance_kick_, double MSTAR_stopping_condition_tolerance_, double MSTAR_output_time_tolerance_, \
    double nbody_analysis_fractional_semimajor_axis_change_parameter_, double nbody_analysis_fractional_integration_time_, double nbody_analysis_minimum_integration_time_, double nbody_analysis_maximum_integration_time_, \
    double nbody_dynamical_instability_direct_integration_time_multiplier_, double nbody_semisecular_direct_integration_time_multiplier_, double nbody_supernovae_direct_integration_time_multiplier_, double nbody_other_direct_integration_time_multiplier_, \
    double chandrasekhar_mass_, double eddington_accretion_factor_, double nova_accretion_factor_, double alpha_wind_accretion_, double beta_wind_accretion_, \
    double triple_mass_transfer_primary_star_accretion_efficiency_no_disk_, double triple_mass_transfer_secondary_star_accretion_efficiency_no_disk_, double triple_mass_transfer_primary_star_accretion_efficiency_disk_, double triple_mass_transfer_secondary_star_accretion_efficiency_disk_, double triple_mass_transfer_inner_binary_alpha_times_lambda_, \
    double effective_radius_multiplication_factor_for_collisions_stars_, double effective_radius_multiplication_factor_for_collisions_compact_objects_, \
    bool MSTAR_include_PN_acc_10_,bool MSTAR_include_PN_acc_20_,bool MSTAR_include_PN_acc_25_,bool MSTAR_include_PN_acc_30_,bool MSTAR_include_PN_acc_35_,bool MSTAR_include_PN_acc_SO_,bool MSTAR_include_PN_acc_SS_,bool MSTAR_include_PN_acc_Q_,bool MSTAR_include_PN_spin_SO_,bool MSTAR_include_PN_spin_SS_,bool MSTAR_include_PN_spin_Q_, \
    bool stop_after_root_found_, \
    double wall_time_max_s_, \
    int NS_model_, int ECSNe_model, \
    int system_index_, \
    int binary_evolution_mass_transfer_model_, int binary_evolution_SNe_Ia_single_degenerate_model_, int binary_evolution_SNe_Ia_double_degenerate_model_, double binary_evolution_SNe_Ia_double_degenerate_model_minimum_eccentricity_for_eccentric_collision_, double binary_evolution_SNe_Ia_double_degenerate_model_minimum_primary_mass_CO_CO_);

int get_random_seed(int *value);
int set_random_seed(int value);

int get_verbose_flag(int *value);
int set_verbose_flag(int value);

/***********
 * Testing *
 * ********/

int unit_tests_interface(int mode);

int determine_compact_object_merger_properties_interface(double m1, double m2, double chi1, double chi2, \
    double spin_vec_1_unit_x, double spin_vec_1_unit_y, double spin_vec_1_unit_z, \
    double spin_vec_2_unit_x, double spin_vec_2_unit_y, double spin_vec_2_unit_z, \
    double h_vec_unit_x, double h_vec_unit_y, double h_vec_unit_z, \
    double e_vec_unit_x, double e_vec_unit_y, double e_vec_unit_z, \
    double *v_recoil_vec_x, double *v_recoil_vec_y, double *v_recoil_vec_z, \
    double *alpha_vec_final_x, double *alpha_vec_final_y, double *alpha_vec_final_z, \
    double *M_final);

int sample_from_3d_maxwellian_distribution_interface(double sigma, double *vx, double *vy, double *vz);
double sample_from_normal_distribution_interface(double mu, double sigma);
double sample_from_kroupa_93_imf_interface();
int sample_spherical_coordinates_unit_vectors_from_isotropic_distribution_interface(double *r_hat_vec_x, double *r_hat_vec_y,double *r_hat_vec_z, \
    double *theta_hat_vec_x, double *theta_hat_vec_y,double *theta_hat_vec_z, \
    double *phi_hat_vec_x, double *phi_hat_vec_y,double *phi_hat_vec_z);



/***********
 * Logging *
 * ********/


int get_size_of_log_data();
int get_internal_index_in_particlesMap_log(int log_index, int absolute_index);
bool get_is_binary_log(int log_index, int particle_index);
int get_log_entry_properties(int log_index, double *time, int *event_flag, int *integration_flag, int *N_particles, int *index1, int *index2, int *binary_index, double *kick_speed_km_s, int *SNe_type, int *SNe_info, int *eccentric_collision);
int get_body_properties_from_log_entry(int log_index, int particle_index, int *parent, double *mass, double *radius, int *stellar_type, double *core_mass, double *sse_initial_mass, double *convective_envelope_mass, \
    double *epoch, double *age, double *core_radius, double *convective_envelope_radius, double *luminosity, double *ospin, double *X, double *Y, double *Z, double *VX, double *VY, double *VZ, int *object_type, double *metallicity, 
    double *spin_vec_x, double *spin_vec_y, double *spin_vec_z);
int get_binary_properties_from_log_entry(int log_index, int particle_index, int *parent, int *child1, int *child2, \
    double *mass, double *a, double *e, double *TA, double *INCL, double *AP, double *LAN, \
    double *h_vec_x, double *h_vec_y, double *h_vec_z, \
    double *e_vec_x, double *e_vec_y, double *e_vec_z);
int write_final_log_entry_interface(double t, int integration_flag);

}
