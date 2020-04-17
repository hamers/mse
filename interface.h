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

int get_number_of_particles(int *N_particles);
int get_is_binary(int index, bool *is_binary);
int get_is_bound(int index, bool *is_bound);

int set_mass(int index, double mass);
int get_mass(int index, double *mass);

int set_mass_transfer_terms(int index, bool include_mass_transfer_terms);
int get_mass_dot(int index, double *mass_dot);

int set_radius(int index, double radius, double radius_dot);
int get_radius(int index, double *radius, double *radius_dot);

int get_level(int index, int *level);

int set_stellar_evolution_properties(int index, int stellar_type, bool evolve_as_star, double sse_initial_mass, double metallicity, double sse_time_step, double epoch, double age, 
    double convective_envelope_mass, double convective_envelope_radius, double core_mass, double core_radius, double luminosity, double apsidal_motion_constant, double gyration_radius, double tides_viscous_time_scale, int tides_viscous_time_scale_prescription);
int get_stellar_evolution_properties(int index, int *stellar_type, bool *evolve_as_star, double *sse_initial_mass, double *metallicity, double *sse_time_step, double *epoch, double *age, 
    double *convective_envelope_mass, double *convective_envelope_radius, double *core_mass, double *core_radius, double *luminosity, double *apsidal_motion_constant, double *gyration_radius, double *tides_viscous_time_scale, double *roche_lobe_radius_pericenter);    
int set_kick_properties(int index, int kick_distribution, double kick_distribution_sigma);
int get_kick_properties(int index, int *kick_distribution, double *kick_distribution_sigma);

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
int set_PN_terms(int index, int include_pairwise_1PN_terms, int include_pairwise_25PN_terms);
int get_PN_terms(int index, int *include_pairwise_1PN_terms, int *include_pairwise_25PN_terms);


/*********
/* tides *
 *********/
int set_tides_terms(int index, bool include_tidal_friction_terms, int tides_method, bool include_tidal_bulges_precession_terms, bool include_rotation_precession_terms, double minimum_eccentricity_for_tidal_precession);
int get_tides_terms(int index, bool *include_tidal_friction_terms, bool *tides_method, bool *include_tidal_bulges_precession_terms, bool *include_rotation_precession_terms, double *minimum_eccentricity_for_tidal_precession);

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
int apply_user_specified_instantaneous_perturbation_interface();
int set_positions_and_velocities_interface();
int clear_internal_particles();
int initialize_code();

/**********************************************
/* orbital element/vector conversion routines *
 **********************************************/
int get_inclination_relative_to_parent(int index, double *inclination_relative_to_parent);
int get_de_dt(int index, double *de_dt);

/************************
/* parameters *
 ************************/


int set_constants(double CONST_G_, double CONST_C_, double CONST_MSUN_, double CONST_R_SUN_, double CONST_L_SUN_, double CONST_KM_PER_S_, double CONST_PER_PC3_);
int set_parameters(double relative_tolerance_, double absolute_tolerance_eccentricity_vectors_, 
     bool include_quadrupole_order_terms_, bool include_octupole_order_binary_pair_terms_, bool include_octupole_order_binary_triplet_terms_,
     bool include_hexadecupole_order_binary_pair_terms_, bool include_dotriacontupole_order_binary_pair_terms_,  bool include_double_averaging_corrections_,
     bool include_flybys_, int flybys_reference_binary_, bool flybys_correct_for_gravitational_focussing_, int flybys_velocity_distribution_, int flybys_mass_distribution_,
     double flybys_mass_distribution_lower_value_, double flybys_mass_distribution_upper_value_, double flybys_encounter_sphere_radius_, 
     double flybys_stellar_density_, double flybys_stellar_relative_velocity_dispersion_,
     int binary_evolution_CE_energy_flag_, int binary_evolution_CE_spin_flag_);

int get_random_seed(int *value);
int set_random_seed(int value);

}
