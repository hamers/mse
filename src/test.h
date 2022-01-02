#include "types.h"
extern "C"
{

int test_tools();
int test_a_P_orb_conversion();
int test_a_h_conversion();
int test_orbital_element_conversion();
int test_kepler_equation_solver();
int test_orbital_vectors_cartesian_conversion();
int test_BH_spin_conversion();
int test_interpolate_1D_data_table_linear();

int test_nbody(int mode);
int test_nbody_two_body_stopping_conditions();
struct RegularizedRegion *generate_binary_ICs(double m1, double m2, double R1, double R2, double a, double e, int stopping_condition_mode1, int stopping_condition_mode2, double gbs_tolerance, double stopping_condition_tolerance);
double test_nbody_compute_time_of_collision(double CONST_G, double M, double r_col, double a, double e);
int test_nbody_compute_elements(double CONST_G, double M, double *r, double *v, double *a, double *e);
int test_nbody_two_body_kick();
struct RegularizedRegion *generate_pythagorean_ICs(double m1, double m2, double m3, int stopping_condition_mode, double gbs_tolerance, double stopping_condition_tolerance);
int test_nbody_pythagorean();
int test_nbody_inspiral(int mode);
int test_nbody_spin_orbit(int mode);
int test_nbody_custom(int mode);
int test_nbody_custom2(int mode);

int test_flybys();
int test_flybys_integrals();
int test_flybys_compute_effects_of_flyby_on_system();
int test_flybys_perturber_sampling(double R_enc, double n_star, double sigma, double M_int, double *M_per, double *b_vec_x, double *b_vec_y, double *b_vec_z, double *V_vec_x, double *V_vec_y, double *V_vec_z);

int test_stellar_evolution(int mode);
int test_binary_evolution_SNe_Ia_single_degenerate_model_1_accumulation_efficiency(double M_WD, double accretion_rate, double luminosity, double *eta, int *WD_accretion_mode);
int test_binary_evolution_SNe_Ia_single_degenerate_model_1_explosion(double M_WD, double accretion_rate, double M_He, double luminosity, bool *explosion);
int test_binary_evolution_SNe_Ia_single_degenerate_model_1_white_dwarf_hydrogen_accretion_boundaries(double m_WD, double *m_dot_lower, double *m_dot_upper);
int test_mass_transfer_with_degenerate_objects_single_degenerate_model_1();
int test_SNe_Ia_single_and_double_degenerate_model_1();
int test_spin_conversion();
int test_apsidal_motion_constant();
int test_sse();
int test_sse_specific_model(double m, double z, int *kw_final, double *m_init_final, double *m_final, double *R_final, double *ospin_final, double *L_final, double *m_core_final, double *m_env_final, double *epoch_final);
int test_sse_specific_model_evolve_function(double m, double z, int *kw_final, double *m_init_final, double *m_final, double *R_final, double *ospin_final, double *L_final, double *m_core_final, double *m_env_final, double *epoch_final);
int test_sse_specific_model_stopping_conditions(double m, double z, int *kw_final, double *m_init_final, double *m_final, double *R_final, double *ospin_final, double *L_final, double *m_core_final, double *m_env_final, double *epoch_final);
int test_sse_custom();
int test_kick_velocity(int kick_distribution, double m, int *kw, double *v_norm);
int test_remove_massless_remnants_from_system();
int test_NS_models(int mode);
int test_determine_sse_compact_object_radius(int mode);

int test_binary_evolution();
int test_compute_Kelvin_Helmholtz_timescale();
int test_compute_Eddington_accretion_rate();
int test_handle_instantaneous_and_adiabatic_mass_changes_in_orbit();
int test_wind_accretion();
int test_mass_accretion_events_with_degenerate_objects();
int test_mass_accretion_events_with_degenerate_objects_single_degenerate_model_1();
int test_compute_bse_mass_transfer_amount_averaged();
int test_binary_common_envelope_evolution();
int test_binary_evolution_emt_model_optimised_functions();

int test_collisions();

int test_collision_stars(double m1, int kw1, double m2, int kw2, int integration_flag);

int test_triple_interactions();
int test_triple_common_envelope_evolution();

}
