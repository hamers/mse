#include "types.h"
extern "C"
{

int test_tools();
int test_a_P_orb_conversion();
int test_a_h_conversion();
int test_orbital_element_conversion();
int test_kepler_equation_solver();
int test_orbital_vectors_cartesian_conversion();

int test_stellar_evolution();
int test_apsidal_motion_constant();
int test_sse();
int test_sse_specific_model(double m, double z, int *kw_final, double *m_init_final, double *m_final, double *R_final, double *ospin_final, double *L_final, double *m_core_final, double *m_env_final, double *epoch_final);
int test_kick_velocity(int kick_distribution, double m, int *kw, double *v_norm);

int test_binary_evolution();
int test_compute_Kelvin_Helmholtz_timescale();
int test_compute_Eddington_accretion_rate();
int test_handle_instantaneous_and_adiabatic_mass_changes_in_orbit();

int test_collisions();
int test_collision_MS_MS();
int test_collision_giant_MS();

int test_collision_stars(double m1, int kw1, double m2, int kw2, int integration_flag);

}
