#include "types.h"
extern "C"
{
int determine_binary_parents_and_levels(ParticlesMap *particlesMap, int *N_bodies, int *N_binaries, int *N_root_finding, int *N_ODE_equations);
void set_binary_masses_from_body_masses(ParticlesMap *particlesMap);
void set_positions_and_velocities(ParticlesMap *particlesMap);
void update_masses_positions_and_velocities_of_all_binaries(ParticlesMap *particlesMap);
void update_orbital_vectors_in_binaries_from_positions_and_velocities(ParticlesMap *particlesMap);

void determine_internal_mass_and_semimajor_axis(ParticlesMap *particlesMap);
double determine_longest_orbital_period_in_system(ParticlesMap *particlesMap);
double determine_shortest_orbital_period_in_system(ParticlesMap *particlesMap);
int determine_number_of_bodies_in_system(ParticlesMap *particlesMap);

void update_structure(ParticlesMap *particlesMap);

void set_up_derived_quantities(ParticlesMap *particlesMap);

double compute_rp_out_crit_MA01(double a_in, double q_out, double e_out, double rel_INCL);
bool check_system_for_dynamical_stability(ParticlesMap *particlesMap, int *integration_flag);
}
