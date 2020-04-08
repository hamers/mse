#include "types.h"

extern "C"
{
void apply_user_specified_instantaneous_perturbation(ParticlesMap *particlesMap);
void apply_user_specified_instantaneous_perturbation_nbody(ParticlesMap *particlesMap);
void reset_instantaneous_perturbation_quantities(ParticlesMap *particlesMap);
void update_masses_positions_and_velocities_of_all_bodies_instantaneous_perturbation(ParticlesMap *particlesMap);

void update_position_vectors_external_particles(ParticlesMap *particlesMap, double time);
void compute_position_vectors_external_particle(ParticlesMap *particlesMap, Particle *perturber, double time, double *r_per, double r_per_vec[3]);

void compute_EOM_Newtonian_external_perturbations(double time, ParticlesMap *particlesMap, Particle *p, double *hamiltonian, double *KS_V, bool compute_hamiltonian_only); 
void compute_EOM_binary_pairs_external_perturbation(ParticlesMap *particlesMap, int binary_index, int perturber_index, double time, double *hamiltonian, double *KS_V, bool compute_hamiltonian_only);

int apply_external_perturbation_assuming_integrated_orbits(ParticlesMap *particlesMap);
void apply_external_perturbation_assuming_integrated_orbits_binary_pair(ParticlesMap *particlesMap, int binary_index, int perturber_index);

double retrieve_D_function(int function_index, double ep, double ep_p2, double ef, double ef_p2, double ef_p4, double one_div_ef_p1, double one_div_ef_p2, double one_div_ef_p3, double ef_p2_minus_one, double sqrt_ef_p2_minus_one, double asec_minus_ef);
double retrieve_D_function_e_derivative(int function_index, double ep, double ep_p2, double ef, double ef_p2, double ef_p4, double one_div_ef_p1, double one_div_ef_p2, double one_div_ef_p3, double ef_p2_minus_one, double sqrt_ef_p2_minus_one, double asec_minus_ef);
}
