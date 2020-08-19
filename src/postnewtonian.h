#include "types.h"
extern "C"
{
void compute_EOM_Post_Newtonian_for_particle(ParticlesMap *particlesMap, Particle *p, double *hamiltonian, bool compute_hamiltonian_only);
double compute_EOM_pairwise_1PN(ParticlesMap *particlesMap, int binary_index, bool compute_hamiltonian_only);
double compute_EOM_pairwise_25PN(ParticlesMap *particlesMap, int binary_index, bool compute_hamiltonian_only);

double compute_spin_parameter_from_spin_frequency(double m, double Omega);
double compute_spin_frequency_from_spin_parameter(double m, double chi);

double compute_1PN_timescale(double a, double M, double e);
}
