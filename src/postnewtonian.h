#include "types.h"
extern "C"
{
double compute_EOM_pairwise_1PN(ParticlesMap *particlesMap, int binary_index, bool compute_hamiltonian_only);
double compute_EOM_pairwise_25PN(ParticlesMap *particlesMap, int binary_index, bool compute_hamiltonian_only);

double compute_spin_parameter_from_spin_frequency(double m, double Omega);
double compute_spin_frequency_from_spin_parameter(double m, double chi);

}
