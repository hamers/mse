#include "types.h"
extern "C"
{

int handle_binary_evolution(ParticlesMap *particlesMap, double *dt_binary_evolution, int *integration_flag);

int compute_mass_transfer_properties(ParticlesMap *particlesMap, double *dt_binary_evolution, int *integration_flag);

double compute_q_crit_for_common_envelope_evolution(int kw, double mass, double core_mass);
int common_envelope_evolution(ParticlesMap *particlesMap, int binary_index, int donor_index, int accretor_index, int *integration_flag);

}
