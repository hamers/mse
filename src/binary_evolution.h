#include "types.h"
extern "C"
{

int handle_binary_evolution(ParticlesMap *particlesMap, double *dt_binary_evolution);
int compute_mass_transfer_rates(ParticlesMap *particlesMap, double *dt_binary_evolution);

}
