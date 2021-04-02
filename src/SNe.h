#include "types.h"
extern "C"
{

int handle_SNe_in_system(ParticlesMap *particlesMap, bool *unbound_orbits, int *integration_flag);

int sample_kick_velocity(Particle *p, int kick_distribution, double *vx, double *vy, double *vz);
bool check_for_unbound_orbits(ParticlesMap *particlesMap);
void remove_massless_remnants_from_system(ParticlesMap *particlesMap, int *integration_flag);
}
