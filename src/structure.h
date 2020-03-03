#include "types.h"
extern "C"
{
int determine_binary_parents_and_levels(ParticlesMap *particlesMap, int *N_bodies, int *N_binaries, int *N_root_finding);
void set_binary_masses_from_body_masses(ParticlesMap *particlesMap);
void determine_internal_mass_and_semimajor_axis(ParticlesMap *particlesMap);
}
