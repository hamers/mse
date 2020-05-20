#include "types.h"
extern "C"
{

int initialize_stars(ParticlesMap *particlesMap);
int evolve_stars(ParticlesMap *particlesMap, double start_time, double end_time, double *stellar_evolution_timestep, bool get_timestep_only, bool *apply_SNe_effects);
double get_new_dt(int kw, double mass, double mt, double age, double dt, double *zpars);

void update_stellar_evolution_properties(Particle *p);

}
