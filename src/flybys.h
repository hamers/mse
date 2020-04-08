#include "types.h"
extern "C"
{

int handle_next_flyby(ParticlesMap *particlesMap, bool initialize, bool *unbound_orbits, int integration_flag);
int sample_next_flyby(ParticlesMap *particlesMap, bool *apply_flyby, double *t_next_encounter, int *N_enc, int *N_not_impulsive, double *M_per, double b_vec[3], double V_vec[3]);
int sample_flyby_position_and_velocity_at_R_enc(ParticlesMap *particlesMap, double R_vec[3], double V_vec[3]);
int compute_effects_of_flyby_on_system(ParticlesMap *particlesMap, double M_per, double b_per_vec[3], double V_per_vec[3], bool *unbound_orbits, int integration_flag);
int compute_total_encounter_rate_and_density_at_R_enc(double *total_encounter_rate, double *stellar_density);
double x_function(double M_tot);
double W_function(double x);
double V_function(double x);
double compute_W_max();
bool correct_mass_function(double M);
double sample_flyby_mass_at_R_enc();
double sample_flyby_mass_at_infinity();


}
