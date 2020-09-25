#include "types.h"
extern "C"
{
double compute_AM_time_scale(Particle *P_p);
void cross_section_function(Particle *p,double *cross_section);

int root_finding_functions(realtype t, N_Vector y, realtype *root_functions, void *data_);
double roche_radius_pericenter_eggleton(double rp, double q);
double roche_radius_pericenter_sepinsky(double rp, double q, double e, double f);

void check_for_roots(ParticlesMap *particlesMap, bool use_root_functions, realtype *root_functions);
int investigate_roots_in_system(ParticlesMap *particlesMap, double t, int integration_flag);

int read_root_finding_data(ParticlesMap *particlesMap, int *roots_found);
int check_for_initial_roots(ParticlesMap *particlesMap);

void handle_roots(ParticlesMap *particlesMap, int root_flag, int *integration_flag, int *CVODE_flag, double t, double *dt_stev, double *dt_binary_evolution);
}
