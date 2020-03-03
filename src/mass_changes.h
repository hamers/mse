#include "types.h"
extern "C"
{

int ODE_handle_stellar_winds(Particle *p);
int ODE_handle_RLOF(Particle *p, Particle *child1, Particle *child2);
int check_if_RLOF_has_occurred(ParticlesMap *particlesMap);

int compute_RLOF_emt_model(Particle *p, Particle *donor, Particle *accretor, double x, double E_0);
int ODE_handle_RLOF_emt(Particle *p, Particle *child1, Particle *child2);

int determine_E_0(double e, double x, double *E_0, bool *in_RLOF);

double fm_function(double e, double x, double E0, double Etau);
double fa_function(double e, double x, double E0, double Etau);
double fe_function(double e, double x, double E0, double Etau);
double fomega_function(double e, double x, double E0, double Etau);
double ga_function(double e, double x, double E0);
double ge_function(double e, double x, double E0);
double ha_function(double e, double x, double E0);
double he_function(double e, double x, double E0);
double XL0_q_function(double q);
}
