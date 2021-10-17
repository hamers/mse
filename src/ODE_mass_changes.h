#include "types.h"
extern "C"
{

int ODE_handle_stellar_winds(Particle *p);
int ODE_handle_RLOF(ParticlesMap *particlesMap, Particle *p);

double compute_Omega_dot_from_J_dot_mass_transfer(int stellar_type, double J_spin_dot, double Omega, double mass, double mass_dot, double core_mass, double radius, double radius_dot, double core_radius, double sse_k2, double sse_k3);
double compute_a_dot_circular_non_conservative_mass_transfer(double a, double m_donor, double m_donor_dot, double m_accretor, double beta, double gamma);

int compute_RLOF_emt_model(Particle *p, Particle *donor, Particle *accretor, double x, double E_0);
int ODE_handle_RLOF_emt(Particle *p, Particle *child1, Particle *child2);
int ODE_handle_RLOF_triple_mass_transfer(ParticlesMap *particlesMap, Particle *outer_orbit, Particle *donor, Particle *inner_binary);

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

double fm_function_old(double e, double x, double E0, double Etau);
double fa_function_old(double e, double x, double E0, double Etau);
double fe_function_old(double e, double x, double E0, double Etau);
double fomega_function_old(double e, double x, double E0, double Etau);
double ga_function_old(double e, double x, double E0);
double ge_function_old(double e, double x, double E0);
double ha_function_old(double e, double x, double E0);
double he_function_old(double e, double x, double E0);
double XL0_q_function_old(double q);
}
