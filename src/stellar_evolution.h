#include "types.h"
extern "C"
{
void check_sse_error_codes();
int initialize_stars(ParticlesMap *particlesMap);
int evolve_stars(ParticlesMap *particlesMap, double start_time, double end_time, double *stellar_evolution_timestep, bool get_timestep_only, bool *apply_SNe_effects, int *integration_flag);
double get_new_dt_sse(int kw, double mass, double mt, double age, double dt, double *zpars);
void update_stellar_evolution_quantities_directly_secular(ParticlesMap *particlesMap, double dt_assumed, double t_old, double t_out);
void update_stellar_evolution_quantities_directly_nbody(ParticlesMap *particlesMap, double t, double dt);

void check_for_critical_rotation(ParticlesMap *particlesMap);
void update_stellar_evolution_properties(Particle *p);
double compute_moment_of_inertia(int stellar_type, double mass, double core_mass, double radius, double core_radius, double k2, double k3);
double compute_moment_of_inertia_dot(int stellar_type, double mass, double core_mass, double radius, double core_radius, double k2, double k3, double mass_dot, double radius_dot);
void get_core_masses_by_composition(int kw, double core_mass, double *He_core_mass, double *CO_core_mass, double *Ne_core_mass);

double compute_spin_angular_momentum_from_spin_frequency(double spin_frequency, int stellar_type, int object_type, double mass, double core_mass, double radius, double core_radius, double k2, double k3);
double compute_spin_frequency_from_spin_angular_momentum(double spin_angular_momentum, int stellar_type, int object_type, double mass, double core_mass, double radius, double core_radius, double k2, double k3);

double compute_breakup_angular_frequency(double mass, double radius);

double determine_sse_compact_object_radius_RSun(int kw, double m);

void compute_NS_formation_properties_Ye19_model(bool merger_event, double *ospin, double *B_G);
bool determine_if_NS_is_MSP(double ospin, double B_G);
double compute_NS_magnetic_field_Ye19_model(double B_0, double T, double t_acc, double Delta_m);
void ODE_handle_NS_properties_Ye19_model(Particle *p, double dt);
void ODE_Ye19_model_update_magnetic_field(ParticlesMap *particlesMap, double dt);
void nbody_handle_NS_properties_Ye19_model(Particle *p, double ospin, double dt);
}
