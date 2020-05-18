#include "types.h"
extern "C"
{

int handle_binary_evolution(ParticlesMap *particlesMap, double t_old, double t, double *dt_binary_evolution, int *integration_flag);

int handle_mass_transfer(ParticlesMap *particlesMap, double t_old, double t, double *dt_binary_evolution, int *integration_flag);
int handle_mass_transfer_cases(ParticlesMap *particlesMap, int parent_index, int donor_index, int accretor_index, int *integration_flag, double t_old, double t, double *dt_binary_evolution);
double compute_q_crit_for_common_envelope_evolution(int kw, double mass, double core_mass);
int determine_mass_transfer_stability(Particle *parent, Particle *donor, Particle *accretor);
double compute_orbit_averaged_mass_transfer_rate_emt_model(double M_donor, double fm, double P_orb);
double compute_stellar_dynamical_timescale(double M, double R);


int dynamical_mass_transfer_low_mass_donor(ParticlesMap *particlesMap, int parent_index, int donor_index, int accretor_index, double t_old, double t, int *integration_flag);
int dynamical_mass_transfer_WD_donor(ParticlesMap *particlesMap, int parent_index, int donor_index, int accretor_index, double t_old, double t, int *integration_flag);
int mass_transfer_NS_donor(ParticlesMap *particlesMap, int parent_index, int donor_index, int accretor_index, double t_old, double t, int *integration_flag);
int mass_transfer_BH_donor(ParticlesMap *particlesMap, int parent_index, int donor_index, int accretor_index, double t_old, double t, int *integration_flag);

int common_envelope_evolution(ParticlesMap *particlesMap, int binary_index, int donor_index, int accretor_index, int *integration_flag);
void set_old_parameters_for_adiabatic_mass_loss(ParticlesMap *particlesMap);
void compute_new_orbits_assuming_adiabatic_mass_loss(ParticlesMap *particlesMap, double mass_loss_timescale);

int stable_mass_transfer_evolution(ParticlesMap *particlesMap, int parent_index, int donor_index, int accretor_index, double t_old, double t, int *integration_flag);

}
