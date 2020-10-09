#include "types.h"
extern "C"
{

int handle_binary_evolution(ParticlesMap *particlesMap, double t_old, double t, double *dt_binary_evolution, int *integration_flag);

int handle_mass_transfer(ParticlesMap *particlesMap, double t_old, double t, double *dt_binary_evolution, int *integration_flag);
int handle_mass_transfer_cases(ParticlesMap *particlesMap, int parent_index, int donor_index, int accretor_index, int *integration_flag, double t_old, double t, double *dt_binary_evolution);
double compute_q_crit_for_common_envelope_evolution(int kw, double mass, double core_mass);
//int determine_mass_transfer_stability(Particle *parent, Particle *donor, Particle *accretor);
double compute_orbit_averaged_mass_transfer_rate_emt_model(double M_donor, double fm, double P_orb);
double compute_stellar_dynamical_timescale(double M, double R);

int dynamical_mass_transfer_low_mass_donor(ParticlesMap *particlesMap, int parent_index, int donor_index, int accretor_index, double t_old, double t, int *integration_flag);
double compute_Kelvin_Helmholtz_timescale(int kw, double mass, double core_mass, double radius, double luminosity);
double compute_Eddington_accretion_rate(double radius, double hydrogen_mass_fraction);

int dynamical_mass_transfer_WD_donor(ParticlesMap *particlesMap, int parent_index, int donor_index, int accretor_index, double t_old, double t, int *integration_flag);
int mass_transfer_NS_BH_donor(ParticlesMap *particlesMap, int parent_index, int donor_index, int accretor_index, double t_old, double t, int *integration_flag);

int handle_wind_accretion(ParticlesMap *particlesMap, double t_old, double t, double *dt_binary_evolution, int *integration_flag);

int stable_mass_transfer_evolution(ParticlesMap *particlesMap, int parent_index, int donor_index, int accretor_index, double t_old, double t, int *integration_flag, double *dt_binary_evolution);
int binary_stable_mass_transfer_evolution(ParticlesMap *particlesMap, int parent_index, int donor_index, int accretor_index, double t_old, double t, int *integration_flag, double *dt_binary_evolution);
int triple_stable_mass_transfer_evolution(ParticlesMap *particlesMap, int parent_index, int donor_index, int accretor_index, double t_old, double t, int *integration_flag, double *dt_binary_evolution);

double compute_bse_mass_transfer_amount(int kw1, double m_donor, double core_mass_donor, double R_donor, double R_RL_av_donor, double dt, double t_dyn_donor, double t_KH_donor);

}
