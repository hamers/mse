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
double compute_Kelvin_Helmholtz_timescale(int kw, double mass, double core_mass, double radius, double luminosity);
double compute_Eddington_accretion_rate(double radius, double hydrogen_mass_fraction);

int dynamical_mass_transfer_WD_donor(ParticlesMap *particlesMap, int parent_index, int donor_index, int accretor_index, double t_old, double t, int *integration_flag);
int mass_transfer_NS_BH_donor(ParticlesMap *particlesMap, int parent_index, int donor_index, int accretor_index, double t_old, double t, int *integration_flag);

int common_envelope_evolution(ParticlesMap *particlesMap, int binary_index, int index1, int index2, double t, int *integration_flag);

void handle_instantaneous_and_adiabatic_mass_changes_in_orbit(ParticlesMap *particlesMap, Particle *star1, Particle *star2, double Delta_m1, double Delta_m2, double mass_loss_timescale, int *integration_flag);
void set_old_parameters_for_adiabatic_mass_loss(ParticlesMap *particlesMap);
void compute_new_orbits_assuming_adiabatic_mass_loss(ParticlesMap *particlesMap, double mass_loss_timescale);
void compute_new_positions_and_velocities_given_new_semimajor_axis_and_eccentricity(double M1_old, double R1_vec_old[3], double V1_vec_old[3], double M2_old, double R2_vec_old[3], double V2_vec_old[3], double M1, double R1_vec[3], double V1_vec[3], double M2, double R2_vec[3], double V2_vec[3], double a, double e);

int handle_wind_accretion(ParticlesMap *particlesMap, double t_old, double t, double *dt_binary_evolution, int *integration_flag);

int stable_mass_transfer_evolution(ParticlesMap *particlesMap, int parent_index, int donor_index, int accretor_index, double t_old, double t, int *integration_flag);

}
