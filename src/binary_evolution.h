#include "types.h"
extern "C"
{

int handle_binary_evolution(ParticlesMap *particlesMap, double t_old, double t, double *dt_binary_evolution, int *integration_flag);

int handle_mass_transfer(ParticlesMap *particlesMap, double t_old, double t, double *dt_binary_evolution, int *integration_flag);
int handle_mass_transfer_cases(ParticlesMap *particlesMap, int parent_index, int donor_index, int accretor_index, int *integration_flag, double t_old, double t, double *dt_binary_evolution);
double compute_q_crit_for_common_envelope_evolution(int kw, double mass, double core_mass);
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

double compute_bse_mass_transfer_amount_averaged(int kw1, double m_donor, double core_mass_donor, double R_donor, double m_accretor, double a, double e, double dt, double t_dyn_donor, double t_KH_donor);
double compute_bse_mass_transfer_amount(int kw1, double m_donor, double core_mass_donor, double R_donor, double R_RL_donor, double dt, double t_dyn_donor, double t_KH_donor);

void handle_mass_accretion_events_with_degenerate_objects(ParticlesMap *particlesMap, double t_old, double t, int *integration_flag, double *dt_binary_evolution);
void reset_ODE_mass_dot_quantities(Particle *p);

void white_dwarf_helium_mass_accumulation_efficiency(double m_WD, double m_dot, double WD_luminosity, double *eta, int *WD_accretion_mode);
void white_dwarf_helium_mass_accumulation_efficiency_KH04(double m_WD, double m_dot, double *eta, double *m_dot_KH_min, double *m_dot_KH_max);
double determine_WK11_max_accretion_rate(double M_WD, double WD_luminosity);
double determine_WK11_max_accretion_rate_for_given_luminosity(double M_WD, double WD_luminosity, const double (*SD_table)[SINGLE_DEGENERATE_WK_DATA_MAX_ACCRETION_RATE_TABLEWIDTH], int SD_table_length);
bool determine_if_He_accreting_WD_explodes(double M_WD, double M_dot_acc, double M_He, double WD_luminosity);
double determine_if_He_accreting_WD_explodes_for_given_luminosity(double M_WD, double M_dot_acc, double M_He, const double (*SD_table)[SINGLE_DEGENERATE_WK_DATA_TABLEWIDTH], int SD_table_length);
double determine_if_He_accreting_WD_explodes_for_given_luminosity_and_WD_mass(double M_He, double M_dot_acc_1e8, double acc_data[SINGLE_DEGENERATE_WK_DATA_LONGEST_TABLELENGTH_FIXED_M_WD], double He_data[SINGLE_DEGENERATE_WK_DATA_LONGEST_TABLELENGTH_FIXED_M_WD], int data_length);

void white_dwarf_helium_mass_accumulation_efficiency_KH04_index0(double log_m_dot, double *eta, double *log_m_dot_min, double *log_m_dot_max);
void white_dwarf_helium_mass_accumulation_efficiency_KH04_index1(double log_m_dot, double *eta, double *log_m_dot_min, double *log_m_dot_max);
void white_dwarf_helium_mass_accumulation_efficiency_KH04_index2(double log_m_dot, double *eta, double *log_m_dot_min, double *log_m_dot_max);
void white_dwarf_helium_mass_accumulation_efficiency_KH04_index3(double log_m_dot, double *eta, double *log_m_dot_min, double *log_m_dot_max);
void white_dwarf_helium_mass_accumulation_efficiency_KH04_index4(double log_m_dot, double *eta, double *log_m_dot_min, double *log_m_dot_max);
void white_dwarf_helium_mass_accumulation_efficiency_KH04_index5(double log_m_dot, double *eta, double *log_m_dot_min, double *log_m_dot_max);

void white_dwarf_hydrogen_accretion_boundaries_WBBP13(double m_WD, double *m_dot_lower, double *m_dot_upper);

}
