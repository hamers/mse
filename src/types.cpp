/* MSE */

#include "types.h"

extern "C"
{

int highest_particle_index = 0;
ParticlesMap particlesMap;

// Default constants //
double CONST_G = 4.0*M_PI*M_PI; 
double CONST_G_P2 = CONST_G*CONST_G;
double CONST_G_P3 = CONST_G_P2*CONST_G;
double CONST_C_LIGHT = 63239.72638679138;
double CONST_C_LIGHT_P2 = CONST_C_LIGHT*CONST_C_LIGHT;
double CONST_C_LIGHT_P3 = CONST_C_LIGHT_P2*CONST_C_LIGHT;
double CONST_C_LIGHT_P4 = CONST_C_LIGHT_P3*CONST_C_LIGHT;
double CONST_C_LIGHT_P5 = CONST_C_LIGHT_P4*CONST_C_LIGHT;
double CONST_C_LIGHT_P6 = CONST_C_LIGHT_P5*CONST_C_LIGHT;
double CONST_MSUN = 1.0;
double CONST_R_SUN = 0.004649130343817401;
double CONST_L_SUN = 0.0002710404109745588;
double CONST_KM_PER_S = 0.210862;
double CONST_PER_PC3 = 1.14059e-16;


// Default parameters //
double relative_tolerance = 1.0e-12;
double absolute_tolerance_eccentricity_vectors = 1.0e-10;
bool include_quadrupole_order_terms = true;
bool include_octupole_order_binary_pair_terms = true;
bool include_octupole_order_binary_triplet_terms = true;
bool include_hexadecupole_order_binary_pair_terms = true;
bool include_dotriacontupole_order_binary_pair_terms = true;
bool include_double_averaging_corrections = false;
int random_seed = 0;

double epsilon = 1.0e-15; /* used for tiny numbers close to machine precision */

bool include_flybys = true;
bool flybys_correct_for_gravitational_focussing = true;
int flybys_velocity_distribution = 0;
int flybys_mass_distribution = 0;
int flybys_N_enc = 0;
int flybys_N_not_impulsive = 0;
int flybys_reference_binary = -1;
double flybys_t_next_encounter = 0.0;
double flybys_mass_distribution_lower_value = 0.1;
double flybys_mass_distribution_upper_value = 100.0;
double flybys_encounter_sphere_radius = 1.0e5;
double flybys_stellar_density = 0.1*CONST_PER_PC3; /* density at infinity */
double flybys_stellar_relative_velocity_dispersion = 30.0*CONST_KM_PER_S;
double flybys_W_max = 0.0;
double flybys_total_encounter_rate_at_R_enc = 0.0;
double flybys_stellar_density_at_R_enc = 0.0; /* density at infinity */
double flybys_internal_mass = 0.0;
double flybys_internal_semimajor_axis = 0.0;

int binary_evolution_CE_energy_flag = 0;
int binary_evolution_CE_spin_flag = 1;

double chandrasekhar_mass = 1.44;
double eddington_accretion_factor = 10.0;
double nova_accretion_factor = 1.0e-3;
double alpha_wind_accretion = 1.5;
double beta_wind_accretion = 0.125;

double mstar_gbs_tolerance = 1.0e-10;
double mstar_collision_tolerance = 1.0e-10;

double triple_mass_transfer_primary_star_accretion_efficiency_no_disk = 0.1;
double triple_mass_transfer_secondary_star_accretion_efficiency_no_disk = 0.1;
double triple_mass_transfer_primary_star_accretion_efficiency_disk = 0.9;
double triple_mass_transfer_secondary_star_accretion_efficiency_disk = 0.9;
double triple_mass_transfer_inner_binary_alpha_times_lambda = 5.0;

}
