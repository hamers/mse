#include "types.h"
extern "C"
{
double compute_orbital_period_from_semimajor_axis(double M, double a);
double compute_semimajor_axis_from_orbital_period(double M, double P);

double compute_spin_angular_frequency_from_spin_period(double P);
double compute_spin_period_from_spin_angular_frequency(double Omega);

double generate_random_number_between_zero_and_unity();
int sample_from_3d_maxwellian_distribution(double sigma, double v[3]);
double sample_from_y_times_maxwellian_distribution(double sigma);
double sample_from_normal_distribution(double mu, double sigma);
int sample_spherical_coordinates_unit_vectors_from_isotropic_distribution(double r_hat_vec[3], double theta_hat_vec[3], double phi_hat_vec[3]);
double sample_from_power_law_distribution(double alpha, double y_lower, double y_upper);
double sample_from_Kroupa_93_imf();

int compute_h_tot_vector(ParticlesMap* particlesMap, double h_tot_vec[3]);
int compute_orbital_vectors_from_orbital_elements(double child1_mass, double child2_mass, double semimajor_axis, double eccentricity, double inclination, double argument_of_pericenter,double longitude_of_ascending_node, double *e_vec_x, double *e_vec_y, double *e_vec_z, double *h_vec_x, double *h_vec_y, double *h_vec_z);
int compute_orbital_vectors_from_orbital_elements_unit(double inclination, double argument_of_pericenter,double longitude_of_ascending_node, double *e_hat_vec_x, double *e_hat_vec_y, double *e_hat_vec_z, double *h_hat_vec_x, double *h_hat_vec_y, double *h_hat_vec_z);
int compute_orbital_elements_from_orbital_vectors(double child1_mass, double child2_mass, double h_tot_vec[3], double e_vec_x, double e_vec_y, double e_vec_z, double h_vec_x, double h_vec_y, double h_vec_z, double *semimajor_axis, double *eccentricity, double *inclination, double *argument_of_pericenter,double *longitude_of_ascending_node);
void compute_eccentric_anomaly_from_mean_anomaly(double mean_anomaly, double eccentricity, double *cos_eccentric_anomaly, double *sin_eccentric_anomaly);
void compute_true_anomaly_from_eccentric_anomaly(double cos_eccentric_anomaly, double sin_eccentric_anomaly, double eccentricity, double *cos_true_anomaly, double *sin_true_anomaly);
double compute_true_anomaly_from_mean_anomaly(double mean_anomaly, double eccentricity);
double compute_mean_anomaly_from_true_anomaly(double true_anomaly, double eccentricity);
void compute_true_anomaly_from_mean_anomaly_hyperbolic_orbit(double mean_anomaly, double eccentricity,double *cos_true_anomaly,double *sin_true_anomaly);
double sample_random_true_anomaly(double eccentricity);
void from_orbital_vectors_to_cartesian(double child1_mass, double child2_mass, double e_vec[3], double h_vec[3], double true_anomaly, double r[3], double v[3]);
void from_cartesian_to_orbital_vectors(double child1_mass, double child2_mass, double r[3], double v[3], double e_vec[3], double h_vec[3], double *true_anomaly);
void compute_semimajor_axis_and_eccentricity_from_orbital_vectors(double m1, double m2, double e_vec[3], double h_vec[3], double *semimajor_axis, double *eccentricity);

bool compute_hyperbolic_orbit_properties_from_position_and_velocity(double M_tot, double R_vec[3], double V_vec[3], double *e_per, double *Q_per, double e_vec_unit_per[3], double h_vec_unit_per[3]);

void get_position_and_velocity_vectors_from_particle(Particle *p, double r[3], double v[3]);
void set_position_and_velocity_vectors_in_particle(Particle *p,  double r[3], double v[3]);

void get_unit_vector(double vec[3], double vec_unit[3]);

void copy_particlesMap(ParticlesMap *source, ParticlesMap *target);
void copy_all_body_properties(Particle *source, Particle *target);

void create_nested_system(ParticlesMap &particlesMap, int N_bodies, double *masses, int *stellar_types, int *object_types, double *smas, double *es, double *TAs, double *INCLs, double *APs, double *LANs);
void create_2p2_quadruple_system(ParticlesMap &particlesMap, double *masses, int *stellar_types, int *object_types, double *smas, double *es, double *TAs, double *INCLs, double *APs, double *LANs);
void print_system(ParticlesMap *particlesMap, int integration_flag);

void get_parallel_and_perpendicular_vectors_and_components(double a_vec[3], double h_vec[3], double *a_par, double *a_perp, double a_perp_vec[3]);
double get_mutual_angle(double a_vec[3], double b_vec[3]);
int get_inclination_relative_to_parent(ParticlesMap *particlesMap, int index, double *inclination_relative_to_parent);

double compute_a_from_h(double m1, double m2, double h, double e);
double compute_h_from_a(double m1, double m2, double a, double e);
double compute_h_dot_div_h(double m1, double m1_dot, double m2, double m2_dot, double a, double a_dot, double e, double e_dot);
double compute_a_dot_div_a(double m1, double m1_dot, double m2, double m2_dot, double h, double h_dot, double e, double e_dot);

bool equal_number(double x1, double x2, double tol);
void check_number(double x, char *source, char *description, bool exit_on_error);

int clear_particles(ParticlesMap *particlesMap);

void MSTAR_compute_center_of_mass_position_and_velocity(struct RegularizedRegion *R, double R_cm[3], double V_cm[3]);
void compute_center_of_mass_position_and_velocity(ParticlesMap *particlesMap, double R_cm[3], double V_cm[3]);

double find_nearest_neighbor_separation(ParticlesMap *particlesMap, int primary_index, double primary_R_vec[3]);
void remove_binaries_from_system(ParticlesMap *particlesMap);

void rescale_vector(double v[3], double factor);

double interpolate_1D_data_table_linear(double x, const double (*data_table)[2], int table_length, bool use_constant_extrapolation);
}
