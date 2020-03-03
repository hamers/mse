#include "types.h"
extern "C"
{

int sample_from_3d_maxwellian_distribution(double sigma, double v[3]);
double sample_from_y_times_maxwellian_distribution(double sigma);
int sample_spherical_coordinates_unit_vectors_from_isotropic_distribution(double r_hat_vec[3], double theta_hat_vec[3], double phi_hat_vec[3]);
double sample_from_power_law_distribution(double alpha, double y_lower, double y_upper);

int compute_h_tot_vector(ParticlesMap* particlesMap, double h_tot_vec[3]);
int compute_orbital_vectors_from_orbital_elements(double child1_mass, double child2_mass, double semimajor_axis, double eccentricity, double inclination, double argument_of_pericenter,double longitude_of_ascending_node, double *e_vec_x, double *e_vec_y, double *e_vec_z, double *h_vec_x, double *h_vec_y, double *h_vec_z);
int compute_orbital_vectors_from_orbital_elements_unit(double inclination, double argument_of_pericenter,double longitude_of_ascending_node, double *e_hat_vec_x, double *e_hat_vec_y, double *e_hat_vec_z, double *h_hat_vec_x, double *h_hat_vec_y, double *h_hat_vec_z);
int compute_orbital_elements_from_orbital_vectors(double child1_mass, double child2_mass, double h_tot_vec[3], double e_vec_x, double e_vec_y, double e_vec_z, double h_vec_x, double h_vec_y, double h_vec_z, double *semimajor_axis, double *eccentricity, double *inclination, double *argument_of_pericenter,double *longitude_of_ascending_node);
void compute_eccentric_anomaly_from_mean_anomaly(double mean_anomaly, double eccentricity, double *cos_eccentric_anomaly, double *sin_eccentric_anomaly);
void compute_true_anomaly_from_eccentric_anomaly(double cos_eccentric_anomaly, double sin_eccentric_anomaly, double eccentricity, double *cos_true_anomaly, double *sin_true_anomaly);
double compute_true_anomaly_from_mean_anomaly(double mean_anomaly, double eccentricity);
double sample_random_true_anomaly(double eccentricity);
void from_orbital_vectors_to_cartesian(double child1_mass, double child2_mass, double e_vec[3], double h_vec[3], double true_anomaly, double r[3], double v[3]);
void from_cartesian_to_orbital_vectors(double child1_mass, double child2_mass, double r[3], double v[3], double e_vec[3], double h_vec[3]);

}
