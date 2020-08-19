#include "types.h"

extern "C"
{
void compute_EOM_Newtonian_for_particle(ParticlesMap *particlesMap, Particle *p, double *hamiltonian, double *KS_V, bool compute_hamiltonian_only);

void compute_EOM_binary_pairs(ParticlesMap *particlesMap, int inner_binary_index, int outer_binary_index, int connecting_child_in_outer_binary, double *hamiltonian, double *KS_V, bool compute_hamiltonian_only);

void compute_EOM_binary_pairs_non_averaged(ParticlesMap *particlesMap, int inner_binary_index, int outer_binary_index, int child_connected_to_p_in_parent_k, double *hamiltonian, double *KS_V, bool compute_hamiltonian_only);
void get_gravity_binary_pairs_order_n(
    Particle *p, Particle *k, \
    int int_n, double n, double C_n, \
    double *hamiltonian, double *KS_V, \
    double mu_p, double mu_k, \
    double *r_p_vec, double *r_k_vec, \
    double r_p_P1, double r_k_P1, double r_p_P2, double r_k_P2, \
    double r_p_div_r_k_P2, double r_p_div_r_k_Pn, double r_p_vec_dot_r_k_vec);

void compute_EOM_binary_pairs_single_averaged(ParticlesMap *particlesMap, int inner_binary_index, int outer_binary_index, int child_connected_to_p_in_parent_k, double *hamiltonian, double *KS_V, bool compute_hamiltonian_only);
void compute_EOM_binary_pairs_double_averaged(ParticlesMap *particlesMap, int inner_binary_index, int outer_binary_index, int connecting_child_in_outer_binary, double *hamiltonian, double *KS_V, bool compute_hamiltonian_only);

void compute_EOM_binary_triplets(ParticlesMap *particlesMap, int binary_A_index, int binary_B_index, int binary_C_index, int connecting_child_in_binary_B_to_binary_A, int connecting_child_in_binary_C_to_binary_B, double *hamiltonian, double *KS_V, bool compute_hamiltonian_only);

double compute_order_of_magnitude_secular_timescale_for_pair(ParticlesMap *particlesMap, int inner_binary_index, int outer_binary_index, int connecting_child_in_outer_binary);

}
