#include "types.h"
extern "C"
{

int determine_merger_type(int kw1, int kw2);

void handle_collisions(ParticlesMap *particlesMap, double t, int *integration_flag);

void handle_destruction_of_binary_in_system(ParticlesMap *particlesMap, Particle *b);
void collision_product(ParticlesMap *particlesMap, int binary_index, int child1_index, int child2_index, double t, int *integration_flag);
void collision_product_star_planet(ParticlesMap *particlesMap, int binary_index, int child1_index, int child2_index, double t, int *integration_flag);
void collision_product_planet_planet(ParticlesMap *particlesMap, int binary_index, int child1_index, int child2_index, double t, int *integration_flag);

void determine_compact_object_merger_properties(double m1_, double m2_, double chi1_, double chi2_, double spin_vec_1_unit_[3], double spin_vec_2_unit_[3], double h_vec_unit[3], double e_vec_unit[3], double v_recoil_vec[3], double alpha_vec_final[3], double *M_final);

double determine_effective_radius_for_collision(double radius, int stellar_type, int integration_flag);
}
