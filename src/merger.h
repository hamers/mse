#include "types.h"
extern "C"
{

int determine_merger_type(int kw1, int kw2);

void handle_collisions(ParticlesMap *particlesMap, int *integration_flag);
void collision_product(ParticlesMap *particlesMap, int binary_index, int child1_index, int child2_index, int *integration_flag);

void apply_instantaneous_mass_changes_and_kicks(ParticlesMap *particlesMap, int *integration_flag);
void determine_compact_object_merger_properties(double m1, double m2, double chi1, double chi2, double spin_vec_1_unit[3], double spin_vec_2_unit[3], double h_vec_unit[3], double e_vec_unit[3], double v_recoil_vec[3]);
}
