#include "types.h"
extern "C"
{

int determine_merger_type(int kw1, int kw2);

void handle_collisions(ParticlesMap *particlesMap, int *integration_flag);
void collision_product(ParticlesMap *particlesMap, int binary_index, int child1_index, int child2_index, int *integration_flag);

}
