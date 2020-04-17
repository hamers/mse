/* MSE */
/* Adrian Hamers November 2019 */

#include "evolve.h"
#include "merger.h"

extern "C"
{

void handle_collisions(ParticlesMap *particlesMap)
{
    /* invoke CE if one star has k in 2,3,4,5,6,8,9 */
}

int determine_merger_type(int kw1, int kw2)
{
    return MERGER_TABLE[kw1][kw2];
}





}
