#include "types.h"
extern "C"
{

double compute_apsidal_motion_constant(Particle *star);

double AMC_data_function(double log_m, Particle *star);

void limit_tau(double *tau);

}
