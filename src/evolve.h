#include <math.h>
#include <cstdlib>
#include <map>
#include <algorithm>
#include <time.h>

#include "types.h"
#include "../interface.h"

#include "cvode/cvode.h"					/* prototypes for CVODE fcts., consts. */
#include "cvode/nvector_serial.h"			/* serial N_Vector types, fcts., macros */
#include "cvode/cvode_dense.h"				/* prototype for CVDense */
#include "cvode/sundials_dense.h"			/* definitions DlsMat DENSE_ELEM */
#include "cvode/sundials_types.h"			/* definition of type realtype */


#include "structure.h"
#include "ODE_system.h"
#include "ODE_root_finding.h"
#include "ODE_newtonian.h"
#include "ODE_postnewtonian.h"
#include "ODE_tides.h"
#include "ODE_mass_changes.h"
#include "ODE_VRR.h"
#include "external.h"
#include "stellar_evolution.h"
#include "SNe.h"
#include "flybys.h"
#include "tools.h"
#include "binary_evolution.h"
#include "nbody_evolution.h"
#include "collision.h"
#include "test.h"
#include "apsidal_motion_constant.h"
#include "common_envelope_evolution.h"
#include "mstar/regularization.h"
#include "logging.h"

extern "C"
{
int initialize_code(ParticlesMap *particlesMap);
void custom_segfault_handler(int s);
int evolve(ParticlesMap *particlesMap, double start_time, double end_time, double *output_time, double *hamiltonian, int *state, int *CVODE_flag, int *CVODE_error_code, int *integration_flag);
}
