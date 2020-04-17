#include <math.h>
#include <cstdlib>
#include <map>

#include "types.h"
#include "../interface.h"

#include "cvode/cvode.h"					/* prototypes for CVODE fcts., consts. */
#include "cvode/nvector_serial.h"			/* serial N_Vector types, fcts., macros */
#include "cvode/cvode_dense.h"				/* prototype for CVDense */
#include "cvode/sundials_dense.h"			/* definitions DlsMat DENSE_ELEM */
#include "cvode/sundials_types.h"			/* definition of type realtype */


#include "structure.h"
#include "ODE_system.h"
#include "root_finding.h"
#include "newtonian.h"
#include "postnewtonian.h"
#include "tides.h"
#include "external.h"
#include "VRR.h"
#include "stellar_evolution.h"
#include "SNe.h"
#include "flybys.h"
#include "tools.h"
#include "mass_changes.h"
#include "binary_evolution.h"
#include "nbody_evolution.h"
#include "merger.h"

extern "C"
{
int evolve(ParticlesMap *particlesMap, double start_time, double end_time, double *output_time, double *hamiltonian, int *state, int *CVODE_flag, int *CVODE_error_code, int *integration_flag);
int integrate_ODE_system(ParticlesMap *particlesMap, double start_time, double end_time, double *output_time, double *hamiltonian, int *output_flag, int *error_code);

static int check_flag(void *flagvalue, char *funcname, int opt);
}
