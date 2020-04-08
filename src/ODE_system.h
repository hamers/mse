#include "types.h"

extern "C"
{
void initialize_direct_integration_quantities(ParticlesMap *particlesMap);
void process_direct_integration_quantities(ParticlesMap *particlesMap, double delta_time);
void compute_KS_EOM(ParticlesMap *particlesMap, double KS_V);

void extract_ODE_variables(ParticlesMap *particlesMap, N_Vector &y, double delta_time);
void write_ODE_variables_dots(ParticlesMap *particlesMap, N_Vector &y_dot);
void reset_ODE_dots(ParticlesMap *particlesMap, N_Vector &y, double delta_time);
void set_up_derived_ODE_quantities(ParticlesMap *particlesMap);

int compute_y_dot(realtype time, N_Vector y, N_Vector y_dot, void *data_);

void set_initial_ODE_variables(ParticlesMap *particlesMap, N_Vector &y, N_Vector &y_abs_tol,double abs_tol_spin_vec, double abs_tol_e_vec, double abs_tol_h_vec);
void extract_final_ODE_variables(ParticlesMap *particlesMap, N_Vector &y_out);
}
