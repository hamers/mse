#include "types.h"
//#include "mstar/regularization.h"
extern "C"
{

void integrate_nbody_system(ParticlesMap *particlesMap, int *integration_flag, double t_old, double t, double *dt_nbody);
double determine_nbody_timestep(ParticlesMap *particlesMap, int integration_flag, double P_orb_min, double P_orb_max);

void update_stellar_evolution_quantities_during_nbody_integration(ParticlesMap *particlesMap, double dt);

int handle_dynamical_instability(ParticlesMap *particlesMap);
//void create_mstar_instance_of_system(ParticlesMap *particlesMap, struct RegularizedRegion *R);
struct RegularizedRegion *create_mstar_instance_of_system(ParticlesMap *particlesMap);

void analyze_mstar_system(struct RegularizedRegion *R, bool *stable_system, ParticlesMap *particlesMap, double *P_orb_min, double *P_orb_max,double dt);
void copy_bodies_from_old_to_new_particlesMap(ParticlesMap *old_particlesMap, ParticlesMap *new_particlesMap);
void find_binaries_in_system(ParticlesMap *particlesMap, double *P_orb_min, double *P_orb_max);
void get_all_semimajor_axes_in_system(ParticlesMap *particlesMap, std::vector<double> *semimajor_axes, std::vector<int> *binary_indices);
//void analyze_mstar_system2(struct RegularizedRegion *R);
void extract_pos_vel_from_mstar_system(struct RegularizedRegion *R, ParticlesMap *particlesMap);
void print_state(struct RegularizedRegion *R);

/* Functions below are from MSTAR that need to be accessible from MSE */
void allocate_armst_structs(struct RegularizedRegion **R, int MaxNumPart);
void initialize_mpi_or_serial(void);
void run_integrator(struct RegularizedRegion *R, double time_interval);
void free_data(struct RegularizedRegion *R);
}
