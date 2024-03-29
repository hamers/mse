#include "types.h"

extern "C"
{

void integrate_nbody_system(ParticlesMap *particlesMap, int *integration_flag, double t_old, double t, double *t_out, double *dt_nbody);

int determine_orbits_in_system_using_nbody(ParticlesMap *particlesMap);
int determine_new_integration_flag_using_nbody(struct RegularizedRegion *R, ParticlesMap *particlesMap, double *P_orb_min, double *P_orb_max);
void handle_collisions_nbody(struct RegularizedRegion *R, ParticlesMap *particlesMap, double t, int *integration_flag);

double determine_nbody_timestep(ParticlesMap *particlesMap, int integration_flag, double P_orb_min, double P_orb_max);

void update_stellar_evolution_quantities_directly(ParticlesMap *particlesMap, double t, double dt);

struct RegularizedRegion *create_mstar_instance_of_system(ParticlesMap *particlesMap, int integration_flag);

void analyze_mstar_system(struct RegularizedRegion *R, bool *stable_system, ParticlesMap *particlesMap, double *P_orb_min, double *P_orb_max,double dt);
void copy_bodies_from_old_to_new_particlesMap(ParticlesMap *old_particlesMap, ParticlesMap *new_particlesMap);
void find_binaries_in_system(ParticlesMap *particlesMap, double *P_orb_min, double *P_orb_max);
void get_all_semimajor_axes_in_system(ParticlesMap *particlesMap, std::vector<double> *semimajor_axes, std::vector<int> *binary_indices);

void extract_pos_vel_from_mstar_system_and_reset_particlesMap(struct RegularizedRegion *R, ParticlesMap *particlesMap);
void update_pos_vel_from_mstar_system(struct RegularizedRegion *R, ParticlesMap *particlesMap);
void mstar_print_state(struct RegularizedRegion *R);
double compute_nbody_total_energy(struct RegularizedRegion *R);

void integrate_nbody_system_with_mass_loss(double end_time, int Nsteps, std::vector<double> &masses, std::vector<double> &delta_masses, std::vector<std::vector<double> > &R_vecs, std::vector<std::vector<double> > &V_vecs);

/* Functions below are from MSTAR that need to be accessible from MSE */
void allocate_armst_structs(struct RegularizedRegion **R, int MaxNumPart);
void initialize_mpi_or_serial(void);
void run_integrator(struct RegularizedRegion *R, double time_interval, double *end_time, int *collision_occurred);
void free_data(struct RegularizedRegion *R);
void into_CoM_frame(struct RegularizedRegion *R);

int check_particlesMap_for_inclusion_in_MSTAR(ParticlesMap *particlesMap);

}
