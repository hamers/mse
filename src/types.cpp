#include "types.h"

extern "C"
{

// Default constants //
double CONST_G = 4.0*M_PI*M_PI; 
double CONST_G_P2 = CONST_G*CONST_G;
double CONST_G_P3 = CONST_G_P2*CONST_G;
double CONST_C_LIGHT = 63239.72638679138;
double CONST_C_LIGHT_P2 = CONST_C_LIGHT*CONST_C_LIGHT;
double CONST_C_LIGHT_P4 = CONST_C_LIGHT_P2*CONST_C_LIGHT_P2;
double CONST_C_LIGHT_P5 = CONST_C_LIGHT_P4*CONST_C_LIGHT;
double CONST_MSUN = 1.0;
double CONST_R_SUN = 0.004649130343817401;
double CONST_L_SUN = 0.0002710404109745588;
double CONST_KM_PER_S = 0.210862;
double CONST_PER_PC3 = 1.14059e-16;

// Default parameters //
double relative_tolerance = 1.0e-12;
double absolute_tolerance_eccentricity_vectors = 1.0e-10;
bool include_quadrupole_order_terms = true;
bool include_octupole_order_binary_pair_terms = true;
bool include_octupole_order_binary_triplet_terms = true;
bool include_hexadecupole_order_binary_pair_terms = true;
bool include_dotriacontupole_order_binary_pair_terms = true;
int random_seed = 0;

bool include_flybys = true;
bool flybys_correct_for_gravitational_focussing = true;
int flybys_velocity_distribution = 0;
int flybys_mass_distribution = 0;
int flybys_N_enc = 0;
int flybys_N_not_impulsive = 0;
int flybys_reference_binary = -1;
double flybys_t_next_encounter = 0.0;
double flybys_mass_distribution_lower_value = 0.1;
double flybys_mass_distribution_upper_value = 100.0;
double flybys_encounter_sphere_radius = 1.0e5;
double flybys_stellar_density = 0.1*CONST_PER_PC3; /* density at infinity */
double flybys_stellar_relative_velocity_dispersion = 30.0*CONST_KM_PER_S;
double flybys_W_max = 0.0;
double flybys_total_encounter_rate_at_R_enc = 0.0;
double flybys_stellar_density_at_R_enc = 0.0; /* density at infinity */
double flybys_internal_mass = 0.0;
double flybys_internal_semimajor_axis = 0.0;

int highest_particle_index = 0;
ParticlesMap particlesMap;

#ifdef IGNORE
void Particle::set_ODE_quantities(double delta_time)
{
    for (int i=0; i<3; i++)
    {        
        dspin_vec_dt[i] = 0.0;
        de_vec_dt[i] = 0.0;
        dh_vec_dt[i] = 0.0;
    }
    
    e = norm3(e_vec);
    h = norm3(h_vec);
    spin_vec_norm = norm3(spin_vec);

    for (int i=0; i<3; i++)
    {        
        e_vec_unit[i] = e_vec[i]/e;
        h_vec_unit[i] = h_vec[i]/h;
    }
    
    e_p2 = e*e;
    j_p2 = 1.0 - e_p2;
    j = sqrt(j_p2);
    j_p3 = j*j_p2;
    j_p4 = j*j_p3;
    j_p5 = j*j_p4;

    if (is_binary == 1)
    {
        a = h*h*child1_mass_plus_child2_mass/( CONST_G*child1_mass_times_child2_mass*child1_mass_times_child2_mass*j_p2 );
    }
    else
    {
        a = 0.0;
    }

    dmass_dt = mass_dot;
    //printf("dmass_dt %g %g\n",dmass_dt,radius_dot);
    dradius_dt = radius_dot + radius_ddot*delta_time;
    //printf("ODE dR/dt %g\n",dradius_dt);
    //printf("types.cpp -- index %d a %g h %g child1_mass_plus_child2_mass %g\n",index,a,h,child1_mass_plus_child2_mass);

}

void Particle::reset_ODE_quantities()
{
    for (int i=0; i<3; i++)
    {        
        dspin_vec_dt[i] = 0.0;
        
        de_vec_dt[i] = 0.0;
        dh_vec_dt[i] = 0.0;
    }
}
#endif

}
