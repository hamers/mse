#include <stdio.h>
#include <math.h>
#include "interface.h"

//#include "src/types.h"
//#include "interface.h"
#include "src/evolve.h"


extern "C"
{
    

/*******************
/* basic interface *
 ******************/
 
 
int add_particle(int *index, int is_binary, int is_external)
{
    *index = highest_particle_index;

    Particle *p = new Particle(highest_particle_index, is_binary);
    particlesMap[highest_particle_index] = p;

    p->is_external = is_external;
    
    highest_particle_index++;
    return 0;
}


int delete_particle(int index)
{
    if (index > highest_particle_index)
    {
        return -1;
    }
  
    particlesMap.erase(index);

    return 0;
}

int set_children(int index, int child1, int child2)
{
    if (index > highest_particle_index)
    {
      return -1;
    }

    Particle *p = particlesMap[index];
    p->child1 = child1;
    p->child2 = child2;
    //printf("c1 %d c2 %d\n",child1,child2);
    return 0;
}
int get_children(int index, int *child1, int *child2)
{
    if (index > highest_particle_index)
    {
      return -1;
    }
  
    Particle * p = particlesMap[index];
    *child1 = p->child1;
    *child2 = p->child2;
    
    return 0;
}

int set_mass(int index, double mass)
{
    if (index > highest_particle_index)
    {
      return -1;
    }

    Particle *p = particlesMap[index];
    p->mass = mass;
    p->sse_initial_mass = mass;
    
    return 0;
}
int get_mass(int index, double *mass)
{
    if (index > highest_particle_index)
    {
      return -1;
    }
    Particle *p = particlesMap[index];
    *mass = p->mass;
    return 0;
}
//int set_mass_dot(int index, double mass_dot)
//{
//    if (index > highest_particle_index)
//    {
//      return -1;
//    }

//    Particle *p = particlesMap[index];
//    p->mass_dot = mass_dot;

//    return 0;
//}
int get_mass_dot(int index, double *mass_dot)
{
    if (index > highest_particle_index)
    {
      return -1;
    }

    Particle *p = particlesMap[index];
    *mass_dot = p->mass_dot_wind + p->mass_dot_RLOF;

    return 0;
}

int set_radius(int index, double radius, double radius_dot)
{
    if (index > highest_particle_index)
    {
      return -1;
    }

    Particle * p = particlesMap[index];
    p->radius = radius;
    p->radius_dot = radius_dot;
    
    return 0;
}
int get_radius(int index, double *radius, double *radius_dot)
{
    if (index > highest_particle_index)
    {
      return -1;
    }
  
    Particle * p = particlesMap[index];
    *radius = p->radius;
    *radius_dot = p->radius_dot;
    
    return 0;
}


/****************************
/* orbital vectors/elements *
 ****************************/

int set_orbital_vectors(int index, double e_vec_x, double e_vec_y, double e_vec_z, \
    double h_vec_x, double h_vec_y, double h_vec_z)
{
    if (index > highest_particle_index)
    {
      return -1;
    }
  
    Particle *p = particlesMap[index];
    p->e_vec_x = e_vec_x;
    p->e_vec_y = e_vec_y;
    p->e_vec_z = e_vec_z;
    p->h_vec_x = h_vec_x;
    p->h_vec_y = h_vec_y;
    p->h_vec_z = h_vec_z;
    //printf("ok %g %g %g %g %g %g\n",e_vec_x,e_vec_y,e_vec_z,h_vec_x,h_vec_y,h_vec_z);
    return 0;
}
int get_orbital_vectors(int index, double *e_vec_x, double *e_vec_y, double *e_vec_z, \
    double *h_vec_x, double *h_vec_y, double *h_vec_z)
{
    if (index > highest_particle_index)
    {
      return -1;
    }
  
    Particle *p = particlesMap[index];
    *e_vec_x = p->e_vec_x;
    *e_vec_y = p->e_vec_y;
    *e_vec_z = p->e_vec_z;
    *h_vec_x = p->h_vec_x;
    *h_vec_y = p->h_vec_y;
    *h_vec_z = p->h_vec_z;

    return 0;
}

int set_orbital_elements(int index, double semimajor_axis, double eccentricity, double true_anomaly, \
    double inclination, double argument_of_pericenter, double longitude_of_ascending_node, int sample_orbital_phase_randomly)
{

    if (index > highest_particle_index)
    {
        return -1;
    }

    Particle * p = particlesMap[index];    

    if (p->is_binary == false)
    {
        return 0;
    }

    /* determine masses in all binaries */
    int N_bodies, N_binaries, N_root_finding;
    determine_binary_parents_and_levels(&particlesMap, &N_bodies, &N_binaries, &N_root_finding);
    set_binary_masses_from_body_masses(&particlesMap);

    compute_orbital_vectors_from_orbital_elements(p->child1_mass, p->child2_mass, semimajor_axis, eccentricity, \
        inclination, argument_of_pericenter, longitude_of_ascending_node, \
        &(p->e_vec_x), &(p->e_vec_y), &(p->e_vec_z), &(p->h_vec_x), &(p->h_vec_y), &(p->h_vec_z) );

    p->true_anomaly = true_anomaly;
    p->sample_orbital_phase_randomly = sample_orbital_phase_randomly;
//    printf("soe a %g e %g TA %g I %g AP %g LAN %g SOPR %d\n",semimajor_axis,eccentricity,true_anomaly,inclination,argument_of_pericenter,longitude_of_ascending_node,sample_orbital_phase_randomly);
//    printf("set_orbital_elements %g %g %g %g %g %g\n",p->e_vec_x,p->e_vec_y,p->e_vec_z,p->h_vec_x,p->h_vec_y,p->h_vec_z);
    
    return 0;
}
int get_orbital_elements(int index, double *semimajor_axis, double *eccentricity, \
    double *inclination, double *argument_of_pericenter, double *longitude_of_ascending_node)
{

    if (index > highest_particle_index)
    {
        return -1;
    }
  
    Particle * p = particlesMap[index];    

    if (p->is_binary == false)
    {
        return 0;
    }

    double h_tot_vec[3];
    compute_h_tot_vector(&particlesMap,h_tot_vec);

    /* determine masses in all binaries */
    int N_bodies, N_binaries, N_root_finding;
    determine_binary_parents_and_levels(&particlesMap, &N_bodies, &N_binaries, &N_root_finding);
    set_binary_masses_from_body_masses(&particlesMap);
    
    compute_orbital_elements_from_orbital_vectors(p->child1_mass, p->child2_mass, h_tot_vec, \
        p->e_vec_x,p->e_vec_y,p->e_vec_z,p->h_vec_x,p->h_vec_y,p->h_vec_z,
        semimajor_axis, eccentricity, inclination, argument_of_pericenter, longitude_of_ascending_node);
    return 0;
}


int get_level(int index, int *value)
{
    if (index > highest_particle_index)
    {
      return -1;
    }
  
    Particle * p = particlesMap[index];
    *value = p->level;
    
    return 0;
}

        
int set_stellar_evolution_properties(int index, int stellar_type, int evolve_as_star, double sse_initial_mass, double metallicity, double sse_time_step, double epoch, double age, 
    double convective_envelope_mass, double convective_envelope_radius, double core_mass, double core_radius, double luminosity, double apsidal_motion_constant, double gyration_radius, double tides_viscous_time_scale)
{
    //printf("set_stellar_evolution_properties %g\n",metallicity);
    if (index > highest_particle_index)
    {
      return -1;
    }

    Particle * p = particlesMap[index];
    p->stellar_type = stellar_type;
    p->evolve_as_star = evolve_as_star;
    p->sse_initial_mass = sse_initial_mass;
    p->metallicity = metallicity;
    p->sse_time_step = sse_time_step;
    p->epoch = epoch;
    p->age = age;
    p->convective_envelope_mass = convective_envelope_mass;
    p->convective_envelope_radius = convective_envelope_radius;
    p->core_mass = core_mass;
    p->core_radius = core_radius;
    p->luminosity = luminosity;
    p->apsidal_motion_constant = apsidal_motion_constant;
    p->gyration_radius = gyration_radius;
    p->tides_viscous_time_scale = tides_viscous_time_scale;
    
    return 0;
}
int get_stellar_evolution_properties(int index, int *stellar_type, int *evolve_as_star, double *sse_initial_mass, double *metallicity, double *sse_time_step, double *epoch, double *age, 
    double *convective_envelope_mass, double *convective_envelope_radius, double *core_mass, double *core_radius, double *luminosity, double *apsidal_motion_constant, double *gyration_radius, double *tides_viscous_time_scale, double *roche_lobe_radius_pericenter)

{
    if (index > highest_particle_index)
    {
      return -1;
    }
  
    Particle * p = particlesMap[index];
    *stellar_type = p->stellar_type;
    *evolve_as_star = p->evolve_as_star;
    *sse_initial_mass = p->sse_initial_mass;
    *metallicity = p->metallicity;
    *sse_time_step = p->sse_time_step;
    *epoch = p->epoch;
    *age = p->age;
    *convective_envelope_mass = p->convective_envelope_mass;
    *convective_envelope_radius = p->convective_envelope_radius;
    *core_mass = p->core_mass;
    *core_radius = p->core_radius;
    *luminosity = p->luminosity;
    *apsidal_motion_constant = p->apsidal_motion_constant;
    *gyration_radius = p->gyration_radius;
    *tides_viscous_time_scale = p->tides_viscous_time_scale;
    
    if (p->is_binary == 0) /* TO DO: streamline with same code in root_finding.cpp */
    {
        Particle *parent = (particlesMap)[p->parent];
        Particle *sibling = (particlesMap)[p->sibling];
                    
        double a = parent->a;
        double e = parent->e;
        double rp = a*(1.0 - e);
        double subject_mass = p->mass;
        double companion_mass = sibling->mass;
    
        double spin_angular_frequency = p->spin_vec_norm;
        double orbital_angular_frequency_periapse = sqrt( CONST_G*(subject_mass + companion_mass)*(1.0 + e)/(rp*rp*rp) );
        double f = spin_angular_frequency/orbital_angular_frequency_periapse;
            
        double roche_radius_pericenter;
        if (p->check_for_RLOF_at_pericentre_use_sepinsky_fit == 0)
        {
            roche_radius_pericenter = roche_radius_pericenter_eggleton(rp, subject_mass/companion_mass);
        }
        else
        {
            roche_radius_pericenter = roche_radius_pericenter_sepinsky(rp, subject_mass/companion_mass, e, f);
        }
        *roche_lobe_radius_pericenter = roche_radius_pericenter;
    }
    
    return 0;
}

int set_kick_properties(int index, int kick_distribution, double kick_distribution_sigma)
{
    //printf("set_stellar_evolution_properties %g\n",metallicity);
    if (index > highest_particle_index)
    {
      return -1;
    }

    Particle * p = particlesMap[index];
    p->kick_distribution = kick_distribution;
    p->kick_distribution_sigma = kick_distribution_sigma;
    
    return 0;
}
int get_kick_properties(int index, int *kick_distribution, double *kick_distribution_sigma)

{
    if (index > highest_particle_index)
    {
      return -1;
    }
  
    Particle * p = particlesMap[index];
    *kick_distribution = p->kick_distribution;
    *kick_distribution_sigma = p->kick_distribution_sigma;
    
    return 0;
}


/*******************************
 * instantaneous perturbations *
 * ****************************/

int set_instantaneous_perturbation_properties(int index, double delta_mass, double delta_x, double delta_y, double delta_z, double delta_vx, double delta_vy, double delta_vz)
{
    if (index > highest_particle_index)
    {
      return -1;
    }
    
    Particle *p = particlesMap[index];

    p->instantaneous_perturbation_delta_mass = delta_mass;
    p->instantaneous_perturbation_delta_x = delta_x;
    p->instantaneous_perturbation_delta_y = delta_y;
    p->instantaneous_perturbation_delta_z = delta_z;
    p->instantaneous_perturbation_delta_vx = delta_vx;
    p->instantaneous_perturbation_delta_vy = delta_vy;
    p->instantaneous_perturbation_delta_vz = delta_vz;
    
    return 0;
}


/************
 * external *
 * *********/

int set_external_particle_properties(int index, double external_t_ref, double e, double external_r_p, double INCL, double AP, double LAN)
{
    if (index > highest_particle_index)
    {
        return -1;
    }

    Particle *p = particlesMap[index];    
    
    /* determine masses in all binaries */
//    int N_bodies, N_binaries, N_root_finding;
//    determine_binary_parents_and_levels(&particlesMap, &N_bodies, &N_binaries, &N_root_finding);
//    set_binary_masses_from_body_masses(&particlesMap);
    
    p->external_t_ref = external_t_ref;
    p->external_e = e;
    p->external_r_p = external_r_p;
    /* e & h vectors for external particles are understood to be unit vectors */
    compute_orbital_vectors_from_orbital_elements_unit(INCL,AP,LAN,&(p->e_vec_x), &(p->e_vec_y), &(p->e_vec_z), &(p->h_vec_x), &(p->h_vec_y), &(p->h_vec_z) ); 
    
    //printf("set_external_particle_properties inputs %g %g %g %g %g\n",external_t_ref, e, external_r_p, INCL, AP, LAN);
    //printf("set_external_particle_properties OE %g %g %g %g %g %g\n",p->e_vec_x,p->e_vec_y,p->e_vec_z,p->h_vec_x,p->h_vec_y,p->h_vec_z);
    
    return 0;
}




/****************
/* spin vectors *
 ****************/

int set_spin_vector(int index, double spin_vec_x, double spin_vec_y, double spin_vec_z)
{
    if (index > highest_particle_index)
    {
      return -1;
    }
  
    Particle * p = particlesMap[index];
    p->spin_vec_x = spin_vec_x;
    p->spin_vec_y = spin_vec_y;
    p->spin_vec_z = spin_vec_z;
    //printf("set spin %g %g %g\n",spin_vec_x,spin_vec_y,spin_vec_z);
    return 0;
}
int get_spin_vector(int index, double *spin_vec_x, double *spin_vec_y, double *spin_vec_z)
{
    if (index > highest_particle_index)
    {
      return -1;
    }
  
    Particle * p = particlesMap[index];
    *spin_vec_x = p->spin_vec_x;
    *spin_vec_y = p->spin_vec_y;
    *spin_vec_z = p->spin_vec_z;
    
    return 0;
}

int set_spin_vector_dot(int index, double spin_vec_x_dot, double spin_vec_y_dot, double spin_vec_z_dot)
{
    if (index > highest_particle_index)
    {
      return -1;
    }
  
    Particle * p = particlesMap[index];
    p->spin_vec_x_dot = spin_vec_x_dot;
    p->spin_vec_y_dot = spin_vec_y_dot;
    p->spin_vec_z_dot = spin_vec_z_dot;
    
    return 0;
}
int get_spin_vector_dot(int index, double *spin_vec_x_dot, double *spin_vec_y_dot, double *spin_vec_z_dot)
{
    if (index > highest_particle_index)
    {
      return -1;
    }
  
    Particle * p = particlesMap[index];
    *spin_vec_x_dot = p->spin_vec_x_dot;
    *spin_vec_y_dot = p->spin_vec_y_dot;
    *spin_vec_z_dot = p->spin_vec_z_dot;
    
    return 0;
}


int set_orbital_vectors_dot(int index, double e_vec_x_dot, double e_vec_y_dot, double e_vec_z_dot, \
    double h_vec_x_dot, double h_vec_y_dot, double h_vec_z_dot)
{
    if (index > highest_particle_index)
    {
      return -1;
    }
  
    Particle * p = particlesMap[index];
    p->e_vec_x_dot = e_vec_x_dot;
    p->e_vec_y_dot = e_vec_y_dot;
    p->e_vec_z_dot = e_vec_z_dot;
    p->h_vec_x_dot = h_vec_x_dot;
    p->h_vec_y_dot = h_vec_y_dot;
    p->h_vec_z_dot = h_vec_z_dot;
    
    return 0;
}
int get_orbital_vectors_dot(int index, double *e_vec_x_dot, double *e_vec_y_dot, double *e_vec_z_dot, \
    double *h_vec_x_dot, double *h_vec_y_dot, double *h_vec_z_dot)
{
    if (index > highest_particle_index)
    {
      return -1;
    }
  
    Particle * p = particlesMap[index];
    *e_vec_x_dot = p->e_vec_x_dot;
    *e_vec_y_dot = p->e_vec_y_dot;
    *e_vec_z_dot = p->e_vec_z_dot;
    *h_vec_x_dot = p->h_vec_x_dot;
    *h_vec_y_dot = p->h_vec_y_dot;
    *h_vec_z_dot = p->h_vec_z_dot;
    
    return 0;
}


int set_position_vector(int index, double x, double y, double z)
{
    if (index > highest_particle_index)
    {
      return -1;
    }
  
    Particle * p = particlesMap[index];
    p->x = x;
    p->x = y;
    p->x = z;
   
    return 0;
}
int get_position_vector(int index, double *x, double *y, double *z)
{
    //printf("get_position_vector\n");
    if (index > highest_particle_index)
    {
      return -1;
    }
  
    set_positions_and_velocities(&particlesMap);
    
    Particle * p = particlesMap[index];
    *x = p->x;
    *y = p->x;
    *z = p->x;
    
    return 0;
}

int set_velocity_vector(int index, double x, double y, double z)
{
    if (index > highest_particle_index)
    {
      return -1;
    }
  
    Particle * p = particlesMap[index];
    p->vx = x;
    p->vy = y;
    p->vz = z;
   
    return 0;
}
int get_velocity_vector(int index, double *x, double *y, double *z)
{
    if (index > highest_particle_index)
    {
      return -1;
    }

    set_positions_and_velocities(&particlesMap);
    
    Particle * p = particlesMap[index];
    *x = p->vx;
    *y = p->vy;
    *z = p->vz;
    
    return 0;
}

/************
/* PN terms *
 ************/

int set_PN_terms(int index, int include_pairwise_1PN_terms, int include_pairwise_25PN_terms)
{
    if (index > highest_particle_index)
    {
        return -1;
    }

    Particle * p = particlesMap[index];
    p->include_pairwise_1PN_terms = include_pairwise_1PN_terms;
    p->include_pairwise_25PN_terms = include_pairwise_25PN_terms;
        
    return 0;
}
int get_PN_terms(int index, int *include_pairwise_1PN_terms, int *include_pairwise_25PN_terms)
{
    if (index > highest_particle_index)
    {
        return -1;
    }
  
    Particle * p = particlesMap[index];
    
    *include_pairwise_1PN_terms = p->include_pairwise_1PN_terms;
    *include_pairwise_25PN_terms = p->include_pairwise_25PN_terms;
        
    return 0;
}


/*********
/* tides *
 *********/
#ifdef IGNORE
int set_tides_terms(int index, int include_tidal_friction_terms, int tides_method, int include_tidal_bulges_precession_terms, int include_rotation_precession_terms, double minimum_eccentricity_for_tidal_precession, 
    double tides_apsidal_motion_constant, double tides_gyration_radius, double tides_viscous_time_scale, int tides_viscous_time_scale_prescription, double convective_envelope_mass, double convective_envelope_radius, double luminosity)
{
    if (index > highest_particle_index)
    {
        return -1;
    }

    Particle *p = particlesMap[index];
    
    p->include_tidal_friction_terms = include_tidal_friction_terms;
    p->tides_method = tides_method;
    p->include_tidal_bulges_precession_terms = include_tidal_bulges_precession_terms;
    p->include_rotation_precession_terms = include_rotation_precession_terms;
    p->minimum_eccentricity_for_tidal_precession = minimum_eccentricity_for_tidal_precession;
    p->tides_apsidal_motion_constant = tides_apsidal_motion_constant;
    p->tides_gyration_radius = tides_gyration_radius;
    p->tides_viscous_time_scale = tides_viscous_time_scale;
    p->tides_viscous_time_scale_prescription = tides_viscous_time_scale_prescription;
    p->convective_envelope_mass = convective_envelope_mass;
    p->convective_envelope_radius = convective_envelope_radius;
    p->luminosity = luminosity;
    //printf("set tides1 %d %d %d %g\n",include_tidal_friction_terms,include_tidal_bulges_precession_terms,include_rotation_precession_terms,minimum_eccentricity_for_tidal_precession);
    //printf("set tides2 %g %g %g %d %g %g %g\n",tides_apsidal_motion_constant,tides_gyration_radius,tides_viscous_time_scale,tides_viscous_time_scale_prescription,convective_envelope_mass,convective_envelope_radius,luminosity);
    return 0;
}
int get_tides_terms(int index, int *include_tidal_friction_terms, int *tides_method, int *include_tidal_bulges_precession_terms, int *include_rotation_precession_terms, double *minimum_eccentricity_for_tidal_precession, 
    double *tides_apsidal_motion_constant, double *tides_gyration_radius, double *tides_viscous_time_scale, int *tides_viscous_time_scale_prescription, double *convective_envelope_mass, double *convective_envelope_radius, double *luminosity)
{
    if (index > highest_particle_index)
    {
        return -1;
    }
  
    Particle *p = particlesMap[index];
    
    *include_tidal_friction_terms = p->include_tidal_friction_terms;
    *tides_method = p->tides_method;
    *include_tidal_bulges_precession_terms = p->include_tidal_bulges_precession_terms;
    *include_rotation_precession_terms = p->include_rotation_precession_terms;
    *minimum_eccentricity_for_tidal_precession = p->minimum_eccentricity_for_tidal_precession;
    *tides_apsidal_motion_constant = p->tides_apsidal_motion_constant;
    *tides_gyration_radius = p->tides_gyration_radius;
    *tides_viscous_time_scale = p->tides_viscous_time_scale;
    *tides_viscous_time_scale_prescription = p->tides_viscous_time_scale_prescription;
    *convective_envelope_mass = p->convective_envelope_mass;
    *convective_envelope_radius = p->convective_envelope_radius;
    *luminosity = p->luminosity;
    return 0;
}
#endif

/****************
 * VRR          *
 ****************/

int set_VRR_properties(int index, int VRR_model, int VRR_include_mass_precession, double VRR_mass_precession_rate, 
    double VRR_Omega_vec_x, double VRR_Omega_vec_y, double VRR_Omega_vec_z, 
    double VRR_eta_20_init, double VRR_eta_a_22_init, double VRR_eta_b_22_init, double VRR_eta_a_21_init, double VRR_eta_b_21_init,
    double VRR_eta_20_final, double VRR_eta_a_22_final, double VRR_eta_b_22_final, double VRR_eta_a_21_final, double VRR_eta_b_21_final,
	double VRR_initial_time, double VRR_final_time)
{
	
	if (index > highest_particle_index)
    {
        return -1;
    }

    Particle *p = particlesMap[index];
    
    p->VRR_model = VRR_model;
    p->VRR_include_mass_precession = VRR_include_mass_precession;
	p->VRR_mass_precession_rate = VRR_mass_precession_rate;
	p->VRR_Omega_vec_x = VRR_Omega_vec_x;
	p->VRR_Omega_vec_y = VRR_Omega_vec_y;
	p->VRR_Omega_vec_z = VRR_Omega_vec_z;
	p->VRR_eta_20_init = VRR_eta_20_init;
	p->VRR_eta_a_22_init = VRR_eta_a_22_init;
	p->VRR_eta_b_22_init = VRR_eta_b_22_init;
	p->VRR_eta_a_21_init = VRR_eta_a_21_init;
	p->VRR_eta_b_21_init = VRR_eta_b_21_init;
	p->VRR_eta_20_final = VRR_eta_20_final;
	p->VRR_eta_a_22_final = VRR_eta_a_22_final;
	p->VRR_eta_b_22_final = VRR_eta_b_22_final;
	p->VRR_eta_a_21_final = VRR_eta_a_21_final;
	p->VRR_eta_b_21_final = VRR_eta_b_21_final;
	p->VRR_initial_time = VRR_initial_time;
	p->VRR_final_time = VRR_final_time;

	return 0;
}

/****************
/* root finding *
 ****************/
int set_root_finding_terms(int index, int check_for_secular_breakdown, int check_for_dynamical_instability, int dynamical_instability_criterion, int dynamical_instability_central_particle, int dynamical_instability_K_parameter,
    int check_for_physical_collision_or_orbit_crossing, int check_for_minimum_periapse_distance, double check_for_minimum_periapse_distance_value, int check_for_RLOF_at_pericentre, int check_for_RLOF_at_pericentre_use_sepinsky_fit, int check_for_GW_condition)
{
    if (index > highest_particle_index)
    {
        return -1;
    }

    Particle *p = particlesMap[index];
    p->check_for_secular_breakdown = check_for_secular_breakdown;
    p->check_for_dynamical_instability = check_for_dynamical_instability;
    p->dynamical_instability_criterion = dynamical_instability_criterion;
    p->dynamical_instability_central_particle = dynamical_instability_central_particle;
    p->dynamical_instability_K_parameter = dynamical_instability_K_parameter;
    p->check_for_physical_collision_or_orbit_crossing = check_for_physical_collision_or_orbit_crossing;
    p->check_for_minimum_periapse_distance = check_for_minimum_periapse_distance;
    p->check_for_minimum_periapse_distance_value = check_for_minimum_periapse_distance_value;
    p->check_for_RLOF_at_pericentre = check_for_RLOF_at_pericentre;
    p->check_for_RLOF_at_pericentre_use_sepinsky_fit = check_for_RLOF_at_pericentre_use_sepinsky_fit;
    p->check_for_GW_condition = check_for_GW_condition;
    return 0;
}
int get_root_finding_terms(int index, int *check_for_secular_breakdown, int *check_for_dynamical_instability, int *dynamical_instability_criterion, int *dynamical_instability_central_particle, int *dynamical_instability_K_parameter,
    int *check_for_physical_collision_or_orbit_crossing, int *check_for_minimum_periapse_distance, double *check_for_minimum_periapse_distance_value, int *check_for_RLOF_at_pericentre, int *check_for_RLOF_at_pericentre_use_sepinsky_fit, int *check_for_GW_condition)
{
    if (index > highest_particle_index)
    {
        return -1;
    }
  
    Particle *p = particlesMap[index];
    *check_for_secular_breakdown = p->check_for_secular_breakdown;
    *check_for_dynamical_instability = p->check_for_dynamical_instability;
    *dynamical_instability_criterion = p->dynamical_instability_criterion;
    *dynamical_instability_central_particle = p->dynamical_instability_central_particle;
    *dynamical_instability_K_parameter = p->dynamical_instability_K_parameter;
    *check_for_physical_collision_or_orbit_crossing = p->check_for_physical_collision_or_orbit_crossing;
    *check_for_minimum_periapse_distance = p->check_for_minimum_periapse_distance;
    *check_for_minimum_periapse_distance_value = p->check_for_minimum_periapse_distance_value;
    *check_for_RLOF_at_pericentre = p->check_for_RLOF_at_pericentre;
    *check_for_RLOF_at_pericentre_use_sepinsky_fit = p->check_for_RLOF_at_pericentre_use_sepinsky_fit;
    *check_for_GW_condition = p->check_for_GW_condition;
    return 0;
}

/* retrieve root finding state */
int set_root_finding_state(int index, int secular_breakdown_has_occurred, int dynamical_instability_has_occurred, int physical_collision_or_orbit_crossing_has_occurred, int minimum_periapse_distance_has_occurred, int RLOF_at_pericentre_has_occurred, int GW_condition_has_occurred)
{
    if (index > highest_particle_index)
    {
        return -1;
    }
  
    Particle * p = particlesMap[index];

    p->secular_breakdown_has_occurred = secular_breakdown_has_occurred;
    p->dynamical_instability_has_occurred = dynamical_instability_has_occurred;
    p->physical_collision_or_orbit_crossing_has_occurred = physical_collision_or_orbit_crossing_has_occurred;
    p->minimum_periapse_distance_has_occurred = minimum_periapse_distance_has_occurred;
    p->RLOF_at_pericentre_has_occurred = RLOF_at_pericentre_has_occurred;
    p->GW_condition_has_occurred = GW_condition_has_occurred;
    
    return 0;
}
int get_root_finding_state(int index, int *secular_breakdown_has_occurred, int *dynamical_instability_has_occurred, int *physical_collision_or_orbit_crossing_has_occurred, int* minimum_periapse_distance_has_occurred, int *RLOF_at_pericentre_has_occurred, int *GW_condition_has_occurred)
{
    if (index > highest_particle_index)
    {
        return -1;
    }
  
    Particle * p = particlesMap[index];

    *secular_breakdown_has_occurred = p->secular_breakdown_has_occurred;
    *dynamical_instability_has_occurred = p->dynamical_instability_has_occurred;
    *physical_collision_or_orbit_crossing_has_occurred = p->physical_collision_or_orbit_crossing_has_occurred;
    *minimum_periapse_distance_has_occurred = p->minimum_periapse_distance_has_occurred;
    *RLOF_at_pericentre_has_occurred = p->RLOF_at_pericentre_has_occurred;
    *GW_condition_has_occurred = p->GW_condition_has_occurred;
    
    return 0;
}


/********************
/* evolve interface *
 ********************/

int evolve_interface(double start_time, double end_time, double *output_time, double *hamiltonian, int *state, int *CVODE_flag, int *CVODE_error_code)
{
    //printf("interface %g %g\n",start_time,time_step);
    srand(random_seed);

    //printf("setting seed %d\n",random_seed);
    int result = evolve(&particlesMap, start_time, end_time, output_time, hamiltonian, state, CVODE_flag, CVODE_error_code);
    
    return result;
}


/* set levels and masses */
int determine_binary_parents_levels_and_masses_interface()
{
    //printf("determine_binary_parents_levels_and_masses_interface\n");
    int N_bodies, N_binaries, N_root_finding;
    determine_binary_parents_and_levels(&particlesMap, &N_bodies, &N_binaries, &N_root_finding);
    set_binary_masses_from_body_masses(&particlesMap);
    
    return 0;
}

int apply_external_perturbation_assuming_integrated_orbits_interface()
{
    //printf("apply_external_perturbation_assuming_integrated_orbits_interface\n");
    apply_external_perturbation_assuming_integrated_orbits(&particlesMap);

    return 0;
}

int apply_user_specified_instantaneous_perturbation_interface()
{
    //printf("apply_user_specified_instantaneous_perturbation\n");
    apply_user_specified_instantaneous_perturbation(&particlesMap);
    
    return 0;
}


int clear_internal_particles()
{
    //printf("clear_internal_particles\n");
    particlesMap.clear();
	highest_particle_index = 0;
    return 0;
}

int initialize_code()
{
    bool unbound_orbits;
    printf("interface.cpp -- set_up_flybys\n");
    handle_next_flyby(&particlesMap,true,&unbound_orbits);
    
    initialize_stars(&particlesMap);
    
    return 0;
}


int set_positions_and_velocities_interface()
{
    set_positions_and_velocities(&particlesMap);
    
    return 0;
}



int get_de_dt(int index, double *de_dt)
{

    if (index > highest_particle_index)
    {
        return -1;
    }
  
    Particle * p = particlesMap[index];
    if (p->is_binary == 0)
    {
        *de_dt = 0.0;
        return 0;
    }

    *de_dt = dot3(p->e_vec_unit,p->de_vec_dt);

    return 0;
}


void get_position_and_velocity_vectors_from_particle(Particle *p, double r[3], double v[3])
{
    r[0] = p->x;
    r[1] = p->y;
    r[2] = p->z;
    v[0] = p->vx;
    v[1] = p->vy;
    v[2] = p->vz;
}
void set_position_and_velocity_vectors_in_particle(Particle *p,  double r[3], double v[3])
{
    p->x = r[0];
    p->y = r[1];
    p->z = r[2];
    p->vx = v[0];
    p->vy = v[1];
    p->vz = v[2];
}
void get_e_and_h_vectors_from_particle(Particle *p, double e_vec[3], double h_vec[3])
{
    e_vec[0] = p->e_vec_x;
    e_vec[1] = p->e_vec_y;
    e_vec[2] = p->e_vec_z;
    h_vec[0] = p->h_vec_x;    
    h_vec[1] = p->h_vec_y;    
    h_vec[2] = p->h_vec_z;    
}
void set_e_and_h_vectors_in_particle(Particle *p, double e_vec[3], double h_vec[3])
{
    p->e_vec_x = e_vec[0];
    p->e_vec_y = e_vec[1];
    p->e_vec_z = e_vec[2];
    p->h_vec_x = h_vec[0];
    p->h_vec_y = h_vec[1];
    p->h_vec_z = h_vec[2];
}


/************************
/* interface parameters *
 ************************/
 
 
int set_constants(double CONST_G_, double CONST_C_, double CONST_MSUN_, double CONST_R_SUN_, double CONST_L_SUN_, double CONST_KM_PER_S_, double CONST_PER_PC3_)
{
    CONST_G = CONST_G_;
    CONST_G_P2 = CONST_G*CONST_G;
    CONST_G_P3 = CONST_G_P2*CONST_G;
    
    CONST_C_LIGHT = CONST_C_;
    CONST_C_LIGHT_P2 = CONST_C_LIGHT*CONST_C_LIGHT;
    CONST_C_LIGHT_P4 = CONST_C_LIGHT_P2*CONST_C_LIGHT_P2;
    CONST_C_LIGHT_P5 = CONST_C_LIGHT_P4*CONST_C_LIGHT;

    CONST_MSUN = CONST_MSUN_;
    CONST_R_SUN = CONST_R_SUN_;
    CONST_L_SUN = CONST_L_SUN_;
    
    CONST_KM_PER_S = CONST_KM_PER_S_;
    CONST_PER_PC3 = CONST_PER_PC3_;
    
    //printf("CONSTS %g %g %g\n",CONST_G,CONST_C_LIGHT,CONST_MSUN);
    
    return 0;
}

int set_parameters(double relative_tolerance_, double absolute_tolerance_eccentricity_vectors_, 
     bool include_quadrupole_order_terms_, bool include_octupole_order_binary_pair_terms_, bool include_octupole_order_binary_triplet_terms_,
     bool include_hexadecupole_order_binary_pair_terms_, bool include_dotriacontupole_order_binary_pair_terms_, 
     bool include_flybys_, int flybys_reference_binary_, bool flybys_correct_for_gravitational_focussing_, int flybys_velocity_distribution_, int flybys_mass_distribution_,
     double flybys_mass_distribution_lower_value_, double flybys_mass_distribution_upper_value_, double flybys_encounter_sphere_radius_, 
     double flybys_stellar_density_, double flybys_stellar_relative_velocity_dispersion_)
{
     relative_tolerance = relative_tolerance_;
     absolute_tolerance_eccentricity_vectors = absolute_tolerance_eccentricity_vectors_;
     include_quadrupole_order_terms = include_quadrupole_order_terms_;
     include_octupole_order_binary_pair_terms = include_octupole_order_binary_pair_terms_;
     include_octupole_order_binary_triplet_terms = include_octupole_order_binary_triplet_terms_;
     include_hexadecupole_order_binary_pair_terms = include_hexadecupole_order_binary_pair_terms_;
     include_dotriacontupole_order_binary_pair_terms = include_dotriacontupole_order_binary_pair_terms_;
     
     include_flybys = include_flybys_;
     flybys_correct_for_gravitational_focussing = flybys_correct_for_gravitational_focussing_;
     flybys_velocity_distribution = flybys_velocity_distribution_;
     flybys_mass_distribution = flybys_mass_distribution_;
     flybys_reference_binary = flybys_reference_binary_;
     flybys_mass_distribution_lower_value = flybys_mass_distribution_lower_value_;
     flybys_mass_distribution_upper_value = flybys_mass_distribution_upper_value_;
     flybys_encounter_sphere_radius = flybys_encounter_sphere_radius_;
     flybys_stellar_density = flybys_stellar_density_;
     flybys_stellar_relative_velocity_dispersion = flybys_stellar_relative_velocity_dispersion_;
     
     //printf("PARAMS %g %g %d %d %d %d %d\n",relative_tolerance,absolute_tolerance_eccentricity_vectors,include_quadrupole_order_terms,include_octupole_order_binary_pair_terms,include_octupole_order_binary_triplet_terms,include_hexadecupole_order_binary_pair_terms,include_dotriacontupole_order_binary_pair_terms);

     return 0;
}
 
int get_relative_tolerance(double *value)
{
    *value = relative_tolerance;
    return 0;
}
int set_relative_tolerance(double value)
{
    relative_tolerance = value;
    return 0;
}
int get_absolute_tolerance_eccentricity_vectors(double *value)
{
    *value = absolute_tolerance_eccentricity_vectors;
    return 0;
}
int set_absolute_tolerance_eccentricity_vectors(double value)
{
    absolute_tolerance_eccentricity_vectors = value;
    return 0;
}

int get_include_quadrupole_order_terms(int *value){
    *value = include_quadrupole_order_terms ? 1 : 0;
    return 0;
}
int set_include_quadrupole_order_terms(int value){
    include_quadrupole_order_terms = value == 1;
    return 0;
}

int get_include_octupole_order_binary_pair_terms(int *value){
    *value = include_octupole_order_binary_pair_terms ? 1 : 0;
    return 0;
}
int set_include_octupole_order_binary_pair_terms(int value){
    include_octupole_order_binary_pair_terms = value == 1;
    return 0;
}

int get_include_octupole_order_binary_triplet_terms(int *value){
    *value = include_octupole_order_binary_triplet_terms ? 1 : 0;
    return 0;
}
int set_include_octupole_order_binary_triplet_terms(int value){
    include_octupole_order_binary_triplet_terms = value == 1;
    return 0;
}

int get_include_hexadecupole_order_binary_pair_terms(int *value){
    *value = include_hexadecupole_order_binary_pair_terms ? 1 : 0;
    return 0;
}
int set_include_hexadecupole_order_binary_pair_terms(int value){
    include_hexadecupole_order_binary_pair_terms = value == 1;
    return 0;
}

int get_include_dotriacontupole_order_binary_pair_terms(int *value){
    *value = include_dotriacontupole_order_binary_pair_terms ? 1 : 0;
    return 0;
}
int set_include_dotriacontupole_order_binary_pair_terms(int value){
    include_dotriacontupole_order_binary_pair_terms = value == 1;
    return 0;
}

int get_random_seed(int *value)
{
    *value = random_seed;
    return 0;
}
int set_random_seed(int value)
{
    random_seed = value;
    return 0;
}

}


