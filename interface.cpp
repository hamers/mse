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
 
int add_particle(int *index, bool is_binary, bool is_external)
{
    *index = particlesMap.size();

    Particle *p = new Particle(*index, is_binary);
    particlesMap[*index] = p;

    p->is_external = is_external;
    //highest_particle_index = *index;
       
    return 0;
}

int delete_particle(int index)
{
    if (index > particlesMap.size())
    {
        return -1;
    }
  
    particlesMap.erase(index);

    return 0;
}

int set_children(int index, int child1, int child2)
{
    if (index > particlesMap.size())
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
    if (index > particlesMap.size())
    {
      return -1;
    }
  
    Particle * p = particlesMap[index];
    *child1 = p->child1;
    *child2 = p->child2;
    
    return 0;
}

int get_number_of_particles(int *N_particles)
{
    *N_particles = particlesMap.size();
    
    return 0;
}

int get_is_binary(int index, bool *is_binary)
{
    if (index > particlesMap.size())
    {
      return -1;
    }
  
    Particle * p = particlesMap[index];
    *is_binary = p->is_binary;
    
    return 0;
}

int get_is_bound(int index, bool *is_bound)
{
    if (index > particlesMap.size())
    {
      return -1;
    }
  
    Particle *p = particlesMap[index];
    //printf("getisb %d %d\n",index,p->is_bound);
    *is_bound = p->is_bound;
    
    return 0;
}    

int set_mass(int index, double mass)
{
    if (index > particlesMap.size())
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
    if (index > particlesMap.size())
    {
      return -1;
    }
    Particle *p = particlesMap[index];
    *mass = p->mass;
    return 0;
}
int set_mass_transfer_terms(int index, bool include_mass_transfer_terms)
{
    if (index > particlesMap.size())
    {
      return -1;
    }

    Particle *p = particlesMap[index];
    p->include_mass_transfer_terms = include_mass_transfer_terms;
    
    return 0;
}
int get_mass_dot(int index, double *mass_dot)
{
    if (index > particlesMap.size())
    {
      return -1;
    }

    Particle *p = particlesMap[index];
    *mass_dot = p->mass_dot_wind + p->mass_dot_RLOF;

    return 0;
}

int set_radius(int index, double radius, double radius_dot)
{
    if (index > particlesMap.size())
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
    if (index > particlesMap.size())
    {
      return -1;
    }
  
    Particle * p = particlesMap[index];
    *radius = p->radius;
    *radius_dot = p->radius_dot;
    
    return 0;
}
int set_integration_method(int index, int integration_method, bool KS_use_perturbing_potential)
{
    if (index > particlesMap.size())
    {
      return -1;
    }

    Particle * p = particlesMap[index];
    p->integration_method = integration_method;
    p->KS_use_perturbing_potential = KS_use_perturbing_potential;
    
    return 0;
}


/****************************
/* orbital vectors/elements *
 ****************************/

int set_orbital_vectors(int index, double e_vec_x, double e_vec_y, double e_vec_z, \
    double h_vec_x, double h_vec_y, double h_vec_z)
{
    if (index > particlesMap.size())
    {
      return -1;
    }
  
    Particle *p = particlesMap[index];

    p->e_vec[0] = e_vec_x;
    p->e_vec[1] = e_vec_y;
    p->e_vec[2] = e_vec_z;
    p->h_vec[0] = h_vec_x;
    p->h_vec[1] = h_vec_y;
    p->h_vec[2] = h_vec_z;
    
    return 0;
}
int get_orbital_vectors(int index, double *e_vec_x, double *e_vec_y, double *e_vec_z, \
    double *h_vec_x, double *h_vec_y, double *h_vec_z)
{
    if (index > particlesMap.size())
    {
      return -1;
    }
  
    Particle *p = particlesMap[index];

    *e_vec_x = p->e_vec[0];
    *e_vec_y = p->e_vec[1];
    *e_vec_z = p->e_vec[2];
    *h_vec_x = p->h_vec[0];
    *h_vec_y = p->h_vec[1];
    *h_vec_z = p->h_vec[2];

    return 0;
}

int set_orbital_elements(int index, double semimajor_axis, double eccentricity, double true_anomaly, \
    double inclination, double argument_of_pericenter, double longitude_of_ascending_node, bool sample_orbital_phase_randomly)
{
    if (index > particlesMap.size())
    {
        return -1;
    }

    Particle * p = particlesMap[index];    
    
    if (p->is_binary == false)
    {
        return 0;
    }

    /* determine masses in all binaries */
    int N_bodies, N_binaries, N_root_finding,N_ODE_equations;

    determine_binary_parents_and_levels(&particlesMap, &N_bodies, &N_binaries, &N_root_finding,&N_ODE_equations);

    set_binary_masses_from_body_masses(&particlesMap);

    compute_orbital_vectors_from_orbital_elements(p->child1_mass, p->child2_mass, semimajor_axis, eccentricity, \
        inclination, argument_of_pericenter, longitude_of_ascending_node, \
        &(p->e_vec[0]), &(p->e_vec[1]), &(p->e_vec[2]), &(p->h_vec[0]), &(p->h_vec[1]), &(p->h_vec[2]) );
    
    p->true_anomaly = true_anomaly;
    p->sample_orbital_phase_randomly = sample_orbital_phase_randomly;
    //printf("soe a %g e %g TA %g I %g AP %g LAN %g SOPR %d\n",semimajor_axis,eccentricity,true_anomaly,inclination,argument_of_pericenter,longitude_of_ascending_node,sample_orbital_phase_randomly);
    //printf("set_orbital_elements %g %g %g %g %g %g\n",p->e_vec[0],p->e_vec[1],p->e_vec[2],p->h_vec[0],p->h_vec[1],p->h_vec[2]);
    
    return 0;
}
int get_orbital_elements(int index, double *semimajor_axis, double *eccentricity, double *true_anomaly, \
    double *inclination, double *argument_of_pericenter, double *longitude_of_ascending_node)
{
    if (index > particlesMap.size())
    {
        return -1;
    }
  
    Particle * p = particlesMap[index];    
    
    if (p->is_binary == false)
    {
        return 0;
    }
    
    *true_anomaly = p->true_anomaly;
    
    double h_tot_vec[3];
    compute_h_tot_vector(&particlesMap,h_tot_vec);

    /* determine masses in all binaries */
    int N_bodies, N_binaries, N_root_finding, N_ODE_equations;
    determine_binary_parents_and_levels(&particlesMap, &N_bodies, &N_binaries, &N_root_finding,&N_ODE_equations);
    set_binary_masses_from_body_masses(&particlesMap);
    
    compute_orbital_elements_from_orbital_vectors(p->child1_mass, p->child2_mass, h_tot_vec, \
        p->e_vec[0],p->e_vec[1],p->e_vec[2],p->h_vec[0],p->h_vec[1],p->h_vec[2],
        semimajor_axis, eccentricity, inclination, argument_of_pericenter, longitude_of_ascending_node);
    return 0;
}

int get_inclination_relative_to_parent_interface(int index, double *inclination_relative_to_parent)
{
    get_inclination_relative_to_parent(&particlesMap,index,inclination_relative_to_parent);
    
    return 0;
}

int get_level(int index, int *value)
{
    if (index > particlesMap.size())
    {
      return -1;
    }
  
    Particle * p = particlesMap[index];
    *value = p->level;
    
    return 0;
}

        
int set_stellar_evolution_properties(int index, int stellar_type, bool evolve_as_star, double sse_initial_mass, double metallicity, double sse_time_step, double epoch, double age, 
    double convective_envelope_mass, double convective_envelope_radius, double core_mass, double core_radius, double luminosity, double apsidal_motion_constant, double gyration_radius, double tides_viscous_time_scale, int tides_viscous_time_scale_prescription)
{
    //printf("set_stellar_evolution_properties %g\n",metallicity);
    if (index > particlesMap.size())
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
    p->tides_viscous_time_scale_prescription = tides_viscous_time_scale_prescription;
    
    return 0;
}
int get_stellar_evolution_properties(int index, int *stellar_type, bool *evolve_as_star, double *sse_initial_mass, double *metallicity, double *sse_time_step, double *epoch, double *age, 
    double *convective_envelope_mass, double *convective_envelope_radius, double *core_mass, double *core_radius, double *luminosity, double *apsidal_motion_constant, double *gyration_radius, double *tides_viscous_time_scale, double *roche_lobe_radius_pericenter)

{
    if (index > particlesMap.size())
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
    
    if (p->is_binary == false) /* TO DO: streamline with same code in root_finding.cpp */
    {
        if (p->is_bound == true) /* only update parameters below for bodies in bound orbits */
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
    }
    
    return 0;
}

int set_kick_properties(int index, int kick_distribution, double kick_distribution_sigma)
{
    //printf("set_stellar_evolution_properties %g\n",metallicity);
    if (index > particlesMap.size())
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
    if (index > particlesMap.size())
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

int set_instantaneous_perturbation_properties(int index, double delta_mass, double delta_X, double delta_Y, double delta_Z, double delta_VX, double delta_VY, double delta_VZ)
{
    if (index > particlesMap.size())
    {
      return -1;
    }
    
    Particle *p = particlesMap[index];

    p->instantaneous_perturbation_delta_mass = delta_mass;
    p->instantaneous_perturbation_delta_X = delta_X;
    p->instantaneous_perturbation_delta_Y = delta_Y;
    p->instantaneous_perturbation_delta_Z = delta_Z;
    p->instantaneous_perturbation_delta_VX = delta_VX;
    p->instantaneous_perturbation_delta_VY = delta_VY;
    p->instantaneous_perturbation_delta_VZ = delta_VZ;
    
    return 0;
}


/************
 * external *
 * *********/

int set_external_particle_properties(int index, double external_t_ref, double e, double external_r_p, double INCL, double AP, double LAN)
{
    if (index > particlesMap.size())
    {
        return -1;
    }

    Particle *p = particlesMap[index];    
    
    p->external_t_ref = external_t_ref;
    p->external_e = e;
    p->external_r_p = external_r_p;
    
    /* e & h vectors for external particles are understood to be unit vectors */
    compute_orbital_vectors_from_orbital_elements_unit(INCL,AP,LAN,&(p->e_vec[0]), &(p->e_vec[1]), &(p->e_vec[2]), &(p->h_vec[0]), &(p->h_vec[1]), &(p->h_vec[2]) ); 
    
    p->evolve_as_star = false; // if enabled, can lead to segfaults in stellar_evolution.cpp
    
    //printf("set_external_particle_properties inputs %g %g %g %g %g\n",external_t_ref, e, external_r_p, INCL, AP, LAN);
    //printf("set_external_particle_properties OE %g %g %g %g %g %g\n",p->e_vec[0],p->e_vec[1],p->e_vec[2],p->h_vec[0],p->h_vec[1],p->h_vec[2]);
    
    return 0;
}




/****************
/* spin vectors *
 ****************/

int set_spin_vector(int index, double spin_vec_x, double spin_vec_y, double spin_vec_z)
{
    if (index > particlesMap.size())
    {
      return -1;
    }
  
    Particle * p = particlesMap[index];
    p->spin_vec[0] = spin_vec_x;
    p->spin_vec[1] = spin_vec_y;
    p->spin_vec[2] = spin_vec_z;

    return 0;
}
int get_spin_vector(int index, double *spin_vec_x, double *spin_vec_y, double *spin_vec_z)
{
    if (index > particlesMap.size())
    {
      return -1;
    }
  
    Particle * p = particlesMap[index];
    *spin_vec_x = p->spin_vec[0];
    *spin_vec_y = p->spin_vec[1];
    *spin_vec_z = p->spin_vec[2];
    
    return 0;
}

int set_spin_vector_dot(int index, double spin_vec_x_dot, double spin_vec_y_dot, double spin_vec_z_dot)
{
    if (index > particlesMap.size())
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
    if (index > particlesMap.size())
    {
      return -1;
    }
  
    Particle * p = particlesMap[index];
    *spin_vec_x_dot = p->spin_vec_x_dot;
    *spin_vec_y_dot = p->spin_vec_y_dot;
    *spin_vec_z_dot = p->spin_vec_z_dot;
    
    return 0;
}


int get_relative_position_and_velocity(int index, double *x, double *y, double *z, double *vx, double *vy, double *vz)
{
    if (index > particlesMap.size())
    {
      return -1;
    }
  
    Particle * p = particlesMap[index];
    *x = p->r_vec[0];
    *y = p->r_vec[1];
    *z = p->r_vec[2];
    *vx = p->v_vec[0];
    *vy = p->v_vec[1];
    *vz = p->v_vec[2];
   
    return 0;
}
int get_absolute_position_and_velocity(int index, double *X, double *Y, double *Z, double *VX, double *VY, double *VZ)
{
    if (index > particlesMap.size())
    {
      return -1;
    }
  
    Particle * p = particlesMap[index];
    
    *X = p->R_vec[0];
    *Y = p->R_vec[1];
    *Z = p->R_vec[2];
    *VX = p->V_vec[0];
    *VY = p->V_vec[1];
    *VZ = p->V_vec[2];
      
    return 0;
}

/************
/* PN terms *
 ************/

int set_PN_terms(int index, int include_pairwise_1PN_terms, int include_pairwise_25PN_terms)
{
    if (index > particlesMap.size())
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
    if (index > particlesMap.size())
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

int set_tides_terms(int index, bool include_tidal_friction_terms, int tides_method, bool include_tidal_bulges_precession_terms, bool include_rotation_precession_terms, double minimum_eccentricity_for_tidal_precession)
{
    if (index > particlesMap.size())
    {
        return -1;
    }

    Particle *p = particlesMap[index];
    
    p->include_tidal_friction_terms = include_tidal_friction_terms;
    p->tides_method = tides_method;
    p->include_tidal_bulges_precession_terms = include_tidal_bulges_precession_terms;
    p->include_rotation_precession_terms = include_rotation_precession_terms;
    p->minimum_eccentricity_for_tidal_precession = minimum_eccentricity_for_tidal_precession;
    //printf("set tides1 %d %d %d %g\n",include_tidal_friction_terms,include_tidal_bulges_precession_terms,include_rotation_precession_terms,minimum_eccentricity_for_tidal_precession);
    //printf("set tides2 %g %g %g %d %g %g %g\n",tides_apsidal_motion_constant,tides_gyration_radius,tides_viscous_time_scale,tides_viscous_time_scale_prescription,convective_envelope_mass,convective_envelope_radius,luminosity);
    return 0;
}
int get_tides_terms(int index, bool *include_tidal_friction_terms, bool *tides_method, bool *include_tidal_bulges_precession_terms, bool *include_rotation_precession_terms, double *minimum_eccentricity_for_tidal_precession)
{
    if (index > particlesMap.size())
    {
        return -1;
    }
  
    Particle *p = particlesMap[index];
    
    *include_tidal_friction_terms = p->include_tidal_friction_terms;
    *tides_method = p->tides_method;
    *include_tidal_bulges_precession_terms = p->include_tidal_bulges_precession_terms;
    *include_rotation_precession_terms = p->include_rotation_precession_terms;
    *minimum_eccentricity_for_tidal_precession = p->minimum_eccentricity_for_tidal_precession;
    return 0;
}


/****************
 * VRR          *
 ****************/

int set_VRR_properties(int index, int VRR_model, int VRR_include_mass_precession, double VRR_mass_precession_rate, 
    double VRR_Omega_vec_x, double VRR_Omega_vec_y, double VRR_Omega_vec_z, 
    double VRR_eta_20_init, double VRR_eta_a_22_init, double VRR_eta_b_22_init, double VRR_eta_a_21_init, double VRR_eta_b_21_init,
    double VRR_eta_20_final, double VRR_eta_a_22_final, double VRR_eta_b_22_final, double VRR_eta_a_21_final, double VRR_eta_b_21_final,
	double VRR_initial_time, double VRR_final_time)
{
	
	if (index > particlesMap.size())
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
int set_root_finding_terms(int index, bool check_for_secular_breakdown, bool check_for_dynamical_instability, int dynamical_instability_criterion, int dynamical_instability_central_particle, double dynamical_instability_K_parameter,
    bool check_for_physical_collision_or_orbit_crossing, bool check_for_minimum_periapse_distance, double check_for_minimum_periapse_distance_value, bool check_for_RLOF_at_pericentre, bool check_for_RLOF_at_pericentre_use_sepinsky_fit, bool check_for_GW_condition)
{
    if (index > particlesMap.size())
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
int get_root_finding_terms(int index, bool *check_for_secular_breakdown, bool *check_for_dynamical_instability, int *dynamical_instability_criterion, int *dynamical_instability_central_particle, double *dynamical_instability_K_parameter,
    bool *check_for_physical_collision_or_orbit_crossing, bool *check_for_minimum_periapse_distance, double *check_for_minimum_periapse_distance_value, bool *check_for_RLOF_at_pericentre, bool *check_for_RLOF_at_pericentre_use_sepinsky_fit, bool *check_for_GW_condition)
{
    if (index > particlesMap.size())
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
int set_root_finding_state(int index, bool secular_breakdown_has_occurred, bool dynamical_instability_has_occurred, bool physical_collision_or_orbit_crossing_has_occurred, bool minimum_periapse_distance_has_occurred, bool RLOF_at_pericentre_has_occurred, bool GW_condition_has_occurred)
{
    if (index > particlesMap.size())
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
int get_root_finding_state(int index, bool *secular_breakdown_has_occurred, bool *dynamical_instability_has_occurred, bool *physical_collision_or_orbit_crossing_has_occurred, bool *minimum_periapse_distance_has_occurred, bool *RLOF_at_pericentre_has_occurred, bool *GW_condition_has_occurred)
{
    if (index > particlesMap.size())
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

int initialize_code_interface()
{
    initialize_code(&particlesMap);

    return 0;
}

int evolve_interface(double start_time, double end_time, double *output_time, double *hamiltonian, int *state, int *CVODE_flag, int *CVODE_error_code, int *integration_flag)
{
    //printf("interface %g %g\n",start_time,time_step);
    srand(random_seed);

    //printf("setting seed %d\n",random_seed);
    int result = evolve(&particlesMap, start_time, end_time, output_time, hamiltonian, state, CVODE_flag, CVODE_error_code, integration_flag);
    
    return result;
}


/* set levels and masses */
int determine_binary_parents_levels_and_masses_interface()
{
    //printf("determine_binary_parents_levels_and_masses_interface\n");
    int N_bodies, N_binaries, N_root_finding, N_ODE_equations;
    determine_binary_parents_and_levels(&particlesMap, &N_bodies, &N_binaries, &N_root_finding,&N_ODE_equations);
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


int reset_interface()
{
    clear_particles(&particlesMap);
    
    return 0;
}


int set_positions_and_velocities_interface()
{
    set_positions_and_velocities(&particlesMap);
    
    return 0;
}



int get_de_dt(int index, double *de_dt)
{

    if (index > particlesMap.size())
    {
        return -1;
    }
  
    Particle * p = particlesMap[index];
    if (p->is_binary == false)
    {
        *de_dt = 0.0;
        return 0;
    }

    *de_dt = dot3(p->e_vec_unit,p->de_vec_dt);

    return 0;
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
     bool include_hexadecupole_order_binary_pair_terms_, bool include_dotriacontupole_order_binary_pair_terms_,  bool include_double_averaging_corrections_,
     bool include_flybys_, int flybys_reference_binary_, bool flybys_correct_for_gravitational_focussing_, int flybys_velocity_distribution_, int flybys_mass_distribution_,
     double flybys_mass_distribution_lower_value_, double flybys_mass_distribution_upper_value_, double flybys_encounter_sphere_radius_, 
     double flybys_stellar_density_, double flybys_stellar_relative_velocity_dispersion_,
     int binary_evolution_CE_energy_flag_, int binary_evolution_CE_spin_flag_, \
     double mstar_gbs_tolerance_, double mstar_collision_tolerance_)
{
     relative_tolerance = relative_tolerance_;
     absolute_tolerance_eccentricity_vectors = absolute_tolerance_eccentricity_vectors_;
     include_quadrupole_order_terms = include_quadrupole_order_terms_;
     include_octupole_order_binary_pair_terms = include_octupole_order_binary_pair_terms_;
     include_octupole_order_binary_triplet_terms = include_octupole_order_binary_triplet_terms_;
     include_hexadecupole_order_binary_pair_terms = include_hexadecupole_order_binary_pair_terms_;
     include_dotriacontupole_order_binary_pair_terms = include_dotriacontupole_order_binary_pair_terms_;
     include_double_averaging_corrections = include_double_averaging_corrections_;
     
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
     
     binary_evolution_CE_energy_flag = binary_evolution_CE_energy_flag_;
     binary_evolution_CE_spin_flag = binary_evolution_CE_spin_flag_;
     
     mstar_gbs_tolerance = mstar_gbs_tolerance_;
     mstar_collision_tolerance = mstar_collision_tolerance_;
     
     
     //printf("PARAMS %g %g %d %d %d %d %d\n",relative_tolerance,absolute_tolerance_eccentricity_vectors,include_quadrupole_order_terms,include_octupole_order_binary_pair_terms,include_octupole_order_binary_triplet_terms,include_hexadecupole_order_binary_pair_terms,include_dotriacontupole_order_binary_pair_terms);

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



/***********
 * Testing *
 * ********/

int unit_tests_interface()
{
    int flag=0;
    flag += test_tools();
    flag += test_binary_evolution();
    flag += test_collisions();
    
    return flag;
}

int determine_compact_object_merger_properties_interface(double m1, double m2, double chi1, double chi2, \
    double spin_vec_1_unit_x, double spin_vec_1_unit_y, double spin_vec_1_unit_z, \
    double spin_vec_2_unit_x, double spin_vec_2_unit_y, double spin_vec_2_unit_z, \
    double h_vec_unit_x, double h_vec_unit_y, double h_vec_unit_z, \
    double e_vec_unit_x, double e_vec_unit_y, double e_vec_unit_z, \
    double *v_recoil_vec_x, double *v_recoil_vec_y, double *v_recoil_vec_z, \
    double *alpha_vec_final_x, double *alpha_vec_final_y, double *alpha_vec_final_z, \
    double *M_final)
{
    
    double spin_vec_1_unit[3] = {spin_vec_1_unit_x,spin_vec_1_unit_y,spin_vec_1_unit_z};
    double spin_vec_2_unit[3] = {spin_vec_2_unit_x,spin_vec_2_unit_y,spin_vec_2_unit_z};
    double h_vec_unit[3] = {h_vec_unit_x,h_vec_unit_y,h_vec_unit_z};
    double e_vec_unit[3] = {e_vec_unit_x,e_vec_unit_y,e_vec_unit_z};
    double v_recoil_vec[3],alpha_vec_final[3];
 
//    printf("determine_compact_object_merger_properties_interface %g %g %g %g\n",m1,m2,chi1,chi2);
//    printf("s1 %g %g %g s2 %g %g %g\n",spin_vec_1_unit_x, spin_vec_1_unit_y, spin_vec_1_unit_z,spin_vec_2_unit_x, spin_vec_2_unit_y, spin_vec_2_unit_z);
    determine_compact_object_merger_properties(m1,m2,chi1,chi2,spin_vec_1_unit,spin_vec_2_unit,h_vec_unit,e_vec_unit,v_recoil_vec,alpha_vec_final,M_final);

    *v_recoil_vec_x = v_recoil_vec[0];
    *v_recoil_vec_y = v_recoil_vec[1];
    *v_recoil_vec_z = v_recoil_vec[2];
    *alpha_vec_final_x = alpha_vec_final[0];
    *alpha_vec_final_y = alpha_vec_final[1];
    *alpha_vec_final_z = alpha_vec_final[2];
    
    return 0;
}

}
