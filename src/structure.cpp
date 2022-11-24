/* MSE */

#include "structure.h"
#include "evolve.h"
#include <stdio.h>

extern "C"
{

int determine_binary_parents_and_levels(ParticlesMap *particlesMap, int *N_bodies, int *N_binaries, int *N_root_finding, int *N_ODE_equations)
{
    *N_bodies = 0;
    *N_binaries = 0;
    *N_root_finding = 0;
    *N_ODE_equations = 0;
    
    /* determine parent for each particle */
    ParticlesMapIterator it_p,it_q;
    for (it_p = particlesMap->begin(); it_p != particlesMap->end(); it_p++)
    {
        Particle *P_p = (*it_p).second;

        #ifdef VERBOSE
        if (verbose_flag > 2)
        {
            printf("structure.cpp -- determine_binary_parents_and_levels -- determine_binary_parents_and_levels reset parent %d index %d is_binary %d C1 %d C2 %d\n",P_p->parent,P_p->index,P_p->is_binary,P_p->child1,P_p->child2);
        }
        #endif
        P_p->parent = -1;
    }

    for (it_p = particlesMap->begin(); it_p != particlesMap->end(); it_p++)
    {
        Particle *P_p = (*it_p).second;
        
        if (P_p->is_binary == true)
        {
            (*N_binaries)++;
            
            /* root finding */
            if (P_p->check_for_secular_breakdown == true)
            {
                (*N_root_finding)++;
            }
            if (P_p->check_for_dynamical_instability == true)
            {
                (*N_root_finding)++;
            }
            if (P_p->check_for_physical_collision_or_orbit_crossing == true)
            {
                (*N_root_finding)++;
            }
            if (P_p->check_for_minimum_periapse_distance == true)
            {
                (*N_root_finding)++;
            }
            if (P_p->check_for_GW_condition == true)
            {
                (*N_root_finding)++;
            }
            if (P_p->check_for_entering_LISA_band == true)
            {
                (*N_root_finding)++;
            }

            /* parents and siblings */
            for (it_q = particlesMap->begin(); it_q != particlesMap->end(); it_q++)
            {
                Particle *P_q = (*it_q).second;
                if (P_q->index == P_p->child1)
                {
                    P_q->parent = P_p->index;
                    P_q->sibling = P_p->child2;
                }
                if (P_q->index == P_p->child2)
                {
                    
                    P_q->parent = P_p->index;
                    P_q->sibling = P_p->child1;
                }
            }
            
            /* count number of ODE eqs */
            if (P_p->integration_method==0)
            {
                (*N_ODE_equations) += 6;
            }
            else if (P_p->integration_method==1)
            {
                (*N_ODE_equations) += 10;
            }
            else if (P_p->integration_method==2)
            {
                (*N_ODE_equations) += 6;
            }

        }
        else
        {
            (*N_bodies)++;
            (*N_ODE_equations) += 5;
            
            /* root finding */
            if (P_p->check_for_RLOF_at_pericentre == true)
            {
                (*N_root_finding)++;
            }
            
        }
    }

    /* determine levels and set of parents for each particle */
    int highest_level = 0;
    for (it_p = particlesMap->begin(); it_p != particlesMap->end(); it_p++)
    {
        Particle *P_p = (*it_p).second;

        P_p->connecting_child_in_parents.clear();
        P_p->parents.clear();
        P_p->level=0;
        
        int child = P_p->index;
        int parent = P_p->parent;

        if (parent != -1) /* if parent == -1, P_p is the `top' binary, for which level=0; note: this loop will never be entered if there are no binaries in the system */
        {
            while (parent != -1) /* search parents until reaching the top binary */
            {
                for (it_q = particlesMap->begin(); it_q != particlesMap->end(); it_q++)
                {
                    Particle *P_q = (*it_q).second;
                    if (P_q->index == parent)
                    {
                        if (child==P_q->child1)
                        {
                            P_p->connecting_child_in_parents.push_back(1);
                        }
                        else if (child==P_q->child2)
                        {
                            P_p->connecting_child_in_parents.push_back(2);
                        }
                        P_p->parents.push_back(parent);
                        P_p->level++;
                        
                        child = P_q->index;
                        parent = P_q->parent;
                        
                        #ifdef VERBOSE
                        if (verbose_flag > 2)
                        {
                            printf("structure.cpp -- determine_binary_parents_and_levels -- p %d q %d %d child %d\n",P_p->index,P_q->index,P_p->level,child);
                        }
                        #endif
                        break;
                    }
                }
            }
        }
        highest_level = CV_max(P_p->level,highest_level);
        
    }

    /* write highest level to all particles -- needed for function set_binary_masses_from_body_masses */
    for (it_p = particlesMap->begin(); it_p != particlesMap->end(); it_p++)
    {
        Particle *p = (*it_p).second;
        
        p->highest_level = highest_level;
        
        if (p->parent == -1)
        {
            p->is_bound = false;
        }
        else
        {
            p->is_bound = true;
        }
    }

    return 0;
}

void set_binary_masses_from_body_masses(ParticlesMap *particlesMap)
{

    /* set binary masses -- to ensure this happens correctly, do this from highest level to lowest level */
    ParticlesMapIterator it_p;

    int highest_level = particlesMap->begin()->second->highest_level;
    int level=highest_level;

    while (level > -1)
    {
        for (it_p = particlesMap->begin(); it_p != particlesMap->end(); it_p++)
        {
            Particle *P_p = (*it_p).second;
            highest_level = P_p->highest_level;
            
            if ((P_p->is_binary == true) && (P_p->level == level))
            {
                Particle *P_p_child1 = (*particlesMap)[P_p->child1];
                Particle *P_p_child2 = (*particlesMap)[P_p->child2];
                
                #ifdef VERBOSE
                if (verbose_flag > 2)
                {
                    printf("structure.cpp -- set_binary_masses_from_body_masses -- parent %d C1 %d C2 %d \n",P_p->index,P_p->child1,P_p->child2);
                }
                #endif
                        
                /* these quantities are used in ODE_system.cpp */
                P_p->child1_mass = P_p_child1->mass;
                P_p->child2_mass = P_p_child2->mass;
                P_p->mass = P_p->child1_mass + P_p->child2_mass;

                P_p->child1_mass_dot_wind = P_p_child1->mass_dot_wind;
                P_p->child2_mass_dot_wind = P_p_child2->mass_dot_wind;
                P_p->mass_dot_wind = P_p->child1_mass_dot_wind + P_p->child2_mass_dot_wind;

                P_p->child1_mass_dot_wind_accretion = P_p_child1->mass_dot_wind_accretion;
                P_p->child2_mass_dot_wind_accretion = P_p_child2->mass_dot_wind_accretion;
                P_p->mass_dot_wind_accretion = P_p->child1_mass_dot_wind_accretion + P_p->child2_mass_dot_wind_accretion;

                P_p->child1_mass_dot_adiabatic_ejection = P_p_child1->mass_dot_adiabatic_ejection;
                P_p->child2_mass_dot_adiabatic_ejection = P_p_child2->mass_dot_adiabatic_ejection;
                P_p->mass_dot_adiabatic_ejection = P_p->child1_mass_dot_adiabatic_ejection + P_p->child2_mass_dot_adiabatic_ejection;

                P_p->child1_mass_plus_child2_mass = P_p->child1_mass + P_p->child2_mass;
                P_p->child1_mass_minus_child2_mass = P_p->child1_mass - P_p->child2_mass;
                P_p->child1_mass_times_child2_mass = P_p->child1_mass*P_p->child2_mass;
                
                P_p->mu = P_p->child1_mass * P_p->child2_mass / P_p->mass;
                
                P_p->delta_child1_mass_adiabatic_mass_loss = P_p_child1->delta_m_adiabatic_mass_loss;
                P_p->delta_child2_mass_adiabatic_mass_loss = P_p_child2->delta_m_adiabatic_mass_loss;
                P_p->delta_m_adiabatic_mass_loss = P_p->delta_child1_mass_adiabatic_mass_loss + P_p->delta_child2_mass_adiabatic_mass_loss;


                #ifdef VERBOSE
                if (verbose_flag > 2)
                {
                    printf("structure.cpp -- set_binary_masses_from_body_masses -- level %d m %g highest_level %d\n",level,P_p->mass,highest_level);
                }
                #endif
            }
        }
        level--;
    }

    /* determine total system mass -- needed for hyperbolic external orbits */
    double total_system_mass;
    for (it_p = particlesMap->begin(); it_p != particlesMap->end(); it_p++)
    {
        Particle *P_p = (*it_p).second;

        if (P_p->level==0) /* lowest-level binary */
        {
            total_system_mass = P_p->child1_mass + P_p->child2_mass;
            break;
        }
    }

    for (it_p = particlesMap->begin(); it_p != particlesMap->end(); it_p++)
    {
        Particle *P_p = (*it_p).second;
        
        P_p->total_system_mass = total_system_mass;
    }
}

void update_structure(ParticlesMap *particlesMap, int integration_flag)
{
    if (integration_flag != 0)
    {
        return;
    }
    
    int N_bodies, N_binaries,N_root_finding,N_ODE_equations;
    determine_binary_parents_and_levels(particlesMap,&N_bodies,&N_binaries,&N_root_finding,&N_ODE_equations);
    
    if (N_binaries>0)
    {
        set_binary_masses_from_body_masses(particlesMap);
    }
    
}

void set_positions_and_velocities(ParticlesMap *particlesMap) /* TO DO: add to name of function: "_of_all_bodies" */
{
    /* Compute and set the positions and velocities of all bodies */
    /* By default, sample orbital phases randomly 
     * if particle.sample_orbital_phases_randomly == False: look for particle.true_anomaly */
    int N_bodies, N_binaries, N_root_finding, N_ODE_equations;
    determine_binary_parents_and_levels(particlesMap, &N_bodies, &N_binaries, &N_root_finding, &N_ODE_equations);
    set_binary_masses_from_body_masses(particlesMap);

    double true_anomaly;

    int index = 0;
    int i;
    
    double r[3],v[3],r_parent[3],v_parent[3],r_child1[3],v_child1[3],r_child2[3],v_child2[3];
    double parent_mass,child1_mass,child2_mass;
    double e;
    
    /*/ Go from the top of the system (level=0) downwards */
    
    ParticlesMapIterator it;

    int highest_level = particlesMap->begin()->second->highest_level;
    int level = 0;
    while (level < highest_level)
    {
        for (it = particlesMap->begin(); it != particlesMap->end(); it++)
        {
            Particle *parent = (*it).second;
            if ((parent->is_binary == true) && (parent->level == level))
            {
                Particle *child1 = (*particlesMap)[parent->child1];
                Particle *child2 = (*particlesMap)[parent->child2];
                
                child1_mass = child1->mass;
                child2_mass = child2->mass;
                parent_mass = parent->mass;
                e = norm3(parent->e_vec);
    
                if (e < 0 or e>= 1.0)
                {
                    continue;
                }

                if (parent->sample_orbital_phase_randomly == false)
                {
                    true_anomaly = parent->true_anomaly;
                }
                else
                {
                    true_anomaly = sample_random_true_anomaly(e);
                }
                
                check_number(true_anomaly,    "structure.cpp -- set_positions_and_velocities","true_anomaly", true);
                
                from_orbital_vectors_to_cartesian(
                    child1_mass,child2_mass,
                    parent->e_vec,parent->h_vec,
                    true_anomaly,
                    r,v);

                check_number(r[0],    "structure.cpp -- set_positions_and_velocities","r[0]", true);
                check_number(r[1],    "structure.cpp -- set_positions_and_velocities","r[1]", true);
                check_number(r[2],    "structure.cpp -- set_positions_and_velocities","r[2]", true);

                if (parent->level == 0)
                {
                    for (i=0; i<3; i++)
                    {
                        r_parent[i] = parent->R_vec[i];
                        v_parent[i] = parent->V_vec[i];
                    }
                }
                else
                {
                    get_position_and_velocity_vectors_from_particle(parent,r_parent,v_parent);
                }
                
                for (i=0; i<3; i++)
                {
                    r_child1[i] = r_parent[i] + (child2_mass/parent_mass)*r[i];
                    v_child1[i] = v_parent[i] + (child2_mass/parent_mass)*v[i];
                    
                    r_child2[i] = r_parent[i] - (child1_mass/parent_mass)*r[i];
                    v_child2[i] = v_parent[i] - (child1_mass/parent_mass)*v[i];
                }
                set_position_and_velocity_vectors_in_particle(child1,r_child1,v_child1);
                set_position_and_velocity_vectors_in_particle(child2,r_child2,v_child2);
                parent->true_anomaly = true_anomaly;
            }
        }
        level++;
    }
}

void update_positions_unbound_bodies(ParticlesMap *particlesMap, double time_step)
{
    int i;
    ParticlesMapIterator it;
    double m0,m_dot_div_m0,factor_V,factor_R;
    for (it = particlesMap->begin(); it != particlesMap->end(); it++)
    {
        Particle *p = (*it).second;
        if ((p->is_binary == false) && (p->is_bound == false))
        {
            m0 = p->mass - p->mass_dot_wind * time_step; // mass at the beginning of the timestep
            m_dot_div_m0 = p->mass_dot_wind/m0;
            
            if (fabs(m_dot_div_m0) <= epsilon)
            {
                for (i=0; i<3; i++)
                {
                    p->R_vec[i] += p->V_vec[i]*time_step;
                }
            }
            else
            {
                #ifdef IGNORE
                factor_V = exp(- m_dot_div_m * time_step);
                factor_R = (1.0 - factor_V)/m_dot_div_m;

                for (i=0; i<3; i++)
                {
                    p->R_vec[i] += p->V_vec[i]*factor_R;
                    p->V_vec[i] *= factor_V;
                }
                #endif
                
                factor_V = 1.0/(1.0 + m_dot_div_m0 * time_step);
                factor_R = log(1.0 + m_dot_div_m0 * time_step)/m_dot_div_m0;
                
                #ifdef VERBOSE
                if (verbose_flag > 1)
                {
                    printf("structure.cpp -- update_positions_unbound_bodies -- id %d mdot %g T %g %g %g %g \n",p->index,p->mass_dot_wind,m_dot_div_m0,factor_V,factor_R,time_step);
                }
                #endif
                
                for (i=0; i<3; i++)
                {
                    p->R_vec[i] += p->V_vec[i]*factor_R; // V_vec is still the old value here (at the beginning of the timestep)
                    p->V_vec[i] *= factor_V;
                }

            }
        }
    }
}


void update_masses_positions_and_velocities_of_all_binaries(ParticlesMap *particlesMap)
{
    set_binary_masses_from_body_masses(particlesMap);

    /* set binary positions and velocities -- to ensure this happens correctly, do this from highest level to lowest level */
 
    int i;
    double child1_mass,child2_mass;
    double r[3],v[3];
    double r_child1[3],v_child1[3];
    double r_child2[3],v_child2[3];
    
    ParticlesMapIterator it_p;
    
    int highest_level = particlesMap->begin()->second->highest_level;
    int level=highest_level;
    while (level > -1)
    {
        for (it_p = particlesMap->begin(); it_p != particlesMap->end(); it_p++)
        {
            Particle *p = (*it_p).second;
            if ((p->is_binary == true) && (p->level == level))
            {
                Particle *child1 = (*particlesMap)[p->child1];
                Particle *child2 = (*particlesMap)[p->child2];
                
                
                get_position_and_velocity_vectors_from_particle(child1,r_child1,v_child1);
                get_position_and_velocity_vectors_from_particle(child2,r_child2,v_child2);
                
                child1_mass = child1->mass;
                child2_mass = child2->mass;

                for (i=0; i<3; i++)
                {
                    r[i] = (r_child1[i]*child1_mass + r_child2[i]*child2_mass)/(child1_mass+child2_mass);
                    v[i] = (v_child1[i]*child1_mass + v_child2[i]*child2_mass)/(child1_mass+child2_mass);
                }
            
                set_position_and_velocity_vectors_in_particle(p,r,v);
            }
        }
        level--;
    }
}


void update_orbital_vectors_in_binaries_from_positions_and_velocities(ParticlesMap *particlesMap)
{
    #ifdef VERBOSE
    if (verbose_flag > 1)
    {
        printf("structure.cpp -- update_orbital_vectors_in_binaries_from_positions_and_velocities\n");
    }
    #endif
                
    int index = 0;
    int i;
    
    double r[3],v[3],r_child1[3],v_child1[3],r_child2[3],v_child2[3];
    double child1_mass,child2_mass;
    double e_vec[3],h_vec[3];
    
    /*/ Go from the top of the system (level=0) downwards */
    
    ParticlesMapIterator it;

    int highest_level = particlesMap->begin()->second->highest_level;
    int level = 0;
    while (level < highest_level)
    {
        for (it = particlesMap->begin(); it != particlesMap->end(); it++)
        {
            Particle *parent = (*it).second;
            if ((parent->is_binary == true) && (parent->level == level))
            {
                Particle *child1 = (*particlesMap)[parent->child1];
                Particle *child2 = (*particlesMap)[parent->child2];
                
                child1_mass = child1->mass;
                child2_mass = child2->mass;
                
                get_position_and_velocity_vectors_from_particle(child1,r_child1,v_child1);
                get_position_and_velocity_vectors_from_particle(child2,r_child2,v_child2);

                for (i=0; i<3; i++)
                {
                    r[i] = r_child1[i] - r_child2[i];
                    v[i] = v_child1[i] - v_child2[i];
                }
                
                from_cartesian_to_orbital_vectors(
                    child1_mass,child2_mass,
                    r,v,
                    parent->e_vec,parent->h_vec,&parent->true_anomaly);
            }
        }
        level++;
    }
}

void determine_internal_mass_and_semimajor_axis(ParticlesMap *particlesMap)
{
    #ifdef VERBOSE
    if (verbose_flag > 1)
    {
        printf("structure.cpp -- determine_internal_mass_and_semimajor_axis\n");
    }
    #endif

    int N_bodies,N_binaries,N_root_finding,N_ODE_equations;
    determine_binary_parents_and_levels(particlesMap,&N_bodies,&N_binaries,&N_root_finding,&N_ODE_equations);
    set_binary_masses_from_body_masses(particlesMap);

    double h_tot_vec[3];
    double semimajor_axis,eccentricity,inclination,argument_of_pericenter,longitude_of_ascending_node;
    
    ParticlesMapIterator it_p;
    for (it_p = particlesMap->begin(); it_p != particlesMap->end(); it_p++)
    {
        Particle *p = (*it_p).second;
        if (p->is_binary == true and p->parent == -1)
        {
            flybys_internal_mass = p->mass;
            
            flybys_reference_binary = p->index;
            
            /* Below: a bit overkill to compute the semimajor axis, but uses consistent functions for orbital elements */
            compute_h_tot_vector(particlesMap,h_tot_vec);
            compute_orbital_elements_from_orbital_vectors(p->child1_mass, p->child2_mass, h_tot_vec, \
                p->e_vec[0],p->e_vec[1],p->e_vec[2],p->h_vec[0],p->h_vec[1],p->h_vec[2],
                &semimajor_axis, &eccentricity, &inclination, &argument_of_pericenter, &longitude_of_ascending_node); 
            flybys_internal_semimajor_axis = semimajor_axis;
            
            #ifdef VERBOSE
            if (verbose_flag > 1)
            {
                printf("structure.cpp -- determine_internal_mass_and_semimajor_axis -- flybys_reference_binary %d M_int %g a_int %g\n",flybys_reference_binary,flybys_internal_mass,flybys_internal_semimajor_axis);
            }
            #endif
        }
    }
}

double determine_longest_orbital_period_in_system(ParticlesMap *particlesMap)
{
    double P_orb;
    ParticlesMapIterator it_p;
    for (it_p = particlesMap->begin(); it_p != particlesMap->end(); it_p++)
    {
        Particle *p = (*it_p).second;
        if (p->level==0) /* lowest-level binary */
        {
            P_orb = compute_orbital_period_from_semimajor_axis(p->mass,p->a);
            break;
        }
    }
    
    return P_orb;
}

double determine_shortest_orbital_period_in_system(ParticlesMap *particlesMap)
{
    double P_orb;
    double P_orb_min = 1.0e100;
    ParticlesMapIterator it_p;
    for (it_p = particlesMap->begin(); it_p != particlesMap->end(); it_p++)
    {
        Particle *p = (*it_p).second;
        if (p->is_binary == true)
        {
            P_orb = compute_orbital_period_from_semimajor_axis(p->mass,p->a);
            if (P_orb < P_orb_min)
            {
                P_orb_min = P_orb;
            }
        }
    }
    
    return P_orb_min;
}

int determine_number_of_bodies_in_system(ParticlesMap *particlesMap)
{
    int N_bodies = 0;
    ParticlesMapIterator it_p;
    for (it_p = particlesMap->begin(); it_p != particlesMap->end(); it_p++)
    {
        Particle *p = (*it_p).second;
        
        if (p->is_binary == false)
        {
            N_bodies++;
        }
    }
    
    return N_bodies;
}

void set_up_derived_quantities(ParticlesMap *particlesMap)
{
    update_structure(particlesMap, 0);
    
    /* These are derived quantities that are often used in the EOM, so they are calculated here once for speed up */
    int i;
    ParticlesMapIterator it_p;
    for (it_p = particlesMap->begin(); it_p != particlesMap->end(); it_p++)
    {
        Particle *p = (*it_p).second;
        if (p->is_binary == false)
        {
            p->spin_vec_norm = norm3(p->spin_vec);
            if (p->spin_vec_norm <= epsilon)
            {
                p->spin_vec_norm = epsilon;
            }

            for (i=0; i<3; i++)
            {
                p->spin_vec_unit[i] = p->spin_vec[i]/p->spin_vec_norm;
            }

            p->chi = compute_spin_parameter_from_spin_frequency(p->mass,p->spin_vec_norm);
        }
        else
        {
            p->e = norm3(p->e_vec);
            p->h = norm3(p->h_vec);

            if (p->e <= epsilon)
            {
                p->e = epsilon;
            }

            for (i=0; i<3; i++)
            {        
                p->e_vec_unit[i] = p->e_vec[i]/p->e;
                p->h_vec_unit[i] = p->h_vec[i]/p->h;
            }
            
            cross3(p->h_vec_unit,p->e_vec_unit,p->q_vec_unit);
            
            p->e_p2 = p->e*p->e;
            p->j_p2 = 1.0 - p->e_p2;
            p->j = sqrt(p->j_p2);
            p->j_p3 = p->j*p->j_p2;
            p->j_p4 = p->j*p->j_p3;
            p->j_p5 = p->j*p->j_p4;

            p->a = p->h*p->h*p->child1_mass_plus_child2_mass/( CONST_G*p->child1_mass_times_child2_mass*p->child1_mass_times_child2_mass*p->j_p2 );

            p->r = norm3(p->r_vec);
            p->r_p2 = p->r*p->r;
            p->r_p3 = p->r*p->r_p2;
            p->r_pm1 = 1.0/p->r;
            p->r_pm2 = p->r_pm1*p->r_pm1;
            p->r_pm3 = p->r_pm1*p->r_pm2;
            
        }
    }
}

double compute_rp_out_crit_MA01(double a_in, double q_out, double e_out, double rel_INCL)
{
    return a_in*2.8*pow( (1.0+q_out)*(1.0+e_out)/sqrt(1.0-e_out), 2.0/5.0) * (1.0 - 0.3*rel_INCL/M_PI);
}


bool check_system_for_dynamical_stability(ParticlesMap *particlesMap, int *integration_flag)
{
    if (*integration_flag > 0)
    {
        return false;
    }

    set_up_derived_quantities(particlesMap); /* for a, e, etc */
    bool stable = true;
   
    double a_out,e_out,a_in,e_in,M_p,rp_out,rp_out_crit,q_out,rel_INCL;
   
    ParticlesMapIterator it_p;
    for (it_p = particlesMap->begin(); it_p != particlesMap->end(); it_p++)
    {
        Particle *p = (*it_p).second;
        if (p->is_binary == true and p->parent != -1)
        {
            Particle *parent = (*particlesMap)[p->parent];
            Particle *sibling = (*particlesMap)[p->sibling];
                            
            a_out = parent->a;
            e_out = parent->e;
            a_in = p->a;
            e_in = p->e;
            M_p = p->mass;
            rp_out = a_out*(1.0-e_out);
                    
            if (p->dynamical_instability_criterion == 0) /* for mass ratios on the order of unity */
            {
                /* Mardling & Aarseth 2001 */
                q_out = sibling->mass/M_p;
                get_inclination_relative_to_parent(particlesMap,p->index,&rel_INCL);
                rp_out_crit = compute_rp_out_crit_MA01(a_in,q_out,e_out,rel_INCL);
                if (rp_out < rp_out_crit)
                {
                    stable = false;
                }

                #ifdef VERBOSE
                if (verbose_flag > 1)
                {
                    printf("structure.cpp -- check_system_for_dynamical_stability -- rp_out %g rp_out_crit %g a_in %g e_in %g M_p %g q_out %g rel_INCL %g\n",rp_out,rp_out_crit,a_in,e_in,M_p,q_out,rel_INCL);
                }
                #endif
                
            }
            else
            {
                printf("structure.cpp -- check_system_for_dynamical_stability -- p %d dynamical_instability_criterion %d not supported!\n",p->index,p->dynamical_instability_criterion);
                //exit(-1);
                error_code = 18;
                longjmp(jump_buf,1);
            }
        }
    }
    
   
    return stable;
}

void handle_instantaneous_and_adiabatic_mass_changes_in_orbit(ParticlesMap *particlesMap, Particle *star1, Particle *star2, double Delta_m1, double Delta_m2, double mass_loss_timescale, int *integration_flag)
{

    set_old_parameters_for_adiabatic_mass_loss(particlesMap); /* copy the old masses and h & e vectors for use of adiabatic mass loss below */

    /* Next, assume instantaneous mass loss for ALL orbits. 
     * Note: `binary' will be overwritten at the very end below. 
     * Other orbits, for which mass loss is expected to be adiabatic,
     * are adjusted below in compute_new_orbits_assuming_adiabatic_mass_loss. */

    star1->instantaneous_perturbation_delta_mass = Delta_m1; // final minus initial mass
    star2->instantaneous_perturbation_delta_mass = Delta_m2; // final minus initial mass        

    apply_instantaneous_mass_changes_and_kicks(particlesMap, integration_flag);
    
    /* Lastly, `correct' the orbits (h & e vectors) for which the effect of mass loss was actually expected to be adiabatic,
     * i.e., time of CE >> t_orb */
    star1->delta_m_adiabatic_mass_loss = Delta_m1;
    star2->delta_m_adiabatic_mass_loss = Delta_m2;
    compute_new_orbits_assuming_adiabatic_mass_loss(particlesMap, mass_loss_timescale);
}

void set_old_parameters_for_adiabatic_mass_loss(ParticlesMap *particlesMap)
{
    int i;
    ParticlesMapIterator it;
    for (it = particlesMap->begin(); it != particlesMap->end(); it++)
    {
        Particle *b = (*it).second;
        b->delta_m_adiabatic_mass_loss = 0.0;
        if (b->is_binary == true)
        {
            b->P_orb_adiabatic_mass_loss = compute_orbital_period_from_semimajor_axis(b->mass,b->a);
            b->child1_mass_adiabatic_mass_loss = b->child1_mass;
            b->child2_mass_adiabatic_mass_loss = b->child2_mass;
            for (i=0; i<3; i++)
            {
                b->h_vec_old_adiabatic_mass_loss[i] = b->h_vec[i];
                b->e_vec_old_adiabatic_mass_loss[i] = b->e_vec[i];
            }
        }
    }
}

void compute_new_orbits_assuming_adiabatic_mass_loss(ParticlesMap *particlesMap, double mass_loss_timescale)
{
    set_binary_masses_from_body_masses(particlesMap); // to set the correct delta_m_adiabatic for each orbit

    int i;    
    ParticlesMapIterator it;
    for (it = particlesMap->begin(); it != particlesMap->end(); it++)
    {
        Particle *b = (*it).second;
        if (b->is_binary == true)
        {
            #ifdef VERBOSE
            if (verbose_flag > 1)
            {
                printf("structure.cpp -- compute_new_orbits_assuming_adiabatic_mass_loss -- %d %g %g\n",b->index,b->P_orb_adiabatic_mass_loss , mass_loss_timescale);
            }
            #endif
            
            if (b->P_orb_adiabatic_mass_loss < mass_loss_timescale) // criterion for adiabatic mass loss
            {
                
                double m1_old = b->child1_mass_adiabatic_mass_loss;
                double m2_old = b->child2_mass_adiabatic_mass_loss;
                double delta_m1 = b->delta_child1_mass_adiabatic_mass_loss;
                double delta_m2 = b->delta_child2_mass_adiabatic_mass_loss;
                double h_old = norm3(b->h_vec_old_adiabatic_mass_loss);
                double m1_new = m1_old + delta_m1;
                double m2_new = m2_old + delta_m2;
                double h_new = h_old * ((m1_old + m2_old)/(m1_new + m2_new)) * (m1_new/m1_old) * (m2_new/m2_old);
                
                #ifdef VERBOSE
                if (verbose_flag > 1)
                {
                    printf("structure.cpp -- compute_new_orbits_assuming_adiabatic_mass_loss -- adiabatic %g %g h_old %g h_new %g delta m1 %g delta m2 %g m1_old %g m2_old %g\n",b->P_orb_adiabatic_mass_loss , mass_loss_timescale,h_old,h_new,delta_m1,delta_m2,m1_old,m2_old);
                }
                #endif
                
                for (i=0; i<3; i++)
                {
                    /* Directions of h and e and magnitude of e do not change for adiabatic mass loss. */
                    b->h_vec[i] = (b->h_vec_old_adiabatic_mass_loss[i] / h_old) * h_new;
                    b->e_vec[i] = b->e_vec_old_adiabatic_mass_loss[i];
                }
            }
        }
    }
}

void compute_new_positions_and_velocities_given_new_semimajor_axis_and_eccentricity(double M1_old, double R1_vec_old[3], double V1_vec_old[3], double M2_old, double R2_vec_old[3], double V2_vec_old[3], double M1, double R1_vec[3], double V1_vec[3], double M2, double R2_vec[3], double V2_vec[3], double a, double e)
{
    /* Assumes the directions of the h & e vectors are unchanged. */
    
    int i;
    double M_old = M1_old + M2_old;
    double M = M1 + M2;
    double R_CM_vec[3],V_CM_vec[3];
    double r_vec_old[3],v_vec_old[3];
    
    for (i=0; i<3; i++)
    {
        R_CM_vec[i] = (R1_vec_old[i]*M1_old + R2_vec_old[i]*M2_old) / M_old;
        V_CM_vec[i] = (V1_vec_old[i]*M1_old + V2_vec_old[i]*M2_old) / M_old;
        r_vec_old[i] = R1_vec_old[i] - R2_vec_old[i];
        v_vec_old[i] = V1_vec_old[i] - V2_vec_old[i];
    }
    
    /* Compute the old e & h vectors (they determine the direction of the new orbit) */
    double e_vec_old[3],h_vec_old[3];
    double true_anomaly_old;
    from_cartesian_to_orbital_vectors(M1_old, M2_old, r_vec_old, v_vec_old, e_vec_old, h_vec_old, &true_anomaly_old);
    double e_old = norm3(e_vec_old);
    double h_old = norm3(h_vec_old);
    
    /* Compute the new e & h vectors */
    double e_vec[3],h_vec[3];
    double h = compute_h_from_a(M1,M2,a,e);
    for (i=0; i<3; i++)
    {
        e_vec[i] = (e_vec_old[i]/e_old) * e;
        h_vec[i] = (h_vec_old[i]/h_old) * h;
    }
    
    /* Compute new absolute positions and velocities */
    double r_vec[3],v_vec[3];
    from_orbital_vectors_to_cartesian(M1, M2, e_vec, h_vec, true_anomaly_old, r_vec, v_vec);
    for (i=0; i<3; i++)
    {
        R1_vec[i] = R_CM_vec[i] + (M2/M) * r_vec[i];
        V1_vec[i] = V_CM_vec[i] + (M2/M) * v_vec[i];

        R2_vec[i] = R_CM_vec[i] - (M1/M) * r_vec[i];
        V2_vec[i] = V_CM_vec[i] - (M1/M) * v_vec[i];
    }
}


void handle_gradual_mass_loss_event_in_system(ParticlesMap *particlesMap, Particle *star1, Particle *star2, double M1, double M1_old, double M2, double M2_old, double mass_loss_timescale, \
    double r_vec[3], double v_vec[3], double initial_R_CM[3], double initial_V_CM[3], double final_R_CM[3], double final_V_CM[3], double final_momentum[3])
{
    /* Take into account effect of mass loss on the rest of the system */

    int i;
    std::vector<double> masses;
    std::vector<double> delta_masses;
    std::vector<std::vector<double>> R_vecs;
    std::vector<std::vector<double>> V_vecs;

    /* Add the mass losing bodies */
    masses.push_back( M1_old + M2_old );
    delta_masses.push_back( (M1 + M2) - (M1_old + M2_old) );
    R_vecs.push_back( {initial_R_CM[0],initial_R_CM[1],initial_R_CM[2]} );
    V_vecs.push_back( {initial_V_CM[0],initial_V_CM[1],initial_V_CM[2]} );

    #ifdef VERBOSE
    if (verbose_flag > 0)
    {
        printf("structure.cpp -- handle_gradual_mass_loss_event_in_system -- m %g dm %g R %g %g %g V %g %g %g\n",masses[0],delta_masses[0],R_vecs[0][0],R_vecs[0][1],R_vecs[0][2],V_vecs[0][0],V_vecs[0][1],V_vecs[0][2]);
    }
    #endif

    /* Add the other bodies (assumed not to lose mass) */
    ParticlesMapIterator it_p;
    for (it_p = particlesMap->begin(); it_p != particlesMap->end(); it_p++)
    {
        Particle *p = (*it_p).second;
        if (p->is_binary == false and p->index != star1->index and p->index != star2->index)
        {
            
            if (separation_between_vectors(initial_R_CM, p->R_vec) < nbody_maximum_separation_for_inclusion)
            {
                masses.push_back(p->mass);
                delta_masses.push_back(0.0);
                R_vecs.push_back( {p->R_vec[0],p->R_vec[1],p->R_vec[2]} );
                V_vecs.push_back( {p->V_vec[0],p->V_vec[1],p->V_vec[2]} );
                
                #ifdef VERBOSE
                if (verbose_flag > 0)
                {
                    printf("structure.cpp -- handle_gradual_mass_loss_event_in_system -- %g %g %g %g %g %g %g\n",p->mass,p->R_vec[0],p->R_vec[1],p->R_vec[2],p->V_vec[0],p->V_vec[1],p->V_vec[2]);
                }
                #endif
            }
        }
    }
    
    /* Integrate non-trivial systems */
    if (masses.size() > 1)
    {
        integrate_nbody_system_with_mass_loss(mass_loss_timescale, binary_evolution_CE_mass_loss_Nsteps, masses, delta_masses, R_vecs, V_vecs);
    }
    else
    {
        for (i=0; i<3; i++)
        {
            final_R_CM[i] = initial_R_CM[i];
            final_V_CM[i] = initial_V_CM[i];
            final_momentum[i] = final_V_CM[i] * (M1_old + M2_old);
        }

        return;
    }
    
    auto it_R = R_vecs.begin();
    auto it_V = V_vecs.begin();

    for (i=0; i<3; i++)
    {
        final_R_CM[i] = (*it_R)[i];
        final_V_CM[i] = (*it_V)[i];
        final_momentum[i] = final_V_CM[i] * (M1_old + M2_old);
    }
    
    #ifdef VERBOSE
    if (verbose_flag > 0)
    {
        printf("structure.cpp -- handle_gradual_mass_loss_event_in_system -- final_R_CM %g %g %g final V_CM %g %g %g\n",final_R_CM[0],final_R_CM[1],final_R_CM[2],final_V_CM[0],final_V_CM[1],final_V_CM[2]);
    }
    #endif
    
    /* Take into account movement of barycenter of the two stars during the mass loss phase */
    for (i=0; i<3; i++)
    {
        star1->R_vec[i] = final_R_CM[i] + (M2_old/(M1_old+M2_old))*r_vec[i];
        star1->V_vec[i] = final_V_CM[i] + (M2_old/(M1_old+M2_old))*v_vec[i];

        star2->R_vec[i] = final_R_CM[i] - (M1_old/(M1_old+M2_old))*r_vec[i];
        star2->V_vec[i] = final_V_CM[i] - (M1_old/(M1_old+M2_old))*v_vec[i];
    }
    
    /* Update positions and velocities of all other bodies */
    for (it_p = particlesMap->begin(); it_p != particlesMap->end(); it_p++)
    {
        Particle *p = (*it_p).second;
        if (p->is_binary == false and p->index != star1->index and p->index != star2->index)
        {
            if (separation_between_vectors(initial_R_CM, p->R_vec) < nbody_maximum_separation_for_inclusion)
            {
                it_R = std::next(it_R,1);
                it_V = std::next(it_V,1);

                for (i=0; i<3; i++)
                {
                    p->R_vec[i] = (*it_R)[i];
                    p->V_vec[i] = (*it_V)[i];
                }
            }
        }
    }

    return;
}


void handle_gradual_mass_loss_event_in_system_triple_CE(ParticlesMap *particlesMap, Particle *star, Particle *companion1, Particle *companion2, double M3, double M3_old, double mass_loss_timescale, \
    double initial_R_CM[3], double initial_V_CM[3], double final_R_CM[3], double final_V_CM[3], double final_momentum[3])
{
    /* Take into account effect of mass loss on the rest of the system */
    
    int i;
    std::vector<double> masses;
    std::vector<double> delta_masses;
    std::vector<std::vector<double>> R_vecs;
    std::vector<std::vector<double>> V_vecs;
    
    masses.push_back( companion1->mass + companion2->mass + M3_old );
    delta_masses.push_back( M3 - M3_old );
    R_vecs.push_back( {initial_R_CM[0],initial_R_CM[1],initial_R_CM[2]} );
    V_vecs.push_back( {initial_V_CM[0],initial_V_CM[1],initial_V_CM[2]} );

    #ifdef VERBOSE
    if (verbose_flag > 0)
    {
        printf("structure.cpp -- handle_gradual_mass_loss_event_in_system_triple_CE -- m %g dm %g R %g %g %g V %g %g %g\n",masses[0],delta_masses[0],R_vecs[0][0],R_vecs[0][1],R_vecs[0][2],V_vecs[0][0],V_vecs[0][1],V_vecs[0][2]); 
    }
    #endif
    
    
    ParticlesMapIterator it_p;
    for (it_p = particlesMap->begin(); it_p != particlesMap->end(); it_p++)
    {
        Particle *p = (*it_p).second;
        if (p->is_binary == false and p->index != star->index and p->index != companion1->index and p->index != companion2->index)
        {

            if (separation_between_vectors(initial_R_CM, p->R_vec) < nbody_maximum_separation_for_inclusion)
            {
                masses.push_back(p->mass);
                delta_masses.push_back(0.0);
                R_vecs.push_back( {p->R_vec[0],p->R_vec[1],p->R_vec[2]} );
                V_vecs.push_back( {p->V_vec[0],p->V_vec[1],p->V_vec[2]} );
                
                #ifdef VERBOSE
                if (verbose_flag > 0)
                {
                    printf("structure.cpp -- handle_gradual_mass_loss_event_in_system_triple_CE -- including index %d m %g R_vec %g %g %g V_vec %g %g %g\n",p->index,p->mass,p->R_vec[0],p->R_vec[1],p->R_vec[2],p->V_vec[0],p->V_vec[1],p->V_vec[2]);
                }
                #endif
            }
        }
    }
    if (masses.size() > 1)
    {
        integrate_nbody_system_with_mass_loss(mass_loss_timescale, binary_evolution_CE_mass_loss_Nsteps, masses, delta_masses, R_vecs, V_vecs);
    }
    else
    {
        #ifdef VERBOSE
        if (verbose_flag > 0)
        {
            printf("structure.cpp -- handle_gradual_mass_loss_event_in_system_triple_CE -- no bodies external to triple interaction subsystem -- skipping\n"); 
        }
        #endif
        
        return;
    }
    
    auto it_R = R_vecs.begin();
    auto it_V = V_vecs.begin();

    for (i=0; i<3; i++)
    {
        final_R_CM[i] = (*it_R)[i];
        final_V_CM[i] = (*it_V)[i];
        final_momentum[i] = final_V_CM[i] * (M3_old);
    }
    
    #ifdef VERBOSE
    if (verbose_flag > 0)
    {
        printf("structure.cpp -- handle_gradual_mass_loss_event_in_system_triple_CE -- final_R_CM %g %g %g final V_CM %g %g %g\n",final_R_CM[0],final_R_CM[1],final_R_CM[2],final_V_CM[0],final_V_CM[1],final_V_CM[2]);
    }
    #endif
   
    /* Update positions and velocities of all other bodies */
    for (it_p = particlesMap->begin(); it_p != particlesMap->end(); it_p++)
    {
        Particle *p = (*it_p).second;
        if (p->is_binary == false and p->index != star->index and p->index != companion1->index and p->index != companion2->index)
        {
            if (separation_between_vectors(initial_R_CM, p->R_vec) < nbody_maximum_separation_for_inclusion)
            {
                it_R = std::next(it_R,1);
                it_V = std::next(it_V,1);

                for (i=0; i<3; i++)
                {
                    p->R_vec[i] = (*it_R)[i];
                    p->V_vec[i] = (*it_V)[i];
                }
            }
        }
    }

    return;
}

void get_initial_binary_orbital_properties_from_position_and_velocity(double R1_vec[3], double V1_vec[3], double R2_vec[3], double V2_vec[3], double M1, double M2, \
        double r_vec[3], double v_vec[3], double initial_momentum[3], double initial_R_CM[3], double initial_V_CM[3], double h_vec[3], double e_vec[3])
{

    for (int i=0; i<3; i++)
    {
        r_vec[i] = R1_vec[i] - R2_vec[i];
        v_vec[i] = V1_vec[i] - V2_vec[i];
        initial_R_CM[i] = (M1 * R1_vec[i] + M2 * R2_vec[i]) / (M1 + M2);
        initial_momentum[i] = M1 * V1_vec[i] + M2 * V2_vec[i];
        initial_V_CM[i] = initial_momentum[i]/(M1 + M2);
    }
    double true_anomaly;
    from_cartesian_to_orbital_vectors(M1,M2,r_vec,v_vec,e_vec,h_vec,&true_anomaly);

    return;
}

void reset_RLOF_flags(ParticlesMap *particlesMap)
{
    ParticlesMapIterator it_p;
    
    for (it_p = particlesMap->begin(); it_p != particlesMap->end(); it_p++)
    {
        Particle *p = (*it_p).second;
        p->RLOF_flag = 0;
    }
}

}
