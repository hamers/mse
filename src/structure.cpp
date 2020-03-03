/*
*/

//#include "types.h"
#include "evolve.h"
#include "structure.h"
//#include <stdio.h>

extern "C"
{

int determine_binary_parents_and_levels(ParticlesMap *particlesMap, int *N_bodies, int *N_binaries, int *N_root_finding)
{

    *N_bodies = 0;
    *N_binaries = 0;
    *N_root_finding = 0;
    
    /* determine parent for each particle */
    ParticlesMapIterator it_p,it_q;
    
    for (it_p = particlesMap->begin(); it_p != particlesMap->end(); it_p++)
    {
        Particle *P_p = (*it_p).second;
        //printf("structure.cpp -- include_tidal_friction_terms %d\n",P_p->include_tidal_friction_terms);
        P_p->parent = -1;
    }
    
    for (it_p = particlesMap->begin(); it_p != particlesMap->end(); it_p++)
    {
        Particle *P_p = (*it_p).second;
        //printf("structure.cpp -- include_tidal_friction_terms %d\n",P_p->include_tidal_friction_terms);
        //P_p->parent = -1;

        if (P_p->is_binary == 1)
        {
            (*N_binaries)++;
            
            /* root finding */
            if (P_p->check_for_secular_breakdown == 1)
            {
                (*N_root_finding)++;
            }
            if (P_p->check_for_dynamical_instability == 1)
            {
                (*N_root_finding)++;
            }
            if (P_p->check_for_physical_collision_or_orbit_crossing == 1)
            {
                (*N_root_finding)++;
            }
            if (P_p->check_for_minimum_periapse_distance == 1)
            {
                (*N_root_finding)++;
            }
            if (P_p->check_for_GW_condition == 1)
            {
                (*N_root_finding)++;
            }

            /* parents and siblings */
            for (it_q = particlesMap->begin(); it_q != particlesMap->end(); it_q++)
            {
                Particle *P_q = (*it_q).second;
                
//                if ((P_q->index == P_p->child1) || (P_q->index == P_p->child2))
//                {
//                    P_q->parent = P_p->index;
//                }
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
        }
        else
        {
            (*N_bodies)++;
            
            /* root finding */
            if (P_p->check_for_RLOF_at_pericentre == 1)
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

        if (parent != -1) /* if parent == -1, P_p is the `top' binary, for which level=0 */
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
//                        printf("p %d q %d %d child %d\n",P_p->index,P_q->index,P_p->level,child);
                        break;
                    }
                }
            }
        }
        highest_level = max(P_p->level,highest_level);
    }
    
    /* write highest level to all particles -- needed for function set_binary_masses_from_body_masses */
    for (it_p = particlesMap->begin(); it_p != particlesMap->end(); it_p++)
    {
        Particle *P_p = (*it_p).second;
        
        P_p->highest_level = highest_level;
    }
    
//    for (it_p = particlesMap->begin(); it_p != particlesMap->end(); it_p++)
//    {
//        Particle *P_p = (*it_p).second;
//        printf("particle %d mass %g e_vec_x %g h_vec_x %g\n",P_p->index,P_p->mass,P_p->e_vec_x,P_p->h_vec_x);
//    }
    return 0;
}

void set_binary_masses_from_body_masses(ParticlesMap *particlesMap)
{
    /* set binary masses -- to ensure this happens correctly, do this from highest level to lowest level */
    ParticlesMapIterator it_p;
    int highest_level = (*particlesMap)[0]->highest_level;
    int level=highest_level;
    while (level > -1)
    {
        for (it_p = particlesMap->begin(); it_p != particlesMap->end(); it_p++)
        {
            Particle *P_p = (*it_p).second;
            if ((P_p->is_binary == 1) && (P_p->level == level))
            {
                Particle *P_p_child1 = (*particlesMap)[P_p->child1];
                Particle *P_p_child2 = (*particlesMap)[P_p->child2];
                
                /* these quantities are used in ODE_system.cpp */
                P_p->child1_mass = P_p_child1->mass;
                P_p->child2_mass = P_p_child2->mass;
                P_p->mass = P_p->child1_mass + P_p->child2_mass;

                P_p->child1_mass_dot_wind = P_p_child1->mass_dot_wind;
                P_p->child2_mass_dot_wind = P_p_child2->mass_dot_wind;
                P_p->mass_dot_wind = P_p->child1_mass_dot_wind + P_p->child2_mass_dot_wind;

                P_p->child1_mass_plus_child2_mass = P_p->child1_mass + P_p->child2_mass;
                P_p->child1_mass_minus_child2_mass = P_p->child1_mass - P_p->child2_mass;
                P_p->child1_mass_times_child2_mass = P_p->child1_mass*P_p->child2_mass;
                
//                printf("level %d m %g hl %d\n",level,P_p->mass,highest_level);
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

void determine_internal_mass_and_semimajor_axis(ParticlesMap *particlesMap)
{
    int N_bodies,N_binaries,N_root_finding;
    determine_binary_parents_and_levels(particlesMap,&N_bodies,&N_binaries,&N_root_finding);
    set_binary_masses_from_body_masses(particlesMap);
    
    double h_tot_vec[3];
    double semimajor_axis,eccentricity,inclination,argument_of_pericenter,longitude_of_ascending_node;
    
    ParticlesMapIterator it_p;
    for (it_p = particlesMap->begin(); it_p != particlesMap->end(); it_p++)
    {
        Particle *p = (*it_p).second;
        if (p->is_binary == 1 and p->level == 0)
        {
            flybys_internal_mass = p->mass;
            
            /* Below: a bit overkill to compute the semimajor axis, but uses consistent functions for orbital elements */
            compute_h_tot_vector(particlesMap,h_tot_vec);
            compute_orbital_elements_from_orbital_vectors(p->child1_mass, p->child2_mass, h_tot_vec, \
                p->e_vec_x,p->e_vec_y,p->e_vec_z,p->h_vec_x,p->h_vec_y,p->h_vec_z,
                &semimajor_axis, &eccentricity, &inclination, &argument_of_pericenter, &longitude_of_ascending_node); 
            flybys_internal_semimajor_axis = semimajor_axis;
            //printf("structure.cpp -- determine_internal_mass_and_semimajor_axis -- M_int %g a_int %g\n",flybys_internal_mass,flybys_internal_semimajor_axis);
        }
    }
}

}
