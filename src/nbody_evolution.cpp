/* MSE */

#include "nbody_evolution.h"
#include "evolve.h"
#include "mstar/regularization.h"

extern "C"
{

void integrate_nbody_system(ParticlesMap *particlesMap, int *integration_flag, double t_old, double t, double *t_out, double *dt_nbody)
{
    /* Integration flags:
     * 0: secular
     * 1: N-body (initiated after dynamical instability; could last for a while)
     * 2: N-body (initiated after semisecular; will last indefinitely)
     * 3: N-body after unbound orbit(s) due to SNe
     * 4: N-body after unbound orbit(s) due to flyby */
     
    /* Here, it is assumed that R_vec and V_vec are already correctly set when integrate_nbody_system() is called */

    struct RegularizedRegion *R = create_mstar_instance_of_system(particlesMap);

    //double tend = determine_nbody_timestep(particlesMap,*integration_flag);
    double dt = t - t_old;
    double dt_reached;
    int collision_occurred;
    
    printf("nbody_evolution.cpp -- integrate_nbody_system -- t %g dt %g\n",t,dt);
    print_system(particlesMap);

    printf("pre... dt %g\n",dt);
    print_state(R);
    run_integrator(R, dt, &dt_reached, &collision_occurred);
    printf("post\n");
    print_state(R);

    *t_out = t_old + dt_reached;

    if (collision_occurred == 1)
    {
        printf("COL\n");
        
        handle_collisions_nbody(R, particlesMap, integration_flag);
        *integration_flag = 1; // continue with direct N-body after the collision, at least initially

        free_data(R);
        return;
    }

    if (*integration_flag == 2) // Continue with direct N-body after semisecular regime; need to make sure particlesMap is updated with pos/vel of bodies
    {
        extract_pos_vel_from_mstar_system(R,particlesMap);
        
        /* Update orbits */
        update_masses_positions_and_velocities_of_all_binaries(particlesMap);
        update_orbital_vectors_in_binaries_from_positions_and_velocities(particlesMap);        
        
        printf("nbody_evolution.cpp -- integrate_nbody_system -- new integration flag %d\n",*integration_flag);
        free_data(R);
        return;
    }

    /* The system is potentially stable.
     * Analyse the system for stability. */

    bool stable_system;
    ParticlesMap new_particlesMap;

    double P_orb_min,P_orb_max;
    analyze_mstar_system(R,&stable_system,&new_particlesMap,&P_orb_min,&P_orb_max,dt); // Will create new particlesMap with new orbits 

    //printf("done an\n");
    //printf("test %g\n",(new_particlesMap)[0]->metallicity);
    //*particlesMap = *new_particlesMap;
    copy_bodies_from_old_to_new_particlesMap(particlesMap,&new_particlesMap); // Copy properties of all bodies from old to new particlesMap
    //printf("ok test 2 %d\n",new_particlesMap[0]->include_mass_transfer_terms);
   
    copy_particlesMap(&new_particlesMap,particlesMap); /* overwrite everything in the particlesMap that was passed onto integrate_nbody_system so the system is updated */

    /* When doing secular integration, some stellar evolution quantities are updated as parts of the ODE solution. In the N-body case, these quantities need to be updated manually */
    update_stellar_evolution_quantities_directly(particlesMap,dt); /* the dt here should be consistent with the dt used to compute the time derivatives in stellar_evolution.cpp */

    if (stable_system == true)
    {
        *integration_flag = 0; // Switch back to secular
    }
    else
    {
        *integration_flag = 1; // Continue running direct N-body in "default" mode
    }

    //int N_bodies, N_binaries,N_root_finding,N_ODE_equations;
    //determine_binary_parents_and_levels(particlesMap,&N_bodies,&N_binaries,&N_root_finding,&N_ODE_equations);
    //set_binary_masses_from_body_masses(particlesMap);

    update_structure(particlesMap);


    *dt_nbody = determine_nbody_timestep(particlesMap,*integration_flag,P_orb_min,P_orb_max);
        
    printf("nbody_evolution.cpp -- integrate_nbody_system -- stable_system %d new integration flag %d dt_nbody %g\n",stable_system,*integration_flag,*dt_nbody);

    //printf("nbody_evolution.cpp -- test %g %g %g\n",(*particlesMap)[0]->R_vec[0],(*particlesMap)[0]->R_vec[1],(*particlesMap)[0]->R_vec[2]);

    print_system(particlesMap);
    free_data(R);
    
    
        
    return; // &new_particlesMap;
}

void handle_collisions_nbody(struct RegularizedRegion *R, ParticlesMap *particlesMap, int *integration_flag)
{

    int i,j,k;
    bool is_binary;
    int index;
    int col_part_i = -1;
    int col_part_j = -1;
    for (i=0; i<R->NumVertex; i++)
    {
        //printf("i %d mass %g \n",i,R->Mass[i]);
        //printf("i %d pos %g %g %g \n",i,R->Pos[3 * i + 0], R->Pos[3 * i + 1],R->Pos[3 * i + 2]);
        //printf("i %d vel %g %g %g \n",i,R->Vel[3 * i + 0], R->Vel[3 * i + 1],R->Vel[3 * i + 2]);
        //printf("R->Collision_Partner[i]; %d\n",R->Collision_Partner[i]);
        
        col_part_i = R->Collision_Partner[i];
        if (col_part_i != -1)
        {
            for (j=0; j<R->NumVertex; j++)
            {
                col_part_j = R->Collision_Partner[j];
                {
                    if (col_part_j != -1 and col_part_j != col_part_i)
                    {
                        break;
                    }
                }
            }
            break;
        }
    }
    if (col_part_i == -1 or col_part_j == -1)
    {
        printf("nbody_evolution.cpp -- error in handle_collisions_nbody: unable to find pair of colliding bodies\n");
        exit(-1);
    }

    (*particlesMap)[col_part_i]->Collision_Partner = col_part_i;
    (*particlesMap)[col_part_j]->Collision_Partner = col_part_j;

    /* Positions & velocities are needed for collision handling */
    for (k=0; k<3; k++)
    {
        (*particlesMap)[col_part_i]->R_vec[k] = R->Pos[3 * col_part_i + k];
        (*particlesMap)[col_part_j]->R_vec[k] = R->Pos[3 * col_part_j + k];

        (*particlesMap)[col_part_i]->V_vec[k] = R->Vel[3 * col_part_i + k];
        (*particlesMap)[col_part_j]->V_vec[k] = R->Vel[3 * col_part_j + k];
        
    }
    //collision_product(particlesMap, binary_index, col_part_i, col_part_j, integration_flag);
    handle_collisions(particlesMap,integration_flag);
    
}

void update_stellar_evolution_quantities_directly(ParticlesMap *particlesMap, double dt)
{
    int i;
    double spin_vec_norm;
    
    ParticlesMapIterator it_p;
    for (it_p = particlesMap->begin(); it_p != particlesMap->end(); it_p++)
    {
        Particle *p = (*it_p).second;
        if (p->is_binary == false and p->evolve_as_star == true)
        {
            p->mass += (p->mass_dot_wind + p->mass_dot_wind_accretion) * dt;
            p->radius += p->radius_dot * dt;
            //printf("dm %.15f %.15f %.15f\n",p->mass_dot_wind,dt,p->mass_dot_wind * dt);
            spin_vec_norm = norm3(p->spin_vec);
            if (spin_vec_norm == 0.0)
            {
                spin_vec_norm = epsilon;
            }
            //printf("spin_vec_norm %g\n",spin_vec_norm);
            for (i=0; i<3; i++)
            {
                p->spin_vec[i] += p->ospin_dot * (p->spin_vec[i]/spin_vec_norm) * dt;
            }
        }
    }
}

double determine_nbody_timestep(ParticlesMap *particlesMap, int integration_flag, double P_orb_min, double P_orb_max)
{
    //double P_orb_max = determine_longest_orbital_period_in_system(particlesMap);
    
    // TO DO: make adjustable
    double dynamical_instability_direct_integration_time_multiplier = 1.5;
    //double dynamical_instability_direct_integration_time_multiplier = 0.01;
    double semisecular_direct_integration_time_multiplier = 1.0e2;
    double supernovae_direct_integration_time_multiplier = 1.5;    
    
    double f;
    if (integration_flag == 1) // dyn. inst.
    {
        f = dynamical_instability_direct_integration_time_multiplier;
    }
    else if (integration_flag == 2) // semisecular
    {
        f = semisecular_direct_integration_time_multiplier;
    }
    else if (integration_flag == 3) // SNe
    {
        f = supernovae_direct_integration_time_multiplier;
    }

    double dt = P_orb_max * f;
    
    return dt;
}

int handle_dynamical_instability(ParticlesMap *particlesMap)
{
    //create_mstar_instance_of_system(particlesMap,R);
    struct RegularizedRegion *R = create_mstar_instance_of_system(particlesMap);
    
//    double P_orb_max = determine_longest_orbital_period_in_system(particlesMap);


//    double tend=P_orb_max;
    
    

    return 0;
}

void copy_bodies_from_old_to_new_particlesMap(ParticlesMap *old_particlesMap, ParticlesMap *new_particlesMap)
{
    int index;
    //double *R_vec,V_vec;
    int i;
    
    ParticlesMapIterator it_p;
    for (it_p = new_particlesMap->begin(); it_p != new_particlesMap->end(); it_p++)
    {
        Particle *p = (*it_p).second;
        
        if (p->is_binary == false)
        {

            index = p->index;

            Particle *p_old = (*old_particlesMap)[index];
            //R_vec = p->R_vec;
            //V_vec = p->V_vec;
            for (i=0; i<3; i++)
            {
                p_old->R_vec[i] = p->R_vec[i];
                p_old->V_vec[i] = p->V_vec[i];
            }

            (*new_particlesMap)[index] = p_old;
//            printf("test copy %g\n",(*new_particlesMap)[index]->metallicity);
        }
    }
}


void analyze_mstar_system(struct RegularizedRegion *R, bool *stable_system, ParticlesMap *particlesMap, double *P_orb_min, double *P_orb_max, double dt)
{
    //printf("ok..... %d\n",particlesMap->size());
//    int highest_new_particle_index = 0;
    //double P_orb_min,P_orb_max;
    std::vector<double> semimajor_axes, future_semimajor_axes;
    std::vector<int> binary_indices, future_binary_indices;
    
    //printf("0\n");
    extract_pos_vel_from_mstar_system(R, particlesMap);
    //printf("1\n");
    find_binaries_in_system(particlesMap,P_orb_min,P_orb_max);
    //printf("2\n");
    get_all_semimajor_axes_in_system(particlesMap, &semimajor_axes, &binary_indices);    
    //printf("3\n");
//    printf("P_orb_max %g\n",P_orb_max);

    /* Integrate system directly for the longest orbital system and compare the semimajor axes,
     * to check if the system is "dynamically stable". */
     
    //double dt_an = *P_orb_min*M_PI*10.0;
    //double dt_an = *P_orb_max * 0.01;
    double dt_an = dt * 0.1;
    
    double maximum_analysis_time = 1.0e5; // TO DO: make user-adjustable
    
    double dt_reached;
    int collision_occurred;
    
    dt_an = min(dt_an, maximum_analysis_time);
    //double dt_an = 0.001 * dt;
    printf("analyze_mstar_system dt %g\n",dt_an);
    run_integrator(R, dt_an, &dt_reached, &collision_occurred);
    printf("done int\n");
    //new_particlesMap.clear();
    //semimajor_axes.clear();

    ParticlesMap future_particlesMap;    
    double future_P_orb_min,future_P_orb_max;
    extract_pos_vel_from_mstar_system(R, &future_particlesMap);
    find_binaries_in_system(&future_particlesMap,&future_P_orb_min,&future_P_orb_max);
//    printf("P_orb_max2 %g\n",P_orb_max);
    get_all_semimajor_axes_in_system(&future_particlesMap, &future_semimajor_axes, &future_binary_indices);    

    std::vector<double>::iterator it;
    double a,a_fut,delta_a;
    int index;
    
    *stable_system = true;
    if (semimajor_axes.size() != future_semimajor_axes.size()) // number of bound orbits should not change!
    {
        *stable_system = false;
    }
    else
    {
        for (it = semimajor_axes.begin(); it != semimajor_axes.end(); it++)
        {
            index = std::distance(semimajor_axes.begin(), it);
            a = semimajor_axes[index];
            a_fut = future_semimajor_axes[index];

            
            delta_a = fabs((a - a_fut))/a;
            printf("AN old %g new %g delta %g\n",a,a_fut,delta_a);

            Particle *p = (*particlesMap)[binary_indices[index]];
            
            if (delta_a > 0.01) // TO DO: make adjustable parameter
            {
                //printf("unstable!\n");
                *stable_system = false;
                p->stable = false;
                
            }
        }
    }
    
    if (collision_occurred == 1)
    {
        *stable_system = false;
    }
    
    printf("nbody_evolution.cpp -- analyze_mstar_system -- stable_system = %d\n",*stable_system);
    //printf("test %g\n",(new_particlesMap)[0]->metallicity);
    
    //return &particlesMap;
}

void extract_pos_vel_from_mstar_system(struct RegularizedRegion *R, ParticlesMap *particlesMap)
{
    int i,j;
    bool is_binary;
    int index;
    
    for (i=0; i<R->NumVertex; i++)
    {
        //printf("i %d mass %g \n",i,R->Mass[i]);
        //printf("i %d pos %g %g %g \n",i,R->Pos[3 * i + 0], R->Pos[3 * i + 1],R->Pos[3 * i + 2]);
        //printf("i %d vel %g %g %g \n",i,R->Vel[3 * i + 0], R->Vel[3 * i + 1],R->Vel[3 * i + 2]);

        is_binary = false;
        //index = particlesMap->size();
        index = R->Index[i];
        Particle *p = new Particle(index, is_binary);
        (*particlesMap)[index] = p;
        
        p->mass = R->Mass[i];
        for (j=0; j<3; j++)
        {
            p->R_vec[j] = R->Pos[3 * i + j];
            p->V_vec[j] = R->Vel[3 * i + j];
            //printf("extract %g \n",p->R_vec[j]);
        }
        
    }
}    

void get_all_semimajor_axes_in_system(ParticlesMap *particlesMap, std::vector<double> *semimajor_axes, std::vector<int> *binary_indices)
{
    ParticlesMapIterator it;
    for (it=particlesMap->begin(); it != particlesMap->end(); it++)
    {
        Particle *p = (*it).second;
        if (p->is_binary == true)
        {
            (*binary_indices).push_back(p->index);
            (*semimajor_axes).push_back(p->a);
            //printf("i %d a %g\n",p->index,p->a);
        }
    }
}
    
void find_binaries_in_system(ParticlesMap *particlesMap, double *P_orb_min, double *P_orb_max)
{
    int i,j;
    bool is_binary;
    double *R1_vec,*V1_vec;
    double *R2_vec,*V2_vec;
    double m1,m2,M;
    
    double r_vec[3],v_vec[3];
    double e_vec[3],h_vec[3];
    double a,e,true_anomaly;
   
    ParticlesMapIterator it,it_p1,it_p2,it_bin;
    //std::map<int, int> binary_child1_indices;
    //std::map<int, int> binary_child2_indices;
    //int binary_child1_indices_i;
    //int binary_child2_indices_i;
    //std::vector<int> binary_child1_indices;
    //std::vector<int> binary_child2_indices;
    
//    ParticlesMapIterator it_p1_begin,it_p1_end;
    //ParticlesMapIterator it_p2_begin,it_p2_end;
    
    ParticlesMap particles_to_be_added;
    //int particles_to_be_added_highest_index;
    int index;

    //double max_a;
    //int N_bodies, N_binaries;
    //int N_root_finding;
    //int N_ODE_equations;
    
    int particles_to_be_added_lower_index = 1 + (*particlesMap->end()).first;
    
    bool binary_already_found;
    double P_orb;
    double delta_a;
    *P_orb_min = 1.0e100;
    *P_orb_max = 0.0;    
    //bool stable_system = false;    
    //int highest_new_particle_index = particlesMap->size();
    
    int new_index=-1;
    for (it = particlesMap->begin(); it != particlesMap->end(); it++)
        //for (it_p1 = it_p1_begin; it_p1 != it_p1_end; it_p1++)
    {
        Particle *p = (*it).second;
        if (p->index > new_index)
        {
            new_index = p->index;
        }
        p->has_found_parent = false;
    }
    new_index += 1;
    //int new_index = (*particlesMap->end()).first;
    //printf("new_index %d\n",new_index);
    
    bool found_new_orbit = true;
    while (found_new_orbit == true)
    {
        //found_new_orbit = false;
        
        //it_p1_begin = new_particlesMap.begin();
        //it_p1_end = new_particlesMap.end();
        //it_p2_begin = new_particlesMap.begin();
        //it_p2_end = new_particlesMap.end();

        particles_to_be_added.clear();
        //particles_to_be_added_highest_index = 0;
        for (it_p1 = particlesMap->begin(); it_p1 != particlesMap->end(); it_p1++)
        //for (it_p1 = it_p1_begin; it_p1 != it_p1_end; it_p1++)
        {
            Particle *p1 = (*it_p1).second;
            
            if (p1->has_found_parent == true)
            {
                continue;
            }
            
            m1 = p1->mass;
            R1_vec = p1->R_vec;
            V1_vec = p1->V_vec;
            
            for (it_p2 = particlesMap->begin(); it_p2 != particlesMap->end(); it_p2++)
            //for (it_p2 = it_p2_begin; it_p2 != it_p2_end; it_p2++)
            {
                Particle *p2 = (*it_p2).second;
                
                if (p2->index <= p1->index) // only consider unique pairs
                {
                    continue;
                }

                if (p2->has_found_parent == true)
                {
                    continue;
                }
                printf("is %d %d\n",p1->index,p2->index);
                m2 = p2->mass;
                M = m1 + m2;
                R2_vec = p2->R_vec;
                V2_vec = p2->V_vec;

                for (j=0; j<3; j++)
                {
                    r_vec[j] = R1_vec[j] - R2_vec[j];
                    v_vec[j] = V1_vec[j] - V2_vec[j];
                }

                from_cartesian_to_orbital_vectors(m1,m2,r_vec,v_vec,e_vec,h_vec,&true_anomaly);
                compute_semimajor_axis_and_eccentricity_from_orbital_vectors(m1,m2,e_vec,h_vec,&a,&e);
                //printf("i1 %d i2 %d r %g %g %g a %g e %g\n",p1->index,p2->index,r_vec[0],r_vec[1],r_vec[2],a,e);
                
                
                if ( a>= 0.0 and e >= 0.0 and e < 1.0 and p1->has_found_parent == false and p2->has_found_parent == false)
                {
                    //if ( (binary_child1_indices.count(p1->index) > 0) and (binary_child2_indices.count(p2->index) > 0) )
                    //if already in new particles
                    

                        //binary_child1_indices[binary_child1_indices_i] = p1->index;
                        //binary_child2_indices[binary_child2_indices_i] = p2->index;
                        //add binary
                        //p1->has_found_parent = true;
                        //p2->has_found_parent = true;

                    p1->has_found_parent = true;
                    p2->has_found_parent = true;
                    
                    is_binary = true;
                    //index = particles_to_be_added.size();
                    //index = particles_to_be_added_lower_index + particles_to_be_added.size();
                    index = new_index;
                    Particle *b = new Particle(index, is_binary);
                    //new_particlesMap[highest_new_particle_index] = b;
                    //++;
                    particles_to_be_added[index] = b;
                    //particles_to_be_added_highest_index++;
                    
                    for (j=0; j<3; j++)
                    {
                        b->e_vec[j] = e_vec[j];
                        b->h_vec[j] = h_vec[j];
                        b->R_vec[j] = (m1 * R1_vec[j] + m2 * R2_vec[j])/M;
                        b->V_vec[j] = (m1 * V1_vec[j] + m2 * V2_vec[j])/M;
                    }
                    b->true_anomaly = true_anomaly;
                    b->mass = M;
                    b->a = a;
                    b->e = e;
                    b->child1 = p1->index;
                    b->child2 = p2->index;
                    
                    //found_new_orbit = true;
                    printf("New binary %d C1 %d C2 %d a %g e %g TA %g M %g\n",b->index,p1->index,p2->index,a,e,true_anomaly,M);
                    
                    b->child1_mass_plus_child2_mass = M;
                    P_orb = compute_orbital_period(b);
                    if (P_orb > *P_orb_max)
                    {
                        *P_orb_max = P_orb;
                    }
                    if (P_orb < *P_orb_min)
                    {
                        *P_orb_min = P_orb;
                    }
                    
                    new_index++;
                }
            }
        }
        
        if (particles_to_be_added.size() == 0)
        {
            found_new_orbit = false;
        }
        else
        {
            //printf("pre add\n");
            for (it = particles_to_be_added.begin(); it != particles_to_be_added.end(); it++)
            {
                Particle *p = (*it).second;
                //p->index = particlesMap->size();
                //p->index = (*particlesMap->end()).first + 1;
                printf("p->index %d\n",p->index);
                
                (*particlesMap)[p->index] = p;
                //highest_new_particle_index++;
                //printf("adding particle %d a %g\n",p->index,p->a);
            }
            //printf("done add\n");
        }
    }
//    printf("test %d\n",new_particlesMap.size());
//    printf("test %g\n",new_particlesMap[0]->R_vec[0]);

//    printf("%d\n",R->NumVertex);
    //printf("find binaries done!\n");
    printf("done an\n");
    print_system(particlesMap);
}


void analyze_mstar_system2(struct RegularizedRegion *R)
{
    int i,j;
    
    /* Copy data from MSTAR to a new particlesMap */
    int highest_new_particle_index = 0;
    ParticlesMap new_particlesMap;
    bool is_binary;

    //int N_bodies = R->NumVertex;
    for (i=0; i<R->NumVertex; i++)
    {
        printf("i %d mass %g \n",i,R->Mass[i]);
        printf("i %d pos %g %g %g \n",i,R->Pos[3 * i + 0], R->Pos[3 * i + 1],R->Pos[3 * i + 2]);
        printf("i %d vel %g %g %g \n",i,R->Vel[3 * i + 0], R->Vel[3 * i + 1],R->Vel[3 * i + 2]);

        is_binary = false;
        Particle *p = new Particle(highest_new_particle_index, is_binary);
        new_particlesMap[highest_new_particle_index] = p;
        highest_new_particle_index++;
        
        p->mass = R->Mass[i];
        for (j=0; j<3; j++)
        {
            p->R_vec[j] = R->Pos[3 * i + j];
            p->V_vec[j] = R->Vel[3 * i + j];
        }
        
    }

    double *R1_vec,*V1_vec;
    double *R2_vec,*V2_vec;
    double m1,m2,M;
    
    double r_vec[3],v_vec[3];
    double e_vec[3],h_vec[3];
    double a,e,true_anomaly;
   
    ParticlesMapIterator it,it_p1,it_p2,it_bin;
    //std::map<int, int> binary_child1_indices;
    //std::map<int, int> binary_child2_indices;
    //int binary_child1_indices_i;
    //int binary_child2_indices_i;
    //std::vector<int> binary_child1_indices;
    //std::vector<int> binary_child2_indices;
    
    ParticlesMapIterator it_p1_begin,it_p1_end;
    ParticlesMapIterator it_p2_begin,it_p2_end;
    
    ParticlesMap particles_to_be_added;
    int particles_to_be_added_highest_index;

    //double max_a;
    //int N_bodies, N_binaries;
    //int N_root_finding;
    //int N_ODE_equations;
    
    bool binary_already_found;
    double P_orb,P_orb_max;
    double delta_a;
    
    bool found_new_orbit = true;
    while (found_new_orbit == true)
    {
        //found_new_orbit = false;
        
        it_p1_begin = new_particlesMap.begin();
        it_p1_end = new_particlesMap.end();
        it_p2_begin = new_particlesMap.begin();
        it_p2_end = new_particlesMap.end();

        P_orb_max = 0.0;
        particles_to_be_added.clear();
        particles_to_be_added_highest_index = 0;
        for (it_p1 = new_particlesMap.begin(); it_p1 != new_particlesMap.end(); it_p1++)
        //for (it_p1 = it_p1_begin; it_p1 != it_p1_end; it_p1++)
        {
            Particle *p1 = (*it_p1).second;
            
            if (p1->has_found_parent == true)
            {
                continue;
            }
            
            m1 = p1->mass;
            R1_vec = p1->R_vec;
            V1_vec = p1->V_vec;
            
            for (it_p2 = new_particlesMap.begin(); it_p2 != new_particlesMap.end(); it_p2++)
            //for (it_p2 = it_p2_begin; it_p2 != it_p2_end; it_p2++)
            {
                Particle *p2 = (*it_p2).second;
                
                if (p2->index <= p1->index) // only consider unique pairs
                {
                    continue;
                }

                if (p2->has_found_parent == true)
                {
                    continue;
                }

                m2 = p2->mass;
                M = m1 + m2;
                R2_vec = p2->R_vec;
                V2_vec = p2->V_vec;

                for (j=0; j<3; j++)
                {
                    r_vec[j] = R1_vec[j] - R2_vec[j];
                    v_vec[j] = V1_vec[j] - V2_vec[j];
                }

                from_cartesian_to_orbital_vectors(m1,m2,r_vec,v_vec,e_vec,h_vec,&true_anomaly);
                compute_semimajor_axis_and_eccentricity_from_orbital_vectors(m1,m2,e_vec,h_vec,&a,&e);
                printf("i1 %d i2 %d r %g %g %g a %g e %g\n",p1->index,p2->index,r_vec[0],r_vec[1],r_vec[2],a,e);
                
                
                if ( a>= 0.0 and e >= 0.0 and e < 1.0)
                {
                    //if ( (binary_child1_indices.count(p1->index) > 0) and (binary_child2_indices.count(p2->index) > 0) )
                    //if already in new particles
                    binary_already_found = false;
                    #ifdef IGNORE
                    for (it_bin = new_particlesMap.begin(); it_bin != new_particlesMap.end(); it_bin++)
                    {
                        Particle *b = (*it_bin).second;
                    
                        if (b->is_binary == true)
                        {
                            if ((b->child1 == p1->index) and (b->child2 == p2->index))
                            {
                                //found_new_orbit = false;
                                binary_already_found = true;
                                printf("Previously found b %d C1 %d C2 %d b->a %g a %g\n",b->index,p1->index,p2->index,b->a,a);
                                delta_a = fabs((b->a - a)/b->a);
                                if (delta_a > 0.2)
                                {
                                    new_particlesMap.erase(b->index);
                                    highest_new_particle_index --;
                                }
                                else
                                {
                                    p1->has_found_parent = true;
                                    p2->has_found_parent = true;
                                }
                            }
                            //check if sma did not change too much
                        }
                    }
                    #endif
                    
                    if (binary_already_found == false)
                    {

                        //binary_child1_indices[binary_child1_indices_i] = p1->index;
                        //binary_child2_indices[binary_child2_indices_i] = p2->index;
                        //add binary
                        p1->has_found_parent = true;
                        p2->has_found_parent = true;
                    
                        is_binary = true;
                        Particle *b = new Particle(particles_to_be_added_highest_index, is_binary);
                        //new_particlesMap[highest_new_particle_index] = b;
                        //++;
                        particles_to_be_added[particles_to_be_added_highest_index] = b;
                        particles_to_be_added_highest_index++;
                        
                        for (j=0; j<3; j++)
                        {
                            b->e_vec[j] = e_vec[j];
                            b->h_vec[j] = h_vec[j];
                            b->R_vec[j] = (m1 * R1_vec[j] + m2 * R2_vec[j])/M;
                            b->V_vec[j] = (m1 * V1_vec[j] + m2 * V2_vec[j])/M;
                        }
                        b->mass = M;
                        b->a = a;
                        b->e = e;
                        b->child1 = p1->index;
                        b->child2 = p2->index;
                        
                        //found_new_orbit = true;
                        printf("New binary %d C1 %d C2 %d a %g e %g\n",b->index,p1->index,p2->index,a,e);
                        
                        b->child1_mass_plus_child2_mass = M;
                        P_orb = compute_orbital_period(b);
                        if (P_orb > P_orb_max)
                        {
                            P_orb_max = P_orb;
                        }
                    }
                }
            }
        }
        
        if (particles_to_be_added.size() == 0)
        {
            found_new_orbit = false;
        }
        else
        {
            printf("pre add\n");
            for (it = particles_to_be_added.begin(); it != particles_to_be_added.end(); it++)
            {
                Particle *p = (*it).second;
                p->index = highest_new_particle_index;
                new_particlesMap[highest_new_particle_index] = p;
                highest_new_particle_index++;
                printf("adding particle %d a %g\n",p->index,p->a);
            }
            printf("done add\n");
            
            //determine_binary_parents_and_levels(&new_particlesMap,&N_bodies,&N_binaries,&N_root_finding,&N_ODE_equations);
            //P_orb_max = determine_longest_orbital_period_in_system(&new_particlesMap);
            printf("P_orb_max %g\n",P_orb_max);
            
            double dt_reached;
            int collision_occurred;
            run_integrator(R, 0.5*P_orb_max, &dt_reached, &collision_occurred);
            printf("done int\n");
            
            for (i=0; i<R->NumVertex; i++)
            {
                printf("i %d mass %g \n",i,R->Mass[i]);
                printf("i %d pos %g %g %g \n",i,R->Pos[3 * i + 0], R->Pos[3 * i + 1],R->Pos[3 * i + 2]);
                printf("i %d vel %g %g %g \n",i,R->Vel[3 * i + 0], R->Vel[3 * i + 1],R->Vel[3 * i + 2]);

                Particle *p = new_particlesMap[i];
                p->mass = R->Mass[i];
                for (j=0; j<3; j++)
                {
                    p->R_vec[j] = R->Pos[3 * i + j];
                    p->V_vec[j] = R->Vel[3 * i + j];
                }
            }
            
            for (it = new_particlesMap.begin(); it != new_particlesMap.end(); it++)
            {
                Particle *b = (*it).second;
                if (b->is_binary == true)
                {
                    Particle *child1 = new_particlesMap[b->child1];
                    Particle *child2 = new_particlesMap[b->child2];
                    
                    m1 = child1->mass;
                    m2 = child2->mass;
                    M = m1 + m2;
                    for (j=0; j<3; j++)
                    {
                        b->R_vec[j] = (m1 * child1->R_vec[j] + m2 * child2->R_vec[j])/M;
                        b->V_vec[j] = (m1 * child1->V_vec[j] + m2 * child2->V_vec[j])/M;
                    }
                }
            }
            
        }
    }
//    printf("test %d\n",new_particlesMap.size());
//    printf("test %g\n",new_particlesMap[0]->R_vec[0]);

//    printf("%d\n",R->NumVertex);
    printf("analyse done!\n");

}


struct RegularizedRegion *create_mstar_instance_of_system(ParticlesMap *particlesMap)
{
     
    initialize_mpi_or_serial(); // This needs to be done to initialize MSTAR
    //printf("?%d \n",NumGbsGroup);
    struct RegularizedRegion *R;

    int i=0;
    int j;
    int N_bodies = determine_number_of_bodies_in_system(particlesMap);
        
    //my_barrier(); // In serial mode, this does nothing (so can be omitted)

    //int N_particles = particlesMap->size();
    //int N_bodies, N_binaries;
    //int N_root_finding;
    //int N_ODE_equations;
    
    //determine_binary_parents_and_levels(particlesMap,&N_bodies,&N_binaries,&N_root_finding,&N_ODE_equations);
    
//    int MAXPART = N_bodies; // (Maximum) number of bodies; should be consistent with number of bodies in input file
    //my_barrier();

    allocate_armst_structs(&R, N_bodies); // Initialize the data structure
    
    double *R_vec, *V_vec;
    ParticlesMapIterator it_p;
    for (it_p = particlesMap->begin(); it_p != particlesMap->end(); it_p++)
    {
        Particle *p = (*it_p).second;
        
        if (p->is_binary == false)
        {
            R_vec = p->R_vec;
            V_vec = p->V_vec;

            R->Vertex[i].type = 0;
            R->Index[i] = p->index;
            R->Mass[i] = p->mass;
            R->Radius[i] = p->radius;
            for (j=0; j<3; j++)
            {
                R->Pos[3 * i + j] = R_vec[j];
                R->Vel[3 * i + j] = V_vec[j];
            }
            
            i++;
        }    
    }   
    
    R->gbs_tolerance = mstar_gbs_tolerance;
    R->collision_tolerance = mstar_collision_tolerance;
    return R;
    
}

void print_state(struct RegularizedRegion *R)
{
    for (int i=0; i<R->NumVertex; i++)
    {
        printf("i %d index %d mass %g radius %g collision partner %d\n",i,R->Index[i],R->Mass[i],R->Radius[i],R->Collision_Partner[i]);
        printf("i %d pos %g %g %g \n",i,R->Pos[3 * i + 0], R->Pos[3 * i + 1],R->Pos[3 * i + 2]);
        printf("i %d vel %g %g %g \n",i,R->Vel[3 * i + 0], R->Vel[3 * i + 1],R->Vel[3 * i + 2]);
    }
    printf("=========\n");
}

}
