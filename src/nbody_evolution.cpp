/* MSE */

#include "nbody_evolution.h"
#include "evolve.h"
//#include "mstar/regularization.h"

extern "C"
{

void integrate_nbody_system(ParticlesMap *particlesMap, int *integration_flag, double t_old, double t, double *t_out, double *dt_nbody)
{
    
    /* Integration flags:
     * 0: secular
     * 1: N-body (initiated after dynamical instability; could last for a while)
     * 2: N-body (initiated after semisecular)
     * 3: N-body after unbound orbit(s) due to SNe
     * 4: N-body after unbound orbit(s) due to flyby */
     
    /* Here, it is assumed that R_vec and V_vec are already correctly set when integrate_nbody_system() is called */

    double dt = t - t_old;
    double dt_reached;
    
    int N_bodies_eff = check_particlesMap_for_inclusion_in_MSTAR(particlesMap);
    
    if (N_bodies_eff < 1)
    {

        #ifdef VERBOSE
        if (verbose_flag > 1)
        {
            printf("nbody_evolution.cpp -- no bodies to integrate N_bodies_eff %d!\n",N_bodies_eff);
            print_system(particlesMap,*integration_flag);
        }
        #endif

        *dt_nbody = 1.0e100;
        *t_out = t_old + dt;
        *integration_flag = 1;
        update_stellar_evolution_quantities_directly(particlesMap,t, dt);
        remove_binaries_from_system(particlesMap);
        return;
    }
    
    struct RegularizedRegion *R = create_mstar_instance_of_system(particlesMap,*integration_flag);

    int collision_occurred;

    double E_init,E_fin;
    MSTAR_verbose = false;
    #ifdef VERBOSE
    if (verbose_flag > 0)
    {
        E_init = compute_nbody_total_energy(R);
        printf("nbody_evolution.cpp -- integrate_nbody_system -- starting integration -- t %g dt %g E_init %g\n",t,dt,E_init);
        print_system(particlesMap,*integration_flag);
        mstar_print_state(R);
        MSTAR_verbose = true;
    }
    #endif

    run_integrator(R, dt, &dt_reached, &collision_occurred);


    #ifdef VERBOSE
    if (verbose_flag > 0)
    {
        E_fin = compute_nbody_total_energy(R);
        printf("nbody_evolution.cpp -- integrate_nbody_system -- finished integration -- t %g dt_reached %g E_fin %g E_err %g\n",t,dt_reached,E_fin,fabs((E_init-E_fin)/E_init));
        mstar_print_state(R);
    }
    #endif

    *t_out = t_old + dt_reached;

    if (collision_occurred == 1)
    {
        #ifdef VERBOSE
        if (verbose_flag > 0)
        {
            printf("nbody_evolution.cpp -- integrate_nbody_system -- finished integration -- COLLISION\n");
            mstar_print_state(R);
        }
        #endif
       
        update_pos_vel_from_mstar_system(R,particlesMap);
        handle_collisions_nbody(R, particlesMap, t, integration_flag);
        *integration_flag = 1; // continue with direct N-body after the collision, at least initially

        return;
    }


    /* The system is potentially stable.
     * Analyse the system for stability. */

    double P_orb_min,P_orb_max;
    *integration_flag = determine_new_integration_flag_using_nbody(R,particlesMap,&P_orb_min,&P_orb_max);

    *dt_nbody = determine_nbody_timestep(particlesMap,*integration_flag,P_orb_min,P_orb_max);

    update_stellar_evolution_quantities_directly(particlesMap,t,dt);

    #ifdef VERBOSE
    if (verbose_flag > 0)
    {
        printf("nbody_evolution.cpp -- integrate_nbody_system -- done -- new integration flag %d dt_nbody %g P_orb_min %g P_orb_max %g\n",*integration_flag,*dt_nbody,P_orb_min,P_orb_max);
        print_system(particlesMap,*integration_flag);
    }
    #endif

    MSTAR_verbose = false;    
    free_data(R);
      
    return;
}

int determine_orbits_in_system_using_nbody(ParticlesMap *particlesMap)
{
    #ifdef VERBOSE
    if (verbose_flag > 0)
    {
        printf("nbody_evolution.cpp -- determine_orbits_in_system_using_nbody -- begin (setting integration_flag = 0 for print_system) \n");
        print_system(particlesMap,1);
    }
    #endif
    int N_bodies_eff = check_particlesMap_for_inclusion_in_MSTAR(particlesMap);

    if (N_bodies_eff < 1)
    {
        return 1;
    }
    int integration_flag = 1;
    struct RegularizedRegion *R = create_mstar_instance_of_system(particlesMap,integration_flag);

    double P_orb_min,P_orb_max;
    integration_flag = determine_new_integration_flag_using_nbody(R,particlesMap,&P_orb_min,&P_orb_max);
    
    free_data(R);

    #ifdef VERBOSE
    if (verbose_flag > 0)
    {
        printf("nbody_evolution.cpp -- determine_orbits_in_system_using_nbody -- end\n");
        print_system(particlesMap,integration_flag);
    }
    #endif
    
    return integration_flag;
}

int determine_new_integration_flag_using_nbody(struct RegularizedRegion *R, ParticlesMap *particlesMap, double *P_orb_min, double *P_orb_max)
{
    
    bool stable_system;
    ParticlesMap new_particlesMap;

    analyze_mstar_system(R,&stable_system,&new_particlesMap,P_orb_min,P_orb_max,nbody_analysis_minimum_integration_time); // Will create new particlesMap with new orbits 

    copy_bodies_from_old_to_new_particlesMap(particlesMap,&new_particlesMap); // Copy properties of all bodies from old to new particlesMap

    copy_particlesMap(&new_particlesMap,particlesMap); /* overwrite everything in the particlesMap that was passed onto integrate_nbody_system so the system is updated */

    int integration_flag;
    int dummy = 0;
    bool stable_MA01 = check_system_for_dynamical_stability(particlesMap, &dummy);

    if (stable_system == true and stable_MA01 == true)
    {
        integration_flag = 0; // Switch back to secular
    }
    else
    {
        integration_flag = 1; // Continue running direct N-body in "default" mode
    }
    
    return integration_flag;
    
}

void handle_collisions_nbody(struct RegularizedRegion *R, ParticlesMap *particlesMap, double t, int *integration_flag)
{

    int i,j,k;
    bool is_binary;
    int index;
    int col_part_i = -1;
    int col_part_j = -1;
    for (i=0; i<R->NumVertex; i++)
    {
        col_part_i = R->Stopping_Condition_Partner[i];
        if (col_part_i != -1)
        {
            for (j=0; j<R->NumVertex; j++)
            {
                col_part_j = R->Stopping_Condition_Partner[j];
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
    free_data(R);
    
    if (col_part_i == -1 or col_part_j == -1)
    {
        printf("nbody_evolution.cpp -- error in handle_collisions_nbody: unable to find pair of colliding bodies\n");
        error_code = 9;
        //exit(-1);
    }

    #ifdef VERBOSE
    if (verbose_flag > 0)
    {
        printf("nbody_evolution.cpp -- handle_collisions_nbody -- col_part_i %d col_part_j %d\n",col_part_i,col_part_j);
        print_system(particlesMap,1);
    }
    #endif

    (*particlesMap)[col_part_i]->Stopping_Condition_Partner = col_part_i;
    (*particlesMap)[col_part_j]->Stopping_Condition_Partner = col_part_j;

    handle_collisions(particlesMap,t,integration_flag);
    
}

void update_stellar_evolution_quantities_directly(ParticlesMap *particlesMap, double t, double dt)
{
    int i;
    double spin_vec_norm;
    
    ParticlesMapIterator it_p;
    for (it_p = particlesMap->begin(); it_p != particlesMap->end(); it_p++)
    {
        Particle *p = (*it_p).second;
        if (p->is_binary == false and p->object_type == 1)
        {
            p->RLOF_flag = 0;
            
            p->mass += (p->mass_dot_wind + p->mass_dot_wind_accretion) * dt;
            p->radius += p->radius_dot * dt;

            #ifdef VERBOSE
            if (verbose_flag > 1)
            {
                printf("nbody_evolution.cpp -- update_stellar_evolution_quantities_directly m %g r %g\n",p->mass,p->radius);
            }
            #endif

            spin_vec_norm = norm3(p->spin_vec);
            if (spin_vec_norm == 0.0)
            {
                spin_vec_norm = epsilon;
            }
            
            /* Spin changes according to SSE */
            for (i=0; i<3; i++)
            {
                p->spin_vec[i] += p->ospin_dot * (p->spin_vec[i]/spin_vec_norm) * dt;
            }
            
            /* Spin changes due to NS spindown (not modelled in SSE) */
            if (NS_model == 1)
            {
                nbody_handle_NS_properties_Ye19_model(p, spin_vec_norm, dt);
            }
        }
    }
}

double determine_nbody_timestep(ParticlesMap *particlesMap, int integration_flag, double P_orb_min, double P_orb_max)
{
    double f;
    if (integration_flag == 1) // dyn. inst.
    {
        f = nbody_dynamical_instability_direct_integration_time_multiplier;
    }
    else if (integration_flag == 2) // semisecular
    {
        f = nbody_semisecular_direct_integration_time_multiplier;
    }
    else if (integration_flag == 3) // SNe
    {
        f = nbody_supernovae_direct_integration_time_multiplier;
    }
    else
    {
        f = nbody_other_direct_integration_time_multiplier;
    }
    
    double dt = P_orb_max * f;
    
    return dt;
}

void copy_bodies_from_old_to_new_particlesMap(ParticlesMap *old_particlesMap, ParticlesMap *new_particlesMap)
{
    int index;
    int i;
    
    /* Copy particles that were included in the MSTAR integration */
    ParticlesMapIterator it_p;
    for (it_p = new_particlesMap->begin(); it_p != new_particlesMap->end(); it_p++)
    {
        Particle *p = (*it_p).second;
        
        if (p->is_binary == false)
        {
            index = p->index;

            Particle *p_old = (*old_particlesMap)[index];
            double S = norm3(p->spin_AM_vec);
            if (S < epsilon)
            {
                S = epsilon;
            }

            double k2;
            if (p_old->object_type == 1)
            {
                k2 = p_old->sse_k2;
            }
            else
            {
                k2 = p_old->gyration_radius;
            }
            
            double Omega = compute_spin_frequency_from_spin_angular_momentum(S, p_old->stellar_type, p_old->object_type, p_old->mass, p_old->core_mass, p_old->radius, p_old->core_radius, k2, p_old->sse_k3);
 
            for (i=0; i<3; i++)
            {
                p_old->R_vec[i] = p->R_vec[i];
                p_old->V_vec[i] = p->V_vec[i];
                p_old->spin_vec[i] = Omega * (p->spin_AM_vec[i]/S);
            }

            (*new_particlesMap)[index] = p_old;
        }
    }
    
    /* Copy particles that were excluded in the MSTAR integration */
    for (it_p = old_particlesMap->begin(); it_p != old_particlesMap->end(); it_p++)
    {
        Particle *p = (*it_p).second;
        if (p->include_in_MSTAR == false)
        {
            (*new_particlesMap)[p->index] = p;
        }
    }
    
}


void analyze_mstar_system(struct RegularizedRegion *R, bool *stable_system, ParticlesMap *particlesMap, double *P_orb_min, double *P_orb_max, double dt)
{
    /* Integrate system directly for the longest orbital system and compare the semimajor axes,
     * to check if the system is "dynamically stable". */

    std::vector<double> semimajor_axes, future_semimajor_axes;
    std::vector<int> binary_indices, future_binary_indices;
    
    extract_pos_vel_from_mstar_system_and_reset_particlesMap(R, particlesMap);
    find_binaries_in_system(particlesMap,P_orb_min,P_orb_max);
    get_all_semimajor_axes_in_system(particlesMap, &semimajor_axes, &binary_indices);    
     
    double dt_an = dt * nbody_analysis_fractional_integration_time;
    
    double maximum_analysis_time = nbody_analysis_maximum_integration_time;
    
    double dt_reached;
    int collision_occurred;
    
    dt_an = CV_min(dt_an, maximum_analysis_time);
    
    MSTAR_verbose = false;
    #ifdef VERBOSE
    if (verbose_flag > 0)
    {
        printf("nbody_evolution.cpp -- analyze_mstar_system -- starting integration\n");
        MSTAR_verbose = true;
    }
    #endif
   
    run_integrator(R, dt_an, &dt_reached, &collision_occurred);

    ParticlesMap future_particlesMap;    
    double future_P_orb_min,future_P_orb_max;
    extract_pos_vel_from_mstar_system_and_reset_particlesMap(R, &future_particlesMap);
    find_binaries_in_system(&future_particlesMap,&future_P_orb_min,&future_P_orb_max);
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

            #ifdef VERBOSE
            if (verbose_flag > 0)
            {
                printf("nbody_evolution.cpp -- analyze_mstar_system -- a old %g a new %g delta %g\n",a,a_fut,delta_a);
            }
            #endif

            Particle *p = (*particlesMap)[binary_indices[index]];
            
            if (delta_a > nbody_analysis_fractional_semimajor_axis_change_parameter)
            {
                *stable_system = false;
                p->stable = false;
                
            }
        }
    }
    
    if (collision_occurred == 1)
    {
        *stable_system = false;
    }
    
    #ifdef VERBOSE
    if (verbose_flag > 0)
    {
        printf("nbody_evolution.cpp -- analyze_mstar_system -- stable_system = %d\n",*stable_system);
    }
    #endif
}

void extract_pos_vel_from_mstar_system_and_reset_particlesMap(struct RegularizedRegion *R, ParticlesMap *particlesMap)
{
    int i,j;
    bool is_binary;
    int index;
    
    for (i=0; i<R->NumVertex; i++)
    {
        is_binary = false;
        index = R->MSEIndex[i];
        Particle *p = new Particle(index, is_binary);
        (*particlesMap)[index] = p;
        
        p->mass = R->Mass[i];
        for (j=0; j<3; j++)
        {
            p->R_vec[j] = R->Pos[3 * i + j];
            p->V_vec[j] = R->Vel[3 * i + j];
            p->spin_AM_vec[j] = R->Spin_S[3 * i + j];
        }
        
    }
}    

void update_pos_vel_from_mstar_system(struct RegularizedRegion *R, ParticlesMap *particlesMap)
{
    int i,j;
    bool is_binary;
    int index;
    
    for (i=0; i<R->NumVertex; i++)
    {

        is_binary = false;
        index = R->MSEIndex[i];
        Particle *p = (*particlesMap)[index];
        
        for (j=0; j<3; j++)
        {
            p->R_vec[j] = R->Pos[3 * i + j];
            p->V_vec[j] = R->Vel[3 * i + j];
            p->spin_AM_vec[j] = R->Spin_S[3 * i + j];
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

    ParticlesMap particles_to_be_added;
    int index;

    int particles_to_be_added_lower_index = 1 + (*particlesMap->end()).first;
    
    bool binary_already_found;
    double P_orb;
    double delta_a;
    *P_orb_min = 1.0e100;
    *P_orb_max = 0.0;    
    
    int new_index=-1;
    for (it = particlesMap->begin(); it != particlesMap->end(); it++)
    {
        Particle *p = (*it).second;
        if (p->index > new_index)
        {
            new_index = p->index;
        }
        p->has_found_parent = false;
    }
    new_index += 1;
    
    bool found_new_orbit = true;
    while (found_new_orbit == true)
    {
        particles_to_be_added.clear();

        for (it_p1 = particlesMap->begin(); it_p1 != particlesMap->end(); it_p1++)
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
                
                #ifdef VERBOSE
                if (verbose_flag > 1)
                {
                    printf("nbody_evolution.cpp -- find_binaries_in_system -- i1 %d i2 %d r %g %g %g a %g e %g\n",p1->index,p2->index,r_vec[0],r_vec[1],r_vec[2],a,e);
                }
                #endif
                
                if ( a>= 0.0 and e >= 0.0 and e < 1.0 and p1->has_found_parent == false and p2->has_found_parent == false)
                {
                    p1->has_found_parent = true;
                    p2->has_found_parent = true;
                    
                    is_binary = true;
                    index = new_index;
                    Particle *b = new Particle(index, is_binary);
                    particles_to_be_added[index] = b;
                    
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
                    
                    #ifdef VERBOSE
                    if (verbose_flag > 1)
                    {
                        printf("nbody_evolution.cpp -- find_binaries_in_system -- New binary %d C1 %d C2 %d a %g e %g TA %g M %g\n",b->index,p1->index,p2->index,a,e,true_anomaly,M);
                    }
                    #endif

                    b->child1_mass_plus_child2_mass = M;
                    P_orb = compute_orbital_period_from_semimajor_axis(b->mass,b->a);
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
            for (it = particles_to_be_added.begin(); it != particles_to_be_added.end(); it++)
            {
                Particle *p = (*it).second;
                
                (*particlesMap)[p->index] = p;
            }
        }
    }
}


struct RegularizedRegion *create_mstar_instance_of_system(ParticlesMap *particlesMap, int integration_flag)
{
    initialize_mpi_or_serial(); // This needs to be done to initialize MSTAR

    struct RegularizedRegion *R;

    int i=0;
    int j;
    int N_bodies = check_particlesMap_for_inclusion_in_MSTAR(particlesMap);
        
    allocate_armst_structs(&R, N_bodies); // Initialize the data structure
    
    double *R_vec, *V_vec;
    double Omega,S; // Spin frequency, spin angular momentum
    
    ParticlesMapIterator it_p;
    for (it_p = particlesMap->begin(); it_p != particlesMap->end(); it_p++)
    {
        Particle *p = (*it_p).second;
        
        if (p->is_binary == false and p->include_in_MSTAR == true)
        {
            R_vec = p->R_vec;
            V_vec = p->V_vec;

            R->Vertex[i].type = 0;

            R->MSEIndex[i] = p->index;

            R->Mass[i] = p->mass;
            R->Radius[i] = determine_effective_radius_for_collision(p->radius, p->stellar_type, 1);
            
            Omega = norm3(p->spin_vec);
            if (Omega < epsilon)
            {
                Omega = epsilon;
            }

            double k2;
            if (p->object_type == 1)
            {
                k2 = p->sse_k2;
            }
            else
            {
                k2 = p->gyration_radius;
            }
            S = compute_spin_angular_momentum_from_spin_frequency(Omega, p->stellar_type, p->object_type, p->mass, p->core_mass, p->radius, p->core_radius, k2, p->sse_k3);

            for (j=0; j<3; j++)
            {
                R->Pos[3 * i + j] = R_vec[j];
                R->Vel[3 * i + j] = V_vec[j];
                R->Spin_S[3 * i + j] = S * p->spin_vec[j]/Omega;
                p->spin_AM_vec[i] = R->Spin_S[3 * i + j];
            }
            
            R->Stopping_Condition_Mode[i] = 0;

            i++;
        }    
    }   
    
    double gbs_tolerance = MSTAR_gbs_tolerance_default;
    if (integration_flag == 3)
    {
        gbs_tolerance = MSTAR_gbs_tolerance_kick;
    }

    R->gbs_tolerance = gbs_tolerance;
    R->stopping_condition_tolerance = MSTAR_stopping_condition_tolerance;
    R->output_time_tolerance = MSTAR_output_time_tolerance;

    #ifdef VERBOSE
    if (verbose_flag > 0)
    {
        printf("nbody_evolution.cpp -- create_mstar_instance_of_system -- R->gbs_tolerance %g R->stopping_condition_tolerance %g R->output_time_tolerance %g\n",R->gbs_tolerance,R->stopping_condition_tolerance,R->output_time_tolerance);
    }
    #endif
    
    return R;
    
}

double compute_nbody_total_energy(struct RegularizedRegion *R)
{
    int i,j,k;
    double m_i,m_j,R_i[3],V_i[3],R_j[3],V_j[3],r[3],v[3];
    double T = 0.0;
    double V = 0.0;
    
    for (i=0; i<R->NumVertex; i++)
    {
        m_i = R->Mass[i];
        for (k=0; k<3; k++)
        {
            R_i[k] = R->Pos[3 * i + k];
            V_i[k] = R->Vel[3 * i + k];
        }
        
        T += 0.5 * m_i * norm3_squared(V_i);
        
        for (j=i+1; j<R->NumVertex; j++)
        {
            m_j = R->Mass[j];
            
            for (k=0; k<3; k++)
            {
                
                R_j[k] = R->Pos[3 * j + k];
                r[k] = R_i[k] - R_j[k];
            }
                
            V += - CONST_G * m_i * m_j/norm3(r);
        }
    }
    
    return T + V;
}

void mstar_print_state(struct RegularizedRegion *R)
{
    for (int i=0; i<R->NumVertex; i++)
    {
        printf("i %d MSEIndex %d mass %g radius %g stopping condition partner %d\n",i,R->MSEIndex[i],R->Mass[i],R->Radius[i],R->Stopping_Condition_Partner[i]);
        printf("i %d pos %g %g %g \n",i,R->Pos[3 * i + 0], R->Pos[3 * i + 1],R->Pos[3 * i + 2]);
        printf("i %d vel %g %g %g \n",i,R->Vel[3 * i + 0], R->Vel[3 * i + 1],R->Vel[3 * i + 2]);
        printf("i %d spin AM %g %g %g \n",i,R->Spin_S[3 * i + 0], R->Spin_S[3 * i + 1],R->Spin_S[3 * i + 2]);
    }
    printf("=========\n");
}

void integrate_nbody_system_with_mass_loss(double end_time, int Nsteps, std::vector<double> &masses, std::vector<double> &delta_masses, std::vector<std::vector<double> > &R_vecs, std::vector<std::vector<double> > &V_vecs)
{
    /* Used to take into account effect of mass loss during events such as CE.
     * No collision detection. */

    #ifdef VERBOSE
    if (verbose_flag > 0)
    {
        printf("nbody_evolution.cpp -- integrate_nbody_system_with_mass_loss end_time %g Nsteps %d \n",end_time,Nsteps);
    }
    #endif
        
    initialize_mpi_or_serial(); // This needs to be done to initialize MSTAR

    struct RegularizedRegion *R;

    int i=0;
    int j;
    int N_bodies = masses.size();
        
    allocate_armst_structs(&R, N_bodies); // Initialize the data structure

    for (auto it_m = masses.begin(); it_m != masses.end(); it_m++)
    {
        R->Vertex[i].type = 0;
        R->MSEIndex[i] = i;

        R->Mass[i] = *it_m;
        R->Radius[i] = 0.0;
        for (j=0; j<3; j++)
        {
            R->Pos[3 * i + j] = R_vecs[i][j];
            R->Vel[3 * i + j] = V_vecs[i][j];
        }
        
        #ifdef VERBOSE
        if (verbose_flag > 1)
        {
            printf("nbody_evolution.cpp -- integrate_nbody_system_with_mass_loss m %g R %g %g %g V %g %g %g\n",R->Mass[i],R->Pos[3 * i + 0],R->Pos[3 * i + 1],R->Pos[3 * i + 2],R->Vel[3 * i + 0],R->Vel[3 * i + 1],R->Vel[3 * i + 2]);
        }
        #endif

        i++;
    }

    R->gbs_tolerance = MSTAR_gbs_tolerance_default;
    R->stopping_condition_tolerance = MSTAR_stopping_condition_tolerance;
    R->output_time_tolerance = MSTAR_output_time_tolerance;

    MSTAR_verbose = false;

    double t=0;
    double t_reached;
    int collision_occurred;
    double dt = end_time/( (double) Nsteps);
    while (t <= end_time)
    {
        run_integrator(R, dt, &t_reached, &collision_occurred);

        t += t_reached;

        /* Increment the masses */
        i=0;
        for (auto it_dm = delta_masses.begin(); it_dm != delta_masses.end(); it_dm++)
        {
            R->Mass[i] += *it_dm/ ((double) Nsteps);
            i++;
        }
    }
    
    i=0;
    for (auto it_m = masses.begin(); it_m != masses.end(); it_m++)
    {
        for (j=0; j<3; j++)
        {
            R_vecs[i][j] = R->Pos[3 * i + j];
            V_vecs[i][j] = R->Vel[3 * i + j];
        }

        #ifdef VERBOSE
        if (verbose_flag > 1)
        {
            printf("nbody_evolution.cpp -- integrate_nbody_system_with_mass_loss -- done m %g R %g %g %g V %g %g %g\n",R->Mass[i],R->Pos[3 * i + 0],R->Pos[3 * i + 1],R->Pos[3 * i + 2],R->Vel[3 * i + 0],R->Vel[3 * i + 1],R->Vel[3 * i + 2]);
        }
        #endif

        i++;
    }
    
    MSTAR_verbose = true;

    return;
}

int check_particlesMap_for_inclusion_in_MSTAR(ParticlesMap *particlesMap)
{
    int N_bodies = 0;

    ParticlesMapIterator it_p;
    for (it_p = particlesMap->begin(); it_p != particlesMap->end(); it_p++)
    {
        Particle *p = (*it_p).second;
        if (p->is_binary == false)
        {
            if (find_nearest_neighbor_separation(particlesMap, p->index, p->R_vec) > nbody_maximum_separation_for_inclusion)
            {
                #ifdef VERBOSE
                if (verbose_flag > 1)
                {
                    printf("nbody_evolution.cpp -- check_particlesMap_for_inclusion_in_MSTAR -- Excluding body index %d norm3(R_vec) %g sep %g\n",p->index,norm3(p->R_vec),find_nearest_neighbor_separation(particlesMap, p->index, p->R_vec));
                }
                #endif
                p->include_in_MSTAR = false;
            }
            else
            {
                p->include_in_MSTAR = true;
                N_bodies += 1;
            }
        }
    }

    return N_bodies;
}

}
