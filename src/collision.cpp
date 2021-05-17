/* MSE */

#include "evolve.h"
#include "collision.h"

extern "C"
{

void handle_collisions(ParticlesMap *particlesMap, double t, int *integration_flag)
{
    /* invoke CE if one star has k in 2,3,4,5,6,8,9 */

    ParticlesMapIterator it_p;
        
    if (*integration_flag == 0) /* secular mode */
    {
    
        std::vector<int> CE_parent_indices;
        std::vector<int> CE_donor_indices;
        std::vector<int> CE_accretor_indices;

        std::vector<int> col_parent_indices;
        std::vector<int> col_C1_indices;
        std::vector<int> col_C2_indices;
    
        for (it_p = particlesMap->begin(); it_p != particlesMap->end(); it_p++)
        {
            Particle *p = (*it_p).second;

            if (p->is_binary == true and p->merged == true)
            {
                Particle *child1 = (*particlesMap)[p->child1];
                Particle *child2 = (*particlesMap)[p->child2];

                int kw1 = child1->stellar_type;
                int kw2 = child2->stellar_type;
                if (kw1 >= 2 and kw1 <= 9 and kw1 != 7) /* CE in case of collision of giant + any other star */
                {
                    //common_envelope_evolution(particlesMap, p->index, child1->index, child2->index, t, integration_flag);
                    CE_parent_indices.push_back(p->index);
                    CE_donor_indices.push_back(child1->index);
                    CE_accretor_indices.push_back(child2->index);
                }
                else if (kw2 >= 2 and kw2 <= 9 and kw2 != 7) /* CE in case of collision of giant + any other star */
                {
                    //common_envelope_evolution(particlesMap, p->index, child2->index, child1->index, t, integration_flag);
                    CE_parent_indices.push_back(p->index);
                    CE_donor_indices.push_back(child2->index);
                    CE_accretor_indices.push_back(child1->index);
                }
                else /* "pure" collision in other cases */
                {
                    //collision_product(particlesMap, p->index, child1->index, child2->index, t, integration_flag);
                    col_parent_indices.push_back(p->index);
                    col_C1_indices.push_back(child1->index);
                    col_C2_indices.push_back(child2->index);
                }
            }
        }
        
        int i;
        for (i=0; i<CE_parent_indices.size(); i++)
        {
            common_envelope_evolution(particlesMap, CE_parent_indices[i], CE_donor_indices[i], CE_accretor_indices[i], t, integration_flag);
        }
        for (i=0; i<col_parent_indices.size(); i++)
        {
            collision_product(particlesMap, col_parent_indices[i], col_C1_indices[i], col_C2_indices[i], t, integration_flag);
        }
        
    }
    else /* N-body mode */
    {
        int col_part_i = -1;
        int col_part_j = -1;
        ParticlesMapIterator it_i,it_j;
        Particle *pi, *pj;
        for (it_i = particlesMap->begin(); it_i != particlesMap->end(); it_i++)
        {
            pi = (*it_i).second;
            col_part_i = pi->Stopping_Condition_Partner;
            if (col_part_i != -1)
            {
                for (it_j = particlesMap->begin(); it_j != particlesMap->end(); it_j++)
                {
                    pj = (*it_j).second;
                    col_part_j = pj->Stopping_Condition_Partner;
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
            printf("merger.cpp -- error in handle_collisions: unable to find pair of colliding bodies; col_part_i %d col_part_j %d\n",col_part_i,col_part_j);
            error_code = 7;
            longjmp(jump_buf,1);
            //exit(-1);
        }
        /* Reset Stopping_Condition_Partner */
        pi->Stopping_Condition_Partner = -1;
        pj->Stopping_Condition_Partner = -1;
        
        //printf("col_part_i %d col_part_j %d pi %d pj %d\n",col_part_i,col_part_j,pi->index,pj->index);
        
        int binary_index = -1; // binary index from particlesMap does not apply in this case        
        int kw1 = pi->stellar_type;
        int kw2 = pj->stellar_type;
        if (kw1 >= 2 and kw1 <= 9 and kw1 != 7) /* CE in case of collision of giant + any other star */
        {
            common_envelope_evolution(particlesMap, binary_index, pi->index, pj->index, t, integration_flag);
        }
        else if (kw2 >= 2 and kw2 <= 9 and kw2 != 7) /* CE in case of collision of giant + any other star */
        {
            common_envelope_evolution(particlesMap, binary_index, pj->index, pi->index, t, integration_flag);
        }
        else /* "pure" collision in other cases */
        {
            collision_product(particlesMap, binary_index, pi->index, pj->index, t, integration_flag);
        }
    }
    
    *integration_flag = determine_orbits_in_system_using_nbody(particlesMap);
  
    #ifdef VERBOSE
    if (verbose_flag > 1)
    {
        printf("collision.cpp -- handle_collisions -- done\n");
        print_system(particlesMap,*integration_flag);
    }
    #endif
}


void collision_product(ParticlesMap *particlesMap, int binary_index, int child1_index, int child2_index, double t, int *integration_flag)
{
    
    #ifdef LOGGING
    Log_info_type log_info;
    log_info.binary_index = binary_index;
    log_info.index1 = child1_index;
    log_info.index2 = child2_index;
    update_log_data(particlesMap, t, *integration_flag, LOG_COL_START, log_info);
    #endif

    #ifdef VERBOSE
    if (verbose_flag > 0)
    {
        printf("collision.cpp -- collision_product -- CP binary %d child1 %d child2 %d t %g int flag %d\n",binary_index,child1_index,child2_index,t,*integration_flag);
        print_system(particlesMap,*integration_flag);
    }
    #endif

    Particle *child1 = (*particlesMap)[child1_index];
    Particle *child2 = (*particlesMap)[child2_index];
    
   
    /* Check for collisions that are NOT of the star-star type */
    if (child1->object_type == 1 and child2->object_type == 2)
    {
        collision_product_star_planet(particlesMap, binary_index, child1_index, child2_index,t,integration_flag);
        return;
    }
    else if (child1->object_type == 2 and child2->object_type == 1)
    {
        collision_product_star_planet(particlesMap, binary_index, child2_index, child1_index,t,integration_flag);
        return;
    }
    else if (child1->object_type == 2 and child2->object_type == 2)
    {
        collision_product_planet_planet(particlesMap, binary_index, child1_index, child2_index,t,integration_flag);
        return;
    }
    
    /* Below: star-star case */
    double m1 = child1->mass;
    double m2 = child2->mass;
    double m = m1 + m2;

    int i;

    double initial_momentum[3],initial_R_CM[3],initial_V_CM[3];
    double h_vec[3],h_vec_unit[3],e_vec[3],e_vec_unit[3],r_vec[3],v_vec[3];

    get_initial_binary_orbital_properties_from_position_and_velocity(child1->R_vec, child1->V_vec, child2->R_vec, child2->V_vec, m1, m2, \
        r_vec, v_vec, initial_momentum, initial_R_CM, initial_V_CM, h_vec, e_vec);
    get_unit_vector(h_vec,h_vec_unit);
    get_unit_vector(e_vec,e_vec_unit);

    double a_old = compute_a_from_h(m1,m2,norm3(h_vec),norm3(e_vec));
    if (a_old <= epsilon) /* Take instantenous separation if semimajor axis not well defined */
    {
        double r[3];
        for (int i=0; i<3; i++)
        {
            r[i] = child1->R_vec[i] - child2->R_vec[i];
        }
        a_old = norm3(r);
    }
    double P_old = compute_orbital_period_from_semimajor_axis(m,a_old);
    double n_old = TWOPI/P_old;
    
    #ifdef VERBOSE
    if (verbose_flag > 0)
    {
        printf("collision.cpp -- collision_product -- a_old %g P_old %g n_old %g\n",a_old,P_old,n_old);
    }
    #endif
    
    int kw = determine_merger_type(child1->stellar_type,child2->stellar_type);

    int kw1 = child1->stellar_type;
    int kw2 = child2->stellar_type;

    double age1 = child1->age*yr_to_Myr;
    double age2 = child2->age*yr_to_Myr;

    double tm1 = child1->sse_main_sequence_timescale*yr_to_Myr;
    double tm2 = child2->sse_main_sequence_timescale*yr_to_Myr;

    double *spin_vec_1 = child1->spin_vec;
    double *spin_vec_2 = child2->spin_vec;

    double spin_vec_1_norm = norm3(spin_vec_1);
    double spin_vec_2_norm = norm3(spin_vec_2);

    /* Spin parameters for NS/BHs */
    double chi1 = compute_spin_parameter_from_spin_frequency(m1, spin_vec_1_norm);
    double chi2 = compute_spin_parameter_from_spin_frequency(m2, spin_vec_2_norm);
    double spin_vec_1_unit[3],spin_vec_2_unit[3];
    get_unit_vector(spin_vec_1, spin_vec_1_unit);
    get_unit_vector(spin_vec_2, spin_vec_2_unit);

    /* Properties of merged object */
    double m0,mc,age,t_MS; /* new initial mass, new core mass, new age, new MS timescale */
    age = -1;
    m0 = -1;
    
    double tm,tn;

    double *GB,*tscls,*lums;
    GB = new double[10];
    tscls = new double[20];
    lums = new double[10];    
    double *zpars;
    zpars = child1->zpars; /* for now: take metallicity of star 1 (arbitrary); should think of Z-mixing in future */
  
    bool destroyed = false;
    bool reset_spin_vec = true;
    
    int kick_distribution = child1->kick_distribution; /* by default, adopt kick distribution of child1 */

    double V_kick_vec[3] = {0.0,0.0,0.0};
    double spin_vec[3];
    double M_final; /* used for compact objects */
    double alpha_vec_final[3]; /* used for compact objects */

    /* Two H/He MS stars forming new MS star */
    if ((kw1 >= 0 and kw1 <= 1 or kw1 == 7) and (kw2 >= 0 and kw2 <= 1 or kw2 == 7)) 
    {
        m0 = m;
        star_(&kw, &m0, &m, &tm, &tn, tscls, lums, GB, zpars);
        
        age = 0.1 * (tm / m) * ( (m1 * age1/tm1) + (m2 * age2/tm2) ); /* HTP02 eq. (80) */
    }
    
    /* Compact star sinks to center of MS star forming giant star with core mass of compact star */
    else if ((kw1 >= 0 and kw1 <= 1) and (kw2 == 7 or kw2 >= 10 and kw2 <= 12)) /* compact star is child2 */
    {
        mc = m2;
        gntage_(&mc,&m,&kw,zpars,&m0,&age);
    }
    else if ((kw2 >= 0 and kw2 <= 1) and (kw1 == 7 or kw1 >= 10 and kw1 <= 12)) /* compact star is child1 */
    {
        mc = m1;
        gntage_(&mc,&m,&kw,zpars,&m0,&age);
    }
    
    /* Rejuvenation of HeMS star by absorbing He WD*/
    else if ((kw1 == 7) and (kw2 == 10))
    {
        m0 = m;
        star_(&kw, &m0, &m, &tm, &tn, tscls, lums, GB, zpars);
        age = tm * m1 * age1 / (m * tm1);
    }
    else if ((kw1 == 10) and (kw2 == 7))
    {
        m0 = m;
        star_(&kw, &m0, &m, &tm, &tn, tscls, lums, GB, zpars);
        age = tm * m2 * age2 / (m * tm2);
    }
    
    /* Making evolved He star by absorbing CO/ONe WD */
    else if ((kw1 == 7) and (kw2 >= 11 and kw2 <= 12))
    {
        mc = m2;
        gntage_(&mc,&m,&kw,zpars,&m0,&age);
    }
    else if ((kw2 == 7) and (kw1 >= 11 and kw1 <= 12))
    {
        mc = m1;
        gntage_(&mc,&m,&kw,zpars,&m0,&age);
    }
        
    /* He WD + He WD */
    else if (kw1 == 10 and kw2 == 10)
    {
        #ifdef LOGGING
        Log_info_type log_info;
        log_info.index1 = child1->index;
        log_info.index2 = child2->index;
        update_log_data(particlesMap, t, *integration_flag, LOG_SNE_START, log_info);
        #endif
        
        destroyed = true;
    }
    
    /* He WD + CO/ONe WD: make new evolved He star */
    else if (kw1 == 10 and kw2 >= 11 and kw2 <= 12)
    {
        mc = m2;
        gntage_(&mc,&m,&kw,zpars,&m0,&age);
    }
    else if (kw2 == 10 and kw1 >= 11 and kw1 <= 12)
    {
        mc = m1;
        gntage_(&mc,&m,&kw,zpars,&m0,&age);
    }
    
    /* CO WD + CO WD */
    else if (kw1 == 11 and kw2 == 11)
    {
        m0 = m;
        age = 0.0;
        if (m >= chandrasekhar_mass)
        {
            destroyed = true; /* SNe Ia */

            #ifdef LOGGING
            Log_info_type log_info;
            log_info.index1 = child1->index;
            log_info.index2 = child2->index;
            update_log_data(particlesMap, t, *integration_flag, LOG_SNE_START, log_info);
            #endif
           
        }
    }
    
    /* ONe WD + CO/ONe WD */
    else if (kw1 == 12 and kw2 >= 11 and kw2 <= 12)
    {
        m0 = m;
        age = 0.0;
        if (m <= chandrasekhar_mass)
        {
            kw = 12;
        }
        else
        {
            kw = 13; /* AIC to NS */
        }
    }
    else if (kw2 == 12 and kw1 >= 11 and kw1 <= 12)
    {
        m0 = m;
        age = 0.0;
        if (m <= chandrasekhar_mass)
        {
            kw = 12;
        }
        else
        {
            kw = 13; /* AIC to NS */
        }
    }
    
    /* Thorne-Zytkow object */
    else if ((kw1 >= 13 and kw1 <= 14) and (kw2 == 0 or kw2 == 1 or kw2 == 7) )
    {
        kw = kw1;
        m = m1;
        m0 = m;
        age = 0.0;
    }
    else if ((kw2 >= 13 and kw2 <= 14) and (kw1 == 0 or kw1 == 1 or kw1 == 7) )
    {
        kw = kw2;
        m = m2;
        m0 = m;
        age = 0.0;
    }
    
    /* Other compact object mergers */
    else if (kw1 >= 13 and kw1 <= 14 and kw2 >= 10 and kw2 <= 14)
    {
        double v_recoil_vec[3];
        determine_compact_object_merger_properties(m1,m2,chi1,chi2,spin_vec_1_unit,spin_vec_2_unit,h_vec_unit,e_vec_unit,v_recoil_vec,alpha_vec_final,&M_final);
        m = M_final;

        reset_spin_vec = false;
        double chi = norm3(alpha_vec_final);
        double Omega = compute_spin_frequency_from_spin_parameter(m,chi);
        for (i=0; i<3; i++)
        {
            V_kick_vec[i] = v_recoil_vec[i];
            spin_vec[i] = Omega * alpha_vec_final[i] / chi;
        }
        
        kw = kw1;
        m0 = m;
        age = 0.0;
        if (kw == 13 and m > 1.8)
        {
            kw = 14;
        }
        
        #ifdef VERBOSE
        if (verbose_flag > 0)
        {
            printf("collision.cpp -- collision product -- kw2 >= 13 and kw2 <= 14 and kw1 >= 10 and kw1 <= 14 m %g chi %g Omega %g v_k %g %g %g spin %g %g %g\n",m,chi,Omega,v_recoil_vec[0],v_recoil_vec[1],v_recoil_vec[2],spin_vec[0],spin_vec[1],spin_vec[2]);
        }
        #endif
        
    }
    else if (kw2 >= 13 and kw2 <= 14 and kw1 >= 10 and kw1 <= 14)
    {
        double v_recoil_vec[3];
        determine_compact_object_merger_properties(m2,m1,chi2,chi1,spin_vec_2_unit,spin_vec_1_unit,h_vec_unit,e_vec_unit,v_recoil_vec,alpha_vec_final,&M_final);
        m = M_final;

        reset_spin_vec = false;
        double chi = norm3(alpha_vec_final);
        double Omega = compute_spin_frequency_from_spin_parameter(m,chi);
        for (i=0; i<3; i++)
        {
            V_kick_vec[i] = v_recoil_vec[i];
            spin_vec[i] = Omega * alpha_vec_final[i] / chi;
        }
        
        kw = kw2;
        m0 = m;
        age = 0.0;
        if (kw == 13 and m > 1.8)
        {
            kw = 14;
        }
    }
    else
    {
        printf("merger.cpp -- collision_product -- unknown outcome! Will ignore the collision. kw1 %d kw2 %d kw %d\n",kw1,kw2,kw);
        error_code = 6;
        longjmp(jump_buf,1);
        return;
    }
    
    if (destroyed == false and (m0 == -1 or age == -1))
    {
        printf("merger.cpp -- collision_product -- destroyed = false -- was not able to determine all properties of merged object! Will ignore the collision. \n");
        error_code = 6;
        longjmp(jump_buf,1);
        return;
    }

    /* stellar properties of merged object */
    double Omega_crit,Omega;
    double r,lum,rc,menv,renv,k2;
    double z_new;
    double *zpars_new;
    zpars_new = new double[20];

    if (destroyed == false)
    {
        star_(&kw, &m0, &m, &tm, &tn, tscls, lums, GB, zpars);

        hrdiag_(&m0,&age,&m,&tm,&tn,tscls,lums,GB,zpars, \
            &r,&lum,&kw,&mc,&rc,&menv,&renv,&k2);

        z_new = child1->metallicity;

        zcnsts_(&z_new,zpars_new);
    }

    *integration_flag = 1;
    if (destroyed == false)
    {
        /* child1 will become the merged object */
        particlesMap->erase(child2->index);
        //it_p = particlesMap->erase(particlesMap->find(child2->index));

        /* Update properties of merged star */
        child1->stellar_type = kw;
        child1->mass = m;
        child1->sse_initial_mass = m0;
        
        child1->age = age*Myr_to_yr;
        child1->epoch = t - child1->age;
        child1->sse_main_sequence_timescale = tm*Myr_to_yr;

        child1->radius = r*CONST_R_SUN;
        child1->luminosity = lum*CONST_L_SUN;
        child1->core_mass = mc;
        child1->core_radius = rc*CONST_R_SUN;
        child1->convective_envelope_mass = menv;
        child1->convective_envelope_radius = renv*CONST_R_SUN;

        child1->metallicity = z_new;
        child1->zpars = zpars_new;

        for (i=0; i<3; i++)
        {
            child1->R_vec[i] = initial_R_CM[i];
            child1->V_vec[i] = initial_momentum[i]/m; /* set new velocity according to linear momentum conservation */
            
            child1->V_vec[i] += V_kick_vec[i]; /* apply possible kick */
        }

        /* Unless calculated above, set the spin equal to the orbital frequency just before collision
         * Assume the direction is equal to the previous orbital orientation. */
        if (reset_spin_vec == true)
        {
            Omega_crit = compute_breakup_angular_frequency(child1->mass,child1->radius);
            Omega = n_old;
            if (n_old >= Omega_crit)
            {
                #ifdef VERBOSE
                if (verbose_flag > 0)
                {
                    printf("collision.cpp -- collision product -- Limiting spin frequency of star %d from %g to breakup rate %g\n",child1->index,n_old,Omega_crit);
                }
                #endif

                Omega = Omega_crit;
            }
            
            for (i=0; i<3; i++)
            {
                child1->spin_vec[i] = Omega * h_vec_unit[i];
            }
        }
        else
        {
            for (i=0; i<3; i++)
            {
                child1->spin_vec[i] = spin_vec[i];
            }
        }
        reset_ODE_mass_dot_quantities(child1);
        //child1->RLOF_flag = 0;
        //child1->apply_kick = false;
        //child1->mass_dot_wind = 0.0;
        //child1->mass_dot_wind_accretion = 0.0;
        //child1->radius_dot = 0.0;
        //child1->ospin_dot = 0.0;
        
        if (NS_model == 1)
        {
            if (kw == 13)
            {
                if (kw1 < 13 and kw2 < 13) /* A NS was formed */
                {
                    double ospin;
                    compute_NS_formation_properties_Ye19_model(false, &ospin, &child1->magnetic_field_strength_gauss);
                    rescale_vector(child1->spin_vec, ospin/norm3(child1->spin_vec));
                    
                    child1->time_of_NS_formation = t;
                    child1->initial_NS_period_s = compute_spin_period_from_spin_angular_frequency(ospin) * yr_to_s;
                    child1->initial_magnetic_field_strength_gauss = child1->magnetic_field_strength_gauss;
                }
                
                /* Check for existing MSPs */
                if ((kw1 == 13 and kw2 < 13) or (kw2 == 13 and kw1 < 13))
                {
                    if ((kw1 == 13 and (determine_if_NS_is_MSP(spin_vec_1_norm,child1->magnetic_field_strength_gauss) == true)) or (kw2 == 13 and (determine_if_NS_is_MSP(spin_vec_2_norm,child2->magnetic_field_strength_gauss) == true)))
                    {
                        /* This block is entered if either of the original objects was a MSP, and the new object is also a NS */
                        double ospin;
                        compute_NS_formation_properties_Ye19_model(true, &ospin, &child1->magnetic_field_strength_gauss);
                        rescale_vector(child1->spin_vec,ospin/norm3(child1->spin_vec));

                        child1->time_of_NS_formation = t;
                        child1->initial_NS_period_s = compute_spin_period_from_spin_angular_frequency(ospin) * yr_to_s;
                        child1->initial_magnetic_field_strength_gauss = child1->magnetic_field_strength_gauss;
                    }
                }
            }
        }
        
    }
    else
    {
        particlesMap->erase(child1->index);
        particlesMap->erase(child2->index);
    }

    reset_RLOF_flags(particlesMap);

    #ifdef LOGGING
    //Log_info_type log_info;
    log_info.binary_index = binary_index;
    log_info.index1 = child1_index;
    log_info.index2 = child2_index;
    update_log_data(particlesMap, t, *integration_flag, LOG_COL_END, log_info);
    #endif

    #ifdef VERBOSE
    if (verbose_flag > 0)
    {
        printf("collision.cpp -- collision product -- done\n");
        print_system(particlesMap,*integration_flag);
    }
    #endif

    delete[] GB;
    delete[] tscls;
    delete[] lums;

    return;
}

int determine_merger_type(int kw1, int kw2)
{
    return MERGER_TABLE[kw1][kw2];
}

void collision_product_star_planet(ParticlesMap *particlesMap, int binary_index, int star_index, int planet_index, double t, int *integration_flag)
{
    #ifdef VERBOSE
    if (verbose_flag > 0)
    {
        printf("collision.cpp -- collision_product_star_planet  binary_index %d  star_index %d  planet_index %d  t %g integration_flag %d\n",binary_index, star_index, planet_index,t,&integration_flag);
        print_system(particlesMap,*integration_flag);
    }
    #endif

    Particle *star = (*particlesMap)[star_index];
    Particle *planet = (*particlesMap)[planet_index];
    
    int i;
    double m1 = star->mass;
    double m2 = planet->mass;
    double m = m1 + m2;

    double initial_momentum[3],initial_R_CM[3],initial_V_CM[3];
    double h_vec[3],h_vec_unit[3],e_vec[3],e_vec_unit[3],r_vec[3],v_vec[3];

    get_initial_binary_orbital_properties_from_position_and_velocity(star->R_vec, star->V_vec, planet->R_vec, planet->V_vec, m1, m2, \
        r_vec, v_vec, initial_momentum, initial_R_CM, initial_V_CM, h_vec, e_vec);
    get_unit_vector(h_vec,h_vec_unit);
    get_unit_vector(e_vec,e_vec_unit);

    double a_old = compute_a_from_h(m1,m2,norm3(h_vec),norm3(e_vec));
    double P_old = compute_orbital_period_from_semimajor_axis(m,a_old);
    double n_old = TWOPI/P_old;
    
    int stellar_type = star->stellar_type;
    double mass = star->mass;
    double sse_initial_mass = star->sse_initial_mass;
    
    double age = star->age;
    double sse_main_sequence_timescale = star->sse_main_sequence_timescale;

    double radius = star->radius;
    double luminosity = star->luminosity;
    double core_mass = star->core_mass;
    double core_radius = star->core_radius;
    double convective_envelope_mass = star->convective_envelope_mass;
    double convective_envelope_radius = star->convective_envelope_radius;
    double metallicity = star->metallicity;
    double *zpars = star->zpars;

    double v_kick_vec[3] = {0.0,0.0,0.0};
    bool reset_spin_vec = false;
    double *spin_vec = star->spin_vec;
    
    *integration_flag = 1;
    
    particlesMap->erase(planet->index);
    star->stellar_type = stellar_type;
    star->mass = mass;
    star->sse_initial_mass = sse_initial_mass;
    
    star->age = age;
    star->epoch = t - star->age;
    star->sse_main_sequence_timescale = sse_main_sequence_timescale;

    star->radius = radius;
    star->luminosity = luminosity;
    star->core_mass = core_mass;
    star->core_radius = core_radius;
    star->convective_envelope_mass = convective_envelope_mass;
    star->convective_envelope_radius = convective_envelope_radius;

    star->metallicity = metallicity;
    star->zpars = zpars;

    for (i=0; i<3; i++)
    {
        star->R_vec[i] = initial_R_CM[i];
        star->V_vec[i] = initial_momentum[i]/m; /* set new velocity according to linear momentum conservation */
        
        star->V_vec[i] += v_kick_vec[i]; /* apply possible kick */
    }


    /* Unless calculated above, set the spin equal to the orbital frequency just before collision
     * Assume the direction is equal to the previous orbital orientation. */
    double Omega_crit, Omega;
    if (reset_spin_vec == true)
    {

        Omega_crit = compute_breakup_angular_frequency(star->mass,star->radius);
        Omega = n_old;
        if (n_old >= Omega_crit)
        {
            #ifdef VERBOSE
            if (verbose_flag > 0)
            {
                printf("collision.cpp -- collision_product_star_planet --  Limiting spin frequency of star %d from %g to breakup rate %g\n",star->index,n_old,Omega_crit);
                print_system(particlesMap,*integration_flag);
            }
            #endif

            Omega = Omega_crit;
        }
        
        for (i=0; i<3; i++)
        {
            star->spin_vec[i] = Omega * h_vec_unit[i];
        }
    }
    else
    {
        for (i=0; i<3; i++)
        {
            star->spin_vec[i] = spin_vec[i];
        }
    }
    star->apply_kick = false; /* Default value for (stellar evolution) bodies */
    star->RLOF_flag = 0;
    
    #ifdef VERBOSE
    if (verbose_flag > 0)
    {
        printf("collision.cpp -- collision_product_star_planet post  binary_index %d  star_index %d  planet_index %d  t %g integration_flag %d\n",binary_index, star_index, planet_index,t,&integration_flag);
        print_system(particlesMap,*integration_flag);
    }
    #endif

    return;
}

void collision_product_planet_planet(ParticlesMap *particlesMap, int binary_index, int planet1_index, int planet2_index, double t, int *integration_flag)
{
    #ifdef VERBOSE
    if (verbose_flag > 0)
    {
        printf("collision.cpp -- collision_product_planet_planet  binary_index %d  planet1_index %d  planet2_index %d  t %g integration_flag %d\n",binary_index, planet1_index, planet2_index,t,&integration_flag);
        print_system(particlesMap,*integration_flag);
    }
    #endif
    
    
    Particle *planet1 = (*particlesMap)[planet1_index];
    Particle *planet2 = (*particlesMap)[planet2_index];
    
    int i;
    double m1 = planet1->mass;
    double m2 = planet2->mass;
    double m = m1 + m2;

    double initial_momentum[3],initial_R_CM[3],initial_V_CM[3];
    double h_vec[3],h_vec_unit[3],e_vec[3],e_vec_unit[3],r_vec[3],v_vec[3];

    get_initial_binary_orbital_properties_from_position_and_velocity(planet1->R_vec, planet1->V_vec, planet2->R_vec, planet2->V_vec, m1, m2, \
        r_vec, v_vec, initial_momentum, initial_R_CM, initial_V_CM, h_vec, e_vec);
    get_unit_vector(h_vec,h_vec_unit);
    get_unit_vector(e_vec,e_vec_unit);

    double a_old = compute_a_from_h(m1,m2,norm3(h_vec),norm3(e_vec));
    double P_old = compute_orbital_period_from_semimajor_axis(m,a_old);
    double n_old = TWOPI/P_old;
    
    double v_kick_vec[3] = {0.0,0.0,0.0};
    bool reset_spin_vec = false;
    double *spin_vec = planet1->spin_vec;


    *integration_flag = 1;
    
    /* Remove planet2 and merge into planet1 */
    particlesMap->erase(planet2->index);
    planet1->mass = m;
    //planet1->radius = ... /* TO DO: what is radius of merged planet?
    
    for (i=0; i<3; i++)
    {
        planet1->R_vec[i] = initial_R_CM[i];
        planet1->V_vec[i] = initial_momentum[i]/m; /* set new velocity according to linear momentum conservation */
        
        planet1->V_vec[i] += v_kick_vec[i]; /* apply possible kick */
    }


    /* Unless calculated above, set the spin equal to the orbital frequency just before collision
     * Assume the direction is equal to the previous orbital orientation. */
    double Omega_crit, Omega;
    if (reset_spin_vec == true)
    {
        Omega_crit = compute_breakup_angular_frequency(planet1->mass,planet1->radius);
        Omega = n_old;
        if (n_old >= Omega_crit)
        {
            #ifdef VERBOSE
            if (verbose_flag > 0)
            {
                printf("collision.cpp -- collision_product_planet_planet -- Limiting spin frequency of star %d from %g to breakup rate %g\n",planet1->index,n_old,Omega_crit);
                print_system(particlesMap,*integration_flag);
            }
            #endif


            Omega = Omega_crit;
        }
        
        for (i=0; i<3; i++)
        {
            planet1->spin_vec[i] = Omega * h_vec_unit[i];
        }
    }
    else
    {
        for (i=0; i<3; i++)
        {
            planet1->spin_vec[i] = spin_vec[i];
        }
    }
    planet1->apply_kick = false; /* Default value for (stellar evolution) bodies */
    planet1->RLOF_flag = 0;
    
    #ifdef VERBOSE
    if (verbose_flag > 0)
    {
        printf("collision.cpp -- collision_product_planet_planet post  binary_index %d  planet1_index %d  planet2_index %d  t %g integration_flag %d\n",binary_index, planet1_index, planet2_index,t,&integration_flag);
        print_system(particlesMap,*integration_flag);
    }
    #endif


    return;
}
    

void determine_compact_object_merger_properties(double m1_, double m2_, double chi1_, double chi2_, double spin_vec_1_unit_[3], double spin_vec_2_unit_[3], double h_vec_unit[3], double e_vec_unit[3], double v_recoil_vec[3], double alpha_vec_final[3], double *M_final)
{
    /* Based on fits from Louto et al. 2010 (https://ui.adsabs.harvard.edu/abs/2010CQGra..27k4006L/abstract) */
    /* The fits assume m1 <= m2; adjust input parameters so this will always be satisfied */

    int i;
    
    double m1,m2,chi1,chi2;
    double spin_vec_1_unit[3],spin_vec_2_unit[3];
    if (m1_ <= m2_)
    {
        m1 = m1_;
        m2 = m2_;
        chi1 = chi1_;
        chi2 = chi2_;
        for (i=0; i<3; i++)
        {
            spin_vec_1_unit[i] = spin_vec_1_unit_[i];
            spin_vec_2_unit[i] = spin_vec_2_unit_[i];
        }
    }
    else
    {
        m1 = m2_;
        m2 = m1_;
        chi1 = chi2_;
        chi2 = chi1_;
        for (i=0; i<3; i++)
        {
            spin_vec_1_unit[i] = spin_vec_2_unit_[i];
            spin_vec_2_unit[i] = spin_vec_1_unit_[i];
        }
    }

    double q = m1/m2; /* q <= 1 */
    double q_p2 = q*q;
    double one_plus_q = 1.0 + q;
    double one_plus_q_p2 = one_plus_q * one_plus_q;
    
    double M = m1 + m2;
    double eta = q/( one_plus_q_p2 );
    double eta_p2 = eta*eta;
    double eta_p3 = eta*eta_p2;

    double alpha1_vec[3],alpha2_vec[3];
    double alpha2_vec_plus_q_times_alpha1_vec[3],alpha2_vec_minus_q_times_alpha1_vec[3];
    double alpha2_perp_vec_plus_q_p2_times_alpha1_perp_vec[3],alpha2_perp_vec_minus_q_times_alpha1_perp_vec[3];
    
    for (i=0; i<3; i++)
    {
        alpha1_vec[i] = chi1 * spin_vec_1_unit[i];
        alpha2_vec[i] = chi2 * spin_vec_2_unit[i];
        alpha2_vec_plus_q_times_alpha1_vec[i] = alpha2_vec[i] + q * alpha1_vec[i];
        alpha2_vec_minus_q_times_alpha1_vec[i] = alpha2_vec[i] - q * alpha1_vec[i];
    }

    double alpha1_par,alpha2_par,alpha1_perp,alpha2_perp;
    double alpha1_perp_vec[3],alpha2_perp_vec[3];
    get_parallel_and_perpendicular_vectors_and_components(alpha1_vec, h_vec_unit, &alpha1_par, &alpha1_perp, alpha1_perp_vec);
    get_parallel_and_perpendicular_vectors_and_components(alpha2_vec, h_vec_unit, &alpha2_par, &alpha2_perp, alpha2_perp_vec);

    for (i=0; i<3; i++)
    {
        alpha2_perp_vec_plus_q_p2_times_alpha1_perp_vec[i] = alpha2_perp_vec[i] + q_p2 * alpha1_perp_vec[i];
        alpha2_perp_vec_minus_q_times_alpha1_perp_vec[i] = alpha2_perp_vec[i] - q * alpha1_perp_vec[i];
    }
   
    double alpha1_par_p2 = alpha1_par * alpha1_par;
    double alpha2_par_p2 = alpha2_par * alpha2_par;
   
    double q_vec_unit[3],n_perp_vec_unit[3];
    cross3(h_vec_unit,e_vec_unit,q_vec_unit);
    
    double e1_vec_unit[3],e2_vec_unit[3];
    double infall_direction_vec[3] = {1.0,0.0,0.0};    
    double Delta_vec[3],Delta_plus_vec[3],Delta_minus_vec[3];
    double S_vec[3];
    for (i=0; i<3; i++)
    {
        e1_vec_unit[i] = e_vec_unit[i];
        e2_vec_unit[i] = q_vec_unit[i]; // duplicate for notational/conceptional reasons
        n_perp_vec_unit[i] = q_vec_unit[i]; // duplicate for notational/conceptional reasons
        Delta_vec[i] = M * ( alpha1_vec[i] * m1 + alpha2_vec[i] * m2 );
        S_vec[i] = alpha1_vec[i] * m1*m1 + alpha2_vec[i] * m2*m2;

        Delta_plus_vec[i] = M * ( alpha2_vec[i] * m2 + alpha1_vec[i] * m1 ); // duplicate for notational/conceptional reasons
        Delta_minus_vec[i] = M * ( alpha2_vec[i] * m2 - alpha1_vec[i] * m1 );
    }
    
    double Theta_Delta = get_mutual_angle(Delta_vec,infall_direction_vec);
    double Theta_S = get_mutual_angle(S_vec,infall_direction_vec);
    double Theta_plus = get_mutual_angle(Delta_plus_vec,infall_direction_vec);
    double Theta_minus = get_mutual_angle(Delta_minus_vec,infall_direction_vec);
    
    double Theta_0 = 0.0;
    double Theta_1 = ((double) rand() / (RAND_MAX)) * 2.0 * M_PI;
    double Theta_2 = 0.0;
    double Theta_3 = ((double) rand() / (RAND_MAX)) * 2.0 * M_PI;
    double Theta_4 = 0.0;
    double Theta_5 = ((double) rand() / (RAND_MAX)) * 2.0 * M_PI;
    
    double zeta = 145.0 * M_PI/180.0;


    /* Recoil velocity */
    /* Unit for constants below where applicable: km/s */
    double A = 1.2e4;
    double B = -0.93;
    double H = 6.9e3;
    double K = 6.2e4;
    double K_S = -0.056;
    double B_K = 0.0;
    double B_H = 0.0;
    double H_S = 0.0;

    double v_m = A * eta_p2 * ((1.0 - q) / (1.0 + q)) * (1.0 + B * eta);
    double v_perp = H * (eta_p2 / (1.0 + q) ) * ( (1.0 + B_H * eta) * (alpha2_par - q * alpha1_par) + H_S * ((1.0 - q)/( one_plus_q_p2)) * (alpha2_par + q_p2 * alpha1_par) );
    double v_par = K * (eta_p2 / (1.0 + q) ) * ( (1.0 + B_K * eta) * fabs(alpha2_perp - q * alpha1_perp) * cos(Theta_Delta - Theta_0) \
        + K_S * ( (1.0-q)/( one_plus_q_p2 ) ) * fabs(alpha2_perp + q_p2 * alpha1_perp) * cos(Theta_S - Theta_1) );
    
    for (i=0; i<3; i++)
    {
        v_recoil_vec[i] = v_m * e1_vec_unit[i] + v_perp * (cos(zeta) * e1_vec_unit[i] + sin(zeta) * e2_vec_unit[i]) + v_par * h_vec_unit[i];
        v_recoil_vec[i] *= CONST_KM_PER_S; /* convert km/s to AU/yr */
    }


    /* New mass */
    double E_2 = 0.341; // ### pm 0.014
    double E_3 = 0.522; // ### pm 0.062
    double E_S = 0.673; // ### pm 0.035
    double E_Delta = -0.3689; // ### pm 0.37
    double E_A = -0.0136; // ### pm 0.021
    double E_B = 0.045; // ### pm 0.010
    double E_C = 0.0;
    double E_D = 0.2611; // ### pm 0.44
    double E_E = 0.09594; // ### pm 0.00045
    double E_F = 0.0;

    double E_ISCO = (1.0 - sqrt(8.0)/3.0) + 0.103803 * eta \
        + ( 1.0/(36.0*sqrt(3.0) * one_plus_q_p2) ) * ( q * (1.0 + 2.0 * q) * alpha1_par + (2.0 + q) * alpha2_par ) \
        - ( 5.0 / (324.0*sqrt(2.0) * one_plus_q_p2) ) * ( dot3(alpha2_vec,alpha2_vec) - 3.0 * alpha2_par_p2 - 2.0 * q * ( dot3(alpha1_vec,alpha2_vec) - 3.0 * alpha1_par * alpha2_par ) + q_p2 * (dot3(alpha1_vec,alpha1_vec) - 3.0 * alpha1_par_p2 ) );

    double delta_M_div_M = eta * E_ISCO + E_2 * eta_p2 + E_3 * eta_p3 \
        + (eta_p2/( one_plus_q_p2) ) * ( E_S * (alpha2_par + q_p2 * alpha1_par) + E_Delta * (1.0 - q)* (alpha2_par - q * alpha1_par) + E_A * dot3( alpha2_vec_plus_q_times_alpha1_vec, alpha2_vec_plus_q_times_alpha1_vec ) \
            + E_B * pow(alpha2_perp + q * alpha1_perp,2.0) * ( pow(cos(Theta_plus - Theta_2),2.0) + E_C) + E_D * dot3(alpha2_vec_minus_q_times_alpha1_vec, alpha2_vec_minus_q_times_alpha1_vec) \
            + E_E * pow(alpha2_perp - q * alpha1_perp,2.0) * ( pow(cos(Theta_minus - Theta_3),2.0) + E_F ) );
    
    *M_final = M - delta_M_div_M * M;


    /* New spin */
    double J_2 = -2.81; // ### pm 0.11
    double J_3 = 1.69; // ### pm 0.51
    double J_A = -2.9667; // ### pm 0.26
    double J_B = -1.7296; // ### pm 0.80
    double J_Delta = 0.0;
    double J_M = 0.0;
    double J_S = 0.0;
    double J_MS = 0.0;
    double J_MDelta = 0.0;

    double J_ISCO_vec[3];
    for (i=0; i<3; i++)
    {
        J_ISCO_vec[i] = ( 2.0*sqrt(3.0) - 1.5255862 * eta - (1.0/(9.0*sqrt(2.0)*one_plus_q_p2)) * ( q*(7.0 + 8.0*q) * alpha1_par + (8.0 + 7.0 * q) * alpha2_par ) \
            + (2.0/(9.0*sqrt(3.0)*one_plus_q_p2)) * ( dot3(alpha2_vec,alpha2_vec) - 3.0*alpha2_par_p2 - 2.0 * q * ( dot3(alpha1_vec,alpha2_vec) - 3.0 * alpha1_par * alpha2_par) + q_p2 * (dot3(alpha1_vec,alpha1_vec) - 3.0 * alpha1_par_p2) ) ) * h_vec_unit[i] \
        - (1.0/(9.0*sqrt(2.0)*one_plus_q_p2)) * ( q*(1.0+4.0*q) * alpha1_vec[i] + (4.0 + q) * alpha2_vec[i] ) + (1.0/eta) * (alpha2_vec[i] + q_p2 * alpha1_vec[i])/( one_plus_q_p2 );
        
        alpha_vec_final[i] = pow( 1.0 - delta_M_div_M,-2.0) * ( eta * J_ISCO_vec[i] + (J_2 * eta_p2 + J_3 * eta_p3) * h_vec_unit[i] \
        + (eta_p2/( one_plus_q_p2 )) * ( (J_A*(alpha2_par + q_p2 * alpha1_par) + J_B*(1.0-q)*(alpha2_par - q*alpha1_par) )*h_vec_unit[i] \
        + (1.0-q)*norm3(alpha2_perp_vec_minus_q_times_alpha1_perp_vec) * sqrt( J_Delta * cos( 2.0*(Theta_Delta - Theta_4) ) + J_MDelta ) * n_perp_vec_unit[i] \
        + norm3(alpha2_perp_vec_plus_q_p2_times_alpha1_perp_vec) * sqrt( J_S * cos( 2.0*(Theta_S - Theta_5) ) + J_MS) * n_perp_vec_unit[i] ) );
    }
}

double determine_effective_radius_for_collision(double radius, int stellar_type, int integration_flag)
{
    double effective_radius = radius;
    if (stellar_type < 10 and integration_flag > 0)
    {
        effective_radius = radius * effective_radius_multiplication_factor_for_collisions_stars;
    }
    else if (stellar_type >= 10)
    {
        effective_radius = radius * effective_radius_multiplication_factor_for_collisions_compact_objects;
    }
    
    #ifdef VERBOSE
    if (verbose_flag > 1)
    {
        printf("collision.cpp -- determine_effective_radius_for_collision R %g st %d int flag %d R_eff %g\n",radius, stellar_type,integration_flag,effective_radius);
    }
    #endif
    
    return effective_radius;
}

}
