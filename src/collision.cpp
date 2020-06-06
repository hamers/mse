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
                    common_envelope_evolution(particlesMap, p->index, child1->index, child2->index, t, integration_flag);
                }
                else if (kw2 >= 2 and kw2 <= 9 and kw2 != 7) /* CE in case of collision of giant + any other star */
                {
                    common_envelope_evolution(particlesMap, p->index, child2->index, child1->index, t, integration_flag);
                }
                else /* "pure" collision in other cases */
                {
                    collision_product(particlesMap, p->index, child1->index, child2->index, t, integration_flag);
                }
            }
        }
        update_structure(particlesMap);
        
        bool stable = check_system_for_dynamical_stability(particlesMap, integration_flag);
        if (stable == false)
        {
            *integration_flag = 1;
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
            col_part_i = pi->Collision_Partner;
            if (col_part_i != -1)
            {
                for (it_j = particlesMap->begin(); it_j != particlesMap->end(); it_j++)
                {
                    pj = (*it_j).second;
                    col_part_j = pj->Collision_Partner;
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
            exit(-1);
        }
        /* Reset Collision_Partner */
        pi->Collision_Partner = -1;
        pj->Collision_Partner = -1;
        
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
        
}


void collision_product(ParticlesMap *particlesMap, int binary_index, int child1_index, int child2_index, double t, int *integration_flag)
{
    /* TO DO: allow for calling this function without specifying the binary_index, i.e., collision during N-body integration (integration_flag>0) */
    
    printf("CP\n");
    print_system(particlesMap,*integration_flag);

    Particle *child1 = (*particlesMap)[child1_index];
    Particle *child2 = (*particlesMap)[child2_index];
    
    double m1 = child1->mass;
    double m2 = child2->mass;
    double m = m1 + m2;

    int i;

    Particle *b;
    double n_old;
    double initial_momentum[3],initial_R_CM[3];
    
    double h_vec_unit[3],e_vec_unit[3];
    if (*integration_flag == 0)
    {
        set_up_derived_quantities(particlesMap); /* for setting a, e, etc. */
        b = (*particlesMap)[binary_index];

        n_old = 2.0*M_PI/compute_orbital_period(b); /* mean motion just prior to collision */
        get_unit_vector(b->h_vec,h_vec_unit);
        get_unit_vector(b->e_vec,e_vec_unit);
    }
    else
    {
        double h_vec[3],e_vec[3],r[3],v[3];
        double true_anomaly;
        
        for (i=0; i<3; i++)
        {
            r[i] = child1->R_vec[i] - child2->R_vec[i];
            v[i] = child1->V_vec[i] - child2->V_vec[i];
            initial_R_CM[i] = (m1 * child1->R_vec[i] + m2 * child2->R_vec[i]) / (m1 + m2);
            initial_momentum[i] = m1 * child1->V_vec[i] + m2 * child2->V_vec[i];
        }
        from_cartesian_to_orbital_vectors(m1,m2,r,v,e_vec,h_vec,&true_anomaly);
        get_unit_vector(h_vec,h_vec_unit);
        get_unit_vector(e_vec,e_vec_unit);
    }
    
    int kw = determine_merger_type(child1->stellar_type,child2->stellar_type);

    int kw1 = child1->stellar_type;
    int kw2 = child2->stellar_type;

    double age1 = child1->age*yr_to_Myr;
    double age2 = child2->age*yr_to_Myr;

    double tm1 = child1->sse_main_sequence_timescale*yr_to_Myr;
    double tm2 = child2->sse_main_sequence_timescale*yr_to_Myr;

    double chi1 = child1->chi;
    double chi2 = child2->chi;
    
    double *spin_vec_1_unit = child1->spin_vec_unit;
    double *spin_vec_2_unit = child2->spin_vec_unit;

    double m0,mc,age,t_MS;
    age = -1;
    m0 = -1;
    
    double tm,tn;

    double *GB,*tscls,*lums;
    GB = new double[10];
    tscls = new double[20];
    lums = new double[10];    
    double *zpars;
    zpars = child1->zpars; /* for now: take metallicity of star 1; should think of Z-mixing in future */
  
    bool destroyed = false;
    //bool apply_kick = false;
    bool reset_spin_vec = true;
    
    int kick_distribution = child1->kick_distribution; /* by default, adopt kick distribution of child1 */


    double v_kick_vec[3] = {0.0,0.0,0.0};
    double spin_vec[3];
    double M_final; /* used for compact objects */
    double alpha_vec_final[3];


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
        determine_compact_object_merger_properties(m1,m2,chi1,chi2,spin_vec_1_unit,spin_vec_2_unit,b->h_vec_unit,b->e_vec_unit,v_recoil_vec,alpha_vec_final,&M_final);
        m = M_final;
        //apply_kick = true;
                
        reset_spin_vec = false;
        double chi = norm3(alpha_vec_final);
        double Omega = compute_spin_frequency_from_spin_parameter(m,chi);
        for (i=0; i<3; i++)
        {
            v_kick_vec[i] = v_recoil_vec[i];
            spin_vec[i] = Omega * alpha_vec_final[i] / chi;
        }
        
        kw = kw1;
        m0 = m;
        age = 0.0;
        if (kw == 13 and m > 1.8)
        {
            kw = 14;
        }
        //if (kw1 == 13 and kw2 == 13)
        //{
            //apply_kick = true;
            //kick_distribution = 1; /* TO DO: implement realistic kick distribution for merged NSs */
        //}
    }
    else if (kw2 >= 13 and kw2 <= 14 and kw1 >= 10 and kw1 <= 14)
    {
        double v_recoil_vec[3];
        determine_compact_object_merger_properties(m2,m1,chi2,chi1,spin_vec_2_unit,spin_vec_1_unit,b->h_vec_unit,b->e_vec_unit,v_recoil_vec,alpha_vec_final,&M_final);
        m = M_final;
        //apply_kick = true;
                
        reset_spin_vec = false;
        double chi = norm3(alpha_vec_final);
        double Omega = compute_spin_frequency_from_spin_parameter(m,chi);
        for (i=0; i<3; i++)
        {
            v_kick_vec[i] = v_recoil_vec[i];
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
        return;
    }
    
    if (destroyed == false and (m0 == -1 or age == -1))
    {
        printf("merger.cpp -- collision_product -- was not able to determine all properties of merged object! Will ignore the collision. \n");
        return;
    }

    
    /* m, m0, age */
    
    
    if (destroyed == true)
    {
        m = 0.0; /* need to set for instantaneous_perturbation_delta_mass below */
        //b->apply_kick = false; /* assume no kicks following destruction events */
    }

    //b->apply_kick = apply_kick;
    //b->kick_distribution = kick_distribution;
    //bool unbound_orbits;
    //handle_SNe_in_system(particlesMap, &unbound_orbits, integration_flag);


    double r,lum,rc,menv,renv,k2;
    double z_new;
    double *zpars_new;
    zpars_new = new double[20];

    if (destroyed == false)
    {

        /* stellar properties of merged object */
        star_(&kw, &m0, &m, &tm, &tn, tscls, lums, GB, zpars);
        hrdiag_(&m0,&age,&m,&tm,&tn,tscls,lums,GB,zpars, \
            &r,&lum,&kw,&mc,&rc,&menv,&renv,&k2);

        z_new = child1->metallicity;

        zcnsts_(&z_new,zpars_new);
    }


    if (*integration_flag == 0) /* secular mode */
    {
        
        b->is_binary = false; /* The binary becomes a body (or destroyed) */
        particlesMap->erase(child1->index);
        particlesMap->erase(child2->index);

        /* Handle effect of fast mass loss and kicks on orbits in the rest of the system */
        b->mass = m1 + m2; // the old total mass
        b->instantaneous_perturbation_delta_mass = m - b->mass; /* new mass minus old one; note: if destroyed, the new mass will be 0 */
        b->instantaneous_perturbation_delta_VX = v_kick_vec[0];
        b->instantaneous_perturbation_delta_VY = v_kick_vec[1];
        b->instantaneous_perturbation_delta_VZ = v_kick_vec[2];

        printf("CP -- delta m %g v_kick_vec %g %g %g\n",b->instantaneous_perturbation_delta_mass,v_kick_vec[0],v_kick_vec[1],v_kick_vec[2]);
        apply_instantaneous_mass_changes_and_kicks(particlesMap, integration_flag); /* this will update the binary's mass */
        

        printf("CP -- destroyed %d\n",destroyed);
        
        /* Update the merged object */
        if (destroyed == false)
        {

            b->stellar_type = kw;
            b->mass = m;
            b->sse_initial_mass = m0;
            
            b->age = age*Myr_to_yr;
            b->epoch = t - b->age;
            b->sse_main_sequence_timescale = tm*Myr_to_yr;

            b->radius = r*CONST_R_SUN;
            b->luminosity = lum*CONST_L_SUN;
            b->core_mass = mc;
            b->core_radius = rc*CONST_R_SUN;
            b->convective_envelope_mass = menv;
            b->convective_envelope_radius = renv*CONST_R_SUN;

            b->metallicity = z_new;
            b->zpars = zpars_new;

            /* Unless calculated above, set the spin equal to the orbital frequency just before collision
             * Assume the direction is equal to the previous orbital orientation. */
            if (reset_spin_vec == true)
            {
                for (i=0; i<3; i++)
                {
                    b->spin_vec[i] = n_old * h_vec_unit[i];
                }
            }
            else
            {
                for (i=0; i<3; i++)
                {
                    b->spin_vec[i] = spin_vec[i];
                }
            }
            b->apply_kick = true; /* Default value for (stellar evolution) bodies */
            b->merged = false;
        }
        else
        {
            /* The binary b has been completely destroyed. 
             * If it has a parent, the latter needs to be replaced
             * from a binary to a body existing of the sibling of b. */

            //particlesMap->erase(child1->index);
            //particlesMap->erase(child2->index);
            
            handle_destruction_of_binary_in_system(particlesMap,b);
        }
    }
    else /* direct N-body mode */
    {
        if (destroyed == false)
        {
            /* child1 will become the merged object */
            particlesMap->erase(child2->index);
            child1->stellar_type = kw;
            child1->mass = m;
            child1->sse_initial_mass = m0;
            
//            child1->epoch = epoch*Myr_to_yr;
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
            }

            /* Unless calculated above, set the spin equal to the orbital frequency just before collision
             * Assume the direction is equal to the previous orbital orientation. */
            if (reset_spin_vec == true)
            {
                for (i=0; i<3; i++)
                {
                    child1->spin_vec[i] = n_old * h_vec_unit[i];
                }
            }
            else
            {
                for (i=0; i<3; i++)
                {
                    child1->spin_vec[i] = spin_vec[i];
                }
            }
            child1->apply_kick = true; /* Default value for (stellar evolution) bodies */
            child1->RLOF_flag = 0;
            
            printf("CP N-body not destroyed post\n");
            print_system(particlesMap,*integration_flag);
        }
        else
        {
            particlesMap->erase(child1->index);
            particlesMap->erase(child2->index);
        }
    }
    
}

void handle_destruction_of_binary_in_system(ParticlesMap *particlesMap, Particle *b)
{
    if (b->parent == -1) /* b was a lone binary; the new system is empty */
    {
        particlesMap->erase(b->index);
    }
    else /* the system survives in some form */
    {
        Particle *p = (*particlesMap)[b->parent];
        Particle *s = (*particlesMap)[b->sibling];
        p->is_binary = false;
        
        if (p->parent == -1) /* the sibling s becomes the top-level binary; p can be safely removed */
        {
            particlesMap->erase(p->index);
        }
        else /* particle p has become obsolete and needs to be replaced by s */
        {
            Particle *pp = (*particlesMap)[p->parent];
            if (p->index == pp->child1)
            {
                pp->child1 = s->index;
            }
            else if (p->index == pp->child2)
            {
                pp->child2 = s->index;
            }
            particlesMap->erase(p->index);
        }
    }
}


int determine_merger_type(int kw1, int kw2)
{
    return MERGER_TABLE[kw1][kw2];
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

}
