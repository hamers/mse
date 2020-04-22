/* MSE */
/* Adrian Hamers November 2019 */

#include "evolve.h"
#include "merger.h"

extern "C"
{

void handle_collisions(ParticlesMap *particlesMap, int *integration_flag)
{
    /* invoke CE if one star has k in 2,3,4,5,6,8,9 */
    
    ParticlesMapIterator it_p;
    for (it_p = particlesMap->begin(); it_p != particlesMap->end(); it_p++)
    {
        Particle *p = (*it_p).second;
        if (p->is_binary == true and p->merged == true)
        {
            Particle *child1 = (*particlesMap)[p->child1];
            Particle *child2 = (*particlesMap)[p->child2];

            int kw = p->stellar_type;
            if (kw >= 2 and kw <= 9 and kw != 7)
            {
                common_envelope_evolution(particlesMap, p->index, child1->index, child2->index, integration_flag);
            }
            else
            {
                collision_product(particlesMap, p->index, child1->index, child2->index, integration_flag);
            }
        }
    }
}


void collision_product(ParticlesMap *particlesMap, int binary_index, int child1_index, int child2_index, int *integration_flag)
{
    Particle *b = (*particlesMap)[binary_index];
    Particle *child1 = (*particlesMap)[child1_index];
    Particle *child2 = (*particlesMap)[child2_index];
    
    double n_old = 2.0*M_PI/compute_orbital_period(b); /* mean motion just prior to collision */
    
    double m1 = child1->mass;
    double m2 = child2->mass;
    double m = m1 + m2;
            
    int kw = determine_merger_type(child1->stellar_type,child2->stellar_type);

    int kw1 = child1->stellar_type;
    int kw2 = child2->stellar_type;

    double age1 = child1->age*yr_to_Myr;
    double age2 = child2->age*yr_to_Myr;

    double tm1 = child1->sse_main_sequence_timescale*yr_to_Myr;
    double tm2 = child2->sse_main_sequence_timescale*yr_to_Myr;
    
    double m0,mc,age,t_MS;
    age = -1;
    m0 = -1;
    double epoch = 0.0; /* Is this OK? */
    
    double tm,tn;

    double *GB,*tscls,*lums;
    GB = new double[10];
    tscls = new double[20];
    lums = new double[10];    
    double *zpars;
    zpars = child1->zpars; /* for now: take metallicity of star 1; should think of Z-mixing in future */
  
    bool destroyed = false;
    bool apply_kick = false;
    int kick_distribution = child1->kick_distribution; /* by default, adopt kick distribution of child1 */
    
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
    
    /* Compact object mergers */
    else if (kw1 >= 13 and kw1 <= 14 and kw2 >= 10 and kw2 <= 14)
    {
        kw = kw1;
        m0 = m;
        age = 0.0;
        if (kw == 13 and m > 1.8)
        {
            kw = 14;
        }
        if (kw1 == 13 and kw2 == 13)
        {
            apply_kick = true;
            kick_distribution = 1; /* TO DO: implement realistic kick distribution for merged NSs */
        }
    }
    else if (kw2 >= 13 and kw2 <= 14 and kw1 >= 10 and kw1 <= 14)
    {
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
        b->apply_kick = false; /* assume no kicks following destruction events */
    }

    /* Handle effect of fast mass loss and kicks on the rest of the system */
    b->mass = m1 + m2; // the old total mass
    b->instantaneous_perturbation_delta_mass = m - b->mass; /* new mass minus old one */
    b->apply_kick = apply_kick;
    b->kick_distribution = kick_distribution;

    bool unbound_orbits;
    handle_SNe_in_system(particlesMap, &unbound_orbits, integration_flag);
    
    /* Update the merged object */
    if (destroyed == false)
    {
        b->is_binary = false; /* The binary becomes a body */
        particlesMap->erase(child1->index);
        particlesMap->erase(child2->index);

        star_(&kw, &m0, &m, &tm, &tn, tscls, lums, GB, zpars);
        double r,lum,rc,menv,renv,k2;        
        hrdiag_(&m0,&age,&m,&tm,&tn,tscls,lums,GB,zpars, \
            &r,&lum,&kw,&mc,&rc,&menv,&renv,&k2);

        b->stellar_type = kw;
        b->mass = m;
        b->sse_initial_mass = m0;
        
        b->epoch = epoch*Myr_to_yr;
        b->age = age*Myr_to_yr;
        b->sse_main_sequence_timescale = tm*Myr_to_yr;

        b->radius = r*CONST_R_SUN;
        b->luminosity = lum*CONST_L_SUN;
        b->core_mass = mc;
        b->core_radius = rc*CONST_R_SUN;
        b->convective_envelope_mass = menv;
        b->convective_envelope_radius = renv*CONST_R_SUN;

        /* Set the spin equal to the orbital frequency just before collision -- TO DO: better solution.
         * Assume the direction is equal to the previous orbital orientation. */
        double spin_vec_unit[3];
        get_unit_vector(b->h_vec,spin_vec_unit);
        for (int i=0; i<3; i++)
        {
            b->spin_vec[i] = n_old * spin_vec_unit[i];
        }
        
        b->apply_kick = true; /* Default value for (stellar evolution) bodies */
    }
    else
    {
        /* The binary b has been completely destroyed. 
         * If it has a parent, the latter needs to be replaced
         * from a binary to a body existing of the sibling of b. */

        particlesMap->erase(child1->index);
        particlesMap->erase(child2->index);
        
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
    
}

int determine_merger_type(int kw1, int kw2)
{
    return MERGER_TABLE[kw1][kw2];
}



}
