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

            int kw1 = child1->stellar_type;
            int kw2 = child2->stellar_type;
            if (kw1 >= 2 and kw1 <= 9 and kw1 != 7) /* CE in case of collision of giant + any other star */
            {
                common_envelope_evolution(particlesMap, p->index, child1->index, child2->index, integration_flag);
            }
            else if (kw2 >= 2 and kw2 <= 9 and kw2 != 7) /* CE in case of collision of giant + any other star */
            {
                common_envelope_evolution(particlesMap, p->index, child2->index, child1->index, integration_flag);
            }
            else /* "true" collision in other cases */
            {
                collision_product(particlesMap, p->index, child1->index, child2->index, integration_flag);
            }
        }
    }

    update_structure(particlesMap);
}


void apply_instantaneous_mass_changes_and_kicks(ParticlesMap *particlesMap, int *integration_flag)
{

    apply_user_specified_instantaneous_perturbation(particlesMap);
    reset_instantaneous_perturbation_quantities(particlesMap);
    
    bool unbound_orbits = check_for_unbound_orbits(particlesMap);
    if (unbound_orbits == true)
    {
        *integration_flag = 5;
        printf("merger.cpp -- apply_instantaneous_mass_changes_and_kicks -- Unbound orbits in system; switching to integration_flag %d\n",*integration_flag);
    }
}

void collision_product(ParticlesMap *particlesMap, int binary_index, int child1_index, int child2_index, int *integration_flag)
{
    printf("CP\n");
    set_up_derived_ODE_quantities(particlesMap); /* for setting a, e, etc. */

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

    //bool unbound_orbits;
    //handle_SNe_in_system(particlesMap, &unbound_orbits, integration_flag);
    apply_instantaneous_mass_changes_and_kicks(particlesMap, integration_flag);

    printf("CP -- destroyed %d\n",destroyed);
    
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


void determine_compact_object_merger_properties(double m1, double m2, double chi1, double chi2, double spin_vec_1_unit[3], double spin_vec_2_unit[3], double h_vec_unit[3], double e_vec_unit[3], double v_recoil_vec[3])
{
    /* Assumes m1 <= m2 */
    
    double alpha1_vec[3],alpha2_vec[3];
    int i;
    for (i=0; i<3; i++)
    {
        alpha1_vec[i] = chi1 * spin_vec_1_unit[i];
        alpha2_vec[i] = chi2 * spin_vec_2_unit[i];
    }

    double alpha1_par,alpha2_par,alpha1_perp,alpha2_perp;
    double alpha1_perp_vec[3],alpha2_perp_vec[3];
    get_parallel_and_perpendicular_vectors_and_components(alpha1_vec, h_vec_unit, &alpha1_par, &alpha1_perp, alpha1_perp_vec);
    get_parallel_and_perpendicular_vectors_and_components(alpha2_vec, h_vec_unit, &alpha2_par, &alpha2_perp, alpha2_perp_vec);
    
    double q = m1/m2; /* q <= 1 */
    double q_p2 = q*q;
    double one_plus_q = 1.0 + q;
    double one_plus_q_p2 = one_plus_q * one_plus_q;
    
    double M = m1 + m2;
    double eta = q/( one_plus_q_p2 );
    double eta_p2 = eta*eta;
    
    double q_vec_unit[3];
    cross3(h_vec_unit,e_vec_unit,q_vec_unit);
    
    double e1_vec_unit[3],e2_vec_unit[3];
    double infall_direction_vec[3] = {1.0,0.0,0.0};    
    double Delta_vec[3];
    double S_vec[3];
    for (i=0; i<3; i++)
    {
        e1_vec_unit[i] = e_vec_unit[i];
        e2_vec_unit[i] = q_vec_unit[i];
        Delta_vec[i] = M * ( alpha1_vec[i] * m1 + alpha2_vec[i] * m2 );
        S_vec[i] = alpha1_vec[i] * m1*m1 + alpha2_vec[i] * m2*m2;
    }
    
    
    double Theta_Delta = get_mutual_angle(Delta_vec,infall_direction_vec);
    double Theta_S = get_mutual_angle(S_vec,infall_direction_vec);
    
    double Theta_0 = 0.0;
    double Theta_1 = ((double) rand() / (RAND_MAX)) * 2.0 * M_PI;
    double zeta = 145.0 * M_PI/180.0;

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
}



#ifdef IGNORE

def determine_recoil_velocity(CONST_G,CONST_C,m1,m2,chi1,chi2,spin_vec_1_unit,spin_vec_2_unit,h_vec_unit,e_vec_unit):
    ### assumes m1 <= m2 ###
    alpha1_vec = chi1 * spin_vec_1_unit
    alpha2_vec = chi2 * spin_vec_2_unit
    alpha1_par, alpha1_perp = get_parallel_and_perpendicular_components(alpha1_vec,h_vec_unit)
    alpha2_par, alpha2_perp = get_parallel_and_perpendicular_components(alpha2_vec,h_vec_unit)
    
    q = m1/m2
    M = m1 + m2
    eta = q/( (1.0+q)**2 )
    
    q_vec_unit = np.cross(h_vec_unit,e_vec_unit)
    
    e1_vec_unit = e_vec_unit
    e2_vec_unit = q_vec_unit

    infall_direction_vec = np.array( [1.0,0.0,0.0] )

    Delta_vec = M * ( alpha1_vec * m1 + alpha2_vec * m2 )
    Theta_Delta = get_mutual_angle(Delta_vec,infall_direction_vec)
    
    S_vec = alpha1_vec * m1**2 + alpha2_vec * m2**2
    Theta_S = get_mutual_angle(S_vec,infall_direction_vec)
    
    Theta_0 = 0.0
    #Theta_1 = 0.0 ### check
    Theta_1 = np.random.random() * 2.0 * np.pi
    zeta = 145.0 * np.pi/180.0

    A = 1.2e4
    B = -0.93
    H = 6.9e3
    K = 6.2e4
    K_S = -0.056
    B_K = 0.0
    B_H = 0.0
    H_S = 0.0

    v_m = A * eta**2 * ((1.0 - q) / (1.0 + q)) * (1.0 + B * eta)
    v_perp = H * (eta**2 / (1.0 + q) ) * ( (1.0 + B_H * eta) * (alpha2_par - q * alpha1_par) + H_S * ((1.0 - q)/( (1.0 + q)**2)) * (alpha2_par + q**2 * alpha1_par) )
    v_par = K * (eta**2 / (1.0 + q) ) * ( (1.0 + B_K * eta) * np.fabs(alpha2_perp - q * alpha1_perp) * np.cos(Theta_Delta - Theta_0) \
        + K_S * ( (1.0-q)/( (1.0 + q)**2 ) ) * np.fabs(alpha2_perp + q**2 * alpha1_perp) * np.cos(Theta_S - Theta_1) )
    
    v_recoil_vec = v_m * e1_vec_unit + v_perp * (np.cos(zeta) * e1_vec_unit + np.sin(zeta) * e2_vec_unit) + v_par * h_vec_unit

    return v_recoil_vec

def determine_fractional_remnant_mass(CONST_G,CONST_C,m1,m2,chi1,chi2,spin_vec_1_unit,spin_vec_2_unit,h_vec_unit,e_vec_unit):
    ### assumes m1 <= m2 ###
    q = m1/m2
    M = m1 + m2
    eta = q/( (1.0+q)**2 )

    alpha1_vec = chi1 * spin_vec_1_unit
    alpha2_vec = chi2 * spin_vec_2_unit
    alpha1_par, alpha1_perp = get_parallel_and_perpendicular_components(alpha1_vec,h_vec_unit)
    alpha2_par, alpha2_perp = get_parallel_and_perpendicular_components(alpha2_vec,h_vec_unit)
    
    E_2 = 0.341 ### pm 0.014
    E_3 = 0.522 ### pm 0.062
    E_S = 0.673 ### pm 0.035
    E_Delta = -0.3689 ### pm 0.37
    E_A = -0.0136 ### pm 0.021
    E_B = 0.045 ### pm 0.010
    E_C = 0.0
    E_D = 0.2611 ### pm 0.44
    E_E = 0.09594 ### pm 0.00045
    E_F = 0.0
    
    Theta_2 = 0.0
    Theta_3 = np.random.random() * 2.0 * np.pi
    
    Delta_plus_vec = M * ( alpha2_vec * m2 + alpha1_vec * m1 )
    Delta_minus_vec = M * ( alpha2_vec * m2 - alpha1_vec * m1 )
    
    infall_direction_vec = np.array( [1.0,0.0,0.0] )
    Theta_plus = get_mutual_angle(Delta_plus_vec,infall_direction_vec)
    Theta_minus = get_mutual_angle(Delta_minus_vec,infall_direction_vec)
    
    E_ISCO = (1.0 - np.sqrt(8.0)/3.0) + 0.103803 * eta \
        + ( 1.0/(36.0*np.sqrt(3.0) * (1.0 + q)**2) ) * ( q * (1.0 + 2.0 * q) * alpha1_par + (2.0 + q) * alpha2_par ) \
        - ( 5.0 / (324.0*np.sqrt(2.0) * (1.0 + q)**2) ) * ( np.dot(alpha2_vec,alpha2_vec) - 3.0 * alpha2_par**2 - 2.0 * q * ( np.dot(alpha1_vec,alpha2_vec) - 3.0 * alpha1_par * alpha2_par ) + q**2 * (np.dot(alpha1_vec,alpha1_vec) - 3.0 * alpha1_par**2 ) )

    delta_M_div_M = eta * E_ISCO + E_2 * eta**2 + E_3 * eta**3 \
        + (eta**2/( (1.0 + q)**2) ) * ( E_S * (alpha2_par + q**2 * alpha1_par) + E_Delta * (1.0 - q)* (alpha2_par - q * alpha1_par) + E_A * np.dot( alpha2_vec + q * alpha1_vec, alpha2_vec + q * alpha1_vec ) \
            + E_B * (alpha2_perp + q * alpha1_perp)**2 * ( np.cos(Theta_plus - Theta_2)**2 + E_C) + E_D * np.dot(alpha2_vec - q * alpha1_vec, alpha2_vec - q * alpha1_vec) \
            + E_E * (alpha2_perp - q * alpha1_perp)**2 * ( np.cos(Theta_minus - Theta_3)**2 + E_F ) )
    
    #print("delta_M_div_M",Theta_plus,Theta_minus)
    #exit(0)
    return delta_M_div_M

def determine_remnant_mass_and_spin(CONST_G,CONST_C,m1,m2,chi1,chi2,spin_vec_1_unit,spin_vec_2_unit,h_vec_unit,e_vec_unit):
    ### assumes m1 <= m2 ###
    delta_M_div_M = determine_fractional_remnant_mass(CONST_G,CONST_C,m1,m2,chi1,chi2,spin_vec_1_unit,spin_vec_2_unit,h_vec_unit,e_vec_unit)

    q = m1/m2
    M = m1 + m2
    eta = q/( (1.0+q)**2 )
        
    alpha1_vec = chi1 * spin_vec_1_unit
    alpha2_vec = chi2 * spin_vec_2_unit
    alpha1_par, alpha1_perp = get_parallel_and_perpendicular_components(alpha1_vec,h_vec_unit)
    alpha2_par, alpha2_perp = get_parallel_and_perpendicular_components(alpha2_vec,h_vec_unit)
    
    alpha1_perp_vec = alpha1_vec - np.dot(alpha1_vec,h_vec_unit) * h_vec_unit
    alpha2_perp_vec = alpha2_vec - np.dot(alpha2_vec,h_vec_unit) * h_vec_unit
    
    q_vec_unit = np.cross(h_vec_unit,e_vec_unit)
    n_perp_vec_unit = q_vec_unit
    
    J_2 = -2.81 ### pm 0.11
    J_3 = 1.69 ### pm 0.51
    J_A = -2.9667 ### pm 0.26
    J_B = -1.7296 ### pm 0.80
    J_Delta = 0.0
    J_M = 0.0
    J_S = 0.0
    J_MS = 0.0
    J_MDelta = 0.0
    Theta_0 = 0.0
    Theta_2 = 0.0
    Theta_4 = 0.0
    Theta_5 = np.random.random() * 2.0 * np.pi

    infall_direction_vec = np.array( [1.0,0.0,0.0] )

    Delta_vec = M * ( alpha1_vec * m1 + alpha2_vec * m2 )
    Theta_Delta = get_mutual_angle(Delta_vec,infall_direction_vec)

    S_vec = alpha1_vec * m1**2 + alpha2_vec * m2**2
    Theta_S = get_mutual_angle(S_vec,infall_direction_vec)
    
    J_ISCO_vec = ( 2.0*np.sqrt(3.0) - 1.5255862 * eta - (1.0/(9.0*np.sqrt(2.0)*(1.0+q)**2)) * ( q*(7.0 + 8.0*q) * alpha1_par + (8.0 + 7.0 * q) * alpha2_par ) \
            + (2.0/(9.0*np.sqrt(3.0)*(1.0 + q)**2)) * ( np.dot(alpha2_vec,alpha2_vec) - 3.0*alpha2_par**2 - 2.0 * q * ( np.dot(alpha1_vec,alpha2_vec) - 3.0 * alpha1_par * alpha2_par) + q**2 * (np.dot(alpha1_vec,alpha1_vec) - 3.0 * alpha1_par**2) ) ) * h_vec_unit \
        - (1.0/(9.0*np.sqrt(2.0)*(1.0+q)**2)) * ( q*(1.0+4.0*q) * alpha1_vec + (4.0 + q) * alpha2_vec ) + (1.0/eta) * (alpha2_vec + q**2 * alpha1_vec)/( (1.0+q)**2 )
    alpha_vec_final = pow( 1.0 - delta_M_div_M,-2.0) * ( eta * J_ISCO_vec + (J_2 * eta**2 + J_3 * eta**3) * h_vec_unit \
        + (eta**2/( (1.0+q)**2 )) * ( (J_A*(alpha2_par + q**2 * alpha1_par) + J_B*(1.0-q)*(alpha2_par - q*alpha1_par) )*h_vec_unit \
        + (1.0-q)*np.linalg.norm(alpha2_perp_vec - q*alpha1_perp_vec) * np.sqrt( J_Delta * np.cos( 2.0*(Theta_Delta - Theta_4) ) + J_MDelta ) * n_perp_vec_unit \
        + np.linalg.norm(alpha2_perp_vec + q**2 * alpha1_perp_vec) * np.sqrt( J_S * np.cos( 2.0*(Theta_S - Theta_5) ) + J_MS) * n_perp_vec_unit ) )
    
    M_final = M - delta_M_div_M * M
    #print("M_init",M,"M_final",M_final,delta_M_div_M)
    return M_final, alpha_vec_final

#endif


}
