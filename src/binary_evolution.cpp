/* MSE */

#include "evolve.h"
#include "binary_evolution.h"

extern "C"
{

int handle_binary_evolution(ParticlesMap *particlesMap, double t_old, double t, double *dt_binary_evolution, int *integration_flag)
{
    *dt_binary_evolution = 1.0e100;
    
    if (*integration_flag != 0) /* Only do binary evolution if we can be sure the system is dynamically stable and we have well-defined orbits. */
    {
        return 0;
    }

    update_structure(particlesMap, *integration_flag);
    set_positions_and_velocities(particlesMap);
    handle_wind_accretion(particlesMap,t_old,t,dt_binary_evolution,integration_flag);
    int binary_flag = handle_mass_transfer(particlesMap,t_old,t,dt_binary_evolution,integration_flag);
    handle_mass_accretion_events_with_degenerate_objects(particlesMap,t_old,t,integration_flag,dt_binary_evolution);

    bool stable = check_system_for_dynamical_stability(particlesMap, integration_flag);
    if (stable == false)
    {
        *integration_flag = 1;
    }

    #ifdef VERBOSE
    if (verbose_flag > 1)
    {
        printf("binary_evolution.cpp -- handle_binary_evolution -- done; stable %d integration_flag %d dt_binary_evolution %g\n",stable,*integration_flag,*dt_binary_evolution);
        print_system(particlesMap,*integration_flag);
    }
    #endif
    
    return binary_flag;
}


int handle_wind_accretion(ParticlesMap *particlesMap, double t_old, double t, double *dt_binary_evolution, int *integration_flag)
{
    #ifdef VERBOSE
    if (verbose_flag > 2)
    {
        printf("binary_evolution.cpp -- handle_wind_accretion\n");
    }
    #endif
    
    set_up_derived_quantities(particlesMap); /* for setting a, e, etc. */
    double v_orb_p2,v_wind_p2,factor;

    /* First, reset mass_dot_wind_accretion for all particles */ 
    ParticlesMapIterator it_p;
    for (it_p = particlesMap->begin(); it_p != particlesMap->end(); it_p++)
    {
        Particle *p = (*it_p).second;
        p->mass_dot_wind_accretion = 0.0;
    }

    for (it_p = particlesMap->begin(); it_p != particlesMap->end(); it_p++)
    {
        
        Particle *p = (*it_p).second;
        
        if (p->is_binary == false and p->is_bound == true)
        {
            /* Wind accretion from particle p onto its companion */
            Particle *parent = (*particlesMap)[p->parent];
            Particle *companion = (*particlesMap)[p->sibling];

            v_orb_p2 = CONST_G * parent->mass / parent->a;
            v_wind_p2 = 2.0 * beta_wind_accretion * CONST_G * p->mass / p->radius; /* squared wind speed from p */
            
            factor = (1.0/parent->j) * pow( CONST_G * companion->mass / v_wind_p2, 2.0) * (alpha_wind_accretion / (2.0 * parent->a*parent->a)) * pow(1.0 + v_orb_p2/v_wind_p2, -1.5);
            companion->mass_dot_wind_accretion = CV_min(1.0, factor) * (- p->mass_dot_wind); /* sanity check (necessary for eccentric orbits): ensure that the companion cannot accrete more than the wind loss from p */

            if (companion->mass_dot_wind_accretion!=companion->mass_dot_wind_accretion)
            {
                printf("binary_evolution.cpp -- handle_wind_accretion -- companion->mass_dot_wind_accretion %g\n",companion->mass_dot_wind_accretion);
                print_system(particlesMap,*integration_flag);
                printf("p->mass %g parent->mass %g p->radius %g v_orb_p2 %g v_wind_p2 %g parent->j %g parent->a %g\n",p->mass,parent->mass,p->radius,v_orb_p2,v_wind_p2,parent->j,parent->h_vec[0]);
                error_code = 2;
                longjmp(jump_buf,1);
            }
        }
    }
    
    return 0;
}


int handle_mass_transfer(ParticlesMap *particlesMap, double t_old, double t, double *dt_binary_evolution, int *integration_flag)
{
    *dt_binary_evolution = 1.0e100;
    double dt = t - t_old; /* Imposed timestep which will be used for secular ODE integration. */

    double epsilon_MT = 0.01;
    double dt_temp;

    #ifdef VERBOSE
    if (verbose_flag > 1)
    {
        printf("binary_evolution.cpp -- handle_mass_transfer -- start\n");
        print_system(particlesMap,*integration_flag);
    }
    #endif

    /* First, determine which stars/orbits are undergoing RLOF */
    std::vector<int> parent_indices,donor_indices,accretor_indices;
    ParticlesMapIterator it_p;
    for (it_p = particlesMap->begin(); it_p != particlesMap->end(); it_p++)
    {
        Particle *donor = (*it_p).second;

        /* Reset some quantities (some of which used in ODE_handle_NS_properties_Ye19_model()) -- could be updated below */
        donor->delta_mass_RLOF = 0.0;
        donor->RLOF_timescale = 0.0;

        donor->mass_dot_RLOF = 0.0;
        donor->mass_dot_RLOF_triple = 0.0;
        donor->mass_dot_adiabatic_ejection = 0.0;
        
        /* Check donor and accretor conditions for RLOF */
        if (donor->is_binary == false and donor->is_bound == true)
        {
            Particle *accretor = (*particlesMap)[donor->sibling];
            if (donor->object_type == 1 and accretor->object_type == 1) /* only allow star-star RLOF */
            {
                if (donor->RLOF_flag == 1 and std::count(parent_indices.begin(), parent_indices.end(), donor->parent) == 0)
                {
                    parent_indices.push_back(donor->parent);
                    donor_indices.push_back(donor->index);
                    accretor_indices.push_back(donor->sibling);
                }
            }
        }
    }

    /* Handle RLOF cases */
    std::vector<int>::iterator it;
    bool reevaluate_integration_flag = false;
    int flag = 0;
    for (int i=0; i<parent_indices.size(); i++)
    {
        flag = handle_mass_transfer_cases(particlesMap, parent_indices[i], donor_indices[i], accretor_indices[i], integration_flag, t_old, t, dt_binary_evolution);
        if (flag != 5)
        {
            reevaluate_integration_flag = true;
        }
    }

    if (reevaluate_integration_flag == true)
    {
        *integration_flag = determine_orbits_in_system_using_nbody(particlesMap);
    }

    update_structure(particlesMap, *integration_flag);

    return flag;
    
}

double compute_q_crit_for_common_envelope_evolution(int kw, double mass, double core_mass)
{
    double q_crit;
    
    if (kw == 2)
    {
        q_crit = 4.0;
    }
    else if (kw == 3 or kw == 5 or kw == 6) // Hjellming & Webbink, 1987, ApJ, 318, 794. 
    {
        q_crit = 0.362 + 1.0/(3.0 * (1.0 - core_mass/mass));
    }
    else if (kw == 8 or kw == 9)
    {
        q_crit = 0.784;
    }
    else
    {
        q_crit = 3.0;
    }
    
    return q_crit;
}


int handle_mass_transfer_cases(ParticlesMap *particlesMap, int parent_index, int donor_index, int accretor_index, int *integration_flag, double t_old, double t, double *dt_binary_evolution)//, ParticlesMapIterator &it_p)
{
    Particle *parent = (*particlesMap)[parent_index];
    Particle *donor = (*particlesMap)[donor_index];
    Particle *accretor = (*particlesMap)[accretor_index];

    #ifdef VERBOSE
    if (verbose_flag > 0)
    {
        printf("binary_evolution.cpp -- handle_mass_transfer_cases -- parent_index %d donor_index %d accretor_index %d donor->RLOF_flag %d accretor->RLOF_flag %d\n",parent_index,donor_index,accretor_index,donor->RLOF_flag,accretor->RLOF_flag);
        print_system(particlesMap,*integration_flag);
    }
    #endif
 
    if (*integration_flag != 0)
    {
        return 0;
    }
    
    int flag = 0;
    
    double P_orb = compute_orbital_period_from_semimajor_axis(parent->mass,parent->a);
    double M_donor = donor->mass;
    double M_accretor = accretor->mass;
    double R_accretor = accretor->radius;
    int kw = donor->stellar_type;
    
    /* MT & dynamical timescales */
    double R_donor = donor->radius;
    double t_dyn_donor = compute_stellar_dynamical_timescale(M_donor, R_donor);

    /* Critical mass ratios */
    double q = M_donor / M_accretor;
    double q_crit = compute_q_crit_for_common_envelope_evolution(kw, M_donor, donor->core_mass);
    double q_crit_low_mass_donor = 0.695;
    double q_crit_high_mass_donor = 3.0;
    double q_crit_WD_donor = 0.628;
    
    /* Prepare some quantities used for the eCAML model from Schreiber et al. (2016, MNRAS, 455, L16) */
    double nudyn = 0.35/M_donor; 
    double zeta_RL = (2.0/3.0) * ( log(1.0+pow(q,1.0/3.0)) - (0.5*pow(q,1.0/3.0))/(1.0+pow(q,1.0/3.0)) )/(0.6*pow(q,2.0/3.0)+log(1.0+pow(q,1.0/3.0))) + 2.0*nudyn + M_donor/(M_accretor+M_donor) - 2.0; // Roche lobe Zeta (Eq. A11 in Belloni et al. 2018, MNRAS, 478, 5639)

    double zeta_AD; // Adiabatic Zeta (Eq. A13 in Belloni et al. 2018, MNRAS, 478, 5639)
    if (M_donor <= 0.38412)
    {
        zeta_AD = -1.0/3.0;
    }
    else
    {
        zeta_AD = 0.782491 - 7.46384*M_donor + 13.9255*pow(M_donor, 2.0) - 5.3892*pow(M_donor, 3.0);
    }

    /* Estimate of mass transfer timescale */
    double a = parent->a;
    double e = parent->e;
    double R_Lc = roche_radius_pericenter_eggleton(a,q); /* with argument "a", actually computes circular Roche lobe radius */
    double x = R_Lc/R_donor;
    double E_0;
    bool in_RLOF;
    determine_E_0(e, x, &E_0, &in_RLOF);
    double fm = fm_function(e,x,E_0,0.0); /* set E_tau = 0 here */

    double m_dot = compute_orbit_averaged_mass_transfer_rate_emt_model(M_donor,fm,P_orb);
    double fabs_m_dot = fabs(m_dot);
    //double m_dot = compute_bse_mass_transfer_amount(donor->stellar_type, M_donor, donor->core_mass, donor->radius, double R_RL_av_donor, double dt, double t_dyn_donor, double t_KH_donor)
    int kw2 = accretor->stellar_type;

    if (accretor->RLOF_flag == 1) /* Contact evolution */
    {
        if ((kw >= 2 and kw <= 9 and kw != 7) or (kw2 >= 2 and kw2 <= 9 and kw2 != 7)) /* CE if either donor or accretor are giants */
        {
            #ifdef VERBOSE
            if (verbose_flag > 0)
            {
                printf("binary_evolution.cpp -- handle_mass_transfer_cases -- Contact evolution -- kw2 %d -- invoking CE \n",kw2);
                print_system(particlesMap,*integration_flag);
            }
            #endif

            flag = 6;
            if (kw >= 2 and kw <= 9 and kw != 7) /* Take donor to be CE donor, irrespective of accretor */
            {
                common_envelope_evolution(particlesMap, parent->index, donor->index, accretor->index, t_old, integration_flag);
            }
            else /* Accretor must be CE donor */
            {
                common_envelope_evolution(particlesMap, parent->index, accretor->index, donor->index, t_old, integration_flag);
            }
            *dt_binary_evolution = ODE_min_dt;
            return flag;
        }
        else /* Otherwise, let the stars merge */
        {
            #ifdef VERBOSE
            if (verbose_flag > 0)
            {
                printf("binary_evolution.cpp -- handle_mass_transfer_cases -- Contact evolution -- kw2 %d -- merging \n",kw2);
                print_system(particlesMap,*integration_flag);
            }
            #endif

            flag = 7;
            collision_product(particlesMap, parent_index, donor_index, accretor_index, t_old, integration_flag);
            *dt_binary_evolution = ODE_min_dt;
            return flag;
        }
    }
    
    double rp = a*(1.0 - e);
    if (rp < (R_donor + R_accretor)) /* Periapsis distance close enough to cause immediate merger rather than RLOF, so invoke coalescence (can happen in some cases, e.g., after SNe kicks) */
    {
        if ((kw >= 2 and kw <= 9 and kw != 7) or (kw2 >= 2 and kw2 <= 9 and kw2 != 7)) /* CE if either donor or accretor are giants */
        {
            #ifdef VERBOSE
            if (verbose_flag > 0)
            {
                printf("binary_evolution.cpp -- handle_mass_transfer_cases -- Immediate merger; CE evolution -- rp %g R_d %g R_a %g \n",rp,R_donor,R_accretor);
                print_system(particlesMap,*integration_flag);
            }
            #endif

            flag = 6;
            common_envelope_evolution(particlesMap, parent->index, donor->index, accretor->index, t_old, integration_flag);
            *dt_binary_evolution = ODE_min_dt;
            return flag;
        }
        else /* Otherwise, let the stars merge */
        {
            #ifdef VERBOSE
            if (verbose_flag > 0)
            {
                printf("binary_evolution.cpp -- handle_mass_transfer_cases -- Immediate merger; collision -- rp %g R_d %g R_a %g \n",rp,R_donor,R_accretor);
                print_system(particlesMap,*integration_flag);
            }
            #endif

            flag = 7;
            collision_product(particlesMap, parent_index, donor_index, accretor_index, t_old, integration_flag);
            *dt_binary_evolution = ODE_min_dt;
            return flag;
        }
    }
    
    double t_MT;
    if (fabs_m_dot <= epsilon)
    {
        t_MT = 1.0e100;
    }
    else
    {
        t_MT = M_donor/fabs_m_dot;
    }

    if (t_MT != t_MT)
    {
        printf("binary_evolution.cpp -- handle_mass_transfer_cases -- ERROR: t_MT = %g; a %g e %g R_Lc %g x %g E_0 %g \n",t_MT,a,e,R_Lc,x,E_0);
        //exit(-1);
        error_code = 3;
        longjmp(jump_buf,1);
    }
    
    #ifdef VERBOSE
    if (verbose_flag > 1)
    {
        printf("binary_evolution.cpp -- handle_mass_transfer_cases -- t_MT %g emt_fm %g\n",t_MT,fm);
        print_system(particlesMap,*integration_flag);
    }
    #endif

    /* Different cases */
    
    /* Low mass donor criterion is user-adjustable -- default or eCAML */
    bool low_mass_donor_criterion;
    if (binary_evolution_use_eCAML_model == false)
    {
        low_mass_donor_criterion = (kw == 0 and q > q_crit_low_mass_donor);
    }
    else
    {
        low_mass_donor_criterion = (kw <= 1 and ( ( M_donor <= 0.8 and zeta_RL > zeta_AD ) or ( M_donor > 0.8 and q > q_crit_high_mass_donor ) ));
    }
    
    if (low_mass_donor_criterion)
    {
        /* `CE'-like evolution from a low-mass MS star to any secondary. */
        flag = 1;
        dynamical_mass_transfer_low_mass_donor(particlesMap, parent->index, donor->index, accretor->index, t_old, t, integration_flag);
        *dt_binary_evolution = ODE_min_dt;
    }
    else if ( (kw >= 2 and kw <= 9 and kw != 7) and ((t_MT <= t_dyn_donor or q > q_crit) or t_MT < P_orb) ) /* include criterion with donor convective envelope mass? */
    {
        /* `Standard CE evolution. */
        flag = 2;
        common_envelope_evolution(particlesMap, parent->index, donor->index, accretor->index, t_old, integration_flag);
        *dt_binary_evolution = ODE_min_dt;
    }
    else if (kw >= 10 and kw <= 12 and q > q_crit_WD_donor)
    {
        /* Dynamical transfer from WD */
        flag = 3;
        dynamical_mass_transfer_WD_donor(particlesMap, parent->index, donor->index, accretor->index, t_old, t, integration_flag);
        *dt_binary_evolution = ODE_min_dt;
    }   
    else if (kw == 13 or kw == 14)
    {
        /* Dynamical transfer from NS */
        flag = 4;
        mass_transfer_NS_BH_donor(particlesMap, parent->index, donor->index, accretor->index, t_old, t, integration_flag);
        *dt_binary_evolution = ODE_min_dt;
    }
    else
    {
        /* Other cases: stable mass transfer. */
        flag = 5;
        stable_mass_transfer_evolution(particlesMap, parent->index, donor->index, accretor->index, t_old, t, integration_flag, dt_binary_evolution);
    }
   

    #ifdef VERBOSE
    if (verbose_flag > 0)
    {
        printf("binary_evolution.cpp -- handle_mass_transfer_cases -- flag %d integration_flag %d dt_binary_evolution %g t %g\n",flag,*integration_flag,*dt_binary_evolution,t);
        print_system(particlesMap,*integration_flag);
    }
    #endif

    
    if (flag == -1)
    {
        printf("binary_evolution.cpp -- fatal error in handle_mass_transfer_cases, parent %d donor %d accretor %d \n",parent->index,donor->index,accretor->index);
        error_code = 3;
        longjmp(jump_buf,1);
        //exit(-1);
    }
    
    return flag;
}

double compute_orbit_averaged_mass_transfer_rate_emt_model(double M_donor, double fm, double P_orb)
{
    return -(M_donor/P_orb) * fm;
}

double compute_stellar_dynamical_timescale(double M, double R)
{
    return sqrt( ( R*R*R )/(CONST_G * M) );
}


int dynamical_mass_transfer_low_mass_donor(ParticlesMap *particlesMap, int parent_index, int donor_index, int accretor_index, double t_old, double t, int *integration_flag)
{
    /* Low-mass dynamical mass transfer case similar to CE (BSE evolv2.f: lines 1150 - 1233).
     * This type of (fast) evolution is not handled by the ODE integrator. 
     * Always leads to the destruction of the low-mass MS donor; the accretor always survives. 
     * Do not allow for this case when the accretor is a binary (this would be unlikely in the first place since low-mass  MS stars have small radii). */

    #ifdef VERBOSE
    if (verbose_flag > 0)
    {
        printf("binary_evolution.cpp -- dynamical_mass_transfer_low_mass_donor\n");
        print_system(particlesMap,*integration_flag);
    }
    #endif
    
    Particle *parent = (*particlesMap)[parent_index];
    Particle *donor = (*particlesMap)[donor_index];
    Particle *accretor = (*particlesMap)[accretor_index];

    if (accretor->is_binary == true)
    {
        printf("binary_evolution.cpp -- dynamical_mass_transfer_low_mass_donor -- accretor is binary; skipping. \n");
        return 0;
    }

    double P_orb = compute_orbital_period_from_semimajor_axis(parent->mass,parent->a);
    double m_donor = donor->mass;
    double m_accretor = accretor->mass;
    double m_accretor_new;
    
    double R_donor = donor->radius;
    double R_accretor = accretor->radius;

    double m_old = m_donor + m_accretor;

    double dm1; /* Mass lost from companion due to RLOF (>0) */
    double dm2; /* Mass accreted by companion (>0). */
    
    double dt = t - t_old;
    
    int kw1 = donor->stellar_type;
    int kw2 = accretor->stellar_type;
    double t_KH_donor = compute_Kelvin_Helmholtz_timescale(kw1,m_donor,donor->core_mass,R_donor,donor->luminosity);
    double t_KH_accretor = compute_Kelvin_Helmholtz_timescale(kw2,m_accretor,accretor->core_mass,R_accretor,accretor->luminosity);

    double t_dyn_donor = compute_stellar_dynamical_timescale(m_donor, R_donor);
    double tau_donor = sqrt(t_KH_donor * t_dyn_donor);
    
    dm1 = m_donor;

    /* Run star_() on accretor */
    double *GB,*tscls,*tscls_new,*lums,*zpars;
    GB = new double[10];
    tscls = new double[20];
    tscls_new = new double[20];
    lums = new double[10];
    zpars = accretor->zpars;
    double tms,tn;
    star_(&kw2,&accretor->sse_initial_mass,&accretor->mass,&tms,&tn,tscls,lums,GB,zpars); /* note: time units from star_ are in Myr */
    
    double tms_new,tn_new;
        
    if (kw2 <= 1)
    {
        /* Restrict accretion to thermal timescale of secondary. */

        dm2 = (tau_donor/t_KH_accretor) * dm1;
        m_accretor_new = accretor->mass + dm2;

        /* Rejuvenate if the star is still on the main sequence. */
        accretor->sse_initial_mass = accretor->mass;
        star_(&kw2,&accretor->sse_initial_mass,&m_accretor_new,&tms_new,&tn_new,tscls_new,lums,GB,zpars); 

        /* If the star has no convective core then the effective age decreases,
         * otherwise it will become younger still. */
         
        if (m_accretor_new < 0.35 or m_accretor_new > 1.25)
        {
            accretor->age = (tms_new/tms) * accretor->age * (m_accretor_new - dm2)/m_accretor_new; 
        }
        else
        {
            accretor->age = (tms_new/tms) * accretor->age;
        }
        accretor->epoch = t - accretor->age;
    }
    else if (kw2 <= 6)
    {
        /* Add all the material to the giant's envelope. */

        dm2 = dm1;
        m_accretor_new = accretor->mass + dm2;
        
        if (kw2 == 2)
        {
            accretor->sse_initial_mass = accretor->mass;
            star_(&kw2,&accretor->sse_initial_mass,&m_accretor_new,&tms_new,&tn_new,tscls_new,lums,GB,zpars); /* will update tscls */
              
            accretor->age = tms_new*Myr_to_yr + (accretor->age - tms*Myr_to_yr) * (tscls_new[0] / tscls[0]);
            accretor->epoch = t - accretor->age;
        }
    }
    else if (kw2 <= 12)
    {
        /* Form a new giant envelope. */

        dm2 = dm1;
        m_accretor_new = accretor->mass + dm2;
        int kw = determine_merger_type(kw1,kw2);

        if (kw == 4)
        {
            //accretor->age = accretor->age / (tms*Myr_to_yr); /* ASH: this line copied from BSE does not make sense to me (unit-wise); I am ignoring it. */
            accretor->core_mass = accretor->mass;
        }
        
        /* Check for planets or low-mass WDs. */
        if ((kw2 == 10 and m_accretor < 0.05) or (kw2 >= 11 and m_accretor < 0.5))
        {
            kw = kw1;
        }
        else
        {
            double new_age;
            gntage_(&accretor->core_mass,&m_accretor_new,&kw,zpars,&accretor->sse_initial_mass,&new_age);
            accretor->age = new_age*Myr_to_yr;
            accretor->stellar_type = kw;
            accretor->epoch = t - accretor->age;
            
        }
        accretor->stellar_type = kw;
    }
    else
    {
        /* The neutron star or black hole simply accretes at (at most) the Eddington rate. */

        double hydrogen_mass_fraction = zpars[10];
        double m_dot_Eddington = compute_Eddington_accretion_rate(R_accretor, hydrogen_mass_fraction);
        double dme = m_dot_Eddington * dt;
       
        dm2 = CV_min(dme*tau_donor/dt, dm1);
        m_accretor_new = accretor->mass + dm2;
    }

    donor->apply_kick = false;
    accretor->apply_kick = false;
    
    double initial_momentum[3],initial_R_CM[3],initial_V_CM[3];
    double h_vec[3],h_vec_unit[3],e_vec[3],e_vec_unit[3],r_vec[3],v_vec[3];

    get_initial_binary_orbital_properties_from_position_and_velocity(donor->R_vec, donor->V_vec, accretor->R_vec, accretor->V_vec, m_donor, m_accretor, \
        r_vec, v_vec, initial_momentum, initial_R_CM, initial_V_CM, h_vec, e_vec);
    
    double final_R_CM[3],final_V_CM[3],final_momentum[3];
    handle_gradual_mass_loss_event_in_system(particlesMap, donor, accretor, 0.0, m_donor, m_accretor + dm2, m_accretor, parent->dynamical_mass_transfer_low_mass_donor_timescale, \
        r_vec, v_vec, initial_R_CM, initial_V_CM, final_R_CM, final_V_CM, final_momentum);

    update_stellar_evolution_properties(accretor);
    //accretor->RLOF_flag = 0;
    reset_ODE_mass_dot_quantities(accretor);
    particlesMap->erase(donor_index);
    
    *integration_flag = 1;


    #ifdef VERBOSE
    if (verbose_flag > 0)
    {
        printf("binary_evolution.cpp -- dynamical_mass_transfer_low_mass_donor -- end\n");
        print_system(particlesMap,*integration_flag);
    }
    #endif
    
    delete[] GB;
    delete[] tscls;
    delete[] tscls_new;
    delete[] lums;

    return 0;
}

double compute_Kelvin_Helmholtz_timescale(int kw, double mass, double core_mass, double radius, double luminosity)
{
    double t_KH = CONST_G * mass / ( 2.0 * radius * luminosity); /* per unit mass */
    if (kw <= 1 or kw == 7 or kw >= 10)
    {
        t_KH *= mass;
    }
    else
    {
        t_KH *= (mass - core_mass);
    }

    return t_KH;
}

double compute_Eddington_accretion_rate(double radius, double hydrogen_mass_fraction)
{
    double kappa = 0.2 * (1.0 + hydrogen_mass_fraction) * 8.88021e6; /* electron scattering opacity in units of AU^2/MSun */
    return eddington_accretion_factor * 4.0 * M_PI * CONST_C_LIGHT * radius/kappa;
}

int dynamical_mass_transfer_WD_donor(ParticlesMap *particlesMap, int parent_index, int donor_index, int accretor_index, double t_old, double t, int *integration_flag)
{
    /* Dynamical mass transfer case WD donor (BSE evolv2.f: lines 1290 - 1342).
     * This type of (fast) evolution is not handled by the ODE integrator. 
     * Always leads to the destruction of the donor; the accretor might survive. 
     * Take into account mass loss from the system, and assume there is no kick imparted on the accretor if it survives.
     * Do not allow for this case when the accretor is a binary (this would be unlikely in the first place given the size of WDs). */

    #ifdef VERBOSE
    if (verbose_flag > 0)
    {
        printf("binary_evolution.cpp -- dynamical_mass_transfer_WD_donor\n");
        print_system(particlesMap,*integration_flag);
    }
    #endif

    Particle *parent = (*particlesMap)[parent_index];
    Particle *donor = (*particlesMap)[donor_index];
    Particle *accretor = (*particlesMap)[accretor_index];

    if (accretor->is_binary == true)
    {
        printf("binary_evolution.cpp -- dynamical_mass_transfer_WD_donor -- accretor is binary; skipping. \n");
        return 0;
    }

    double P_orb = compute_orbital_period_from_semimajor_axis(parent->mass,parent->a);
    double m_donor = donor->mass;
    double m_accretor = accretor->mass;

    double m_old = m_donor + m_accretor;

    double dm1; /* Mass lost from companion due to RLOF (<0) */
    double dm2; /* Mass accreted by companion (>0).  */

    double dt = t - t_old;
    dm1 = m_donor;
    dm2 = dm1;

    bool destroyed;
    double v_kick_vec[3] = {0.0,0.0,0.0};
    
    int kw1 = donor->stellar_type;
    int kw2 = accretor->stellar_type;
    
    if (kw1 == 10 and kw2 == 10)
    {
        /* Assume the energy released by ignition of the triple-alpha reaction is enough to destroy both stars. */
        destroyed = true;
    }
    else if (kw2 >= 10 and kw2 <= 11 and (m_accretor + dm2) > chandrasekhar_mass)
    {
        /* Potentially SNe Ia that destroys the system. */
        destroyed = true;
    }
    else if (kw1 == 10 or kw2 == 10)
    {
        /* Should be helium overflowing onto a CO or ONe core in which case the 
        *  helium swells up to form a giant envelope so a HeGB star is formed. 
        *  Allowance for the rare case of CO or ONe flowing onto He is made. */

        int kst = 9;
        double accretor_core_mass = accretor->core_mass; 
        if (kw2 == 10)
        {
            accretor_core_mass = dm2;
        }
        double accretor_age = accretor->age*yr_to_Myr;
        gntage_(&accretor_core_mass,&m_accretor,&kst,accretor->zpars,&accretor->sse_initial_mass,&accretor_age);
        accretor->age = accretor_age * Myr_to_yr;
        accretor->stellar_type = kst;
        accretor->epoch = t - accretor_age * Myr_to_yr;
        
        destroyed = false;
        //m_new = m_accretor;
    }
    else if (kw2 <= 12)
    {
        accretor->sse_initial_mass = m_accretor;
        if (kw1 == 12 and kw2 == 11) /* Mixture of ONe and CO will result in an ONe product. */
        {
            accretor->stellar_type = 12;
        }
        destroyed = false;
        //m_new = m_accretor;        
    }
    else
    {
        //printf("binary_evolution.cpp -- dynamical_mass_transfer_WD_donor -- unknown case -- parent %d donor %d accretor %d\n",parent_index,donor_index,accretor_index);
        //error_code = 4;
        //longjmp(jump_buf,1);
        //exit(-1);
        collision_product(particlesMap, parent_index, donor_index, accretor_index, t_old, integration_flag);

        *integration_flag = determine_orbits_in_system_using_nbody(particlesMap);

        return 0;
    }

    double initial_momentum[3],initial_R_CM[3],initial_V_CM[3];
    double h_vec[3],h_vec_unit[3],e_vec[3],e_vec_unit[3],r_vec[3],v_vec[3];

    get_initial_binary_orbital_properties_from_position_and_velocity(donor->R_vec, donor->V_vec, accretor->R_vec, accretor->V_vec, m_donor, m_accretor, \
        r_vec, v_vec, initial_momentum, initial_R_CM, initial_V_CM, h_vec, e_vec);
    
    double final_R_CM[3],final_V_CM[3],final_momentum[3];
    handle_gradual_mass_loss_event_in_system(particlesMap, donor, accretor, 0.0, m_donor, m_accretor + dm2, m_accretor, parent->dynamical_mass_transfer_low_mass_donor_timescale, \
        r_vec, v_vec, initial_R_CM, initial_V_CM, final_R_CM, final_V_CM, final_momentum);

    if (destroyed == false)
    {
        update_stellar_evolution_properties(accretor);
        reset_ODE_mass_dot_quantities(accretor);
        
        particlesMap->erase(donor_index);
    }
    else
    {
        particlesMap->erase(donor_index);
        particlesMap->erase(accretor_index);
    }

    *integration_flag = 1;

    return 0;
    
}

int mass_transfer_NS_BH_donor(ParticlesMap *particlesMap, int parent_index, int donor_index, int accretor_index, double t_old, double t, int *integration_flag)
{
    /* "Mass transfer" case NS/BH donor (BSE evolv2.f: lines 1343 - 1365 ).
     * Treat this as a pure collision. 
     * Do not allow for this case when the accretor is a binary (this would be unlikely in the first place given the size of NNS/BHs). */

    #ifdef VERBOSE
    if (verbose_flag > 0)
    {
        printf("binary_evolution.cpp -- mass_transfer_NS_BH_donor\n");
        print_system(particlesMap,*integration_flag);
    }
    #endif

    Particle *parent = (*particlesMap)[parent_index];
    Particle *donor = (*particlesMap)[donor_index];
    Particle *accretor = (*particlesMap)[accretor_index];

    if (accretor->is_binary == true)
    {
        printf("binary_evolution.cpp -- mass_transfer_NS_BH_donor -- accretor is binary; skipping. \n");
        return 0;
    }

    collision_product(particlesMap, parent_index, donor_index, accretor_index, t_old, integration_flag);

    *integration_flag = determine_orbits_in_system_using_nbody(particlesMap);
    
    return 0;
}


int stable_mass_transfer_evolution(ParticlesMap *particlesMap, int parent_index, int donor_index, int accretor_index, double t_old, double t, int *integration_flag, double *dt_binary_evolution)
{
    Particle *donor = (*particlesMap)[donor_index];
    Particle *accretor = (*particlesMap)[accretor_index];
    
    if (accretor->is_binary == false)
    {
        /* Both objects are physical stars. Apply the "standard" mass transfer scheme from HTP02. */
        binary_stable_mass_transfer_evolution(particlesMap, parent_index, donor_index, accretor_index, t_old, t, integration_flag, dt_binary_evolution);
    }
    else
    {
        /* A tertiary star is overflowing onto a companion binary. */
        triple_stable_mass_transfer_evolution(particlesMap, parent_index, donor_index, accretor_index, t_old, t, integration_flag, dt_binary_evolution);
    }
    
    return 0;
}


int binary_stable_mass_transfer_evolution(ParticlesMap *particlesMap, int parent_index, int donor_index, int accretor_index, double t_old, double t, int *integration_flag, double *dt_binary_evolution)
{
    /* Stable mass transfer case (BSE evolv2.f: lines 1370 - 1905) */

    #ifdef VERBOSE
    if (verbose_flag > 0)
    {
        printf("binary_evolution.cpp -- binary_stable_mass_transfer_evolution parent_index %d donor_index %d accretor_index %d\n",parent_index,donor_index,accretor_index);
        print_system(particlesMap,*integration_flag);
    }
    #endif


    Particle *parent = (*particlesMap)[parent_index];
    Particle *donor = (*particlesMap)[donor_index];
    Particle *accretor = (*particlesMap)[accretor_index];

    double P_orb = compute_orbital_period_from_semimajor_axis(parent->mass,parent->a);
    double m_donor = donor->mass;
    double m_accretor = accretor->mass;

    double q = m_donor/m_accretor;

    double R_donor = donor->radius;
    double R_accretor = accretor->radius;

    double m_old = m_donor + m_accretor;

    int kw1 = donor->stellar_type;
    int kw2 = accretor->stellar_type;

    double t_KH_donor = compute_Kelvin_Helmholtz_timescale(kw1,m_donor,donor->core_mass,R_donor,donor->luminosity);
    double t_KH_accretor = compute_Kelvin_Helmholtz_timescale(kw2,m_accretor,accretor->core_mass,R_accretor,accretor->luminosity);

    double t_dyn_donor = compute_stellar_dynamical_timescale(m_donor, R_donor);
    double a = parent->a;
    double e = parent->e;
    double rp = a*(1.0-e);
    
    double dt = t - t_old;

    /* Default value of dm1 -- transferred mass during this time-step */
    double R_RL_av_donor = roche_radius_pericenter_eggleton(rp,q);
    //double dm1 = compute_bse_mass_transfer_amount(kw1,m_donor,donor->core_mass,R_donor,R_RL_av_donor,dt,t_dyn_donor,t_KH_donor);
    double dm1 = compute_bse_mass_transfer_amount_averaged(kw1, m_donor, donor->core_mass, R_donor, m_accretor, a, e, dt, t_dyn_donor, t_KH_donor);
    //dm1 *= (1.0 - e);
    double dm2; /* amount of mass gained by accretor during dt (defined dm2>0) */

    double dms_donor = fabs(donor->mass_dot_wind * dt); /* absolute value of net mass change due to winds */
    double dms_accretor = fabs(accretor->mass_dot_wind * dt); /* absolute value of net mass change due to winds */
     
    double tms,tn;
    double *GB,*tscls,*lums;
    GB = new double[10];
    tscls = new double[20];
    lums = new double[10];  

    /* These quantities are used below to determine aging */
    star_(&kw1,&donor->sse_initial_mass,&m_donor,&tms,&tn,tscls,lums,GB,donor->zpars);
    double tms_donor_old = tms;
    double tbgb_donor_old = tscls[0];
    
    star_(&kw2,&accretor->sse_initial_mass,&m_accretor,&tms,&tn,tscls,lums,GB,accretor->zpars);
    double tms_accretor_old = tms;
    double tbgb_accretor_old = tscls[0];

    bool nova = false;

    /* Determine dm2. */
    /* Eddington accretion rate for compact accretors. */
    double hydrogen_mass_fraction = donor->zpars[10];
    double m_dot_Eddington = compute_Eddington_accretion_rate(R_accretor, hydrogen_mass_fraction);
    double dme = m_dot_Eddington * dt;

    double taum = (m_accretor/dm1) * dt;
    bool gntage_was_called = false;
    
    if (kw2 <= 2 or kw2 == 4)
    {
        /* Limit according to the thermal timescale of the secondary. */

        dm2 = dm1 * CV_min(1.0, 10.0 * taum/t_KH_accretor);
    }
    else if (kw2 >= 7 and kw2 <= 9)
    {
        /* Naked helium star secondary swells up to a core helium burning star
         * or SAGB star unless the primary is also a helium star. */

        if (kw1 >= 7)
        {
            dm2 = dm1 * CV_min(1.0, 10.0 * taum/t_KH_accretor);
        }
        else
        {
            dm2 = dm1;
            
            double dmchk = dm2 - 1.05 * dms_accretor;
            if (dmchk > 0.0 and dm2/m_accretor > 1.0e-4)
            {
                int kst = CV_min(6, 2*kw2 - 10);
                double mcx;
                double accretor_age_Myr = accretor->age*yr_to_Myr;
                if (kst == 4)
                {
                    mcx = m_accretor;
                    accretor_age_Myr /= tms_accretor_old;
                }
                else
                {
                    mcx = accretor->core_mass;
                }
        
                double mt2 = m_accretor + (dm2 - dms_accretor);
                
                gntage_(&mcx, &mt2, &kst, accretor->zpars, &accretor->sse_initial_mass, &accretor_age_Myr); /* The age passed on to gntage should actually be a fractional age, fage, of the CHeB lifetime -- see gntage.f lines 101-102 */
                gntage_was_called = true;
                
                star_(&kst, &accretor->sse_initial_mass, &mt2, &tms, &tn, tscls, lums, GB, accretor->zpars);
                
                double r,lum,rc,menv,renv,k2;
                //printf("BIN 1 tms %g tn %g tscls[0] %g accretor_age_Myr %g accretor->sse_initial_mass %g kst %d\n",tms,tn,tscls[0],accretor_age_Myr,accretor->sse_initial_mass,kst);
                hrdiag_(&accretor->sse_initial_mass,&accretor_age_Myr,&mt2,&tms,&tn,tscls,lums,GB,accretor->zpars, \
                    &r,&lum,&kst,&mcx,&rc,&menv,&renv,&k2);

                accretor->stellar_type = kst; /* accretor becomes CHeB or AGB star; need to update stellar type */
                //printf("BIN 2 tms %g tn %g tscls[0] %g accretor_age_Myr %g accretor->sse_initial_mass %g kst %d\n",tms,tn,tscls[0],accretor_age_Myr,accretor->sse_initial_mass,kst);
                accretor->core_mass = mcx;

                accretor->luminosity = lum*CONST_L_SUN;
                accretor->radius = r*CONST_R_SUN;
                accretor->core_radius = rc*CONST_R_SUN;
                accretor->convective_envelope_mass = menv;
                accretor->convective_envelope_radius = renv*CONST_R_SUN;
                accretor->sse_k2 = k2;                

                accretor->sse_main_sequence_timescale = tms*Myr_to_yr;

                accretor->epoch = t - accretor_age_Myr*Myr_to_yr;
                accretor->age = accretor_age_Myr*Myr_to_yr;
                
                accretor->mass_dot_wind = 0.0;
                accretor->radius_dot = 0.0;
                accretor->ospin_dot = 0.0;
                
                accretor->age_dot = 0.0;
                accretor->sse_initial_mass_dot = 0.0;
                accretor->core_mass_dot = 0.0;
                accretor->core_radius_dot = 0.0;

                accretor->luminosity_dot = 0.0;
                accretor->convective_envelope_mass_dot = 0.0;
                accretor->convective_envelope_radius_dot = 0.0;
                accretor->sse_k2_dot = 0.0;
                
                check_number(accretor->radius,                   "binary_evolution.cpp -- binary_stable_mass_transfer_evolution -- Naked helium star secondary swells up to a core helium burning star","accretor->radius", true);
                check_number(accretor->core_radius,                   "binary_evolution.cpp -- binary_stable_mass_transfer_evolution -- Naked helium star secondary swells up to a core helium burning star","accretor->core_radius", true);
                check_number(accretor->age,                   "binary_evolution.cpp -- binary_stable_mass_transfer_evolution -- Naked helium star secondary swells up to a core helium burning star","accretor->age", true);

            }
        }
    }
    else if (kw1 <= 6 and kw2 >= 10 and kw2 <= 12)
    {
        /* White dwarf secondary. */
        /* Hydrogen-rich donor */

        /* The default boundary values are adopted from HTP02 */
        double m_dot_upper = 2.71e-7; /* Above this accretion rate, the WD will swell up to a giant star */
        double m_dot_lower = 1.03e-7; /* Below this accretion rate, nova outbursts will occur that blow away most of the accreted material */
        
        if (binary_evolution_SNe_Ia_single_degenerate_model == 1)
        {
            
            white_dwarf_hydrogen_accretion_boundaries_WBBP13(m_accretor, &m_dot_lower, &m_dot_upper);
        }

        if (dm1/dt < m_dot_upper) 
        {
            if (dm1/dt < m_dot_lower)
            {
                /* Accrete until a nova explosion blows away most of the accreted material. */
                nova = true;

                dm2 = nova_accretion_factor * CV_min(dm1, dme);

                #ifdef VERBOSE
                if (verbose_flag > 0)
                {
                    printf("binary_evolution.cpp -- binary_stable_mass_transfer_evolution -- WD accretion -- nova m_dot %g dm1 %g dm2 %g m_dot_lower %g m_dot_upper %g\n",dm1/dt,dm1,dm2,m_dot_lower,m_dot_upper);
                }
                #endif

            }
            else
            {
                if (binary_evolution_SNe_Ia_single_degenerate_model == 0)
                {
                    /* Steady burning at the surface (to C/O) */
                    dm2 = dm1;
                    
                }
                else if (binary_evolution_SNe_Ia_single_degenerate_model == 1)
                {
                    /* Assume accreted hydrogen is quickly burned to helium; the latter accumulates as an outer layer on to the WD (100% efficiency) */
                    /* Since SSE does not model any He layer on WDs, treat this separately; dm2 should be zero to avoid the mass being accumulated as C/O. */
                    /* WD_He_layer_mass and m_dot_accretion_SD will be used in handle_mass_accretion_events_with_degenerate_objects() to check for SNe */
                    accretor->WD_He_layer_mass += dm1;
                    accretor->m_dot_accretion_SD = dm1/dt;
                    dm2 = 0.0;
                    
                    #ifdef VERBOSE
                    if (verbose_flag > 0)
                    {
                        printf("binary_evolution.cpp -- binary_stable_mass_transfer_evolution -- WD accretion -- accumulation onto He layer H donor m_dot %g He layer mass %g m_dot_lower %g m_dot_upper %g m_dot_accretion_SD %g\n",dm1/dt,accretor->WD_He_layer_mass,m_dot_lower,m_dot_upper,accretor->m_dot_accretion_SD);
                    }
                    #endif
                    
                }

            }
                
        }
        else
        {
            /* Make a new giant envelope. */

            dm2 = dm1;
            
            #ifdef VERBOSE
            if (verbose_flag > 0)
            {
                printf("binary_evolution.cpp -- binary_stable_mass_transfer_evolution -- WD accretion -- mass transfer rate %g MSun/yr is high enough to make a giant star -- m_dot_lower %g m_dot_upper %g\n",dm1/dt,m_dot_lower,m_dot_upper);
            }
            #endif
            
            int kw;
            if ( (kw2 == 10 and m_accretor < 0.05) or (kw2 >= 11 and m_accretor < 0.5) ) /* Check for planets or low-mass WDs. */
            {
                kw = kw2;
            }
            else
            {
                kw = CV_min(6, 3*kw2 - 27);

                double accretor_age_Myr = accretor->age*yr_to_Myr;
                double mt2 = m_accretor + (dm2 - dms_accretor);
                gntage_(&accretor->core_mass,&mt2,&kw,accretor->zpars,&accretor->sse_initial_mass,&accretor_age_Myr);
                gntage_was_called = true;
                
                star_(&kw, &accretor->sse_initial_mass, &mt2, &tms, &tn, tscls, lums, GB, accretor->zpars);
                
                double r,lum,rc,menv,renv,k2;
                hrdiag_(&accretor->sse_initial_mass,&accretor_age_Myr,&mt2,&tms,&tn,tscls,lums,GB,accretor->zpars, \
                    &r,&lum,&kw,&accretor->core_mass,&rc,&menv,&renv,&k2);

                accretor->stellar_type = kw;

                accretor->luminosity = lum*CONST_L_SUN;
                accretor->radius = r*CONST_R_SUN;
                accretor->core_radius = rc*CONST_R_SUN;
                accretor->convective_envelope_mass = menv;
                accretor->convective_envelope_radius = renv*CONST_R_SUN;
                accretor->sse_k2 = k2;                

                accretor->sse_main_sequence_timescale = tms*Myr_to_yr;

                accretor->epoch = t - accretor_age_Myr*Myr_to_yr;
                accretor->age = accretor_age_Myr*Myr_to_yr;
                
                accretor->mass_dot_wind = 0.0;
                accretor->radius_dot = 0.0;
                accretor->ospin_dot = 0.0;
                
                accretor->age_dot = 0.0;
                accretor->sse_initial_mass_dot = 0.0;
                accretor->core_mass_dot = 0.0;
                accretor->core_radius_dot = 0.0;

                accretor->luminosity_dot = 0.0;
                accretor->convective_envelope_mass_dot = 0.0;
                accretor->convective_envelope_radius_dot = 0.0;
                accretor->sse_k2_dot = 0.0;

                
                check_number(accretor->radius,                   "binary_evolution.cpp -- binary_stable_mass_transfer_evolution -- WD secondary -- make new giant envelope","accretor->radius", true);
                check_number(accretor->core_radius,                   "binary_evolution.cpp -- binary_stable_mass_transfer_evolution -- WD secondary -- make new giant envelope","accretor->core_radius", true);
                check_number(accretor->age,                   "binary_evolution.cpp -- binary_stable_mass_transfer_evolution -- WD secondary -- make new giant envelope","accretor->age", true);
            }
        }
    }
    else if (kw2 >= 10)
    {
        dm2 = CV_min(dm1, dme);
    }
    else
    {
        /* We have a giant whose envelope can absorb any transferred material. */
        dm2 = dm1;
    }
    
    if (binary_evolution_SNe_Ia_single_degenerate_model == 1)
    {
        /* Model for He accretion onto CO WDs *
         * Based on Neunteufel+ 2016, 2016A&A...589A..43N */
        if (kw1 >= 7 and kw1 <= 10 and kw2 == 11)
        {
            double eta;
            int WD_accretion_mode = -1;
            double m_dot_mass_transfer_SD = dm1/dt; // always positive
            white_dwarf_helium_mass_accumulation_efficiency(m_accretor, m_dot_mass_transfer_SD, accretor->luminosity, &eta, &WD_accretion_mode);
            double m_dot_accretion_SD = eta * m_dot_mass_transfer_SD; // always positive
            accretor->m_dot_accretion_SD = m_dot_accretion_SD;
            
            if (WD_accretion_mode == 1) // mass gets added in the form of a He layer on top of the CO core (tracked with the parameter "WD_He_layer_mass"); do not change the SSE CO WD mass
            {
                accretor->WD_He_layer_mass += eta * dm1;
                dm2 = 0.0;
            }
            if (WD_accretion_mode == 2) // no retention (novae)
            {
                dm2 = 0.0;
            }
            if (WD_accretion_mode == 3) // burning of accreted He to C; add to SSE CO WD mass (only represents CO core mass)
            {
                dm2 = eta * dm1;
            }
        }
    }
   
    
    
    #ifdef VERBOSE
    if (verbose_flag > 0)
    {
        printf("binary_evolution.cpp -- binary_stable_mass_transfer_evolution parent_index %d donor_index %d accretor_index %d dm1 %g dm2 %g\n",parent_index,donor_index,accretor_index,dm1,dm2);
        print_system(particlesMap,*integration_flag);
    }
    #endif
    
    /* "New" masses used for aging. Note: the actual masses (and radii) will be updated during the ODE integration.
     * Also, by definition, dm1 > 0, dm2 > 0 */
    double m_donor_new = m_donor - dm1 + donor->mass_dot_wind * dt;
    double m_accretor_new = m_accretor + dm2 + accretor->mass_dot_wind * dt;

    if (kw1 <= 1 or kw1 == 7)
    {
        donor->sse_initial_mass = m_donor_new;
    }
    if (kw2 <= 1 or kw2 == 7 and gntage_was_called == false)
    {
        accretor->sse_initial_mass = m_accretor_new;
    }
    
    /* For a HG star check if the initial mass can be reduced. */    
    if (kw1 == 2 and donor->sse_initial_mass <= donor->zpars[2])
    {
        double m0 = donor->sse_initial_mass;
        donor->sse_initial_mass = m_donor;

        star_(&kw1,&donor->sse_initial_mass,&donor->mass,&tms,&tn,tscls,lums,GB,donor->zpars);
        if (GB[8] < donor->core_mass)
        {
            donor->sse_initial_mass = m0;
        }
    }
    if (kw2 == 2 and accretor->sse_initial_mass <= accretor->zpars[2]) /* the same but with donor and accretor swapped */
    {
        double m0 = accretor->sse_initial_mass;
        accretor->sse_initial_mass = m_accretor;

        star_(&kw2,&accretor->sse_initial_mass,&accretor->mass,&tms,&tn,tscls,lums,GB,accretor->zpars);
        if (GB[8] < accretor->core_mass)
        {
            accretor->sse_initial_mass = m0;
        }
    }

    /* Always rejuvenate the secondary and age the primary if they are on the main sequence. */
    if (kw1 <= 2 or kw1 == 7)
    {
        star_(&kw1,&donor->sse_initial_mass,&m_donor_new,&tms,&tn,tscls,lums,GB,donor->zpars);

        if (kw1 == 2)
        {
            double age = tms + (donor->age*yr_to_Myr - tms_donor_old) * (tscls[0] - tms)/(tbgb_donor_old - tms_donor_old);
            donor->age = age*Myr_to_yr;
        }
        else
        {
            donor->age *= (tms/tms_donor_old);
            //printf("age %g %g %.15f m_donor_new %g m_donor %g\n",tms,tms_donor_old,tms/tms_donor_old,m_donor_new,m_donor);
        }
        donor->epoch = t - donor->age;
    }
    if (kw2 <= 2 or kw2 == 7)
    {
        star_(&kw2,&accretor->sse_initial_mass,&m_accretor_new,&tms,&tn,tscls,lums,GB,accretor->zpars);
        if (kw2 == 2)
        {
            double age = tms + (accretor->age*yr_to_Myr - tms_accretor_old) * (tscls[0] - tms)/(tbgb_accretor_old - tms_accretor_old);
            accretor->age = age*Myr_to_yr;
        }
        else if ((m_accretor_new < 0.35 or m_accretor_new > 1.25) and (kw2 != 7))
        {
            accretor->age *= (tms/tms_accretor_old) * ((m_accretor_new - dm2)/m_accretor_new);
        }
        else
        {
            if (gntage_was_called == false)
            {
                accretor->age *= (tms/tms_accretor_old);
            }
        }
        accretor->epoch = t - accretor->age;
    }

    /* Determine whether or not an accretion disk forms around the secondary
     * Follow prescription of https://ui.adsabs.harvard.edu/abs/1976ApJ...206..509U/abstract */
    accretor->emt_accretion_radius = R_accretor;

    accretor->accretion_disk_is_present = false;
    double q_accretor = m_accretor/m_donor;
    accretor->accretion_disk_r_min = 0.0425 * rp * pow(q_accretor*(1.0 + q_accretor), 0.25); /* ASH: taking rp instead of separation */

    if (accretor->accretion_disk_r_min < R_accretor)
    {
        accretor->accretion_disk_is_present = true;
        accretor->emt_accretion_radius = 1.7*accretor->accretion_disk_r_min;
    }

    /* Determine ejection mode */
    double temp = 1.0 - e;
    double n_peri = (TWOPI/P_orb) * sqrt( (1.0 + e)/(temp*temp*temp) );
    double ospin_hat = norm3(donor->spin_vec) / n_peri;
    double q_donor = m_donor/m_accretor;
    if (ospin_hat < 0.1) /* "Small" donor spin */
    {
        donor->emt_ejection_radius_mode = 1;
    }
    else if (q_donor > 10) /* "Large" mass ratio */
    {
        donor->emt_ejection_radius_mode = 2;
    }
    else /* Other cases: simply to not take into account finite size of donor */
    {
        donor->emt_ejection_radius_mode = 0;
    }
    
    
    /* Set mass transfer rates */
    donor->mass_dot_RLOF = -dm1/dt;
    accretor->mass_dot_RLOF = dm2/dt;

    donor->delta_mass_RLOF = -dm1;
    accretor->delta_mass_RLOF = dm2;
    donor->RLOF_timescale = dt;
    accretor->RLOF_timescale = dt;

    /* Assume that any mass not accreted is ejected from the accretor in an isotropic wind. -- TO DO: could the effect be instantaneous?
     * This somewhat mimics non-conservative mass transfer. */
    donor->mass_dot_adiabatic_ejection = 0.0;
    accretor->mass_dot_adiabatic_ejection = - (dm1 - dm2)/dt;

    #ifdef VERBOSE
    if (verbose_flag > 0)
    {
        printf("binary_evolution.cpp -- stable_mass_transfer_evolution -- donor %d accretor %d donor->mass_dot_RLOF %g accretor->mass_dot_RLOF %g accretor->mass_dot_adiabatic_ejection %g dt_binary_evolution old %g\n",donor->index,accretor->index,donor->mass_dot_RLOF,accretor->mass_dot_RLOF,accretor->mass_dot_adiabatic_ejection,*dt_binary_evolution);
        //print_system(particlesMap,*integration_flag);
    }
    #endif

    
    double fabs_m_dot_donor = fabs(dm1 / dt);
    double fabs_m_dot_accretor = fabs(dm2 / dt);
    
    if (dm1 <= epsilon)
    {
        dm1 = epsilon;
    }
    if (dm2 <= epsilon)
    {
        dm2 = epsilon;
    }
    *dt_binary_evolution = dt * binary_evolution_mass_transfer_timestep_parameter * CV_min( (m_donor/dm1), (m_accretor/dm2) );
    *dt_binary_evolution = CV_max(ODE_min_dt, *dt_binary_evolution);

    #ifdef VERBOSE
    if (verbose_flag > 1)
    {
        printf("binary_evolution.cpp -- stable_mass_transfer_evolution -- dt_binary_evolution new %g dm1 %g dm2 %g\n",*dt_binary_evolution,dm1,dm2);
    }
    #endif

    delete[] GB;
    delete[] tscls;
    delete[] lums;

    return 0;
}


void white_dwarf_hydrogen_accretion_boundaries_WBBP13(double m_WD, double *m_dot_lower, double *m_dot_upper)
{
    /* Based on the simulations of Wolf, Bildsten, Brooks, and Paxton (2013; https://ui.adsabs.harvard.edu/abs/2013ApJ...777..136W/abstract) 
     * Extrapolates the two boundaries for giant star formation (m_dot > m_dot_upper), and for nova flashes (m_dot < m_dot_lower) */
    
    /* When extrapolation is called for, use linear extrapolation (last argument given to interpolation function is "false") */
    *m_dot_lower = interpolate_1D_data_table_linear(m_WD, SINGLE_DEGENERATE_WBBP_DATA_LOWER_ACCRETION_RATE, SINGLE_DEGENERATE_WBBP_DATA_LOWER_ACCRETION_RATE_TABLELENGTH,false);
    *m_dot_upper = interpolate_1D_data_table_linear(m_WD, SINGLE_DEGENERATE_WBBP_DATA_UPPER_ACCRETION_RATE, SINGLE_DEGENERATE_WBBP_DATA_UPPER_ACCRETION_RATE_TABLELENGTH,false);
    
    /* Make sure that extrapolation does not cause unphysical results (i.e., negative accretion rates) */
    if (*m_dot_lower <= 0.0)
    {
        *m_dot_lower = 0.0;
    }
    if (*m_dot_upper <= 0.0)
    {
        *m_dot_upper = 0.0;
    }
}


void white_dwarf_helium_mass_accumulation_efficiency(double m_WD, double m_dot, double WD_luminosity, double *eta, int *WD_accretion_mode)
{
    double m_dot_WK_max = determine_WK11_max_accretion_rate(m_WD, WD_luminosity);
    double m_dot_KH_min = 0.0;
    double m_dot_KH_max = 0.0;
    double eta_ret = -1;
    white_dwarf_helium_mass_accumulation_efficiency_KH04(m_WD,m_dot,&eta_ret,&m_dot_KH_min,&m_dot_KH_max);

    if (m_dot > 0.0 and m_dot < m_dot_WK_max)
    {
        eta_ret = 1.0;
        *WD_accretion_mode = 1; // accumulation of He onto the WD
    }
    else if (m_dot >= m_dot_WK_max and m_dot <= m_dot_KH_min)
    {
        eta_ret = 0.0;
        *WD_accretion_mode = 2; // no retention
    }
    else if (m_dot > m_dot_KH_min and m_dot < m_dot_KH_max)
    {
        *WD_accretion_mode = 3; // burning of accreted He to C; simply add to SSE WD mass; non-trivial efficiency
    }
    else if (m_dot >= m_dot_KH_max)
    {
        eta_ret = 1.0;
        *WD_accretion_mode = 3; // burning of accreted He to C; simply add to SSE WD mass
    }
    
    if (eta_ret == -1)// or m_dot_WK_max > m_dot_KH_min or m_dot_WK_max > m_dot_KH_max)
    {
        printf("binary_evolution.cpp -- white_dwarf_helium_mass_accumulation_efficiency() -- error in determining eta=%g; m_WD=%g; m_dot=%g m_dot_WK_max=%g m_dot_KH_min=%g m_dot_KH_max=%g\n",eta_ret,m_WD,m_dot,m_dot_WK_max,m_dot_KH_min,m_dot_KH_max);
        error_code = 42;
        longjmp(jump_buf,1);
    }
    
    //printf("binary_evolution.cpp -- white_dwarf_helium_mass_accumulation_efficiency -- eta %g m_dot_WK_max %g m_dot_KH_min %g m_dot_KH_max %g\n",eta_ret,m_dot_WK_max,m_dot_KH_min,m_dot_KH_max);
    
    *eta = eta_ret;
    
}


void white_dwarf_helium_mass_accumulation_efficiency_KH04(double m_WD, double m_dot, double *eta, double *m_dot_KH_min, double *m_dot_KH_max)
{
    /* Based on Kato & Hachisu (2004; https://ui.adsabs.harvard.edu/abs/2004ApJ...613L.129K/abstract) */
    double log_m_dot = log10(m_dot); // Log 10 of Helium accretion rate
    
    int KH04_WD_masses_length = 6;
    double KH04_WD_masses[KH04_WD_masses_length] = {0.8,0.9,1.0,1.2,1.3,1.35};

    int i,i_low,i_up;
    double log_m_dot_KH_min,log_m_dot_KH_max;

    /* First check if m_WD exceeds boundary values of KH04 */
    if (m_WD < KH04_WD_masses[0])
    {
        white_dwarf_helium_mass_accumulation_efficiency_KH04_index0(log_m_dot,eta,&log_m_dot_KH_min,&log_m_dot_KH_max);
        *m_dot_KH_min = pow(10.0,log_m_dot_KH_min);
        *m_dot_KH_max = pow(10.0,log_m_dot_KH_max);
        return;
    }
    else if (m_WD > KH04_WD_masses[KH04_WD_masses_length-1])
    {
        white_dwarf_helium_mass_accumulation_efficiency_KH04_index5(log_m_dot,eta,&log_m_dot_KH_min,&log_m_dot_KH_max);
        *m_dot_KH_min = pow(10.0,log_m_dot_KH_min);
        *m_dot_KH_max = pow(10.0,log_m_dot_KH_max);
        return;
    }

    typedef void (*functionarray) (double log_m_dot, double *eta, double *log_m_dot_min, double *log_m_dot_max);
    functionarray KH04_functions[] = 
    {
        white_dwarf_helium_mass_accumulation_efficiency_KH04_index0, \
        white_dwarf_helium_mass_accumulation_efficiency_KH04_index1, \
        white_dwarf_helium_mass_accumulation_efficiency_KH04_index2, \
        white_dwarf_helium_mass_accumulation_efficiency_KH04_index3, \
        white_dwarf_helium_mass_accumulation_efficiency_KH04_index4, \
        white_dwarf_helium_mass_accumulation_efficiency_KH04_index5
    };

    /* Determine which mass bin applies */
    double m_low,m_up;
    for (i=0; i<KH04_WD_masses_length-1; i++)
    {
        m_low = KH04_WD_masses[i];
        m_up = KH04_WD_masses[i+1];
        if (m_WD >= m_low and m_WD < m_up)
        {
            i_low = i;
            i_up = i+1;
            break;
        }
    }
    
    /* Determine efficiencies and critical m_dots for the boundaries in the mass bins, and interpolate linearly */
    double eta_low,log_m_dot_min_low,log_m_dot_max_low;
    double eta_up,log_m_dot_min_up,log_m_dot_max_up;
    KH04_functions[i_low](log_m_dot, &eta_low, &log_m_dot_min_low, &log_m_dot_max_low);
    KH04_functions[i_up](log_m_dot, &eta_up, &log_m_dot_min_up, &log_m_dot_max_up);
    
    *eta = (eta_up - eta_low) * (m_WD - m_low) / (m_up - m_low) + eta_low;
    *m_dot_KH_min = pow(10.0, (log_m_dot_min_up - log_m_dot_min_low) * (m_WD - m_low) / (m_up - m_low) + log_m_dot_min_low);
    *m_dot_KH_max = pow(10.0, (log_m_dot_max_up - log_m_dot_max_low) * (m_WD - m_low) / (m_up - m_low) + log_m_dot_max_low);
    
    //printf("KH04 log_m_dot %g i_low %d i_up %d m_low %g m_up %g m_WD %g eta %g m_dot_KH_min %g m_dot_KH_max %g eta_low %g eta_up %g\n",log_m_dot,i_low,i_up,m_low,m_up,m_WD,*eta,*m_dot_KH_min,*m_dot_KH_max,eta_low,eta_up);
}

void white_dwarf_helium_mass_accumulation_efficiency_KH04_index0(double log_m_dot, double *eta, double *log_m_dot_min, double *log_m_dot_max)
{
    *log_m_dot_min = -6.5;
    *log_m_dot_max = -6.34;
    
    if (log_m_dot >= *log_m_dot_max)
    {
        *eta = 1.0;
    }
    else if (log_m_dot < *log_m_dot_max and log_m_dot > *log_m_dot_min)
    {
        *eta = -0.35 * (log_m_dot + 6.1)*(log_m_dot + 6.1) + 1.02;
    }
    else
    {
        *eta = 0.0;
    }
}
void white_dwarf_helium_mass_accumulation_efficiency_KH04_index1(double log_m_dot, double *eta, double *log_m_dot_min, double *log_m_dot_max)
{
    *log_m_dot_min = -6.88;
    *log_m_dot_max = -6.05;
    if (log_m_dot >= *log_m_dot_max)
    {
        *eta = 1.0;
    }
    else if (log_m_dot < *log_m_dot_max and log_m_dot > *log_m_dot_min)
    {
        *eta = -0.35 * (log_m_dot + 5.6)*(log_m_dot + 5.6) + 1.07;
    }
    else
    {
        *eta = 0.0;
    }
}
void white_dwarf_helium_mass_accumulation_efficiency_KH04_index2(double log_m_dot, double *eta, double *log_m_dot_min, double *log_m_dot_max)
{
    *log_m_dot_min = -6.92;
    *log_m_dot_max = -5.93;
    if (log_m_dot >= *log_m_dot_max)
    {
        *eta = 1.0;
    }
    else if (log_m_dot < *log_m_dot_max and log_m_dot > *log_m_dot_min)
    {
        *eta = -0.35 * (log_m_dot + 5.6)*(log_m_dot + 5.6) + 1.01;
    }
    else
    {
        *eta = 0.0;
    }
}
void white_dwarf_helium_mass_accumulation_efficiency_KH04_index3(double log_m_dot, double *eta, double *log_m_dot_min, double *log_m_dot_max)
{
    *log_m_dot_min = -7.06;
    *log_m_dot_max = -5.76;
    if (log_m_dot >= *log_m_dot_max)
    {
        *eta = 1.0;
    }
    else if (log_m_dot < *log_m_dot_max and log_m_dot >= -5.95)
    {
        *eta = -0.54 * (log_m_dot + 5.6)*(log_m_dot + 5.6) + 1.01;
    }
    else if (log_m_dot < -5.95 and log_m_dot > *log_m_dot_min)
    {
        *eta = 0.54 * log_m_dot + 4.16;
    }
    else
    {
        *eta = 0.0;
    }
}
void white_dwarf_helium_mass_accumulation_efficiency_KH04_index4(double log_m_dot, double *eta, double *log_m_dot_min, double *log_m_dot_max)
{
    *log_m_dot_min = -7.35;
    *log_m_dot_max = -5.83;
    if (log_m_dot >= *log_m_dot_max)
    {
        *eta = 1.0;
    }
    else if (log_m_dot < *log_m_dot_max and log_m_dot > *log_m_dot_min)
    {
        *eta = -0.175 * (log_m_dot + 5.35)*(log_m_dot + 5.35) + 1.03;
    }
    else
    {
        *eta = 0.0;
    }
}
void white_dwarf_helium_mass_accumulation_efficiency_KH04_index5(double log_m_dot, double *eta, double *log_m_dot_min, double *log_m_dot_max)
{
    *log_m_dot_min = -7.4;
    *log_m_dot_max = -6.05;
    if (log_m_dot >= *log_m_dot_max)
    {
        *eta = 1.0;
    }
    else if (log_m_dot < *log_m_dot_max and log_m_dot > *log_m_dot_min)
    {
        *eta = -0.115 * (log_m_dot + 5.7)*(log_m_dot + 5.7) + 1.01;
    }
    else
    {
        *eta = 0.0;
    }
}

double determine_WK11_max_accretion_rate(double M_WD, double WD_luminosity)
{
    #ifdef VERBOSE
    if (verbose_flag > 1)
    {
        printf("binary_evolution.cpp -- determine_WK11_max_accretion_rate -- M_WD %g WD_luminosity %g\n",M_WD,WD_luminosity/CONST_L_SUN);
    }
    #endif
    
   
    double L_LSun = WD_luminosity/CONST_L_SUN;
   
    /* First check the luminosity for boundary excession, and interpolate if necessary over luminosity */
    if (L_LSun < SINGLE_DEGENERATE_WK_DATA_COLD_LUMINOSITY_LSUN)
    {
        return determine_WK11_max_accretion_rate_for_given_luminosity(M_WD, WD_luminosity, SINGLE_DEGENERATE_WK_DATA_MAX_ACCRETION_RATE_COLD, SINGLE_DEGENERATE_WK_DATA_MAX_ACCRETION_RATE_TABLELENGTH);
    }
    else if (L_LSun > SINGLE_DEGENERATE_WK_DATA_HOT_LUMINOSITY_LSUN)
    {
        return determine_WK11_max_accretion_rate_for_given_luminosity(M_WD, WD_luminosity, SINGLE_DEGENERATE_WK_DATA_MAX_ACCRETION_RATE_HOT, SINGLE_DEGENERATE_WK_DATA_MAX_ACCRETION_RATE_TABLELENGTH);
    }
    else
    {
        double m_dot1 = determine_WK11_max_accretion_rate_for_given_luminosity(M_WD, WD_luminosity, SINGLE_DEGENERATE_WK_DATA_MAX_ACCRETION_RATE_COLD, SINGLE_DEGENERATE_WK_DATA_MAX_ACCRETION_RATE_TABLELENGTH);
        double m_dot2 = determine_WK11_max_accretion_rate_for_given_luminosity(M_WD, WD_luminosity, SINGLE_DEGENERATE_WK_DATA_MAX_ACCRETION_RATE_HOT, SINGLE_DEGENERATE_WK_DATA_MAX_ACCRETION_RATE_TABLELENGTH);

        double m_dot = (m_dot2 - m_dot1) * (L_LSun - SINGLE_DEGENERATE_WK_DATA_COLD_LUMINOSITY_LSUN) / (SINGLE_DEGENERATE_WK_DATA_HOT_LUMINOSITY_LSUN-SINGLE_DEGENERATE_WK_DATA_COLD_LUMINOSITY_LSUN) + m_dot1;

        return m_dot;
    }
}

double determine_WK11_max_accretion_rate_for_given_luminosity(double M_WD, double WD_luminosity, const double (*SD_table)[SINGLE_DEGENERATE_WK_DATA_MAX_ACCRETION_RATE_TABLEWIDTH], int SD_table_length)
{

    double tab_m_max = SD_table[0][0];
    double tab_m_min = SD_table[SD_table_length-1][0];
    
    double ret_val;
    
    /* If the WD mass is beyond the range in the table, assume the boundary accretion rate value */
    if (M_WD > tab_m_max)
    {
        ret_val = SD_table[0][1] * 1.0e-8;
    }
    else if (M_WD < tab_m_min)
    {
        ret_val = SD_table[SD_table_length-1][1] * 1.0e-8;
    }
    else
    {
        double m_low,m_up,acc_low,acc_up;
        int i;
        
        for (i=0; i<SD_table_length-1; i++)
        {
            m_low = SD_table[i+1][0];
            acc_low = SD_table[i+1][1];
            m_up = SD_table[i][0];
            acc_up = SD_table[i][1];
            
            if (M_WD >= m_low and M_WD < m_up)
            {
                break;
            }
        }
        double acc = (acc_up - acc_low) * (M_WD - m_low) / (m_up - m_low) + acc_low;
        
        ret_val = acc * 1.0e-8;
    }

    //printf("binary_evolution.cpp -- determine_WK11_max_accretion_rate -- ret_val %g\n",ret_val);
    
    return ret_val;
}

bool determine_if_He_accreting_WD_explodes(double M_WD, double M_dot_acc, double M_He, double WD_luminosity)
{
    #ifdef VERBOSE
    if (verbose_flag > 1)
    {
        printf("binary_evolution.cpp -- determine_if_He_accreting_WD_explodes -- M_WD %g M_dot_acc %g  M_He %g WD_luminosity %g\n",M_WD,M_dot_acc,M_He,WD_luminosity/CONST_L_SUN);
    }
    #endif
   
    double L_LSun = WD_luminosity/CONST_L_SUN;

    double ret_val = -1;
    
    /* First check the luminosity for boundary excession, and interpolate if necessary over luminosity */
    if (L_LSun < SINGLE_DEGENERATE_WK_DATA_COLD_LUMINOSITY_LSUN)
    {
        ret_val = determine_if_He_accreting_WD_explodes_for_given_luminosity(M_WD, M_dot_acc, M_He, SINGLE_DEGENERATE_WK_DATA_COLD, SINGLE_DEGENERATE_WK_DATA_COLD_TABLELENGTH);
    }
    else if (L_LSun > SINGLE_DEGENERATE_WK_DATA_HOT_LUMINOSITY_LSUN)
    {
        ret_val = determine_if_He_accreting_WD_explodes_for_given_luminosity(M_WD, M_dot_acc, M_He, SINGLE_DEGENERATE_WK_DATA_HOT, SINGLE_DEGENERATE_WK_DATA_HOT_TABLELENGTH);
    }
    else
    {
        double ret_val1 = determine_if_He_accreting_WD_explodes_for_given_luminosity(M_WD, M_dot_acc, M_He, SINGLE_DEGENERATE_WK_DATA_COLD, SINGLE_DEGENERATE_WK_DATA_COLD_TABLELENGTH);
        double ret_val2 = determine_if_He_accreting_WD_explodes_for_given_luminosity(M_WD, M_dot_acc, M_He, SINGLE_DEGENERATE_WK_DATA_HOT, SINGLE_DEGENERATE_WK_DATA_HOT_TABLELENGTH);
        
        ret_val = (ret_val2 - ret_val1) * (L_LSun - SINGLE_DEGENERATE_WK_DATA_COLD_LUMINOSITY_LSUN) / (SINGLE_DEGENERATE_WK_DATA_HOT_LUMINOSITY_LSUN-SINGLE_DEGENERATE_WK_DATA_COLD_LUMINOSITY_LSUN) + ret_val1;
    }

    bool explosion = false;

    if (ret_val == -1 or ret_val != ret_val)
    {
        printf("binary_evolution.cpp -- determine_if_He_accreting_WD_explodes() -- ERROR: ret_val=%g\n",ret_val);
        error_code = 43;
        longjmp(jump_buf,1);
    }
    else
    {
        /* ret_val interpolates between 0 and 1, but ultimately a binary decision needs to be made; make the "split" at ret_val = 0.5 */
        if (ret_val < 0.5)
        {
            explosion = false;
        }
        else
        {
            explosion = true;
        }
    }
    
    return explosion;
}

double determine_if_He_accreting_WD_explodes_for_given_luminosity(double M_WD, double M_dot_acc, double M_He, const double (*SD_table)[SINGLE_DEGENERATE_WK_DATA_TABLEWIDTH], int SD_table_length)
{

    /* Minimum and maximum WD masses in the data table */
    double tab_m_max = SD_table[0][0];
    double tab_m_min = SD_table[SD_table_length-1][0];
    double M_dot_acc_1e8 = M_dot_acc*1.0e8;
        
    int i;
    double tab_m, tab_acc, tab_He;

    double ret_val;
    
    if (M_WD > tab_m_max or M_WD < tab_m_min) /* If the WD mass is beyond the range in the table, assume no explosion occurs */
    {
        ret_val = 0.0;
    }
    else
    {
        /* Extract (M_dot_acc,M_He) data for the two WD masses nearest to M_WD in the data table (M_WD1 & M_WD2); this will later be interpolated over both M_WD and M_He. */
        
        double acc_data1[SINGLE_DEGENERATE_WK_DATA_LONGEST_TABLELENGTH_FIXED_M_WD];
        double He_data1[SINGLE_DEGENERATE_WK_DATA_LONGEST_TABLELENGTH_FIXED_M_WD];
        int data_length1 = 0;

        double acc_data2[SINGLE_DEGENERATE_WK_DATA_LONGEST_TABLELENGTH_FIXED_M_WD];
        double He_data2[SINGLE_DEGENERATE_WK_DATA_LONGEST_TABLELENGTH_FIXED_M_WD];
        int data_length2 = 0;

        double M_WD1 = -1; // Lower WD mass from data table
        double M_WD2 = -1; // Upper WD mass from data table
        
        for (i=0; i<SD_table_length; i++)
        {
            tab_m = SD_table[i][0]; // WD mass
            tab_acc = SD_table[i][1]; // unit: 1e-8 MSun/yr
            tab_He = SD_table[i][2]; // He layer mass

            if (M_WD == tab_m) // WD mass happens to match exactly one of the tabulated values; no need to interpolate over WD mass and hence M_WD1 and M_WD2 do not apply
            {
                He_data1[data_length1] = tab_He;
                acc_data1[data_length1] = tab_acc;
                data_length1++;

                M_WD1 = -2;
                M_WD2 = -2;
            }
            else if (M_WD < tab_m and M_WD > tab_m - SINGLE_DEGENERATE_WK_DATA_DELTA_M_WD)
            {
                He_data1[data_length1] = tab_He;
                acc_data1[data_length1] = tab_acc;
                data_length1++;
                if (M_WD1 == -1)
                {
                    M_WD1 = tab_m;
                }
            }
            else if (M_WD > tab_m and M_WD < tab_m + SINGLE_DEGENERATE_WK_DATA_DELTA_M_WD)
            {
                He_data2[data_length2] = tab_He;
                acc_data2[data_length2] = tab_acc;
                data_length2++;
                if (M_WD2 == -1)
                {
                    M_WD2 = tab_m;
                }
            }
            else
            {
                continue;
            }
        }

        if (M_WD1 == -1 or M_WD2 == -1)
        {
            /* Upper and lower WD masses were not found in the table; this is an error */
            printf("binary_evolution.cpp -- determine_if_He_accreting_WD_explodes_for_given_luminosity() -- ERROR: M_WD=%g M_WD1 = %g; M_WD2 = %g\n",M_WD,M_WD1,M_WD2);
            error_code = 43;
            longjmp(jump_buf,1);
        }
        
        if (M_WD1 == -2 and M_WD2 == -2)
        {
            /* WD mass happens to match exactly one of the tabulated values; no need to interpolate over WD mass */
            ret_val = determine_if_He_accreting_WD_explodes_for_given_luminosity_and_WD_mass(M_He,M_dot_acc_1e8,acc_data1,He_data1,data_length1);
        }
        else
        {
            /* Interpolate over WD mass */
            double ret_val1 = determine_if_He_accreting_WD_explodes_for_given_luminosity_and_WD_mass(M_He,M_dot_acc_1e8,acc_data1,He_data1,data_length1);
            double ret_val2 = determine_if_He_accreting_WD_explodes_for_given_luminosity_and_WD_mass(M_He,M_dot_acc_1e8,acc_data2,He_data2,data_length2);
            
            ret_val = (ret_val2 - ret_val1) * (M_WD - M_WD1) / (M_WD2 - M_WD1) + ret_val1;
        }
    }

    return ret_val;
}

double determine_if_He_accreting_WD_explodes_for_given_luminosity_and_WD_mass(double M_He, double M_dot_acc_1e8, double acc_data[SINGLE_DEGENERATE_WK_DATA_LONGEST_TABLELENGTH_FIXED_M_WD], double He_data[SINGLE_DEGENERATE_WK_DATA_LONGEST_TABLELENGTH_FIXED_M_WD], int data_length)
{
    /* If the He layer mass lies within available data range, interpolate linearly. Otherwise, assume an explosion never occurs. */
    int i;
    double He_low,He_up,acc_low,acc_up;
    double acc_crit = -1;
    double ret_val;
    
    for (i=0; i<data_length-1; i++)
    {
        He_low = He_data[i];
        He_up = He_data[i+1];
        acc_low = acc_data[i];
        acc_up = acc_data[i+1];

        if (M_He > He_low and M_He <= He_up)
        {
            //printf("M_He %g He_low %g He_up %g\n",M_He,He_low,He_up);
            acc_crit = (acc_up - acc_low) * (M_He - He_low) / (He_up - He_low) + acc_low;
            break;
        }
    }
    
    if (acc_crit != -1 and M_dot_acc_1e8 < acc_crit)
    {
        /* He layer mass lies within valid range and the corresponding accretion rate is below the critical value for an explosion. */
        ret_val = 1.0;
    }
    else
    {
        ret_val = 0.0;
    }

    return ret_val;
}

int triple_stable_mass_transfer_evolution(ParticlesMap *particlesMap, int parent_index, int donor_index, int accretor_index, double t_old, double t, int *integration_flag, double *dt_binary_evolution)
{
    #ifdef VERBOSE
    if (verbose_flag > 0)
    {
        printf("binary_evolution.cpp -- triple_stable_mass_transfer_evolution -- start\n");
        print_system(particlesMap,*integration_flag);
    }
    #endif
    

    Particle *outer_binary = (*particlesMap)[parent_index];
    Particle *donor = (*particlesMap)[donor_index];
    Particle *inner_binary = (*particlesMap)[accretor_index];

    /* Check whether all prerequisites are met. */
    if (donor->is_binary == true or inner_binary->is_binary == false)
    {
        printf("binary_evolution.cpp -- triple_stable_mass_transfer_evolution() -- ERROR: donor should be a star, and accretor a binary\n");
        error_code = 5;
        longjmp(jump_buf,1);
        //exit(-1);
    }

    Particle *c1 = (*particlesMap)[inner_binary->child1];
    Particle *c2 = (*particlesMap)[inner_binary->child2];

    if (c1->is_binary == true or c2->is_binary == true)
    {
        printf("binary_evolution.cpp -- triple_stable_mass_transfer_evolution() -- at least of the components in the inner binary is itself a binary -- skipping. \n");
        return 0;
    }
    
    if (*integration_flag != 0)
    {
        printf("binary_evolution.cpp -- triple_stable_mass_transfer_evolution() -- integration flag is nonzero; skipping. \n");
        return 0;
    }

    /* Define primary (star1) as the most massive star in the inner binary */
    Particle *star1, *star2;
    if (c1->mass > c2->mass)
    {
        star1 = (*particlesMap)[inner_binary->child1];
        star2 = (*particlesMap)[inner_binary->child2];
    }
    else
    {
        star1 = (*particlesMap)[inner_binary->child2];
        star2 = (*particlesMap)[inner_binary->child1];
    }

    set_up_derived_quantities(particlesMap); /* for setting a, e, etc. for all orbits */
    
    double m_donor = donor->mass;
    double m_inner_binary = inner_binary->mass;
    double m1 = star1->mass;
    double m2 = star2->mass;

    double q = m_donor/m_inner_binary;

    double R_donor = donor->radius;

    int kw1 = donor->stellar_type;

    double t_KH_donor = compute_Kelvin_Helmholtz_timescale(kw1,m_donor,donor->core_mass,R_donor,donor->luminosity);
    double t_dyn_donor = compute_stellar_dynamical_timescale(m_donor, R_donor);

    /* These quantities are used below to determine aging */
    double tms,tn;
    double *GB,*tscls,*lums;
    GB = new double[10];
    tscls = new double[20];
    lums = new double[10];  

    star_(&kw1,&donor->sse_initial_mass,&m_donor,&tms,&tn,tscls,lums,GB,donor->zpars);
    double tms_donor_old = tms;
    double tbgb_donor_old = tscls[0];
    
    star_(&star1->stellar_type,&star1->sse_initial_mass,&m1,&tms,&tn,tscls,lums,GB,star1->zpars);
    double tms_star1_old = tms;
    double tbgb_star1_old = tscls[0];
    
    star_(&star2->stellar_type,&star2->sse_initial_mass,&m2,&tms,&tn,tscls,lums,GB,star2->zpars);
    double tms_star2_old = tms;
    double tbgb_star2_old = tscls[0];

    double a_out = outer_binary->a;
    double e_out = outer_binary->e;
    double rp_out = a_out*(1.0-e_out);

    double a_in = inner_binary->a;
    double e_in = inner_binary->e;
    double rp_in = a_in*(1.0-e_in);
    
    double dt = t - t_old;

    /* Default value of dm3 -- transferred mass during this time-step */
    double R_RL_av_donor = roche_radius_pericenter_eggleton(rp_out,q);
    //double dm3 = compute_bse_mass_transfer_amount(kw1,m_donor,donor->core_mass,R_donor,R_RL_av_donor,dt,t_dyn_donor,t_KH_donor);
    double dm3 = compute_bse_mass_transfer_amount_averaged(kw1, m_donor, donor->core_mass, R_donor, m_inner_binary, a_out, e_out, dt, t_dyn_donor, t_KH_donor);
    
    donor->mass_dot_RLOF_triple = -dm3/dt;
    
    /* Determine the amount of mass accreted in the inner binary. 
     * Depends on whether or not a circumbinary accretion disk can form. */
    double q_accretor = m_inner_binary/m_donor;
    double ra_in = a_in * (1.0 + e_in);
    double circum_inner_binary_accretion_disk_r_min = 0.0425 * rp_out * pow(q_accretor*(1.0 + q_accretor), 0.25); /* https://ui.adsabs.harvard.edu/abs/1976ApJ...206..509U/abstract */

    bool circum_inner_binary_accretion_disk_is_present = false;
    if (circum_inner_binary_accretion_disk_r_min < ra_in)
    {
        circum_inner_binary_accretion_disk_is_present = true;
    }
    
    double dm1,dm2;
    if (circum_inner_binary_accretion_disk_is_present == false)
    {
        /* Expect that most of the mass from the tertiary will not be accreted by the inner binary */
        dm1 = triple_mass_transfer_primary_star_accretion_efficiency_no_disk * dm3;
        dm2 = triple_mass_transfer_secondary_star_accretion_efficiency_no_disk * dm3;
    }
    else
    {
        /* Expect that more mass is accreted when a circumbinary disk is present. */
        dm1 = triple_mass_transfer_primary_star_accretion_efficiency_disk * dm3;
        dm2 = triple_mass_transfer_secondary_star_accretion_efficiency_disk * dm3;
    }

    star1->mass_dot_RLOF_triple = dm1/dt;
    star2->mass_dot_RLOF_triple = dm2/dt;
    
    
    /****************************
     *  Aging and rejuvenation *
     ***************************/

    double m1_new = m1 + dm1;
    double m2_new = m2 + dm2;
    double m_donor_new = m_donor + dm3;
    
    /* MS stars: simply set m_init = m_new */
    if (kw1 <= 1 or kw1 == 7)
    {
        donor->sse_initial_mass = m_donor_new;
    }
    if (star1->stellar_type <= 1 or star1->stellar_type == 7)
    {
        star1->sse_initial_mass = m1_new;
    }
    if (star2->stellar_type <= 1 or star2->stellar_type == 7)
    {
        star2->sse_initial_mass = m2_new;
    }
    
    /* For HG stars check if the initial mass can be reduced. */    
    if (kw1 == 2 and donor->sse_initial_mass <= donor->zpars[2])
    {
        double m0 = donor->sse_initial_mass;
        donor->sse_initial_mass = m_donor;

        star_(&kw1,&donor->sse_initial_mass,&donor->mass,&tms,&tn,tscls,lums,GB,donor->zpars);
        if (GB[8] < donor->core_mass)
        {
            donor->sse_initial_mass = m0;
        }
    }
    if (star1->stellar_type == 2 and star1->sse_initial_mass <= star1->zpars[2])
    {
        double m0 = star1->sse_initial_mass;
        star1->sse_initial_mass = m1;

        star_(&star1->stellar_type,&star1->sse_initial_mass,&star1->mass,&tms,&tn,tscls,lums,GB,star1->zpars);
        if (GB[8] < star1->core_mass)
        {
            star1->sse_initial_mass = m0;
        }
    }
    if (star2->stellar_type == 2 and star2->sse_initial_mass <= star2->zpars[2])
    {
        double m0 = star2->sse_initial_mass;
        star2->sse_initial_mass = m2;

        star_(&star2->stellar_type,&star2->sse_initial_mass,&star2->mass,&tms,&tn,tscls,lums,GB,star2->zpars);
        if (GB[8] < star2->core_mass)
        {
            star2->sse_initial_mass = m0;
        }
    }

    /* Always rejuvenate the accretors and age the donor if they are on the main sequence. */
    if (kw1 <= 2 or kw1 == 7)
    {
        star_(&kw1,&donor->sse_initial_mass,&m_donor_new,&tms,&tn,tscls,lums,GB,donor->zpars);

        if (kw1 == 2)
        {
            double age = tms + (donor->age*yr_to_Myr - tms_donor_old) * (tscls[0] - tms)/(tbgb_donor_old - tms_donor_old);
            donor->age = age*Myr_to_yr;
        }
        else
        {
            donor->age *= (tms/tms_donor_old);
            //printf("age %g %g %.15f m_donor_new %g m_donor %g\n",tms,tms_donor_old,tms/tms_donor_old,m_donor_new,m_donor);
        }
        donor->epoch = t - donor->age;
    }
    if (star1->stellar_type <= 2 or star1->stellar_type == 7)
    {
        star_(&star1->stellar_type,&star1->sse_initial_mass,&m1_new,&tms,&tn,tscls,lums,GB,star1->zpars);
        if (star1->stellar_type == 2)
        {
            double age = tms + (star1->age*yr_to_Myr - tms_star1_old) * (tscls[0] - tms)/(tbgb_star1_old - tms_star1_old);
            star1->age = age*Myr_to_yr;
        }
        else if ((m1_new < 0.35 or m1_new > 1.25) and (star1->stellar_type != 7))
        {
            star1->age *= (tms/tms_star1_old) * ((m1_new - dm1)/m1_new);
        }
        star1->epoch = t - star1->age;
    }
    if (star2->stellar_type <= 2 or star2->stellar_type == 7)
    {
        star_(&star2->stellar_type,&star2->sse_initial_mass,&m2_new,&tms,&tn,tscls,lums,GB,star2->zpars);
        if (star2->stellar_type == 2)
        {
            double age = tms + (star2->age*yr_to_Myr - tms_star2_old) * (tscls[0] - tms)/(tbgb_star2_old - tms_star2_old);
            star2->age = age*Myr_to_yr;
        }
        else if ((m2_new < 0.35 or m2_new > 1.25) and (star2->stellar_type != 7))
        {
            star2->age *= (tms/tms_star2_old) * ((m2_new - dm2)/m2_new);
        }
        star2->epoch = t - star2->age;
    }
    
    double a_in_f = a_in * ( (m1 + dm1)*(m2 + dm2) / (m1*m2 + 2.0*(m1+m2)*dm3/triple_mass_transfer_inner_binary_alpha_times_lambda) );
    inner_binary->triple_mass_transfer_a_in_dot = (a_in_f - a_in)/dt;
    
    /* Effect of mass loss of triple subsystem on the rest of the system. */
    /* Assume the mass not accreted by the inner binary is lost in an adiabatic wind. */
    inner_binary->mass_dot_wind = -( dm3 - (dm1 + dm2) )/dt;
 
    check_number(donor->mass_dot_RLOF_triple,                  "binary_evolution.cpp -- triple_stable_mass_transfer_evolution()","donor->mass_dot_RLOF_triple", true);
    check_number(star1->mass_dot_RLOF_triple,                  "binary_evolution.cpp -- triple_stable_mass_transfer_evolution()","star1->mass_dot_RLOF_triple", true);
    check_number(star2->mass_dot_RLOF_triple,                  "binary_evolution.cpp -- triple_stable_mass_transfer_evolution()","star2->mass_dot_RLOF_triple", true);
    check_number(inner_binary->triple_mass_transfer_a_in_dot,                  "binary_evolution.cpp -- triple_stable_mass_transfer_evolution()","inner_binary->triple_mass_transfer_a_in_dot", true);
    check_number(inner_binary->mass_dot_wind,                  "binary_evolution.cpp -- triple_stable_mass_transfer_evolution()","inner_binary->mass_dot_wind", true);
 
    delete[] GB;
    delete[] tscls;
    delete[] lums;

    #ifdef VERBOSE
    if (verbose_flag > 0)
    {
        printf("binary_evolution.cpp -- triple_stable_mass_transfer_evolution -- end\n");
        print_system(particlesMap,*integration_flag);
    }
    #endif
 
    return 0;
}

double compute_bse_mass_transfer_amount_averaged(int kw1, double m_donor, double core_mass_donor, double R_donor, double m_accretor, double a, double e, double dt, double t_dyn_donor, double t_KH_donor)
{
    double q = m_donor/m_accretor;
    double R_Lc = roche_radius_pericenter_eggleton(a,q); /* with argument "a", actually computes circular Roche lobe radius */
    double x = R_Lc/R_donor;
    double E_0;
    bool in_RLOF;
    determine_E_0(e, x, &E_0, &in_RLOF);
    
    if (in_RLOF == false)
    {
        return 0.0;
    }
    if (E_0 <= epsilon)
    {
        E_0 = epsilon;
    }
    
    int N_points = binary_evolution_numerical_mass_transfer_rate_number_of_points;
    double dE = 2.0 * E_0/(double(N_points));
    double integral = 0.0;
    double E = -E_0;
    double r_div_a,R_L,dm;
    for (int i=0; i<N_points; i++)
    {
        r_div_a = 1.0 - e*cos(E);
        R_L = R_Lc * r_div_a;
        dm = compute_bse_mass_transfer_amount(kw1, m_donor, core_mass_donor, R_donor, R_L, dt, t_dyn_donor, t_KH_donor);
        integral += dm * r_div_a * dE;
        E += dE;
        //printf("E %g dE %g E_0 %g R_L %g r_div_a %g dm %g integral %g\n",E,dE,E_0,R_L,r_div_a,dm,integral);
    }
    double dm_av = integral/(2.0*E_0);

    #ifdef IGNORE
    /* Same calculation but based on true anomaly (not used) */
    N_points = 100;
    double nu,dnu,cos_nu_0,sin_nu_0;
    compute_true_anomaly_from_eccentric_anomaly(cos(E_0), sin(E_0), e, &cos_nu_0, &sin_nu_0);
    double nu_0 = atan2(sin_nu_0,cos_nu_0);
    dnu = 2.0 * nu_0/(double(N_points));
    
    nu = -nu_0;
    integral = 0.0;
    double omesq = 1.0 - e*e;
    double one_div_sqomesq = 1.0/sqrt(omesq);
    double jac;
    for (int i=0; i<N_points; i++)
    {
        r_div_a = omesq/(1.0 + e*cos(nu));
        jac = r_div_a*r_div_a*one_div_sqomesq;
        R_L = R_Lc * r_div_a;
        dm = compute_bse_mass_transfer_amount(kw1, m_donor, core_mass_donor, R_donor, R_L, dt, t_dyn_donor, t_KH_donor);
        integral += jac * dm * dnu;
        nu += dnu;
        //printf("nu %g nu_0 %g R_L %g r_div_a %g dm %g integral %g\n",nu,nu_0,R_L,r_div_a,dm,integral);
    }
    dm_av = integral/(2.0*nu_0);
    //printf("3 dm_av %g %g\n",dm_av, compute_bse_mass_transfer_amount(kw1, m_donor, core_mass_donor, R_donor, roche_radius_pericenter_eggleton(a*(1.0-e),q), dt, t_dyn_donor, t_KH_donor));
    #endif
    
    check_number(dm_av,  "binary_evolution.cpp -- compute_bse_mass_transfer_amount_averaged","dm_av", true);

    #ifdef VERBOSE
    if (verbose_flag > 1)
    {
        printf("compute_bse_mass_transfer_amount_averaged kw1 %d  m_donor %g,  core_mass_donor %g,  R_donor %g,  m_accretor %g,  a %g,  e %g,  dt %g,  t_dyn_donor %g,  t_KH_donor %g E_0 %g dm_av %g dm_av_rp %g\n", kw1, m_donor, core_mass_donor, R_donor, m_accretor, a, e, dt, t_dyn_donor,  t_KH_donor, E_0, dm_av, compute_bse_mass_transfer_amount(kw1, m_donor, core_mass_donor, R_donor, roche_radius_pericenter_eggleton(a*(1.0-e),q), dt, t_dyn_donor, t_KH_donor));
    }
    #endif

    return dm_av;
}
   

double compute_bse_mass_transfer_amount(int kw1, double m_donor, double core_mass_donor, double R_donor, double R_RL_donor, double dt, double t_dyn_donor, double t_KH_donor)
{
    double m_fac = CV_min(m_donor, 5.0);
    m_fac *= m_fac;

    if (R_donor < R_RL_donor)
    {
        
        #ifdef VERBOSE
        if (verbose_flag > 0)
        {
            printf("binary_evolution.cpp -- compute_bse_mass_transfer_amount -- skipping since R_donor < R_RL_av -- R_donor %g RL_av_donor %g\n",R_donor,R_RL_donor);
        }
        #endif
        
        return 0.0;
    }
    
    double log_fac = log(R_donor/R_RL_donor);
    double dm1 = dt * fabs(3.0e-6 * m_fac * log_fac*log_fac*log_fac); /* HPT eq. 58-59 */

    if (kw1 == 2)
    {
        double mew = (m_donor - core_mass_donor)/m_donor;
        dm1 = CV_max(mew,0.01) * dm1;
    }
    else if (kw1 == 10)
    {
        dm1 = 1.0e3* dm1 * m_donor/CV_max(R_donor/CONST_R_SUN,1.0e-4);
    }
    
    if (kw1 >= 2 and kw1 <= 9 and kw1 != 7) /* Limit mass transfer to the thermal rate for giant-like stars */
    {
        dm1 = CV_min(dm1, m_donor * dt/t_KH_donor);
    }
    else /* Limit to dynamical rate for other cases. NOTE ASH: ignore for now the case "rad(j1).gt.10.d0*rol(j1).or.(kstar(j1).le.1.and.kstar(j2).le.1.and.q(j1).gt.qc" */
    {
        dm1 = CV_min(dm1, m_donor * dt/t_dyn_donor);
    }

    #ifdef VERBOSE
    if (verbose_flag > 2)
    {
        printf("binary_evolution.cpp -- compute_bse_mass_transfer_amount -- dm1 %g R_donor %g R_RL_donor %g\n",dm1,R_donor,R_RL_donor);
    }
    #endif

    return dm1;
}

void handle_mass_accretion_events_with_degenerate_objects(ParticlesMap *particlesMap, double t_old, double t, int *integration_flag, double *dt_binary_evolution)
{
    /* Special cases (in particular, when WDs become "too" massive) -- will potentially return before reaching the end of the current function. */
    if (*integration_flag != 0)
    {
        return;
    }

    #ifdef VERBOSE
    if (verbose_flag > 2)
    {
        printf("binary_evolution.cpp -- handle_mass_accretion_events_with_degenerate_objects -- begin \n");
        print_system(particlesMap,*integration_flag);
    }
    #endif
    
    double dt = t - t_old;
    
    double m1,m2,mdot1,mdot2,m1_new,m2_new,Delta_m1,Delta_m2;
    int kw1,kw2;
    
    double initial_momentum[3],initial_R_CM[3],initial_V_CM[3];
    double h_vec[3],h_vec_unit[3],e_vec[3],e_vec_unit[3],r_vec[3],v_vec[3];

    double final_R_CM[3],final_V_CM[3],final_momentum[3];

    ParticlesMapIterator it_p;
    bool found_new_event = true;
    bool structure_change = false;
    
    while (found_new_event == true)
    {

        found_new_event = false;
        for (it_p = particlesMap->begin(); it_p != particlesMap->end(); it_p++)
        {
            Particle *star2 = (*it_p).second;
            kw2 = star2->stellar_type;

            if (star2->is_binary == false and star2->is_bound == true and star2->object_type == 1 and kw2 >= 10 and kw2 <= 12)
            {
                Particle *star1 = (*particlesMap)[star2->sibling];
                Particle *parent = (*particlesMap)[star2->parent];
                if (star1->is_binary == false and star1->object_type == 1)
                {
                    m1 = star1->mass;
                    m2 = star2->mass;
                    
                    mdot1 = star1->mass_dot_wind + star1->mass_dot_wind_accretion + star1->mass_dot_adiabatic_ejection + star1->mass_dot_RLOF + star1->mass_dot_RLOF_triple;
                    mdot2 = star2->mass_dot_wind + star2->mass_dot_wind_accretion + star2->mass_dot_adiabatic_ejection + star2->mass_dot_RLOF + star2->mass_dot_RLOF_triple;
                    m1_new = m1 + mdot1 * dt; // new mass donor/accretor would have after this time-step (taking account wind mass losses)
                    m2_new = m2 + mdot2 * dt;
                    
                    kw1 = star1->stellar_type;

                    if (kw1 <= 10 and kw2 == 10 and m2_new >= 0.7)
                    {
                        if (binary_evolution_SNe_Ia_single_degenerate_model == 0 or binary_evolution_SNe_Ia_single_degenerate_model == 1)
                        {
                             /* HeWD can only accrete helium-rich material up to a mass of 0.7 when it is destroyed in a possible Type 1a SN. */
                             
                            found_new_event = true;
                            structure_change = true;
                        
                            Delta_m1 = m1_new - m1;
                            Delta_m2 = -m2;

                            #ifdef VERBOSE
                            if (verbose_flag > 0)
                            {
                                printf("binary_evolution.cpp -- handle_mass_accretion_events_with_degenerate_objects --  HeWD can only accrete helium-rich material up to a mass of 0.7 when it is destroyed in a possible Type 1a SN  -- star1 %d star2 %d - Delta_m1 %g Delta_m2 %g m1 %g m2 %g m1dot %g m2dot %g m1_new %g m2_new %g\n",star1->index,star2->index,Delta_m1,Delta_m2,m1,m2,mdot1,mdot2,m1_new,m2_new);
                                print_system(particlesMap,*integration_flag);
                            }
                            #endif

                            #ifdef LOGGING
                            Log_info_type log_info;
                            log_info.index1 = star2->index;
                            log_info.SNe_type = 1;
                            log_info.SNe_info = 1;
                            update_log_data(particlesMap, t, *integration_flag, LOG_SNE_START, log_info);
                            #endif
                
                            get_initial_binary_orbital_properties_from_position_and_velocity(star1->R_vec, star1->V_vec, star2->R_vec, star2->V_vec, m1, m2, \
                                r_vec, v_vec, initial_momentum, initial_R_CM, initial_V_CM, h_vec, e_vec);
                
                            handle_gradual_mass_loss_event_in_system(particlesMap, star1, star2, m1 + Delta_m1, m1, m2 + Delta_m2, m2, parent->dynamical_mass_transfer_low_mass_donor_timescale, \
                                r_vec, v_vec, initial_R_CM, initial_V_CM, final_R_CM, final_V_CM, final_momentum);

                            star1->mass += Delta_m1;
                            update_stellar_evolution_properties(star1);
                            reset_ODE_mass_dot_quantities(star1);

                            particlesMap->erase(star2->index);

                            *integration_flag = 1;
                            break;
                        }
                    }
                    else if (kw1 <= 10 and kw2 == 11)
                    {
                        if (binary_evolution_SNe_Ia_single_degenerate_model == 0)
                        {
                            if (m2_new - star2->sse_initial_mass >= 0.15)
                            {
                                /* CO and ONeWDs accrete helium-rich material until the accumulated 
                                * material exceeds a mass of 0.15 when it ignites. For a COWD with 
                                * mass less than 0.95 the system will be destroyed as an ELD in a 
                                * possible Type 1a SN. COWDs with mass greater than 0.95 and ONeWDs 
                                * will survive with all the material converted to ONe (JH 30/09/99). 
                                * Now changed to an ELD for all COWDs when 0.15 accreted (JH 11/01/00).  */

                                found_new_event = true;
                                structure_change = true;
                                    
                                Delta_m1 = m1_new - m1;
                                Delta_m2 = -m2;

                                #ifdef VERBOSE
                                if (verbose_flag > 0)
                                {
                                    printf("binary_evolution.cpp -- handle_mass_accretion_events_with_degenerate_objects -- SNe Ia -- star1 %d star2 %d - Delta_m1 %g Delta_m2 %g m1 %g m2 %g m1dot %g m2dot %g m1_new %g m2_new %g\n",star1->index,star2->index,Delta_m1,Delta_m2,m1,m2,mdot1,mdot2,m1_new,m2_new);
                                    print_system(particlesMap,*integration_flag);
                                }
                                #endif

                                #ifdef LOGGING
                                Log_info_type log_info;
                                log_info.index1 = star2->index;
                                log_info.SNe_type = 1;
                                log_info.SNe_info = 1;
                                update_log_data(particlesMap, t, *integration_flag, LOG_SNE_START, log_info);
                                #endif

                                get_initial_binary_orbital_properties_from_position_and_velocity(star1->R_vec, star1->V_vec, star2->R_vec, star2->V_vec, m1, m2, \
                                    r_vec, v_vec, initial_momentum, initial_R_CM, initial_V_CM, h_vec, e_vec);
                            
                                handle_gradual_mass_loss_event_in_system(particlesMap, star1, star2, m1 + Delta_m1, m1, m2 + Delta_m2, m2, parent->compact_object_disruption_mass_loss_timescale, \
                                    r_vec, v_vec, initial_R_CM, initial_V_CM, final_R_CM, final_V_CM, final_momentum);
            
                                star1->mass += Delta_m1;
                                update_stellar_evolution_properties(star1);
                                reset_ODE_mass_dot_quantities(star1);
                                
                                particlesMap->erase(star2->index);

                                *integration_flag = 1;
                                break;
                            }
                        }
                        else if (binary_evolution_SNe_Ia_single_degenerate_model == 1)
                        {
                            /* Model for H or He accretion onto CO WDs */
                            /* Currently, do not distinguish between accumulation of helium, or accumulation of hydrogen which quickly burns to helium */
                            //if (kw1 >= 7 and kw1 <= 9 and kw2 == 11)
                            if (kw1 >= 0 and kw1 <= 10 and kw2 == 11) /* Donor can have hydrogen or helium envelope */
                            {

                                double m_dot_accretion_SD = star2->m_dot_accretion_SD;

                                bool SNe_explosion = determine_if_He_accreting_WD_explodes(m2, m_dot_accretion_SD, star2->WD_He_layer_mass, star2->luminosity);

                                #ifdef VERBOSE
                                if (verbose_flag > 1)
                                {
                                    printf("binary_evolution.cpp -- handle_mass_accretion_events_with_degenerate_objects -- Model for H or He accretion onto CO WDs -- M_WD %g m_dot %g He Mass %g Lum %g\n",m2,m_dot_accretion_SD,star2->WD_He_layer_mass,star2->luminosity,star1->index,star2->index,Delta_m1,Delta_m2,m1,m2,mdot1,mdot2,m1_new,m2_new);
                                    print_system(particlesMap,*integration_flag);
                                }
                                #endif
                                    
                                if (SNe_explosion == true)
                                {
                                    found_new_event = true;
                                    structure_change = true;
                                        
                                    Delta_m1 = m1_new - m1;
                                    Delta_m2 = -m2;

                                    #ifdef VERBOSE
                                    if (verbose_flag > 0)
                                    {
                                        printf("binary_evolution.cpp -- handle_mass_accretion_events_with_degenerate_objects -- Model for H or He accretion onto CO WDs -- SNe Ia -- M_WD %g m_dot %g He Mass %g Lum %g  star1 %d star2 %d - Delta_m1 %g Delta_m2 %g m1 %g m2 %g m1dot %g m2dot %g m1_new %g m2_new %g\n",m2,m_dot_accretion_SD,star2->WD_He_layer_mass,star2->luminosity);
                                        print_system(particlesMap,*integration_flag);
                                    }
                                    #endif

                                    #ifdef LOGGING
                                    Log_info_type log_info;
                                    log_info.index1 = star2->index;
                                    log_info.SNe_type = 1;
                                    log_info.SNe_info = 3;
                                    update_log_data(particlesMap, t, *integration_flag, LOG_SNE_START, log_info);
                                    #endif

                                    get_initial_binary_orbital_properties_from_position_and_velocity(star1->R_vec, star1->V_vec, star2->R_vec, star2->V_vec, m1, m2, \
                                        r_vec, v_vec, initial_momentum, initial_R_CM, initial_V_CM, h_vec, e_vec);
                                
                                    handle_gradual_mass_loss_event_in_system(particlesMap, star1, star2, m1 + Delta_m1, m1, m2 + Delta_m2, m2, parent->compact_object_disruption_mass_loss_timescale, \
                                        r_vec, v_vec, initial_R_CM, initial_V_CM, final_R_CM, final_V_CM, final_momentum);
                
                                    star1->mass += Delta_m1;
                                    update_stellar_evolution_properties(star1);
                                    reset_ODE_mass_dot_quantities(star1);
                                    
                                    particlesMap->erase(star2->index);

                                    *integration_flag = 1;
                                    break;
                                }
                            }
                        }
                    }
                    else if ((kw2 == 10 or kw2 == 11) and m2_new >= chandrasekhar_mass)
                    {
                        
                        /* If the Chandrasekhar limit is exceeded for a white dwarf then destroy
                        * the white dwarf in a supernova. If the WD is ONe then a neutron star
                        * will survive the supernova and we let HRDIAG take care of this when
                        * the stars are next updated. */

                        found_new_event = true;
                        structure_change = true;
                        
                        //dm1 = chandrasekhar_mass - m_accretor + dms_accretor;
                        Delta_m1 = m1_new - m1;
                        Delta_m2 = -m2;

                        #ifdef VERBOSE
                        if (verbose_flag > 0)
                        {
                            printf("binary_evolution.cpp -- handle_mass_accretion_events_with_degenerate_objects -- Chandrasekhar limit of accretor exceeded -- star1 %d star2 %d - Delta_m1 %g Delta_m2 %g m1 %g m2 %g m1dot %g m2dot %g m1_new %g m2_new %g\n",star1->index,star2->index,Delta_m1,Delta_m2,m1,m2,mdot1,mdot2,m1_new,m2_new);
                            print_system(particlesMap,*integration_flag);
                        }
                        #endif

                        #ifdef LOGGING
                        Log_info_type log_info;
                        log_info.index1 = star2->index;
                        log_info.SNe_type = 1;
                        log_info.SNe_info = 1;
                        update_log_data(particlesMap, t, *integration_flag, LOG_SNE_START, log_info);
                        #endif

                        get_initial_binary_orbital_properties_from_position_and_velocity(star1->R_vec, star1->V_vec, star2->R_vec, star2->V_vec, m1, m2, \
                            r_vec, v_vec, initial_momentum, initial_R_CM, initial_V_CM, h_vec, e_vec);
                    
                        handle_gradual_mass_loss_event_in_system(particlesMap, star1, star2, m1 + Delta_m1, m1, m2 + Delta_m2, m2, parent->compact_object_disruption_mass_loss_timescale, \
                            r_vec, v_vec, initial_R_CM, initial_V_CM, final_R_CM, final_V_CM, final_momentum);
    
                        star1->mass += Delta_m1;
                        update_stellar_evolution_properties(star1);
                        reset_ODE_mass_dot_quantities(star1);

                        particlesMap->erase(star2->index);

                        *integration_flag = 1;
                        break;
                    }
                    else if (kw2 == 12 and m2_new > 1.38)
                    {

                        if (ECSNe_model == 1)
                        {
                            /* If an ONeMg WD exceeds 1.38 MSun, assume an electron-capture SNe (ECSN) occurs. */
                        
                            found_new_event = true;
                            structure_change = true;

                            Delta_m1 = m1_new - m1;
                            Delta_m2 = m2_new - m2;

                            #ifdef VERBOSE
                            if (verbose_flag > 0)
                            {
                                printf("binary_evolution.cpp -- handle_mass_accretion_events_with_degenerate_objects -- ONeMg WD exceeds 1.38 MSun - electron-capture SNe -- star1 %d star2 %d - Delta_m1 %g Delta_m2 %g m1 %g m2 %g m1dot %g m2dot %g m1_new %g m2_new %g\n",star1->index,star2->index,Delta_m1,Delta_m2,m1,m2,mdot1,mdot2,m1_new,m2_new);
                                print_system(particlesMap,*integration_flag);
                            }
                            #endif
                           
                            #ifdef LOGGING
                            Log_info_type log_info;
                            log_info.index1 = star2->index;
                            log_info.SNe_type = 3;
                            log_info.SNe_info = 1;
                            update_log_data(particlesMap, t, *integration_flag, LOG_SNE_START, log_info);
                            #endif

                            star2->stellar_type = 13;
                            star2->apply_ECSN_kick = true;

                            star1->instantaneous_perturbation_delta_mass = Delta_m1;
                            star1->instantaneous_perturbation_delta_mass = Delta_m2;
                            star2->apply_kick = true;
                            star2->apply_ECSN_kick == true;

                            bool unbound_orbits;
                            handle_SNe_in_system(particlesMap, &unbound_orbits, integration_flag);

                            update_stellar_evolution_properties(star1);
                            update_stellar_evolution_properties(star2);

                            reset_ODE_mass_dot_quantities(star1);
                            reset_ODE_mass_dot_quantities(star2);
                            
                            if (NS_model == 1)
                            {
                                double ospin;
                                compute_NS_formation_properties_Ye19_model(false, &ospin, &star2->magnetic_field_strength_gauss);
                                rescale_vector(star2->spin_vec, ospin/norm3(star2->spin_vec));
                                star2->time_of_NS_formation = t;
                                star2->initial_NS_period_s = compute_spin_period_from_spin_angular_frequency(ospin) * yr_to_s;
                                star2->initial_magnetic_field_strength_gauss = star2->magnetic_field_strength_gauss;
                            }
                            
                            *integration_flag = 1;
                            break;
                        }
                    }
                }
            }
        }
    }
    
    if (structure_change == true)
    {
        *integration_flag = determine_orbits_in_system_using_nbody(particlesMap);
    }

    update_structure(particlesMap, *integration_flag);

    #ifdef VERBOSE
    if (verbose_flag > 2)
    {
        printf("binary_evolution.cpp -- handle_mass_accretion_events_with_degenerate_objects -- end \n");
        print_system(particlesMap,*integration_flag);
    }
    #endif
    
}

void reset_ODE_mass_dot_quantities(Particle *p)
{
    p->RLOF_flag = 0;
    p->apply_kick = false;
    p->mass_dot_wind = 0.0;
    p->mass_dot_wind_accretion = 0.0;
    p->mass_dot_adiabatic_ejection = 0.0;
    p->mass_dot_RLOF = 0.0;
    p->mass_dot_RLOF_triple = 0.0;
    p->radius_dot = 0.0;
    p->ospin_dot = 0.0;
    p->instantaneous_perturbation_delta_mass = 0.0;
    
    p->age_dot = 0.0;
    p->sse_initial_mass_dot = 0.0;
    p->core_mass_dot = 0.0;
    p->core_radius_dot = 0.0;

    p->luminosity_dot = 0.0;
    p->core_radius_dot = 0.0;
    p->convective_envelope_mass_dot = 0.0;
    p->convective_envelope_radius_dot = 0.0;
    p->sse_k2_dot = 0.0;
    
    p->WD_He_layer_mass = 0.0;
}

}
