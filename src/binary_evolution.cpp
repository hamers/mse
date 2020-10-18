/* MSE */

#include "evolve.h"
#include "binary_evolution.h"

extern "C"
{

//struct ktype__ ktype_;

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

    bool stable = check_system_for_dynamical_stability(particlesMap, integration_flag);
    if (stable == false)
    {
        *integration_flag = 1;
    }
    
//    printf("binary_evolution.cpp -- handle_binary_evolution -- done; stable %d integration_flag %d dt_binary_evolution %g\n",stable,*integration_flag,*dt_binary_evolution);
//    print_system(particlesMap,*integration_flag);
    return binary_flag;
}


int handle_wind_accretion(ParticlesMap *particlesMap, double t_old, double t, double *dt_binary_evolution, int *integration_flag)
{
 //   printf("binary_evolution.cpp -- handle_wind_accretion\n");
   
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
        
        //printf("i %d %d %d\n",donor->index,donor->is_binary,donor->is_bound);
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
                exit(-1);
            }
            //printf("binary_evolution.cpp -- handle_wind_accretion %g\n",companion->mass_dot_wind_accretion);
            //printf("binary_evolution.cpp -- handle_wind_accretion %g\n",companion->mass_dot_wind_accretion);
        }
    }
    
    return 0;
}


int handle_mass_transfer(ParticlesMap *particlesMap, double t_old, double t, double *dt_binary_evolution, int *integration_flag)
{
    *dt_binary_evolution = 1.0e100;
    double min_dt = 1.0e0;
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


    std::vector<int> parent_indices,donor_indices,accretor_indices;
    ParticlesMapIterator it_p;
    for (it_p = particlesMap->begin(); it_p != particlesMap->end(); it_p++)
    //for (it = particlesMap->begin(), next_it = it; it != particlesMap->end(); it = next_it)
    {
        //++next_it;
        Particle *donor = (*it_p).second;
        
        //Particle *accretor = (*particlesMap)[donor->sibling];
        //Particle *donor = (*it).second;

        //printf("handle_mass_transfer donor %d is_binary %d is_bound %d donor->RLOF_flag %d\n",donor->index,donor->is_binary,donor->is_bound,donor->RLOF_flag);
        //if (donor->is_binary == false and donor->is_bound == true and donor->evolve_as_star == true and accretor->evolve_as_star == true)
        if (donor->is_binary == false and donor->is_bound == true)
        {
            Particle *accretor = (*particlesMap)[donor->sibling];
            if (donor->object_type == 1 and accretor->object_type == 1)
            {
                donor->mass_dot_RLOF = 0.0; /* zero by default; could be updated below */
                donor->mass_dot_RLOF_triple = 0.0;
                
                if (donor->RLOF_flag == 1 and std::count(parent_indices.begin(), parent_indices.end(), donor->parent) == 0)
                {
                    parent_indices.push_back(donor->parent);
                    donor_indices.push_back(donor->index);
                    accretor_indices.push_back(donor->sibling);
                    //int flag = handle_mass_transfer_cases(particlesMap, donor->parent, donor->index, donor->sibling, integration_flag, t_old, t, dt_binary_evolution, it_p);
                }
            }
        }
    }

    update_structure(particlesMap, *integration_flag);
    
    std::vector<int>::iterator it;
    //for (it = parent_indices.begin(); it != parent_indices.end(); it++)
    int flag = 0;
    for (int i=0; i<parent_indices.size(); i++)
    {
        //printf("HMT %d %d %d\n",parent_indices[i], donor_indices[i], accretor_indices[i]);
        flag = handle_mass_transfer_cases(particlesMap, parent_indices[i], donor_indices[i], accretor_indices[i], integration_flag, t_old, t, dt_binary_evolution);
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
    #ifdef VERBOSE
    if (verbose_flag > 1)
    {
        printf("binary_evolution.cpp -- handle_mass_transfer_cases -- parent_index %d donor_index %d accretor_index %d\n",parent_index,donor_index,accretor_index);
        print_system(particlesMap,*integration_flag);
    }
    #endif
    
    
    Particle *parent = (*particlesMap)[parent_index];
    Particle *donor = (*particlesMap)[donor_index];
    Particle *accretor = (*particlesMap)[accretor_index];
    
    int flag = 0;
    
    double P_orb = compute_orbital_period_from_semimajor_axis(parent->mass,parent->a);
    double M_donor = donor->mass;
    double M_accretor = accretor->mass;

    /* MT & dynamical timescales */
    //double fm = donor->emt_fm;
    
    
    double R_donor = donor->radius;
    double t_dyn_donor = compute_stellar_dynamical_timescale(M_donor, R_donor);

    /* Critical mass ratios. */
    double q = M_donor / M_accretor;
    double q_crit = compute_q_crit_for_common_envelope_evolution(donor->stellar_type, M_donor, donor->core_mass);
    double q_crit_low_mass_donor = 0.695;
    double q_crit_WD_donor = 0.628;

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
        exit(-1);
    }
    
    int kw = donor->stellar_type;

    #ifdef VERBOSE
    if (verbose_flag > 1)
    {
        printf("binary_evolution.cpp -- handle_mass_transfer_cases -- t_MT %g emt_fm %g\n",t_MT,fm);
        print_system(particlesMap,*integration_flag);
    }
    #endif


    /* Different cases */
    if (kw == 0 and q > q_crit_low_mass_donor)
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
        common_envelope_evolution(particlesMap, parent->index, donor->index, accretor->index, t, integration_flag);//, it_p);
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
    /* TO DO: contact cases! -- applies if accretor->in_RLOF=1 */
    

    #ifdef VERBOSE
    if (verbose_flag > 0)
    {
        printf("binary_evolution.cpp -- handle_mass_transfer_cases -- flag %d integration_flag %d\n",flag,*integration_flag);
        print_system(particlesMap,*integration_flag);
    }
    #endif

    
    if (flag == -1)
    {
        printf("binary_evolution.cpp -- fatal error in handle_mass_transfer_cases, parent %d donor %d accretor %d \n",parent->index,donor->index,accretor->index);
        exit(-1);
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

    //double fm = donor->fm;
    //double m_dot = compute_orbit_averaged_mass_transfer_rate_emt_model(m_donor,fm,P_orb);
    
    double dm1; /* Mass lost from companion due to RLOF (>0) */
    double dm2; /* Mass accreted by companion (>0). */
    
    double dt = t - t_old;
    //dm1 = m_donor;
    
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
        //tms_new *= Myr_to_yr;

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
//            double tbgb = tscls[0] * Myr_to_yr; /* old base of GB timescale */
            star_(&kw2,&accretor->sse_initial_mass,&m_accretor_new,&tms_new,&tn_new,tscls_new,lums,GB,zpars); /* will update tscls */
//            tms_new *= Myr_to_yr;
  //          double tbgb_new = tscls[0] * Myr_to_yr;
            
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
            accretor->epoch = t - accretor->age;
        }
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
    
//    double Delta_m1 = -m_donor;
//    double Delta_m2 = dm2;

    double initial_momentum[3],initial_R_CM[3],initial_V_CM[3];
    double h_vec[3],h_vec_unit[3],e_vec[3],e_vec_unit[3],r_vec[3],v_vec[3];

    get_initial_binary_orbital_properties_from_position_and_velocity(donor->R_vec, donor->V_vec, accretor->R_vec, accretor->V_vec, m_donor, m_accretor, \
        r_vec, v_vec, initial_momentum, initial_R_CM, initial_V_CM, h_vec, e_vec);
    
    double final_R_CM[3],final_V_CM[3],final_momentum[3];
    handle_gradual_mass_loss_event_in_system(particlesMap, donor, accretor, 0.0, m_donor, m_accretor + dm2, m_accretor, parent->dynamical_mass_transfer_low_mass_donor_timescale, \
        r_vec, v_vec, initial_R_CM, initial_V_CM, final_R_CM, final_V_CM, final_momentum);

    update_stellar_evolution_properties(accretor);
    accretor->RLOF_flag = 0;
    particlesMap->erase(donor_index);
    
    *integration_flag = 1;

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

    //double fm = donor->fm;
    //double m_dot = compute_orbit_averaged_mass_transfer_rate_emt_model(m_donor,fm,P_orb);
    
    double dm1; /* Mass lost from companion due to RLOF (<0) */
    double dm2; /* Mass accreted by companion (>0).  */

    double dt = t - t_old;
    //dm1 = m_dot * dt;
    dm1 = m_donor;
    
    //double dme = 2.08d-03*eddfac*(1.d0/(1.d0 + zpars(11)))*rad(j2)*tb // TO DO: include case for eddfac<10
    
    dm2 = dm1;

    //double m_new; /* New mass of remnant (if there is one). */
    bool destroyed;;
    double v_kick_vec[3] = {0.0,0.0,0.0};
    
    int kw1 = donor->stellar_type;
    int kw2 = accretor->stellar_type;
    
    if (kw1 == 10 and kw2 == 10)
    {
        /* Assume the energy released by ignition of the triple-alpha reaction is enough to destroy both stars. */
        destroyed = true;
        //m_new = 0.0;
    }
    else if (kw2 >= 10 and kw2 <= 11 and m_accretor > chandrasekhar_mass)
    {
        /* Potentially SNe Ia that destroys the system. */
        destroyed = true;
        //m_new = 0.0;
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
        printf("binary_evolution.cpp -- dynamical_mass_transfer_WD_donor -- unknown case -- parent %d donor %d accretor %d\n",parent_index,donor_index,accretor_index);
        exit(-1);
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
        accretor->RLOF_flag = 0;
        
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

    collision_product(particlesMap, parent_index, donor_index, accretor_index, t, integration_flag);

    return 0;
}


int stable_mass_transfer_evolution(ParticlesMap *particlesMap, int parent_index, int donor_index, int accretor_index, double t_old, double t, int *integration_flag, double *dt_binary_evolution)
{
    //Particle *parent = (*particlesMap)[parent_index];
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
    
//    double fm = donor->fm;
//    double m_dot = compute_orbit_averaged_mass_transfer_rate_emt_model(m_donor,fm,P_orb); /* NOTE: ignoring mass transfer rates as calculated in BSE */
    
    double dt = t - t_old;
    //double dm1 = -m_dot * dt; /* amount of mass lost by donor during dt; NOTE: defined here dm1>0 */


    /* Default value of dm1 -- transferred mass during this time-step */

    double R_RL_av_donor = roche_radius_pericenter_eggleton(rp,q);
    double dm1 = compute_bse_mass_transfer_amount(kw1,m_donor,donor->core_mass,R_donor,R_RL_av_donor,dt,t_dyn_donor,t_KH_donor);
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
                if (kst == 4)
                {
                    /* ignoring the age adjustment -- units seem to be wrong? */
                    mcx = m_accretor;
                }
                else
                {
                    mcx = accretor->core_mass;
                }
        
                double mt2 = m_accretor + (dm2 - dms_accretor);
                double accretor_age;
                gntage_(&mcx, &mt2, &kst, accretor->zpars, &accretor->sse_initial_mass, &accretor_age);
                
                accretor->epoch = t - accretor_age*Myr_to_yr;
            }
        }
    }
    else if (kw1 <= 6 and kw2 >= 10 and kw2 <= 12)
    {
        /* White dwarf secondary. */

        if (dm1/dt < 2.71e-7) 
        {
            if (dm1/dt < 1.03e-7)
            {
                /* Accrete until a nova explosion blows away most of the accreted material. */
                nova = true;

                dm2 = nova_accretion_factor * CV_min(dm1, dme);
            }
            else
            {
                /* Steady burning at the surface */
                dm2 = dm1;
            }
                
        }
        else
        {
            /* Make a new giant envelope. */

            dm2 = dm1;

            int kw;
            if ( (kw2 == 10 and m_accretor < 0.05) or (kw2 >= 11 and m_accretor < 0.5) ) /* Check for planets or low-mass WDs. */
            {
                kw = kw2;
            }
            else
            {
                kw = CV_min(6, 3*kw2 - 27);

                double new_age;
                double mt2 = m_accretor + (dm2 - dms_accretor);
                gntage_(&accretor->core_mass,&mt2,&kw,accretor->zpars,&accretor->sse_initial_mass,&new_age);
                accretor->age = new_age*Myr_to_yr;
                accretor->epoch = t - accretor->age;
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
    
    if (kw2 >= 10 and kw2 <= 12)
    {
        double mt2 = m_accretor + (dm2 - dms_accretor);
        if (kw1 <= 10 and kw2 == 10 and mt2 >= 0.7)
        {
            /* HeWD can only accrete helium-rich material up to a mass of 0.7 when it is destroyed in a possible Type 1a SN. */
            double Delta_m1 = dm1 + dms_donor;
            double Delta_m2 = -m_accretor;
            
            double initial_momentum[3],initial_R_CM[3],initial_V_CM[3];
            double h_vec[3],h_vec_unit[3],e_vec[3],e_vec_unit[3],r_vec[3],v_vec[3];

            get_initial_binary_orbital_properties_from_position_and_velocity(donor->R_vec, donor->V_vec, accretor->R_vec, accretor->V_vec, m_donor, m_accretor, \
                r_vec, v_vec, initial_momentum, initial_R_CM, initial_V_CM, h_vec, e_vec);
            
            double final_R_CM[3],final_V_CM[3],final_momentum[3];
            handle_gradual_mass_loss_event_in_system(particlesMap, donor, accretor, m_donor + Delta_m1, m_donor, m_accretor + Delta_m2, m_accretor, parent->dynamical_mass_transfer_low_mass_donor_timescale, \
                r_vec, v_vec, initial_R_CM, initial_V_CM, final_R_CM, final_V_CM, final_momentum);

            update_stellar_evolution_properties(donor);
            donor->RLOF_flag = 0;
                
            particlesMap->erase(accretor_index);

            *integration_flag = 1;

            #ifdef IGNORE
            handle_instantaneous_and_adiabatic_mass_changes_in_orbit(particlesMap, donor, accretor, Delta_m1, Delta_m2, parent->compact_object_disruption_mass_loss_timescale, integration_flag); /* TO DO: should use different ML timescale here? */

            /* The binary becomes a body with the donor's properties. */
            update_stellar_evolution_properties(donor);
            parent->is_binary = false; 
            copy_all_body_properties(donor, parent);

            //particlesMap->erase(donor_index);
            particlesMap->erase(accretor_index);
            #endif
        }
        else if (kw1 <= 10 and kw2 >= 11)
        {
            /* CO and ONeWDs accrete helium-rich material until the accumulated 
            * material exceeds a mass of 0.15 when it ignites. For a COWD with 
            * mass less than 0.95 the system will be destroyed as an ELD in a 
            * possible Type 1a SN. COWDs with mass greater than 0.95 and ONeWDs 
            * will survive with all the material converted to ONe (JH 30/09/99). 
            * Now changed to an ELD for all COWDs when 0.15 accreted (JH 11/01/00).  */
    
            if (mt2 - accretor->sse_initial_mass >= 0.15)
            {
                if (kw2 == 11)
                {
                    double Delta_m1 = dm1 + dms_donor;
                    double Delta_m2 = -m_accretor;
                    
                    double initial_momentum[3],initial_R_CM[3],initial_V_CM[3];
                    double h_vec[3],h_vec_unit[3],e_vec[3],e_vec_unit[3],r_vec[3],v_vec[3];

                    get_initial_binary_orbital_properties_from_position_and_velocity(donor->R_vec, donor->V_vec, accretor->R_vec, accretor->V_vec, m_donor, m_accretor, \
                        r_vec, v_vec, initial_momentum, initial_R_CM, initial_V_CM, h_vec, e_vec);
                    
                    double final_R_CM[3],final_V_CM[3],final_momentum[3];
                    handle_gradual_mass_loss_event_in_system(particlesMap, donor, accretor, m_donor + Delta_m1, m_donor, m_accretor + Delta_m2, m_accretor, parent->compact_object_disruption_mass_loss_timescale, \
                        r_vec, v_vec, initial_R_CM, initial_V_CM, final_R_CM, final_V_CM, final_momentum);

                    update_stellar_evolution_properties(donor);
                    donor->RLOF_flag = 0;
                    particlesMap->erase(accretor_index);

                    *integration_flag = 1;

                    #ifdef IGNORE
                    handle_instantaneous_and_adiabatic_mass_changes_in_orbit(particlesMap, donor, accretor, Delta_m1, Delta_m2, parent->compact_object_disruption_mass_loss_timescale, integration_flag);

                    /* The binary becomes a body with the donor's properties. */
                    update_stellar_evolution_properties(donor);
                    parent->is_binary = false; 
                    copy_all_body_properties(donor, parent);

                    particlesMap->erase(donor_index);
                    particlesMap->erase(accretor_index);
                    #endif
                }
                else
                {
                    accretor->sse_initial_mass = mt2;
                }
            }
            else
            {
                accretor->sse_initial_mass = mt2;
            }
        }
        if (kw2 == 10 or kw2 == 11)
        {
            if (mt2 >= chandrasekhar_mass)
            {
                /* If the Chandrasekhar limit is exceeded for a white dwarf then destroy
                * the white dwarf in a supernova. If the WD is ONe then a neutron star
                * will survive the supernova and we let HRDIAG take care of this when
                * the stars are next updated. */
                
                dm1 = chandrasekhar_mass - m_accretor + dms_accretor;
                double Delta_m1 = dm1 + dms_donor;
                double Delta_m2 = -m_accretor;
                
                double initial_momentum[3],initial_R_CM[3],initial_V_CM[3];
                double h_vec[3],h_vec_unit[3],e_vec[3],e_vec_unit[3],r_vec[3],v_vec[3];

                get_initial_binary_orbital_properties_from_position_and_velocity(donor->R_vec, donor->V_vec, accretor->R_vec, accretor->V_vec, m_donor, m_accretor, \
                    r_vec, v_vec, initial_momentum, initial_R_CM, initial_V_CM, h_vec, e_vec);
                
                double final_R_CM[3],final_V_CM[3],final_momentum[3];
                handle_gradual_mass_loss_event_in_system(particlesMap, donor, accretor, m_donor + Delta_m1, m_donor, m_accretor + Delta_m2, m_accretor, parent->compact_object_disruption_mass_loss_timescale, \
                    r_vec, v_vec, initial_R_CM, initial_V_CM, final_R_CM, final_V_CM, final_momentum);

                update_stellar_evolution_properties(donor);
                donor->RLOF_flag = 0;
                particlesMap->erase(accretor_index);
        
                #ifdef IGNORE
                handle_instantaneous_and_adiabatic_mass_changes_in_orbit(particlesMap, donor, accretor, Delta_m1, Delta_m2, parent->compact_object_disruption_mass_loss_timescale, integration_flag);

                /* The binary becomes a body with the donor's properties. */
                update_stellar_evolution_properties(donor);
                parent->is_binary = false; 
                copy_all_body_properties(donor, parent);

                particlesMap->erase(donor_index);
                particlesMap->erase(accretor_index);
                #endif
            }
        }
    }

    /* "New" masses used for aging. Note: the actual masses (and radii) will be updated during the ODE integration.
     * Also, by definition, dm1 > 0, dm2 > 0 */
    double m_donor_new = m_donor - dm1 + donor->mass_dot_wind * dt;
    double m_accretor_new = m_accretor + dm2 + accretor->mass_dot_wind * dt;

    if (kw1 <= 1 or kw1 == 7)
    {
        donor->sse_initial_mass = m_donor_new;
    }
    if (kw2 <= 1 or kw2 == 7)
    {
        accretor->sse_initial_mass = m_donor_new;
    }
    
    /* For a HG star check if the initial mass can be reduced. */    
    if (kw1 == 2 and donor->sse_initial_mass <= donor->zpars[2])
    {
        double m0 = m_donor_new;
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
            accretor->age *= (tms/tms_donor_old);
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
        printf("binary_evolution.cpp -- stable_mass_transfer_evolution -- dt_binary_evolution new %g\n",*dt_binary_evolution);
    }
    #endif

    return 0;
}



int triple_stable_mass_transfer_evolution(ParticlesMap *particlesMap, int parent_index, int donor_index, int accretor_index, double t_old, double t, int *integration_flag, double *dt_binary_evolution)
{
    #ifdef VERBOSE
    if (verbose_flag > 1)
    {
        printf("binary_evolution.cpp -- triple_stable_mass_transfer_evolution\n");
    }
    #endif
    

    Particle *outer_binary = (*particlesMap)[parent_index];
    Particle *donor = (*particlesMap)[donor_index];
    Particle *inner_binary = (*particlesMap)[accretor_index];

    /* Check whether all prerequisites are met. */
    if (donor->is_binary == true or inner_binary->is_binary == false)
    {
        printf("binary_evolution.cpp -- triple_stable_mass_transfer_evolution() -- ERROR: donor should be a star, and accretor a binary\n");
        exit(-1);
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
    

//    double P_orb = compute_orbital_period_from_semimajor_axis(parent->mass,parent->a);
    double m_donor = donor->mass;
    double m_inner_binary = inner_binary->mass;
    double m1 = star1->mass;
    double m2 = star2->mass;


    double q = m_donor/m_inner_binary;

    double R_donor = donor->radius;
    //double R_accretor = accretor->radius;

    //double m_old = m_donor + m_inner_binary;

    int kw1 = donor->stellar_type;
    //int kw2 = accretor->stellar_type;

    double t_KH_donor = compute_Kelvin_Helmholtz_timescale(kw1,m_donor,donor->core_mass,R_donor,donor->luminosity);
    //double t_KH_accretor = compute_Kelvin_Helmholtz_timescale(kw2,m_accretor,accretor->core_mass,R_accretor,accretor->luminosity);

    double t_dyn_donor = compute_stellar_dynamical_timescale(m_donor, R_donor);
    double a_out = outer_binary->a;
    double e_out = outer_binary->e;
    double rp_out = a_out*(1.0-e_out);

    double a_in = inner_binary->a;
    double e_in = inner_binary->e;
    double rp_in = a_in*(1.0-e_in);
    
//    double fm = donor->fm;
//    double m_dot = compute_orbit_averaged_mass_transfer_rate_emt_model(m_donor,fm,P_orb); /* NOTE: ignoring mass transfer rates as calculated in BSE */
    
    double dt = t - t_old;
    //double dm1 = -m_dot * dt; /* amount of mass lost by donor during dt; NOTE: defined here dm1>0 */


    /* Default value of dm3 -- transferred mass during this time-step */

    double R_RL_av_donor = roche_radius_pericenter_eggleton(rp_out,q);
    double dm3 = compute_bse_mass_transfer_amount(kw1,m_donor,donor->core_mass,R_donor,R_RL_av_donor,dt,t_dyn_donor,t_KH_donor);
    
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
    
    double a_in_f = a_in * ( (m1 + dm1)*(m2 + dm2) / (m1*m2 + 2.0*(m1+m2)*dm3/triple_mass_transfer_inner_binary_alpha_times_lambda) );
    inner_binary->triple_mass_transfer_a_in_dot = (a_in_f - a_in)/dt;
    
    /* Effect of mass loss of triple subsystem on the rest of the system. */
    /* Assume the mass not accreted by the inner binary is lost in an adiabatic wind. */
    inner_binary->mass_dot_wind = -( dm3 - (dm1 + dm2) )/dt;
 
    return 0;
}
 
double compute_bse_mass_transfer_amount(int kw1, double m_donor, double core_mass_donor, double R_donor, double R_RL_av_donor, double dt, double t_dyn_donor, double t_KH_donor)
{
    double m_fac = CV_min(m_donor, 5.0);
    m_fac *= m_fac;

    if (R_donor < R_RL_av_donor)
    {
        
        #ifdef VERBOSE
        if (verbose_flag > 0)
        {
            printf("binary_evolution.cpp -- compute_bse_mass_transfer_amount -- skipping since R_donor < R_RL_av -- R_donor %g RL_av_donor %g\n",R_donor,R_RL_av_donor);
        }
        #endif
        
        return 0.0;
    }
    
    double log_fac = log(R_donor/R_RL_av_donor);
    double dm1 = dt * fabs(3.0e-6 * m_fac * log_fac*log_fac*log_fac); /* HPT eq. 58-59 */

    if (kw1 == 2)
    {
        double mew = (m_donor - core_mass_donor)/m_donor;
        dm1 = CV_max(mew,0.01) * dm1;
    }
    else if (kw1 == 10)
    {
        //dm1 = 1.0e3*dm1 / CV_max(R_donor/CONST_R_SUN, 1.0e-4);
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
    if (verbose_flag > 1)
    {
        printf("binary_evolution.cpp -- compute_bse_mass_transfer_amount -- dm1 %g R_donor %g RL_av_donor %g\n",dm1,R_donor,R_RL_av_donor);
    }
    #endif

    return dm1;
}


}
