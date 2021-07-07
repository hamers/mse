/* MSE */

#include "evolve.h"
#include "stellar_evolution.h"

extern "C"
{

struct value1__ value1_;
struct value2__ value2_;
struct value3__ value3_;
struct value4__ value4_;
struct value5__ value5_;
struct flags__ flags_;
struct points__ points_;
struct sse_error_output__ sse_error_output_;

void check_sse_error_codes()
{
    if (sse_error_output_.sse_error_code > 0)
    {
        error_code = sse_error_output_.sse_error_code;
        printf("stellar_evolution.cpp -- check_sse_error_codes() -- nonzero SSE error code %d detected\n",error_code);
        longjmp(jump_buf,1);
    }
}

int initialize_stars(ParticlesMap *particlesMap)
{

    value1_.neta = 0.5;
    value1_.bwind = 0.0;
    value1_.hewind = 0.5;
    value1_.mxns = 3.0;
    value2_.alpha1 = 1.0;
    value2_.lambda = 1.0;
    value3_.idum = 0;
    value4_.sigma = 190.0;
    value4_.bhflag = 1;
    //        value5_.beta = 0.0;
    //        value5_.xi = 0.0;
    //        value5_.acc2 = 0.0;
    //        value5_.epsnov = 0.0;
    //        value5_.eddfac = 0.0;
    //        value5_.gamma = 0.0;
    flags_.ceflag = binary_evolution_CE_energy_flag;
    flags_.tflag = 0;
    flags_.ifflag = 0;
    flags_.nsflag = 1;
    flags_.wdflag = 1;
    points_.pts1 = 0.05;
    points_.pts2 = 0.01;
    points_.pts3 = 0.02;
    sse_error_output_.sse_error_code = 0;
    
    double sse_initial_mass,mass,mt;
    double z;

    set_up_derived_quantities(particlesMap); /* compute a & e for all orbits; used at the end below to adjust h vectors such that possible initial stellar evolution would not change the orbital elements  */

    int i;
    int kw,kw_desired;
    double tm,tn;
    double age;
    double *GB,*tscls,*lums;
    double r,lum,mc,rc,menv,renv,k2;
            
           
    ParticlesMapIterator it_p;
    for (it_p = particlesMap->begin(); it_p != particlesMap->end(); it_p++)
    {
        Particle *p = (*it_p).second;
        if (p->is_binary == false and p->object_type == 1)
        {
            z = p->metallicity;
            
            if (z < 0.0001 or z > 0.03)
            {
                printf("stellar_evolution.cpp -- ERROR: metallicity (given: z = %g for star with index %d and initial mass %g MSun) should be in the range 0.0001 < z < 0.03; terminating program\n",z,p->index,p->mass);
                //exit(-1);
                error_code = 17;
                longjmp(jump_buf,1);
            }

            double *zpars;
            zpars = new double[20];
            //double zpars[20];
            zcnsts_(&z,zpars);
            p->zpars = zpars;
            
            kw_desired = p->stellar_type;
            sse_initial_mass = p->sse_initial_mass;
            mt = p->mass;
            age = 0.0;

            GB = new double[10];
            tscls = new double[20];
            lums = new double[10];    
            
//            star_(&kw, &mass, &mt, &tm, &tn, tscls, lums, GB, zpars);
//            hrdiag_(&mass,&age,&mt,&tm,&tn,tscls,lums,GB,zpars,
//                &r,&lum,&kw,&mc,&rc,&menv,&renv,&k2);
            
            double dtp=0.0;
            double ospin=0.0;
            double epoch = 0.0;
            double tms;
            double tphys = 0.0;
            double tphysf = 0.0;;
            double dt;

            if (kw_desired == 0 or kw_desired == 1) /* ZAMS star; no aging */
            {
                kw = kw_desired;
                evolv1_(&kw,&sse_initial_mass,&mt,&r,&lum,&mc,&rc,&menv,&renv,&ospin,&epoch,&tms,&tphys,&tphysf,&dtp,&z,zpars,&k2);
            }
            else /* age star until reaching desired stellar type */
            {
                double dtm,dtr;
                kw = 1;
                int j=0;
                while (kw < kw_desired)
                {
                    #ifdef VERBOSE
                    if (verbose_flag > 1)
                    {
                        printf("stellar_evolution.cpp -- initialise stars -- input %d kw minit %g mt %g tphys %g tphysf %g\n",kw,p->sse_initial_mass,mt,tphys,tphysf);
                    }
                    #endif

                    evolv1_(&kw,&sse_initial_mass,&mt,&r,&lum,&mc,&rc,&menv,&renv,&ospin,&epoch,&tms,&tphys,&tphysf,&dtp,&z,zpars,&k2);
                    check_sse_error_codes();
                    age = tphysf - epoch;

                    dt = get_new_dt_sse(kw,sse_initial_mass,mt,age,dt,zpars);
    
                    tphys = tphysf;
                    tphysf += dt;
                    j++;
                    
                    if (j>10000)
                    {
                        #ifdef VERBOSE
                        if (verbose_flag > 1)
                        {
                            printf("stellar_evolution.cpp -- unable to evolve star with index %d to desired stellar type %d; current stellar type %d; stopping initial evolution lum %g\n",p->index,kw_desired,kw,lum);
                        }
                        #endif

                        break;
                    }
                }
            }

            if (kw == 13)
            {
                if (NS_model == 1)
                {
                    compute_NS_formation_properties_Ye19_model(false, &ospin, &p->magnetic_field_strength_gauss);
                    rescale_vector(p->spin_vec, ospin/norm3(p->spin_vec));
                    
                    p->time_of_NS_formation = 0.0;
                    p->initial_NS_period_s = compute_spin_period_from_spin_angular_frequency(ospin) * yr_to_s;
                    p->initial_magnetic_field_strength_gauss = p->magnetic_field_strength_gauss;
                }
            }

            #ifdef VERBOSE
            if (verbose_flag > 1)
            {
                printf("stellar_evolution.cpp -- initialize_star -- done evolve; kw %d sse_initial_mass %g mt %g\n",kw,sse_initial_mass,mt);
            }
            #endif

            p->stellar_type = kw;

            p->sse_initial_mass = sse_initial_mass;
            p->mass = mt;
            p->radius = r*CONST_R_SUN;

            p->age = age*Myr_to_yr;
            p->sse_main_sequence_timescale = tms*Myr_to_yr;
            p->epoch = epoch*Myr_to_yr;

            p->luminosity = lum*CONST_L_SUN;
            p->core_mass = mc;
            p->core_radius = rc*CONST_R_SUN;
            p->convective_envelope_mass = menv;
            p->convective_envelope_radius = renv*CONST_R_SUN;

            p->sse_k2 = k2;
            p->sse_k3 = 0.21;

            /* Set up spins parallel with parent orbit (if the star has a parent) */
            if (p->parent != -1)
            {
                Particle *parent = (*particlesMap)[p->parent];
                double *h_vec = parent->h_vec;
                double h = norm3(h_vec);
                
                for (i=0; i<3; i++)
                {
                    p->spin_vec[i] = ospin*h_vec[i]/h;
                }
            }
            else
            {
                p->spin_vec[0] = 0.0;
                p->spin_vec[1] = 0.0;
                p->spin_vec[2] = ospin;
            }

            p->check_for_RLOF_at_pericentre = true;
            
            check_number(p->stellar_type,                  "stellar_evolution.cpp -- initialize_stars","stellar_type", true);
            check_number(p->sse_initial_mass,              "stellar_evolution.cpp -- initialize_stars","sse_initial_mass", true);
            check_number(p->mass,                          "stellar_evolution.cpp -- initialize_stars","mass", true);
            check_number(p->radius,                        "stellar_evolution.cpp -- initialize_stars","radius", true);
            check_number(p->age,                           "stellar_evolution.cpp -- initialize_stars","age", true);
            check_number(p->sse_main_sequence_timescale,   "stellar_evolution.cpp -- initialize_stars","sse_main_sequence_timescale", true);
            check_number(p->epoch,                         "stellar_evolution.cpp -- initialize_stars","epoch", true);
            check_number(p->luminosity,                    "stellar_evolution.cpp -- initialize_stars","luminosity", true);
            check_number(p->core_mass,                     "stellar_evolution.cpp -- initialize_stars","core_mass", true);
            check_number(p->core_radius,                   "stellar_evolution.cpp -- initialize_stars","core_radius", true);
            check_number(p->convective_envelope_mass,      "stellar_evolution.cpp -- initialize_stars","convective_envelope_mass", true);
            check_number(p->convective_envelope_radius,    "stellar_evolution.cpp -- initialize_stars","convective_envelope_radius", true);
            check_number(p->sse_k2,                        "stellar_evolution.cpp -- initialize_stars","sse_k2", true);
            check_number(p->sse_k3,                        "stellar_evolution.cpp -- initialize_stars","sse_k3", true);
            check_number(p->spin_vec[0],                   "stellar_evolution.cpp -- initialize_stars","spin_vec[0]", true);
            check_number(p->spin_vec[1],                   "stellar_evolution.cpp -- initialize_stars","spin_vec[1]", true);
            check_number(p->spin_vec[2],                   "stellar_evolution.cpp -- initialize_stars","spin_vec[2]", true);

            delete[] GB;
            delete[] tscls;
            delete[] lums;

        }
    }



    update_structure(particlesMap, 0); /* set new binary masses */

    /* Adjust h vectors such that possible initial stellar evolution would not change the original orbital elements  */
    for (it_p = particlesMap->begin(); it_p != particlesMap->end(); it_p++)
    {
        Particle *p = (*it_p).second;
        if (p->is_binary == true)
        {
            double h = (p->child1_mass*p->child2_mass*sqrt(CONST_G*p->a/(p->child1_mass + p->child2_mass)))*sqrt(1.0 - p->e * p->e); /* new h based on new masses, but old a & e */
            for (i=0; i<3; i++)
            {
                p->h_vec[i] = p->h_vec_unit[i] * h;
            }
        }
    }

    return 0;
}

int evolve_stars(ParticlesMap *particlesMap, double start_time, double end_time, double *stellar_evolution_timestep, bool get_timestep_only, bool *apply_SNe_effects, int *integration_flag)
{
    /* Go through all stars in the system and evolve them with SSE */
        
    /* SSE global parameters that need to be set */
    value1_.neta = 0.5;
    value1_.bwind = 0.0;
    value1_.hewind = 0.5;
    value1_.mxns = 3.0;
    value2_.alpha1 = 1.0;
    value2_.lambda = 1.0;
    value3_.idum = 0;
    value4_.sigma = 190.0;
    value4_.bhflag = 1;
    //        value5_.beta = 0.0;
    //        value5_.xi = 0.0;
    //        value5_.acc2 = 0.0;
    //        value5_.epsnov = 0.0;
    //        value5_.eddfac = 0.0;
    //        value5_.gamma = 0.0;
    flags_.ceflag = binary_evolution_CE_energy_flag;
    flags_.tflag = 0;
    flags_.ifflag = 0;
    flags_.nsflag = 1;
    flags_.wdflag = 1;
    points_.pts1 = 0.05;
    points_.pts2 = 0.01;
    points_.pts3 = 0.02;
    sse_error_output_.sse_error_code = 0;


    ParticlesMapIterator it_p;
    std::vector<int>::iterator it_parent_p,it_parent_q;

    double sse_initial_mass,sse_initial_mass_old,sse_time_step;
    double mt,tphys,tphysf,dtp,epoch,epoch_old,ospin,ospin_old,age,age_old;
    double dt,dt_min,one_div_dt;
    double r,lum,lum_old,mc,rc,rc_old,menv,menv_old,renv,renv_old,tms,k2,k2_old;
    double r_old,mt_old;
    double rzams,fac,menv_fraction;
    int kw,kw_old;
    double spin_vec[3];
    bool update_other_quantities_immediately;
    dt_min = 1.0e100;
    
    double z;
    double *zpars;
    
    *apply_SNe_effects = false;
            
    for (it_p = particlesMap->begin(); it_p != particlesMap->end(); it_p++)
    {
        Particle *p = (*it_p).second;
        if (p->is_binary == false and p->object_type == 1)
        {
            /* Extract stellar evolution parameters */
            /* Note the different units used in SSE vs. MSE */
            mt = mt_old = p->mass;
            sse_initial_mass = sse_initial_mass_old = p->sse_initial_mass;
            kw = kw_old = p->stellar_type;
            r = r_old = p->radius/CONST_R_SUN;
            
            z = p->metallicity;
            zpars = p->zpars;
            ospin_old = ospin = norm3(p->spin_vec);

            sse_time_step = p->sse_time_step*yr_to_Myr;
            age = age_old = p->age*yr_to_Myr;

            p->core_mass_old = p->core_mass;

            tphys = start_time*yr_to_Myr;
            tphysf = end_time*yr_to_Myr;
            double desired_tphysf = tphysf;

            epoch = epoch_old = tphys - age;

            lum_old = p->luminosity/CONST_L_SUN;
            rc_old = p->core_radius/CONST_R_SUN;
            renv_old = p->convective_envelope_radius/CONST_R_SUN;
            menv_old = p->convective_envelope_mass;
            k2_old = p->sse_k2;
           
            dtp = 0.0;
            
            #ifdef VERBOSE
            if (verbose_flag > 1)
            {
                printf("stellar_evolution.cpp -- evolve_stars -- index %d -- prior to evolv1_() call --  kw %d mt %.15f sse_initial_mass %.15f r %.15f sse_time_step %.15f epoch %.15f age %.15f tphys %.15f tphysf %.15f ospin %.15f tms %.15f epoch %.15f k2 %.15f rc %.15f mc %.15f menv %.15f renv %.15f \n",p->index,kw,mt,sse_initial_mass,r,sse_time_step,epoch,age,tphys,tphysf,ospin,tms,epoch,k2_old,rc_old,p->core_mass,menv_old,renv_old);
            }
            #endif

            if (get_timestep_only == true)
            {
                age = tphys - epoch;
                sse_time_step = get_new_dt_sse(kw,sse_initial_mass,mt,age,sse_time_step,zpars);
            }
            else
            {
                evolv1_(&kw,&sse_initial_mass,&mt,&r,&lum,&mc,&rc,&menv,&renv,&ospin,&epoch,&tms,&tphys,&tphysf,&dtp,&z,zpars,&k2);
                check_sse_error_codes();

                if ( fabs(tphysf - desired_tphysf)/desired_tphysf > epsilon)
                {
                    printf("stellar_evolution.cpp -- Warning tphysf != desired_tphysf -- kw %d mt %g r %g lum %g sse_time_step %g epoch %g age %g  tphys %g tphysf %g ospin %g\n",kw,mt,r,lum,sse_time_step,epoch,age,tphys,tphysf,ospin);
                }

                age = tphysf - epoch;
                sse_time_step = get_new_dt_sse(kw,sse_initial_mass,mt,age,sse_time_step,zpars);

                #ifdef VERBOSE
                if (verbose_flag > 1)
                {
                    printf("stellar_evolution.cpp -- evolve_stars -- index %d -- post evolv1_() call --  kw %d mt %.15f sse_initial_mass %.15f r %.15f sse_time_step %.15f epoch %.15f age %.15f tphys %.15f tphysf %.15f ospin %.15f tms %.15f epoch %.15f k2 %.15f rc %.15f mc %.15f menv %.15f renv %.15f \n",p->index,kw,mt,sse_initial_mass,r,sse_time_step,epoch,age,tphys,tphysf,ospin,tms,epoch,k2,rc,mc,menv,renv);
                }
                #endif
                
                if (ospin <= epsilon)
                {
                    ospin = epsilon;
                }

                /* Write new stellar evolution parameters */
                dt = end_time - start_time;                
                one_div_dt = 1.0/dt;

                if (kw_old==kw) 
                /* No stellar type change; let mass, radius and spin change smoothly during ODE integration *
                 * Other parameters are updated only after ODE integration (which could involve a stopping condition, i.e., 
                 * the actual timestep taken could be shorter than that assumed during stellar evolution, hence the need to use the
                 * (first-order) time derivatives, instead of updating the quantities directly here. */
                {
                    p->mass_dot_wind = (mt - mt_old)*one_div_dt;
                    p->radius_dot = CONST_R_SUN*(r - r_old)*one_div_dt;
                    p->ospin_dot = (ospin - ospin_old)*one_div_dt;

                    p->instantaneous_perturbation_delta_mass = 0.0;
                    p->apply_kick = false;
                    
                    if (kw == 13 or kw == 14)
                    {
                        /* Ignore SSE changes in spins for NS/BHs */
                        p->ospin_dot = 0.0;
                    }

                    update_other_quantities_immediately = false;

                    p->age_dot = Myr_to_yr * (age - age_old)*one_div_dt;
                    //p->epoch_dot = Myr_to_yr * (epoch - epoch_old)*one_div_dt;
                    p->sse_initial_mass_dot = (sse_initial_mass - sse_initial_mass_old)*one_div_dt;
                    p->core_mass_dot = (mc - p->core_mass_old)*one_div_dt;
                    p->core_radius_dot = CONST_R_SUN*(rc - rc_old)*one_div_dt;
                    
                    p->luminosity_dot = CONST_L_SUN * (lum - lum_old)*one_div_dt;
                    p->convective_envelope_mass_dot = (menv - menv_old)*one_div_dt;
                    p->convective_envelope_radius_dot = CONST_R_SUN*(renv - renv_old)*one_div_dt;
                    p->sse_k2_dot = (k2 - k2_old)*one_div_dt;

                }
                else if (kw < 13)
                /* kw change not involving NS/BH *
                 * Let mass, radius and spin change smoothly (mass change will be handled assuming
                 * adiabatic mass loss), but update other parameters immediately to avoid potential problems
                 * with binary_evolution.cpp, common_envelope_evolution.cpp, collision.cpp, ODE_tides.cpp */
                {
                    p->mass_dot_wind = (mt - mt_old)*one_div_dt;
                    p->radius_dot = CONST_R_SUN*(r - r_old)*one_div_dt;
                    p->ospin_dot = (ospin - ospin_old)*one_div_dt;
                    
                    p->instantaneous_perturbation_delta_mass = 0.0;
                    p->apply_kick = false;
                    
                    update_other_quantities_immediately = true;
                }
                
                else 
                /* kw change with new kw=13 or kw=14, i.e., NS/BH formation *
                 * Update all parameters immediately, except for the mass, which will be handled
                 * immediately after this function in SNe.cpp */
                {
                    p->mass_dot_wind = 0.0;
                    p->radius_dot = 0.0;
                    p->ospin_dot = 0.0;

                    /* New mass will be handled by handle_SNe_in_system(), but radius and spin have to be updated here */
                    p->radius = r*CONST_R_SUN;
                    if (kw == 15) // SSE may return invalid spins for massless remnants; override with something arbitrary (the body will be removed later anyway in SNe.cpp -- remove_massless_remnants_from_system()
                    {
                        ospin = 1.0;
                    }
                    rescale_vector(p->spin_vec, ospin/ospin_old); /* update spin (scalar change only) */
                    
                    *apply_SNe_effects = true;
                    
                    p->apply_kick = true;
                    p->instantaneous_perturbation_delta_mass = mt - mt_old; /* Will be used by handle_SNe_in_system() */
                    p->RLOF_flag = 0;

                    update_other_quantities_immediately = true;
                }

                if (update_other_quantities_immediately == true)
                /* Update the "other" (everything except mass, radius and spin) stellar evolution quantities immediately *
                 * The corresponding "dots" should be zero to avoid overupdating */
                {
                    p->age_dot = 0.0;
                    p->sse_initial_mass_dot = 0.0;
                    p->core_mass_dot = 0.0;
                    p->core_radius_dot = 0.0;

                    p->luminosity_dot = 0.0;
                    p->convective_envelope_mass_dot = 0.0;
                    p->convective_envelope_radius_dot = 0.0;
                    p->sse_k2_dot = 0.0;
                    
                    p->stellar_type = kw;

                    p->age = age*Myr_to_yr;
                    p->epoch = epoch*Myr_to_yr;
                    p->sse_initial_mass = sse_initial_mass;
                    p->core_mass = mc;
                    p->stellar_type = kw;
                        
                    p->luminosity = lum*CONST_L_SUN;
                    p->core_radius = rc*CONST_R_SUN;
                    p->convective_envelope_mass = menv;
                    p->convective_envelope_radius = renv*CONST_R_SUN;
                    p->sse_k2 = k2;
                }
                p->mass_dot_wind_accretion = 0.0; /* need to reset this since it is only determined in binary_evolution.cpp in secular integration mode, so errorenous results could occur when the integration flag switches */
                
                /* WD kicks */
                if (kw_old!=kw and kw >= 10 and kw <= 12 and p->include_WD_kicks == true)
                {
                    p->apply_kick = true;
                    *apply_SNe_effects = true;
                    p->RLOF_flag = 0;
                    
                    #ifdef VERBOSE
                    if (verbose_flag > 0)
                    {
                        printf("stellar_evolution.cpp -- will apply WD kick to star %d with old stellar type %d and new stellar type %d\n",p->index,kw_old,kw);
                    }
                    #endif
                }
                
                /* Special treatment for NS formation */
                if (kw_old!=kw and kw == 13)
                {
                    if (NS_model == 1)
                    {
                        compute_NS_formation_properties_Ye19_model(false, &ospin, &p->magnetic_field_strength_gauss);
                        rescale_vector(p->spin_vec, ospin/norm3(p->spin_vec));
                        p->time_of_NS_formation = end_time;
                        p->initial_NS_period_s = compute_spin_period_from_spin_angular_frequency(ospin) * yr_to_s;
                        p->initial_magnetic_field_strength_gauss = p->magnetic_field_strength_gauss;
                    }
                }

                p->sse_main_sequence_timescale = tms*Myr_to_yr;
                p->sse_k3 = 0.21;


                check_number(p->stellar_type,                  "stellar_evolution.cpp -- evolve_stars","stellar_type", true);
                check_number(p->sse_initial_mass,              "stellar_evolution.cpp -- evolve_stars","sse_initial_mass", true);
                check_number(p->mass,                          "stellar_evolution.cpp -- evolve_stars","mass", true);
                check_number(p->radius,                        "stellar_evolution.cpp -- evolve_stars","radius", true);
                check_number(p->age,                           "stellar_evolution.cpp -- evolve_stars","age", true);
                check_number(p->sse_main_sequence_timescale,   "stellar_evolution.cpp -- evolve_stars","sse_main_sequence_timescale", true);
                check_number(p->epoch,                         "stellar_evolution.cpp -- evolve_stars","epoch", true);
                check_number(p->luminosity,                    "stellar_evolution.cpp -- evolve_stars","luminosity", true);
                check_number(p->core_mass,                     "stellar_evolution.cpp -- evolve_stars","core_mass", true);
                check_number(p->core_radius,                   "stellar_evolution.cpp -- evolve_stars","core_radius", true);
                check_number(p->convective_envelope_mass,      "stellar_evolution.cpp -- evolve_stars","convective_envelope_mass", true);
                check_number(p->convective_envelope_radius,    "stellar_evolution.cpp -- evolve_stars","convective_envelope_radius", true);
                check_number(p->sse_k2,                        "stellar_evolution.cpp -- evolve_stars","sse_k2", true);
                check_number(p->sse_k3,                        "stellar_evolution.cpp -- evolve_stars","sse_k3", true);
                check_number(p->spin_vec[0],                   "stellar_evolution.cpp -- evolve_stars","spin_vec[0]", true);
                check_number(p->spin_vec[1],                   "stellar_evolution.cpp -- evolve_stars","spin_vec[1]", true);
                check_number(p->spin_vec[2],                   "stellar_evolution.cpp -- evolve_stars","spin_vec[2]", true);
                check_number(p->mass_dot_wind,                   "stellar_evolution.cpp -- evolve_stars","p->mass_dot_wind", true);
                check_number(p->radius_dot,                   "stellar_evolution.cpp -- evolve_stars","p->radius_dot", true);
                check_number(p->ospin_dot,                   "stellar_evolution.cpp -- evolve_stars","p->ospin_dot", true);
                
                /* Logging */
                #ifdef LOGGING
                if (kw != kw_old)
                {
                    Log_info_type log_info;
                    log_info.index1 = p->index;

                    if (kw < 13)
                    {
                        update_log_data(particlesMap, start_time, *integration_flag, LOG_ST_CHANGE, log_info);
                        
                        if (kw >= 10 and kw <= 12 and p->include_WD_kicks == true)
                        {
                            update_log_data(particlesMap, start_time, *integration_flag, LOG_WD_KICK_START, log_info);
                        }                        
                    }
                    else
                    {
                        
                        update_log_data(particlesMap, start_time, *integration_flag, LOG_SNE_START, log_info);
                    }
                }
                if (NS_model == 1 and p->stellar_type == 13)
                {
                    if (determine_if_NS_is_MSP(norm3(p->spin_vec),p->magnetic_field_strength_gauss) == true)
                    {
                        if (p->has_formed_MSP == false) /* Only consider first-time MSP formation */
                        {
                            Log_info_type log_info2;
                            log_info2.index1 = p->index;

                            update_log_data(particlesMap, start_time, *integration_flag, LOG_MSP_FORMATION, log_info2);
                            p->has_formed_MSP = true;
                        }
                    }
                }
                #endif

            }
            p->sse_time_step = sse_time_step*Myr_to_yr;

            /* Common time step */
            if (p->sse_time_step < dt_min)
            {
                dt_min = p->sse_time_step;
            }

        }
    }
    
    *stellar_evolution_timestep = dt_min;
    
    check_for_critical_rotation(particlesMap);
    
    return 0;
}

void update_stellar_evolution_quantities_directly_secular(ParticlesMap *particlesMap, double dt_assumed, double t_old, double t_out)
{
    ParticlesMapIterator it_p;
    double true_time_step = t_out - t_old;
    double dt_dif = dt_assumed - true_time_step;

    #ifdef VERBOSE
    if (verbose_flag > 2)
    {
        printf("stellar_evolution.cpp -- update_stellar_evolution_quantities_directly_secular -- start -- dt_assumed %g t_old %g t_out %g true_time_step %g dt_dif %g\n",dt_assumed,t_old,t_out,true_time_step,dt_dif);
        print_system(particlesMap,0);
    }
    #endif
        
    for (it_p = particlesMap->begin(); it_p != particlesMap->end(); it_p++)
    {
        Particle *p = (*it_p).second;
        if (p->is_binary == false and p->object_type == 1)
        {
            p->sse_initial_mass += p->sse_initial_mass_dot * true_time_step;
            p->age += p->age_dot * true_time_step;
            p->epoch = t_out - p->age;

            p->core_mass += p->core_mass_dot * true_time_step;
            p->core_radius += p->core_radius_dot * true_time_step;            

            p->luminosity += p->luminosity_dot * true_time_step;
            p->convective_envelope_mass += p->convective_envelope_mass_dot * true_time_step;
            p->convective_envelope_radius += p->convective_envelope_radius_dot * true_time_step;
            p->sse_k2 += p->sse_k2_dot * true_time_step;
        }
    }

    #ifdef VERBOSE
    if (verbose_flag > 2)
    {
        printf("stellar_evolution.cpp -- update_stellar_evolution_quantities_directly_secular -- end -- dt_assumed %g t_old %g t_out %g true_time_step %g dt_dif %g\n",dt_assumed,t_old,t_out,true_time_step,dt_dif);
        print_system(particlesMap,0);
    }
    #endif

}

void update_stellar_evolution_quantities_directly_nbody(ParticlesMap *particlesMap, double t, double dt)
{
    int i;
    double spin_vec_norm;
    
    #ifdef VERBOSE
    if (verbose_flag > 2)
    {
        printf("stellar_evolution.cpp -- update_stellar_evolution_quantities_directly_nbody -- start t %g dt %g\n",t,dt);
        print_system(particlesMap,1);
    }
    #endif
    
    ParticlesMapIterator it_p;
    for (it_p = particlesMap->begin(); it_p != particlesMap->end(); it_p++)
    {
        Particle *p = (*it_p).second;
        if (p->is_binary == false and p->object_type == 1)
        {
            p->RLOF_flag = 0;
            
            p->mass += p->mass_dot_wind * dt;
            p->radius += p->radius_dot * dt;

            p->sse_initial_mass += p->sse_initial_mass_dot * dt;
            p->age += p->age_dot * dt;
            p->epoch = t - p->age;

            p->core_mass += p->core_mass_dot * dt;
            p->core_radius += p->core_radius_dot * dt;
            
            p->luminosity += p->luminosity_dot * dt;
            p->convective_envelope_mass += p->convective_envelope_mass_dot * dt;
            p->convective_envelope_radius += p->convective_envelope_radius_dot * dt;
            p->sse_k2 += p->sse_k2_dot * dt;

            spin_vec_norm = norm3(p->spin_vec);
            if (spin_vec_norm <= epsilon)
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
    
    #ifdef VERBOSE
    if (verbose_flag > 2)
    {
        printf("stellar_evolution.cpp -- update_stellar_evolution_quantities_directly_nbody -- end %g dt %g\n",t,dt);
        print_system(particlesMap,1);
    }
    #endif

    
}

void check_for_critical_rotation(ParticlesMap *particlesMap)
{
    int i;
    ParticlesMapIterator it_p;
    for (it_p = particlesMap->begin(); it_p != particlesMap->end(); it_p++)
    {
        
        Particle *p = (*it_p).second;
        if (p->is_binary == false and p->is_bound == true and p->object_type == 1)
        {
            double spin_vec_norm = norm3(p->spin_vec);
            if (spin_vec_norm < epsilon)
            {
                spin_vec_norm = epsilon;
            }
            double Omega_crit = compute_breakup_angular_frequency(p->mass,p->radius);

            if (spin_vec_norm > Omega_crit)
            {

                #ifdef VERBOSE
                if (verbose_flag > 0)
                {
                    printf("stellar_evolution.cpp -- check_for_critical_rotation -- limiting id %d Omega %g Omega_crit %g\n",p->index,spin_vec_norm,Omega_crit);
                }
                #endif
                
                for (i=0; i<3; i++)
                {
                    p->spin_vec[i] *= Omega_crit/spin_vec_norm;
                }
            }
        }
    }
}

double get_new_dt_sse(int kw, double mass, double mt, double age, double dt, double *zpars)
{
    #ifdef VERBOSE
    if (verbose_flag > 3)
    {
        printf("stellar_evolution.cpp -- get_new_dt_sse -- kw %d mass0 %g mt %g age %g dt %g\n",kw,mass,mt,age,dt);
    }
    #endif

    double dtm,dtr;
    double tm,tn;

    double *GB,*tscls,*lums;
    GB = new double[10];
    tscls = new double[20];
    lums = new double[10];    

    star_(&kw, &mass, &mt, &tm, &tn, tscls, lums, GB, zpars);
    deltat_(&kw,&age,&tm,&tn,tscls,&dt,&dtr);

    dtm = CV_min(dtr, dt);
    dtm = CV_max(dtm,1.0d-07*age);

    dtm = CV_max(dtm,1.0e-6);
    
    delete[] GB;
    delete[] tscls;
    delete[] lums;
    
    return dtm;
}


void update_stellar_evolution_properties(Particle *p)
{
    /* hrdiag: 
     * Computes the new mass, luminosity, radius & stellar type.
    *  Input (MASS, AJ, TM, TN, LUMS & TSCLS) supplied by routine STAR. */
    
    double tm,tn;

    double *GB,*tscls,*lums;
    GB = new double[10];
    tscls = new double[20];
    lums = new double[10];  
    
    double age = p->age*yr_to_Myr;
    double r,rc,lum,menv,renv,k2;
    star_(&p->stellar_type, &p->sse_initial_mass, &p->mass, &tm, &tn, tscls, lums, GB, p->zpars);
    hrdiag_(&p->sse_initial_mass,&age,&p->mass,&tm,&tn,tscls,lums,GB,p->zpars, \
            &r,&lum,&p->stellar_type,&p->core_mass,&rc,&menv,&renv,&k2);

    p->radius = r*CONST_R_SUN;

    p->age = age*Myr_to_yr;

    p->luminosity = lum*CONST_L_SUN;
    p->core_radius = rc*CONST_R_SUN;
    p->convective_envelope_mass = menv;
    p->convective_envelope_radius = renv*CONST_R_SUN;
    
    delete[] GB;
    delete[] tscls;
    delete[] lums;

}

double compute_moment_of_inertia(int stellar_type, double mass, double core_mass, double radius, double core_radius, double k2, double k3)
{
    //printf("compute_moment_of_inertia k %d m %g mc %g r %g rc %g kw %g k3 %g",stellar_type,mass,core_mass,radius,core_radius,k2,k3);
    if (core_mass > mass) /* This can happen in edge cases when the core mass is updated immediately after stellar type change, but the mass is changed linearly during ODE integration; circumvent by setting mc=c */
    {
        core_mass = mass;
    }
    
    if (stellar_type < 10) // stars
    {
        return k2*(mass - core_mass)*radius*radius + k3*core_mass*core_radius*core_radius;
    }
    else if (stellar_type >= 10 and stellar_type <= 13)
    {
        return 0.4 * mass * radius * radius; // moment of inertia of solid sphere
    }
    else
    {
        printf("stellar_evolution.cpp -- compute_moment_of_inertia -- ERROR: stellar_type should not be %d\n",stellar_type);
        error_code = 38;
        longjmp(jump_buf,1);
        return 0.0;
    }
}

double compute_moment_of_inertia_dot(int stellar_type, double mass, double core_mass, double radius, double core_radius, double k2, double k3, double mass_dot, double radius_dot)
{
    if (stellar_type < 10) // stars
    {
        return k2*mass_dot*radius*radius + 2.0*k2*(mass - core_mass)*radius*radius_dot; /* Approximate; neglects changes in core masses and radii and k2 and k3  */
    }
    else if (stellar_type >= 10 and stellar_type <= 13) // WDs/NSs
    {
        return 0.4*radius * (mass_dot*radius + 2.0*mass*radius_dot); // moment of inertia change of solid sphere
    }
    else
    {
        printf("stellar_evolution.cpp -- compute_moment_of_inertia_dot -- ERROR: stellar_type should not be %d\n",stellar_type);
        error_code = 38;
        longjmp(jump_buf,1);
        return 0.0;
    }
}


void get_core_masses_by_composition(int kw, double core_mass, double *He_core_mass, double *CO_core_mass, double *Ne_core_mass)
{
    *He_core_mass = 0.0;
    *CO_core_mass = 0.0;
    *Ne_core_mass = 0.0;
    
    if (kw <= 3 or kw == 10)
    {
        *He_core_mass += core_mass;
    }
    else if (kw == 12)
    {
        *Ne_core_mass += core_mass;
    }
    else
    {
        *CO_core_mass += core_mass;
    }
}

double compute_spin_angular_momentum_from_spin_frequency(double spin_frequency, int stellar_type, int object_type, double mass, double core_mass, double radius, double core_radius, double k2, double k3)
{
    double S;
    if (object_type != 1)
    {
        #ifdef VERBOSE
        if (verbose_flag > 0)
        {
            printf("stellar_evolution.cpp -- compute_spin_angular_momentum_from_spin_frequency -- ST %d OT %d m %g R %g k2 %g\n",stellar_type,  object_type, mass,radius,k2);
        }
        #endif

        double I = k2 * mass * radius * radius;
        S = I * spin_frequency;
    }
    else
    {
        if (stellar_type < 14) // Stars/WDs/NSs
        {
            double I = compute_moment_of_inertia(stellar_type, mass, core_mass, radius, core_radius, k2, k3);
            S = I * spin_frequency;
        }
        else // BHs
        {
            double chi = compute_spin_parameter_from_spin_frequency(mass, spin_frequency);
            S = chi * CONST_G * mass * mass / CONST_C_LIGHT;
        }
    }
    
    check_number(S,"stellar_evolution.cpp -- compute_spin_angular_momentum_from_spin_frequency","S", true);
    return S;
}

double compute_spin_frequency_from_spin_angular_momentum(double spin_angular_momentum, int stellar_type, int object_type, double mass, double core_mass, double radius, double core_radius, double k2, double k3)
{
    double Omega;
    if (object_type != 1)
    {
        double I = k2 * mass * radius * radius;

        Omega = spin_angular_momentum/I;
    }
    else
    {
        if (stellar_type < 14) // Stars/WDs/NSs
        {
            double I = compute_moment_of_inertia(stellar_type, mass, core_mass, radius, core_radius, k2, k3);
            Omega = spin_angular_momentum/I;
        }
        else // BHs
        {
            double chi = spin_angular_momentum * CONST_C_LIGHT / ( CONST_G * mass * mass);
            Omega = compute_spin_frequency_from_spin_parameter(mass, chi);
        }
    }
    check_number(Omega,"stellar_evolution.cpp -- compute_spin_frequency_from_spin_angular_momentum","Omega", true);
    return Omega;
}

double compute_breakup_angular_frequency(double mass, double radius)
{
    return sqrt(CONST_G*mass/(radius*radius*radius));
}

double determine_sse_compact_object_radius_RSun(int kw, double m)
{
    double radius;
    if (kw >= 10 and kw <= 12)
    {
        radius = 0.0115 * sqrt( CV_max( 1.48204e-6, pow(chandrasekhar_mass/m,c_2div3) - pow(m/chandrasekhar_mass,c_2div3) ) );
        radius = CV_min(0.1, radius);
        if (m < 0.0005)
        {
            radius = 0.09;
        }
        if (m < 0.000005)
        {
            radius = 0.009;
        }
    }
    else if (kw == 13)
    {
        radius = 1.4e-5;
    }
    else if (kw == 14)
    {
        radius = 4.24e-6*m;
    }
    else
    {
        printf("stellar_evolution.cpp -- determine_sse_compact_object_radius_RSun -- ERROR: this function was called with invalid kw=%d and m=%g\n",kw,m);
        error_code = 37;
        longjmp(jump_buf,1);
    }
    return radius;
}


void compute_NS_formation_properties_Ye19_model(bool merger_event, double *ospin, double *B_G)
{
    /* https://ui.adsabs.harvard.edu/abs/2019ApJ...877..122Y/abstract */
    
    double P_s;
    if (merger_event == false)
    {
        P_s = NS_Ye19_model_NS_formation_single_P_s_lower + generate_random_number_between_zero_and_unity() * (NS_Ye19_model_NS_formation_single_P_s_upper - NS_Ye19_model_NS_formation_single_P_s_lower);
        *B_G = NS_Ye19_model_NS_formation_single_B_G_lower + generate_random_number_between_zero_and_unity() * (NS_Ye19_model_NS_formation_single_B_G_upper - NS_Ye19_model_NS_formation_single_B_G_lower);
    }
    else
    {
        P_s = NS_Ye19_model_NS_formation_merger_P_s_lower + generate_random_number_between_zero_and_unity() * (NS_Ye19_model_NS_formation_merger_P_s_upper - NS_Ye19_model_NS_formation_merger_P_s_lower);
        *B_G = NS_Ye19_model_NS_formation_merger_B_G_lower + generate_random_number_between_zero_and_unity() * (NS_Ye19_model_NS_formation_merger_B_G_upper - NS_Ye19_model_NS_formation_merger_B_G_lower);
    }
    double P_yr = P_s * s_to_yr;
    *ospin = compute_spin_angular_frequency_from_spin_period(P_yr);
}

bool determine_if_NS_is_MSP(double ospin, double B_G)
{
    bool is_MSP = false;
    
    double P_s = compute_spin_period_from_spin_angular_frequency(ospin) * yr_to_s;
    double B_crit_G = pulsar_death_line_B_crit_const_G * P_s * P_s;
    //printf("determine_if_NS_is_MSP %g %g %g %g\n",B_G,B_crit_G,P_s,MSP_defining_period_s);
    if (B_G > B_crit_G and P_s < MSP_defining_period_s)
    {
        is_MSP = true;
    }
    return is_MSP;
}

double compute_NS_magnetic_field_Ye19_model(double B_0, double T, double t_acc, double Delta_m)
{
    //printf("%g %g %g %g\n",B_0,T,t_acc,Delta_m);
    double B_G = B_0 * exp( - (T - t_acc)/NS_Ye19_model_NS_B_field_decay_timescale) / (1.0 + Delta_m/NS_Ye19_model_NS_RLOF_threshold_accreted_mass);
    if (B_G < NS_Ye19_model_NS_minimum_B_G)
    {
        B_G = NS_Ye19_model_NS_minimum_B_G;
    }
    return B_G;
}

void ODE_handle_NS_properties_Ye19_model(Particle *p, double dt)
{
    /* https://ui.adsabs.harvard.edu/abs/2019ApJ...877..122Y/abstract */
    
    if (p->stellar_type == 13)
    {
        double ospin = p->spin_vec_norm; /* should be up to date when called from compute_y_dot() */
        double B_G = compute_NS_magnetic_field_Ye19_model(p->magnetic_field_strength_gauss, dt, p->RLOF_timescale, p->delta_mass_RLOF);
        double ospin_dot = -ONE_DIV_FOURPISQ * NS_Ye19_model_NS_spin_down_constant_K * ospin * ospin * ospin * B_G * B_G;

        for (int i=0; i<3; i++)
        {
            p->dspin_vec_dt[i] += ospin_dot * p->spin_vec_unit[i];
        }
        #ifdef VERBOSE
        if (verbose_flag > 1)
        {
            printf("stellar_evolution.cpp -- ODE_handle_NS_properties_Ye19_model -- dt %g ospin %g ospin_dot %g B_G %g P_s %g\n",dt,ospin,ospin_dot,B_G,compute_spin_period_from_spin_angular_frequency(ospin) * yr_to_s);
        }
        #endif


    }
}

void ODE_Ye19_model_update_magnetic_field(ParticlesMap *particlesMap, double dt)
{
    /* https://ui.adsabs.harvard.edu/abs/2019ApJ...877..122Y/abstract */
    /* Updates magnetic fields of NSs in particlesMap after ODE evolution */
    
    ParticlesMapIterator it_p;
    for (it_p = particlesMap->begin(); it_p != particlesMap->end(); it_p++)
    {
        Particle *p = (*it_p).second;
        if (p->stellar_type == 13)
        {
            p->magnetic_field_strength_gauss = compute_NS_magnetic_field_Ye19_model(p->magnetic_field_strength_gauss, dt, p->RLOF_timescale, p->delta_mass_RLOF);
        }
    }
}

void nbody_handle_NS_properties_Ye19_model(Particle *p, double ospin, double dt)
{
    if (p->stellar_type == 13)
    {
        double ospin_new = ospin/sqrt( 1.0 - ospin*ospin * NS_Ye19_model_NS_spin_down_constant_K * ONE_DIV_FOURPISQ * p->magnetic_field_strength_gauss*p->magnetic_field_strength_gauss * NS_Ye19_model_NS_B_field_decay_timescale * ( exp(-2.0 * dt/NS_Ye19_model_NS_B_field_decay_timescale) - 1.0) );

        rescale_vector(p->spin_vec,ospin_new/ospin);
        p->magnetic_field_strength_gauss = compute_NS_magnetic_field_Ye19_model(p->magnetic_field_strength_gauss, dt, 0.0, 0.0); /* Do not take into account mass accretion effects in direct N-body mode */

        #ifdef VERBOSE
        if (verbose_flag > 1)
        {
            printf("stellar_evolution.cpp -- nbody_handle_NS_properties_Ye19_model -- dt %g ospin %g ospin_new %g p->magnetic_field_strength_gauss %g P_s %g\n",dt,ospin,ospin_new,p->magnetic_field_strength_gauss,compute_spin_period_from_spin_angular_frequency(ospin_new) * yr_to_s);
        }
        #endif
    }
}

}
