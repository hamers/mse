/* MSE */

#include "evolve.h"
#include "flybys.h"

extern "C"
{

bool flyby_criterion(ParticlesMap *particlesMap, int *integration_flag)
{
    bool enable_flybys = true;

    int N_bound_subsystems = 0;
    ParticlesMapIterator it_p;
    for (it_p = particlesMap->begin(); it_p != particlesMap->end(); it_p++)
    {
        Particle *p = (*it_p).second;
        if (p->is_binary == true and p->parent == -1)
        {
            N_bound_subsystems++;
        }
    }

    if (*integration_flag > 0 or N_bound_subsystems != 1)
    {
        #ifdef VERBOSE
        if (verbose_flag > 1)
        {
            printf("flybys.cpp -- flyby_criterion -- integration_flag %d N_bound_subsystems %d -- not handling flyby \n",*integration_flag,N_bound_subsystems);
        }
        #endif
        
        enable_flybys = false;
    }
    
    return enable_flybys;
}    

int handle_next_flyby(ParticlesMap *particlesMap, bool initialize, int *integration_flag)
{
    #ifdef VERBOSE
    if (verbose_flag > 1)
    {
        printf("flybys.cpp -- sample_next_flyby -- include_flybys %d flybys_correct_for_gravitational_focussing %d flybys_velocity_distribution %d flybys_mass_distribution %d flybys_mass_distribution_lower_value %g flybys_mass_distribution_upper_value %g flybys_encounter_sphere_radius %g flybys_stellar_density %g flybys_stellar_relative_velocity_dispersion %g\n",include_flybys,flybys_correct_for_gravitational_focussing,flybys_velocity_distribution,flybys_mass_distribution,flybys_mass_distribution_lower_value,flybys_mass_distribution_upper_value,flybys_encounter_sphere_radius,flybys_stellar_density,flybys_stellar_relative_velocity_dispersion);
    }
    #endif

    determine_internal_mass_and_semimajor_axis(particlesMap);
    
    if (initialize == true)
    {
        compute_total_encounter_rate_and_density_at_R_enc(&flybys_total_encounter_rate_at_R_enc, &flybys_stellar_density_at_R_enc);
        flybys_W_max = compute_W_max();
    }

    double b_vec[3], V_vec[3];
    double M_per;
    bool apply_flyby;
    
    sample_next_flyby(particlesMap, &apply_flyby, &flybys_t_next_encounter, &flybys_N_enc, &flybys_N_not_impulsive, &M_per, b_vec, V_vec);
    
    bool unbound_orbits = false;
    if (initialize == false and apply_flyby == true)
    {
        compute_effects_of_flyby_on_system(particlesMap, M_per, b_vec, V_vec, &unbound_orbits, true, integration_flag);
    }

    if (unbound_orbits == true)
    {
        #ifdef VERBOSE
        if (verbose_flag > 0)
        {
            printf("flybys.cpp -- handle_next_flyby -- unbound orbits in system due to flyby!\n");
        }
        #endif

        *integration_flag = 4;
    }

    return 0;
}

int sample_next_flyby(ParticlesMap *particlesMap, bool *apply_flyby, double *t_next_encounter, int *N_enc, int *N_not_impulsive, double *M_per, double b_vec[3], double V_vec[3])
{

    double minimum_angular_speed_ratio = 1.0;

    bool satisfied_criteria = false;
    double R_vec[3];
    double M_tot,mu;
    int k;
    double V_vec_unit[3];
    double V,R_vec_dot_V_vec_unit;
    double Q,V_infty,E;
    double theta_dot_peri,n_internal,angular_speed_ratio;
    
    double delta_time_encounter,u;
    
    int i = 0;
    while (satisfied_criteria == false)
    {
        i += 1;
        
        /* Sample perturber mass at R_enc */
        *M_per = sample_flyby_mass_at_R_enc();
        
        #ifdef VERBOSE
        if (verbose_flag > 1)
        {
            printf("flybys.cpp -- sample_next_flyby.cpp -- M_per %g S %d\n",*M_per,random_seed);
        }
        #endif
        
        M_tot = flybys_internal_mass + *M_per;
        mu = CONST_G*M_tot;
        
        /* Sample perturber velocity at R_enc */
        sample_flyby_position_and_velocity_at_R_enc(particlesMap,R_vec,V_vec);
        
        #ifdef VERBOSE
        if (verbose_flag > 1)
        {
            printf("flybys.cpp -- sample_next_flyby.cpp -- R %g %g %g V %g %g %g\n",R_vec[0],R_vec[1],R_vec[2],V_vec[0],V_vec[1],V_vec[2]);
        }
        #endif
        
        /* Compute impact parameter for perturber orbit */
        V = norm3(V_vec);

        for (k=0; k<3; k++)
        {
            V_vec_unit[k] = V_vec[k]/V;
        }
        
        double R_vec_dot_V_vec_unit = dot3(R_vec,V_vec_unit);
        for (k=0; k<3; k++)
        {
            
            b_vec[k] = R_vec[k] - V_vec_unit[k]*R_vec_dot_V_vec_unit;
        }

        /* Compute angular speed ratio and reject unsuitable cases */
        Q = norm3(b_vec); /* effective Q */
        V_infty = V; /* effective V_infty */
        E = 1.0 + Q*V_infty*V_infty/mu;
        
        theta_dot_peri = sqrt((1.0 + E)*mu/(Q*Q*Q));
        n_internal = sqrt( CONST_G*flybys_internal_mass/(flybys_internal_semimajor_axis*flybys_internal_semimajor_axis*flybys_internal_semimajor_axis));
        angular_speed_ratio = theta_dot_peri/n_internal;
        
        /* Compute time of next encounter */
        u = generate_random_number_between_zero_and_unity();
        delta_time_encounter = -log(u) / flybys_total_encounter_rate_at_R_enc;
        
        //printf("delta_time_encounter %g\n",delta_time_encounter);
        if (delta_time_encounter <= 0.0)
        {
            printf("flybys.cpp -- sample_next_flyby -- FATAL ERROR delta_time_encounter <= 0\n");
            //exit(-1);
            error_code = 14;
        }
        *t_next_encounter += delta_time_encounter;

        if (angular_speed_ratio < minimum_angular_speed_ratio)
        {
            *N_not_impulsive += 1;
            *apply_flyby = false;
            
            #ifdef VERBOSE
            if (verbose_flag > 1)
            {
                printf("flybys.cpp -- sample_next_flyby -- not impulsive; angular_speed_ratio = %g; i=%d\n",angular_speed_ratio,i);
            }
            #endif            
            
            break;
        }

        #ifdef VERBOSE
        if (verbose_flag > 1)
        {
            printf("flybys.cpp -- sample_next_flyby -- V %g M_per %g Q %g t_next_encounter %g\n",V,*M_per,Q,*t_next_encounter);
            print_system(particlesMap,0);
        }
        #endif
        

        *apply_flyby = true;
       
        satisfied_criteria = true;
    }
    
    *N_enc += 1;
    
    return 0;
}

int sample_flyby_position_and_velocity_at_R_enc(ParticlesMap *particlesMap, double R_vec[3], double V_vec[3])
{
    if (flybys_velocity_distribution == 0) /* Maxwellian */
    {

        /* Sample velocity magnitude at R_enc */
        double v_prime[3];
        sample_from_3d_maxwellian_distribution(flybys_stellar_relative_velocity_dispersion, v_prime);
        v_prime[2] = sample_from_y_times_maxwellian_distribution(flybys_stellar_relative_velocity_dispersion); /* Need to sample the z-component from a distribution \propto y*Exp[-y^2/(2 sigma^2)] */
        
        /* Sample direction of perturber */
        double r_hat_vec[3],theta_hat_vec[3],phi_hat_vec[3];
        sample_spherical_coordinates_unit_vectors_from_isotropic_distribution(r_hat_vec,theta_hat_vec,phi_hat_vec);

        int k;
        for (k=0; k<3; k++)
        {
            R_vec[k] = r_hat_vec[k]*flybys_encounter_sphere_radius;
            V_vec[k] = v_prime[0]*theta_hat_vec[k] + v_prime[1]*phi_hat_vec[k] - fabs(v_prime[2])*r_hat_vec[k];
        }

    }
    else
    {
        printf("flybys.cpp -- flybys_velocity_distribution = %d is not supported -- exiting\n",flybys_velocity_distribution);
        //exit(-1);
        error_code = 15;
    }
        
    return 0;
}


int compute_effects_of_flyby_on_system(ParticlesMap *particlesMap, double M_per, double b_per_vec[3], double V_per_vec[3], bool *unbound_orbits, bool reset, int *integration_flag)
{
    ParticlesMapIterator it_p;
    int k;
    
    double V_per = norm3(V_per_vec);
    double V_per_vec_unit[3];
    for (k=0; k<3; k++)
    {
        V_per_vec_unit[k] = V_per_vec[k]/V_per;
    }
    
    double b_i_vec[3],Delta_V_vec[3];
    double b_per_vec_minus_R_vec[3];
    double b_per_vec_minus_R_vec_dot_V_per_vec_unit,b_i_temp;
    
    double R_ref[3] = {0.0,0.0,0.0};
    if (flybys_reference_binary != -1) /* Default value -1: use center of mass (center R_vec and V_vec around origin) */
    {
        
        Particle *p = (*particlesMap)[flybys_reference_binary];
        for (k=0; k<3; k++)
        {
            R_ref[k] = p->R_vec[k];
            
        }
    }
    
    #ifdef VERBOSE
    if (verbose_flag > 1)
    {
        printf("flybys.cpp -- compute_effects_of_flyby_on_system -- R_ref %g %g %g\n",R_ref[0],R_ref[1],R_ref[2]);
        print_system(particlesMap,0);
    }
    #endif
    
    for (it_p = particlesMap->begin(); it_p != particlesMap->end(); it_p++)
    {
        Particle *p = (*it_p).second;
        if (p->is_binary == false and p->is_bound == true)
        {
            p->instantaneous_perturbation_delta_mass = 0.0;
            p->instantaneous_perturbation_delta_X = 0.0;
            p->instantaneous_perturbation_delta_Y = 0.0;
            p->instantaneous_perturbation_delta_Z = 0.0;

            if (*integration_flag == 0)
            {
                set_positions_and_velocities(particlesMap); /* To make sure the latest positions are used in secular case */
            }

            for (k=0; k<3; k++)
            {
                b_per_vec_minus_R_vec[k] = b_per_vec[k] - p->R_vec[k] - R_ref[k];
            }
            b_per_vec_minus_R_vec_dot_V_per_vec_unit = dot3(b_per_vec_minus_R_vec,V_per_vec_unit);
            
            for (k=0; k<3; k++)
            {
                b_i_vec[k] = b_per_vec_minus_R_vec[k] - V_per_vec_unit[k]*b_per_vec_minus_R_vec_dot_V_per_vec_unit;
            }
            
            b_i_temp = 2.0*(CONST_G*M_per/V_per)/norm3_squared(b_i_vec);
            for (k=0; k<3; k++)
            {
                Delta_V_vec[k] = b_i_temp*b_i_vec[k];
            }
            
            p->instantaneous_perturbation_delta_VX = Delta_V_vec[0];
            p->instantaneous_perturbation_delta_VY = Delta_V_vec[1];
            p->instantaneous_perturbation_delta_VZ = Delta_V_vec[2];
            
        }
    }
    
    if (*integration_flag == 0)
    {
        apply_user_specified_instantaneous_perturbation(particlesMap);
    }
    else
    {
        apply_user_specified_instantaneous_perturbation_nbody(particlesMap);
    }

    *unbound_orbits = check_for_unbound_orbits(particlesMap);
            
    if (reset == true)
    {
        reset_instantaneous_perturbation_quantities(particlesMap);
    }

    return 0;
}

int compute_total_encounter_rate_and_density_at_R_enc(double *total_encounter_rate, double *stellar_density)
{
    double R_enc = flybys_encounter_sphere_radius;
    double n_star = flybys_stellar_density;
    double sigma_rel = flybys_stellar_relative_velocity_dispersion;
    double Gamma_0 = 2.0*sqrt(2.0*M_PI)*R_enc*R_enc*n_star*sigma_rel;

    if (flybys_correct_for_gravitational_focussing == false) /* do not correct for gravitational focussing */
    {
        *total_encounter_rate = Gamma_0;
        *stellar_density = flybys_stellar_density;
    }
    else
    {
        int N_points = 1000;

        double points[N_points];
        double Ms[N_points];
        int k;
        double log10_m_lower = log10(flybys_mass_distribution_lower_value);
        double log10_m_upper = log10(flybys_mass_distribution_upper_value);
        double r;
        for (k=0; k<N_points; k++)
        {
            r = ((double) k)/( (double) N_points ); 
            points[k] = pow(10.0, log10_m_lower + (log10_m_upper - log10_m_lower)*r );
            Ms[k] = sample_flyby_mass_at_infinity();
        }
        
        double M_lower,M_upper,delta_M,mean_M;
        double x,W,V,mass_fraction;
        double integral_rate = 0.0;
        double integral_density = 0.0;
        int l;
        
        for (k=0; k<N_points; k++)
        {
            if (k==0)
            {
                continue;
            }
            M_lower = points[k-1];
            M_upper = points[k];
            delta_M = M_upper - M_lower;
            mean_M = (1.0/2.0)*(M_lower + M_upper);

            x = x_function(flybys_internal_mass + mean_M);
            V = V_function(x);
            W = W_function(x);

            mass_fraction = 0.0;
            for (l=0; l<N_points; l++)
            {
                if (Ms[l] > M_lower and Ms[l] <= M_upper)
                {
                    mass_fraction += 1.0/N_points;
                }
            }

            integral_density += mass_fraction*W;
            integral_rate += mass_fraction*V;
        }
        
        #ifdef VERBOSE
        if (verbose_flag > 1)
        {
            printf("flybys.cpp -- compute_total_encounter_rate_and_density_at_R_enc -- integral_density %g integral_rate %g \n",integral_density,integral_rate);
        }
        #endif
        
        *total_encounter_rate = Gamma_0*integral_rate;
        *stellar_density = integral_density*flybys_stellar_density;
    }

    return 0;
}

double x_function(double M_tot)
{
    return CONST_G*(M_tot)/(flybys_encounter_sphere_radius*flybys_stellar_relative_velocity_dispersion*flybys_stellar_relative_velocity_dispersion);
}
double W_function(double x)
{
    return 2.0*sqrt(x/M_PI) + exp(x)*erfc(sqrt(x));
}

double V_function(double x)
{
    return 1.0 + x;
}
    
double compute_W_max()
{
    double x_max = x_function(flybys_internal_mass + flybys_mass_distribution_upper_value);

    return W_function(x_max);
}

bool correct_mass_function(double M)
{
    double x = x_function(flybys_internal_mass + M);
    double W = W_function(x);
    bool resample;
    
    double u = generate_random_number_between_zero_and_unity();
    if (u <= (W/flybys_W_max) )
    {
        resample = false;
    }
    else
    {
        resample = true;
    }

    #ifdef VERBOSE
    if (verbose_flag > 1)
    {
        printf("flybys.cpp -- correct_mass_function -- W %g resample %d\n",W,resample);
    }
    #endif

    return resample;
}

double sample_flyby_mass_at_R_enc()
{
    bool sampled_mass = false;
    bool resample;
    double M;
    
    while (sampled_mass == false)
    {
        M = sample_flyby_mass_at_infinity();
        
        if (flybys_correct_for_gravitational_focussing == true)
        {
            resample = correct_mass_function(M);

            if (resample == false)
            {
                sampled_mass = true;
            }
        }
        else
        {
            sampled_mass = true;
        }
    }

    return M;
}

double sample_flyby_mass_at_infinity()
{
    double M;
    if (flybys_mass_distribution == 0) /* Salpeter */
    {
        M = sample_from_power_law_distribution(-2.35, flybys_mass_distribution_lower_value, flybys_mass_distribution_upper_value);
    }
    else if (flybys_mass_distribution == 1)
    {
        M = sample_from_Kroupa_93_imf();
    }
    else
    {
        printf("flybys.cpp -- flybys_mass_distribution = %d is not supported; exiting\n",flybys_mass_distribution);
        //exit(-1);
        error_code = 16;
    }
    
    return M;
}

}
