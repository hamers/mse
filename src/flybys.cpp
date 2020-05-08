/* SecularMultiple */
/* Adrian Hamers November 2019 */

#include "evolve.h"
#include "flybys.h"

extern "C"
{

int handle_next_flyby(ParticlesMap *particlesMap, bool initialize, bool *unbound_orbits, int *integration_flag)
{
   
    //printf("flybys.cpp -- sample_next_flyby -- include_flybys %d flybys_correct_for_gravitational_focussing %d flybys_velocity_distribution %d flybys_mass_distribution %d flybys_mass_distribution_lower_value %g flybys_mass_distribution_upper_value %g flybys_encounter_sphere_radius %g flybys_stellar_density %g flybys_stellar_relative_velocity_dispersion %g\n",include_flybys,flybys_correct_for_gravitational_focussing,flybys_velocity_distribution,flybys_mass_distribution,flybys_mass_distribution_lower_value,flybys_mass_distribution_upper_value,flybys_encounter_sphere_radius,flybys_stellar_density,flybys_stellar_relative_velocity_dispersion);

    //printf("handle_next_flyby %d\n",initialize);

    determine_internal_mass_and_semimajor_axis(particlesMap);
    
    if (initialize == true)
    {
        compute_total_encounter_rate_and_density_at_R_enc(&flybys_total_encounter_rate_at_R_enc, &flybys_stellar_density_at_R_enc);
        flybys_W_max = compute_W_max();
        printf("flybys_W_max %g\n",flybys_W_max);
    }

    double b_vec[3], V_vec[3];
    double M_per;
    bool apply_flyby;
    
    sample_next_flyby(particlesMap, &apply_flyby, &flybys_t_next_encounter, &flybys_N_enc, &flybys_N_not_impulsive, &M_per, b_vec, V_vec);
    //printf("flybys.cpp -- sample_next_flyby -- apply_flyby %d N_enc %d N_non_im %d flybys_t_next_encounter %g\n",apply_flyby,flybys_N_enc,flybys_N_not_impulsive,flybys_t_next_encounter);
    
    if (initialize == false and apply_flyby == true)
    {
        compute_effects_of_flyby_on_system(particlesMap, M_per, b_vec, V_vec, unbound_orbits, integration_flag);
        //printf("flybys.cpp -- compute_effects_of_flyby_on_system\n");
    }

    if (*unbound_orbits == true)
    {
        printf("flybys.cpp -- handle_next_flyby -- unbound orbits in system due to flyby!\n");
//        *state = 4; /* TO DO: make general state macros */
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
        M_tot = flybys_internal_mass + *M_per;
        mu = CONST_G*M_tot;
        
        /* Sample perturber velocity at R_enc */
        sample_flyby_position_and_velocity_at_R_enc(particlesMap,R_vec,V_vec);

        /* Compute properties of perturber orbit */
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

        //#Delta_t = -2.0*R_vec.dot(V_per_vec_unit)/V_per

        /* Compute angular speed ratio and reject unsuitable cases */
        Q = norm3(b_vec); /* effective Q */
        V_infty = V; /* effective V_infty */
        E = 1.0 + Q*V_infty*V_infty/mu;
        
        theta_dot_peri = sqrt((1.0 + E)*mu/(Q*Q*Q));
        n_internal = sqrt( CONST_G*flybys_internal_mass/(flybys_internal_semimajor_axis*flybys_internal_semimajor_axis*flybys_internal_semimajor_axis));
        angular_speed_ratio = theta_dot_peri/n_internal;
        
//        if cmd_options["verbose"] == 1:
  //          print 'sample_next_impulsive_encounter; ', 'AS_ratio',angular_speed_ratio
        
       // #t_ref = time_next_encounter - Delta_t_orbit ### time of pericenter passage; note that Delta_t_orbit < 0
        //#t_from_R_enc_to_periapse = - Delta_t_orbit ### time needed for perturber to move from R_enc to Q
        //#t_passed = t_ref - Delta_t_orbit ### time when again hitting the encounter sphere
        //#delta_time_in_encounter_sphere = -2.0*Delta_t_orbit
        
        /* Compute time of next encounter */
        u = ((double) rand() / (RAND_MAX));
        delta_time_encounter = -log(u) / flybys_total_encounter_rate_at_R_enc;
        //delta_time_encounter *= 0.01;
        
        //printf("delta_time_encounter %g\n",delta_time_encounter);
        if (delta_time_encounter <= 0.0)
        {
            //printf("flybys.cpp -- sample_next_flyby -- FATAL ERROR delta_time_encounter <= 0\n");
            exit(-1);
        }
        *t_next_encounter += delta_time_encounter;

        if (angular_speed_ratio < minimum_angular_speed_ratio)
        {
            *N_not_impulsive += 1;
            *apply_flyby = false;
            //printf("flybys.cpp -- sample_next_flyby -- not impulsive; angular_speed_ratio = %g; i=%d\n",angular_speed_ratio,i);
            break;
        }

        *apply_flyby = true;
        //#next_external_time = time + dt
        
        //#print 'Delta_t',dt.value_in(units.Myr)
        
//#        if (t + dt > tend and i < minimum_number_of_encounters): ### ensure that the minimum number of encounters is reached 
//#            continue

        satisfied_criteria = true;

        //if store_encounter_data == True:
            //output_arrays["external_e_arrays"].append(E | units.none)
            //output_arrays["external_q_arrays"].append(Q)
            //output_arrays["external_angular_speed_ratio_arrays"].append(angular_speed_ratio | units.none)            
            //#output_arrays["N_not_hyperbolic"] = N_not_hyperbolic | units.none
            //output_arrays["N_not_impulsive"] = N_not_impulsive | units.none
            
            //print 'N_not_impulsive',N_not_impulsive
            

        //if cmd_options["verbose"] == -1:
        //#if cmd_options["verbose"] == True:
        //#if 1==1:
          //  #print 'time',time.value_in(units.Myr),'dt',dt.value_in(units.Myr)
            //print 'delta_time_encounter',delta_time_encounter.value_in(units.Myr)
            //print 'time_next_encounter',time_next_encounter.value_in(units.Myr)
            //print 'm_per',m_per.value_in(units.MSun),'M_m',M_m.value_in(units.MSun)
            //print 'mu',mu,'R',R.value_in(units.AU),'R_vec',R_vec.value_in(units.AU),'V_vec',V_vec.value_in(units.km/units.s),'V',np.sqrt(V_vec[0]**2+V_vec[1]**2+V_vec[2]**2).value_in(units.km/units.s)
            //print 'q',q.value_in(units.AU),'e',e,'angular_speed_ratio',angular_speed_ratio
            //print 'N_not_hyperbolic',N_not_hyperbolic,'N_not_impulsive',N_not_impulsive
            
           //# print 'external_particle',external_particle
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
        
        //double vx_prime = helper_functions.sample_distribution(None, None, stellar_relative_velocity_distribution_x)
        //double vy_prime = helper_functions.sample_distribution(None, None, stellar_relative_velocity_distribution_y)
        //double vz_prime = helper_functions.sample_distribution(None, None, stellar_relative_velocity_distribution_z)
        //v_prime_norm = np.sqrt(vx_prime**2+vy_prime**2+vz_prime**2)
        //#print 'v',vx_prime.value_in(units.km/units.s),vy_prime.value_in(units.km/units.s),vz_prime.value_in(units.km/units.s),v.value_in(units.km/units.s)

        /* Sample direction of perturber */
        double r_hat_vec[3],theta_hat_vec[3],phi_hat_vec[3];
        sample_spherical_coordinates_unit_vectors_from_isotropic_distribution(r_hat_vec,theta_hat_vec,phi_hat_vec);

        int k;
        for (k=0; k<3; k++)
        {
            R_vec[k] = r_hat_vec[k]*flybys_encounter_sphere_radius;
            V_vec[k] = v_prime[0]*theta_hat_vec[k] + v_prime[1]*phi_hat_vec[k] - fabs(v_prime[2])*r_hat_vec[k];
        }


        /* */
        if (flybys_reference_binary != -1) /* Default value -1: use center of mass (keep R_vec and V_vec fixed) */
        {
            Particle *p = (*particlesMap)[flybys_reference_binary];
            double *r = p->R_vec;
            double *v = p->V_vec;
            #ifdef IGNORE
            double r[3],v[3];
            r[0] = p->X;
            r[1] = p->Y;
            r[2] = p->Z;
            v[0] = p->VX;
            v[1] = p->VY;
            v[2] = p->VZ;
            #endif
            for (k=0; k<3; k++)
            {
                //printf("0 %g %g\n",r[k],v[k]);
                //printf("1 %g %g \n",R_vec[k],V_vec[k]);
                R_vec[k] -= r[k];
                V_vec[k] -= v[k];
                //printf("2 %g %g \n",R_vec[k],V_vec[k]);
            }
            
        }
    }
    else
    {
        printf("flybys.cpp -- flybys_velocity_distribution = %d is not supported -- exiting\n",flybys_velocity_distribution);
        exit(-1);
    }
        
    return 0;
}


int compute_effects_of_flyby_on_system(ParticlesMap *particlesMap, double M_per, double b_per_vec[3], double V_per_vec[3], bool *unbound_orbits, int *integration_flag)
{
    int flag;
    double vx,vy,vz;
    ParticlesMapIterator it_p;
    //std::vector<int>::iterator it_parent_p,it_parent_q;

//    int seed = orbital_phases_random_seed;
    int k;
    
    double V_per = norm3(V_per_vec);
    double V_per_vec_unit[3];
    for (k=0; k<3; k++)
    {
        V_per_vec_unit[k] = V_per_vec[k]/V_per;
    }
    
    double b_i_vec[3],Delta_V_vec[3];
    double *R_vec;
    double b_per_vec_minus_R_vec[3];
    double b_per_vec_minus_R_vec_dot_V_per_vec_unit,b_i_temp;
    
    for (it_p = particlesMap->begin(); it_p != particlesMap->end(); it_p++)
    {
        Particle *p = (*it_p).second;
        if (p->is_binary == false and p->is_bound == true)
        {
            /* Careful: is it OK to set these to zero here? */
            p->instantaneous_perturbation_delta_mass = 0.0;
            p->instantaneous_perturbation_delta_X = 0.0;
            p->instantaneous_perturbation_delta_Y = 0.0;
            p->instantaneous_perturbation_delta_Z = 0.0;

            if (*integration_flag == 0)
            {
                set_positions_and_velocities(particlesMap); /* To make sure the latest positions are used in secular case */
            }
            R_vec = p->R_vec;

            for (k=0; k<3; k++)
            {
                b_per_vec_minus_R_vec[k] = b_per_vec[k] - R_vec[k];
            }
            b_per_vec_minus_R_vec_dot_V_per_vec_unit = dot3(b_per_vec_minus_R_vec,V_per_vec_unit);
            
            for (k=0; k<3; k++)
            {
                b_i_vec[k] = b_per_vec[k] - R_vec[k] - V_per_vec_unit[k]*b_per_vec_minus_R_vec_dot_V_per_vec_unit;
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
            
    reset_instantaneous_perturbation_quantities(particlesMap);

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

    //        #print 'i',index,mean_M,'x',x
                
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
    //        #print 'f',fraction,'V',V,'x',x
        }
        printf("flybys.cpp -- compute_total_encounter_rate_and_density_at_R_enc -- integral_density %g integral_rate %g \n",integral_density,integral_rate);
        *total_encounter_rate = Gamma_0*integral_rate;
        *stellar_density = integral_density*flybys_stellar_density;
    }

    //print 'compute_total_encounter_rate_function -- integral: ',integral,' N_enc,est',total_rate*simulation_parameters_particle.tend

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
    
    double u = ((double) rand() / (RAND_MAX));
    if (u <= (W/flybys_W_max) )
    {
        resample = false;
    }
    else
    {
        resample = true;
    }

    //#print 'W/W_max',W/W_max,random_number,resample
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
    else
    {
        printf("flybys.cpp -- flybys_mass_distribution = %d is not supported; exiting\n",flybys_mass_distribution);
        exit(-1);
    }
    
    return M;
}


}
