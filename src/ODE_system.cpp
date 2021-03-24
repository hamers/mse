/* MSE */

#include "types.h"
#include "evolve.h"
#include "ODE_system.h"

extern "C"
{
    
int integrate_ODE_system(ParticlesMap *particlesMap, double start_time, double end_time, double *output_time, double *hamiltonian, int *output_flag, int *error_code)
{
    int N_particles = particlesMap->size();
    int N_bodies, N_binaries;
    int N_root_finding;
    int N_ODE_equations;
    
    determine_binary_parents_and_levels(particlesMap,&N_bodies,&N_binaries,&N_root_finding,&N_ODE_equations);
    set_binary_masses_from_body_masses(particlesMap);

    initialize_direct_integration_quantities(particlesMap);

    check_for_integration_exclusion_orbits(particlesMap);

    #ifdef VERBOSE
    if (verbose_flag > 2)
    {
        printf("ODE_system.cpp -- evolve -- N_bodies %d N_binaries %d N_particles %d N_root_finding %d N_ODE_equations %d\n",N_bodies,N_binaries,N_particles,N_root_finding,N_ODE_equations);
        print_system(particlesMap,0);

    }
    #endif
   
    /*********************
     * setup of UserData *
     ********************/
     
	UserData data;
	data = NULL;
	data = (UserData) malloc(sizeof *data);
	data->particlesMap = particlesMap;
    data->N_root_finding = N_root_finding;
    data->start_time = start_time;

    /********************************
     * set ODE tolerances   *
     ********************************/

    #ifdef VERBOSE
    if (verbose_flag > 2)
    {
        for (int i=1; i<=N_ODE_equations; i++)
        {
            printf("ODE_evolve.cpp -- evolve --  relative_tolerance %g absolute_tolerance_spin_vectors %g absolute_tolerance_eccentricity_vectors %g absolute_tolerance_angular_momentum_vectors %g\n",relative_tolerance,absolute_tolerance_spin_vectors,absolute_tolerance_eccentricity_vectors,absolute_tolerance_angular_momentum_vectors);
        }
    }
    #endif
        
    double abs_tol_spin_vec = absolute_tolerance_spin_vectors;
    double abs_tol_e_vec = absolute_tolerance_eccentricity_vectors;
    double abs_tol_h_vec = absolute_tolerance_angular_momentum_vectors;
    double initial_ODE_timestep = ODE_min_dt; /* one year by default */
    int maximum_number_of_internal_ODE_steps = 5e8;
    int maximum_number_of_convergence_failures = 100;    
    //double maximum_ODE_integration_time = 13.8e10;


    /***************************
     * setup of ODE variables  *
     **************************/    
	N_Vector y, y_out, y_abs_tol;
	void *cvode_mem;
	int flag;

	y = y_out = y_abs_tol = NULL;
	cvode_mem = NULL;

    int number_of_ODE_variables = N_ODE_equations;
    y = N_VNew_Serial(number_of_ODE_variables);
	if (check_flag((void *)y, "N_VNew_Serial", 0)) return 1;
    y_out = N_VNew_Serial(number_of_ODE_variables);
	if (check_flag((void *)y_out, "N_VNew_Serial", 0)) return 1;
    y_abs_tol = N_VNew_Serial(number_of_ODE_variables); 
	if (check_flag((void *)y_abs_tol, "N_VNew_Serial", 0)) return 1;         

    
    set_initial_ODE_variables(particlesMap, y, y_abs_tol,abs_tol_spin_vec,abs_tol_e_vec,abs_tol_h_vec);

    #ifdef VERBOSE
    if (verbose_flag > 2)
    {
        for (int i=1; i<=N_ODE_equations; i++)
        {
            printf("ODE_evolve.cpp -- evolve -- i %d Ith(y,i) %g Ith(y_abs,i) %g\n",i,Ith(y,i),Ith(y_abs_tol,i));
        }
    }
    #endif
    
    /***************************
     * setup of ODE integrator *
     **************************/    

    /* use Backward Differentiation Formulas (BDF)
        scheme in conjunction with Newton iteration --
        these choices are recommended for stiff ODEs
        in the CVODE manual                          
    */

    cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);
	if (check_flag((void *)cvode_mem, "CVodeCreate", 0)) return 1;    
    
    /* essential initializations */
    flag = CVodeInit(cvode_mem, compute_y_dot, start_time, y);
	if (check_flag(&flag, "CVodeInit", 1)) return 1;    

	flag = CVodeSetUserData(cvode_mem, data);
	if (check_flag(&flag, "CVodeSetUsetData", 1)) return 1;

    flag = CVodeSVtolerances(cvode_mem, relative_tolerance, y_abs_tol);
	if (check_flag(&flag, "CVodeSVtolerances", 1)) return 1;

	flag = CVDense(cvode_mem, number_of_ODE_variables);
	if (check_flag(&flag, "CVDense", 1)) return 1;

	flag = CVodeSetInitStep(cvode_mem, initial_ODE_timestep);
	if (check_flag(&flag, "CVodeSetInitStep", 1)) return 1;

    /* optional initializations */
//	flag = CVodeSetErrHandlerFn(cvode_mem, ehfun, eh_data); // error handling function
//	if (check_flag(&flag, "CVodeSetErrHandlerFn", 1)) return;
	  		
	flag = CVodeSetMaxNumSteps(cvode_mem, maximum_number_of_internal_ODE_steps);
	if (check_flag(&flag, "CVodeSetMaxNumSteps", 1)) return 1;

//	flag = CVodeSetMinStep(cvode_mem, 1.0e-14); // minimum step size
//	if (check_flag(&flag, "CVodeSetMinStep", 1)) return 1;

	flag = CVodeSetMaxHnilWarns(cvode_mem, 1);
	if (check_flag(&flag, "CVodeSetMaxHnilWarns", 1)) return 1;
			
//	flag = CVodeSetStopTime(cvode_mem, MAXTIME); // maximum time
//	if (check_flag(&flag, "CVodeSetStopTime", 1)) return 1;

	flag = CVodeSetMaxConvFails(cvode_mem, maximum_number_of_convergence_failures);
	if (check_flag(&flag, "CVodeSetMaxConvFails", 1)) return 1;

    /* initialization of root finding */
    int roots_found[N_root_finding];
	flag = CVodeRootInit(cvode_mem, N_root_finding, root_finding_functions);
	if (check_flag(&flag, "CVodeRootInit", 1)) return 1;	



    /***************************
     * ODE integration         *
     **************************/ 
    
    double user_end_time = end_time;
	realtype integrator_end_time;

    if (check_for_initial_roots(particlesMap) > 0) /* do not evolve ODE system when a root was found initially, EXCEPT for RLOF */
    {
        flag = CV_ROOT_RETURN;
        *output_flag = CV_ROOT_RETURN;
        *output_time = start_time;
        
        return 0;
    }
    else
    {
        flag = CVode(cvode_mem, user_end_time, y_out, &integrator_end_time, CV_NORMAL);	
    }
    
	if (flag == CV_SUCCESS)
	{
		*output_flag = CV_SUCCESS;
		*error_code = 0;
        *output_time = integrator_end_time;
	}
	else if (flag == CV_ROOT_RETURN) // a root was found during the integration
	{
		CVodeGetRootInfo(cvode_mem,roots_found);
        read_root_finding_data(particlesMap,roots_found);
        *output_flag = CV_ROOT_RETURN;
        *output_time = integrator_end_time;
    }
    else if (flag == CV_WARNING) // a warning has occurred during the integration
    {
		*output_flag = 99;
		*error_code = flag;
    }
	else // an error has occurred during the integration
    {
		*output_flag = flag;
		*error_code = flag;
    }

    /***************************
     * y_out -> particlesMap   *
     * ************************/

    double actual_time_step = integrator_end_time - start_time;

    extract_final_ODE_variables(particlesMap,y_out);
    set_binary_masses_from_body_masses(particlesMap); /* update masses in binaries, needed for direct integration update below */
    process_direct_integration_quantities(particlesMap,actual_time_step); /* update relative pos & vel, as well as new e, h vectors & TA */

    set_positions_and_velocities(particlesMap); /* update positions and velocities of all bound bodies based on updated binary e, h, and TA */
    update_positions_unbound_bodies(particlesMap, actual_time_step); /* update positions of unbound bodies, which are assumed (integration_flag = 0) to move unaccelerated */

    *hamiltonian = data->hamiltonian;
    
    N_VDestroy_Serial(y);
    N_VDestroy_Serial(y_out);
    N_VDestroy_Serial(y_abs_tol);
    CVodeFree(&cvode_mem);

	return 0;
}

    
    
int compute_y_dot(realtype t, N_Vector y, N_Vector y_dot, void *data_)
{
    
    time_t wall_time_current;
    time(&wall_time_current);
    double wall_time_diff_s = difftime(wall_time_current,wall_time_start);
    if (wall_time_diff_s > wall_time_max_s)
    {
        //printf("ODE_system.cpp -- wall time of %g s has exceeded allowed maximum of %g s -- initiating termination procedure \n",wall_time_diff_s,wall_time_max_s);
        error_code = 35;
        return 0;
    }
    
	UserData data;
	data = (UserData) data_;
    ParticlesMap *particlesMap = data->particlesMap;
    
    double start_time = data->start_time;
    double delta_time = t - start_time;

    #ifdef VERBOSE
    if (verbose_flag > 2)
    {
        printf("ODE_system.cpp -- compute_y_dot1 t=%g start_time=%g delta_time = %g\n",t,start_time,delta_time);
        print_system(particlesMap,0);

    }
    #endif
    
    extract_ODE_variables(particlesMap, y, delta_time);
    reset_ODE_dots(particlesMap, y, delta_time);

    //printf("ODE_system.cpp -- compute_y_dot2 t=%g start_time=%g delta_time = %g\n",time,start_time,delta_time);
    //print_system(particlesMap,0);
    
    /****************************
     * compute right-hand sides *
     * **************************/

    double hamiltonian = 0.0;
    double KS_V = 0.0;
    bool compute_hamiltonian_only = false;
    ParticlesMapIterator it_p,it_f;
    std::vector<int>::iterator it_parent_p,it_parent_q;
    
    //printf("compute_y_dot t=%g\n",time);
    //print_system(particlesMap);
    
    for (it_p = particlesMap->begin(); it_p != particlesMap->end(); it_p++)
    {
        Particle *p = (*it_p).second;
        if (p->is_binary == true)
        {
            /* Newtonian gravitational point mass dynamics (internal system) */
            compute_EOM_Newtonian_for_particle(particlesMap,p,&hamiltonian,&KS_V,compute_hamiltonian_only);
            
            /* Perturbations by flybys (external perturbations) */
            compute_EOM_Newtonian_external_perturbations(t,particlesMap,p,&hamiltonian,&KS_V,compute_hamiltonian_only); 

            /* VRR-related perturbations */

            if (p->VRR_model > 0) /* 0: no VRR perturbation */
            {
                hamiltonian += compute_VRR_perturbations(particlesMap,p->index,t);
            }

            /* PN corrections */
            compute_EOM_Post_Newtonian_for_particle(particlesMap, p, &hamiltonian, compute_hamiltonian_only);

            /* Tidal friction */
            handle_secular_tidal_evolution(particlesMap,p);

            /* Orbital changes due to mass changes */
            ODE_handle_stellar_winds(p);
            ODE_handle_RLOF(particlesMap,p);
            
        }
            
    }
    compute_KS_EOM(particlesMap,KS_V);
    
    write_ODE_variables_dots(particlesMap,y_dot);

    data->hamiltonian = hamiltonian;

    return 0;
}

void initialize_direct_integration_quantities(ParticlesMap *particlesMap)
{
    
    bool KS_setup_required = false;
    
    ParticlesMapIterator it;
    int highest_level = particlesMap->begin()->second->highest_level;
    int level = 0;
    while (level < highest_level)
    {
        for (it = particlesMap->begin(); it != particlesMap->end(); it++)
        {
            Particle *p = (*it).second;
            if ((p->is_binary == true) && (p->level == level))
            {
                if (p->integration_method==0)
                {
                    /* compute initial mean anomaly */
                    p->initial_mean_anomaly = compute_mean_anomaly_from_true_anomaly(p->true_anomaly,norm3(p->e_vec));
                }
                else
                {
                    from_orbital_vectors_to_cartesian(
                        p->child1_mass,p->child2_mass,
                        p->e_vec,p->h_vec,
                        p->true_anomaly,
                        p->r_vec,p->v_vec);
                
                    p->r = norm3(p->r_vec);    
                }

                if (p->integration_method==1)
                {
                    KS_setup_required = true;
                }

                
            }
        }
        level++;
    }
    set_up_derived_quantities(particlesMap);

    if (KS_setup_required==true)
    {

        double alpha_vec[4],beta_vec[4];
    
        double hamiltonian = 0.0;
        double KS_V = 0.0;
        
        for (it = particlesMap->begin(); it != particlesMap->end(); it++)
        {
            Particle *p = (*it).second;
            if (p->is_binary == true)// && (p->integration_method==1))
            {
                compute_EOM_Newtonian_for_particle(particlesMap,p,&hamiltonian,&KS_V,true);
            }
        }

        for (it = particlesMap->begin(); it != particlesMap->end(); it++)
        {
            Particle *p = (*it).second;
            if ((p->is_binary == true) && (p->integration_method==1))
            {
                double V = KS_V/p->mu;

                p->KS_omega = sqrt( c_1div2*( CONST_G*p->mass/p->r - c_1div2*norm3_squared(p->v_vec) - V) );
                p->KS_E = 0.0;

                transform_r_to_u_d0(p->r_vec,p->KS_u_vec);
                transform_v_to_u_star(p->KS_u_vec,p->v_vec,p->KS_omega,p->KS_u_star_vec);

                for (int i=0; i<4; i++)
                {
                    p->KS_alpha_vec[i] = p->KS_u_vec[i];
                    p->KS_beta_vec[i] = 2.0*p->KS_u_star_vec[i];
                }
            }
        }
    }
}

void process_direct_integration_quantities(ParticlesMap *particlesMap, double delta_time)
{
    double r_child1[3],r_child2[3];
    double v_child1[3],v_child2[3];
    
    double u_vec[4],u_star_vec[4];
    
    ParticlesMapIterator it;

    for (it = particlesMap->begin(); it != particlesMap->end(); it++)
    {
        Particle *p = (*it).second;
        if (p->is_binary == true)
        {
            if (p->integration_method==0) /* advance true anomalies of averaged binaries (not done by ODE solver) */
            {
                double n = sqrt(CONST_G*p->mass/(p->a*p->a*p->a));
                double new_MA = p->initial_mean_anomaly + delta_time*n;
                new_MA = remainder(new_MA,TWOPI);

                p->true_anomaly = compute_true_anomaly_from_mean_anomaly(new_MA,norm3(p->e_vec));
                from_orbital_vectors_to_cartesian(p->child1_mass,p->child2_mass,p->e_vec,p->h_vec,p->true_anomaly,p->r_vec,p->v_vec);

                #ifdef VERBOSE
                if (verbose_flag > 2)
                {
                    printf("ODE_system.cpp -- process_direct_integration_quantities -- n %g new_MA %g old_MA %g\n",n,new_MA,p->initial_mean_anomaly);
                }
                #endif
                
            }
            else if (p->integration_method>0)
            {

                from_cartesian_to_orbital_vectors(
                    p->child1_mass,p->child2_mass,
                    p->r_vec,p->v_vec,
                    p->e_vec,p->h_vec,&p->true_anomaly);

                if (p->integration_method==1)
                {
                    transform_alpha_beta_to_u_u_star(p->KS_alpha_vec,p->KS_beta_vec,p->KS_E,u_vec,u_star_vec);
                    transform_u_d0_to_r(u_vec,p->r_vec);
                    transform_u_star_to_v(u_vec,u_star_vec,p->KS_omega,p->v_vec);
                }
            }
        }
    }
    
}

void compute_KS_EOM(ParticlesMap *particlesMap, double KS_V)
{
    ParticlesMapIterator it_p;
    int k=1;
    int i;

    double domega_dE,domega_dt;
    double *u_vec;
    double *u_star_vec;
    double LT_P[4];
    double omega;
    double omega_inv;
    double r,r_inv;
    double E,dE_dt;
    double sin_E_div_2;
    double cos_E_div_2;
    double temp;
    
    for (it_p = particlesMap->begin(); it_p != particlesMap->end(); it_p++)
    {
        Particle *p = (*it_p).second;
        if (p->is_binary == true && p->integration_method==1)
        {
            E = p->KS_E;
            sin_E_div_2 = sin(c_1div2*E);
            cos_E_div_2 = cos(c_1div2*E);
            
            omega = p->KS_omega;
            omega_inv = 1.0/omega;
        
            r = p->r;
            r_inv = 1.0/r;
            dE_dt = 2.0*omega*r_inv;
        
            u_star_vec = p->KS_u_star_vec;
            LT_u_on_vec3(p->KS_u_vec,p->a_vec,LT_P);
            
            if (p->KS_use_perturbing_potential==false)
            {
                domega_dE = -c_1div2 * omega_inv * dot4(u_star_vec,LT_P);
                domega_dt = domega_dE*dE_dt;
                
                for (i=0; i<4; i++)
                {
                    temp = ( -c_1div2 * omega_inv *  LT_P[i] + 4.0 * r_inv * domega_dE * u_star_vec[i] );
                    p->KS_dalpha_vec_dt[i] = temp * sin_E_div_2;
                    p->KS_dbeta_vec_dt[i] = -temp * cos_E_div_2;
                }

                p->KS_domega_dt = domega_dt;
                p->KS_dE_dt = dE_dt;
            }
            else
            {
                u_vec = p->KS_u_vec;                
                for (i=0; i<4; i++)
                {
                    temp = c_1div2 * r_inv * omega_inv * ( (1.0/p->mu) * KS_V * u_vec[i] -  r * LT_P[i] );
                    p->KS_dalpha_vec_dt[i] = temp * sin_E_div_2;
                    p->KS_dbeta_vec_dt[i] = -temp * cos_E_div_2;
                }

                p->KS_domega_dt = 0.0;
                p->KS_dE_dt = dE_dt;
            }
        }
    }
}

void extract_ODE_variables(ParticlesMap *particlesMap, N_Vector &y, double delta_time)
{
    ParticlesMapIterator it_p;
    int k=1;
    int k_component;

    for (it_p = particlesMap->begin(); it_p != particlesMap->end(); it_p++)
    {
        Particle *p = (*it_p).second;
        if (p->is_binary == false)
        {
            for (k_component=0; k_component<3; k_component++)
            {
                p->spin_vec[k_component] = Ith(y,k + k_component);
            }
            p->mass = Ith(y,k + 2 + 1);
            p->radius = Ith(y,k + 2 + 2);
            
            k=k+5;
        }
        if (p->is_binary == true)
        {
            if (p->integration_method==0)
            {
                for (k_component=0; k_component<3; k_component++)
                {
                    p->e_vec[k_component] = Ith(y,k + k_component);
                    p->h_vec[k_component] = Ith(y,k + k_component + 3);
                }
                if (norm3(p->e_vec) > 1.0)
                {
                   
                    #ifdef VERBOSE
                    if (verbose_flag > 0)
                    {
                        printf("ODE_system.cpp -- extract_ODE_variables -- p %d e = %g > 1; original e_vec %g %g %g\n",p->index,norm3(p->e_vec),p->e_vec[0],p->e_vec[1],p->e_vec[2]); 
                    }
                    #endif

                    double e_norm = norm3(p->e_vec);
                    double e_overshoot = e_norm - 1.0;
                    for (k_component=0; k_component<3; k_component++)
                    {
                        p->e_vec[k_component] = fmod(1.0 - e_overshoot,1.0) * p->e_vec[k_component]/e_norm;
                    }

                    #ifdef VERBOSE
                    if (verbose_flag > 0)
                    {
                        printf("ODE_system.cpp -- extract_ODE_variables -- p %d; new adjusted e = %g; e_vec %g %g %g\n",p->index,norm3(p->e_vec),p->e_vec[0],p->e_vec[1],p->e_vec[2]); 
                    }
                    #endif


                }
               
                k=k+6;
            }
            else if (p->integration_method==1)
            {
                for (k_component=0; k_component<4; k_component++)
                {
                    p->KS_alpha_vec[k_component] = Ith(y,k + k_component);
                    p->KS_beta_vec[k_component] = Ith(y,k + k_component + 4);
                }
                
                p->KS_omega = Ith(y,k + 7 + 1);
                p->KS_E = Ith(y,k + 7 + 2);
                
                transform_alpha_beta_to_u_u_star(p->KS_alpha_vec,p->KS_beta_vec,p->KS_E,p->KS_u_vec,p->KS_u_star_vec);
                transform_u_d0_to_r(p->KS_u_vec,p->r_vec);
                transform_u_star_to_v(p->KS_u_vec,p->KS_u_star_vec,p->KS_omega,p->v_vec);
               
                p->r = norm3(p->r_vec);
                
                k=k+10;
            }
            else if (p->integration_method==2)
            {
                for (k_component=0; k_component<3; k_component++)
                {
                    p->r_vec[k_component] = Ith(y,k + k_component);
                    p->v_vec[k_component] = Ith(y,k + k_component + 3);
                }
                
                k=k+6;
            }
        }
    }
    
    /* update children masses etc. based on new body masses */
    set_binary_masses_from_body_masses(particlesMap);
    set_up_derived_quantities(particlesMap);
}
    

void reset_ODE_dots(ParticlesMap *particlesMap, N_Vector &y, double delta_time)
{
    ParticlesMapIterator it_p;
    int i;
    
    for (it_p = particlesMap->begin(); it_p != particlesMap->end(); it_p++)
    {
        Particle *p = (*it_p).second;
        if (p->is_binary == false)
        {

            /* Externally imposed */
            double spin_vec_norm = norm3(p->spin_vec);
            if (spin_vec_norm == 0.0)
            {
                spin_vec_norm = epsilon;
            }
            //printf("spin_vec_norm %g\n",spin_vec_norm);
            for (i=0; i<3; i++)
            {
                p->dspin_vec_dt[i] = p->ospin_dot * (p->spin_vec[i]/spin_vec_norm);
                //p->dspin_vec_dt[i] = 0.0;
            }
            p->dmass_dt = p->mass_dot_wind + p->mass_dot_wind_accretion;
            p->dradius_dt = p->radius_dot + p->radius_ddot*delta_time;
            
            #ifdef VERBOSE
            if (verbose_flag > 2)
            {
                printf("ODE_system.cpp -- reset_ODE_dots -- body i %d dmass_dt %g dradius_dt %g dspin_vec_dt %g %g %g p->ospin_dot %g p->spin_vec %g %g %g\n",p->index,p->dmass_dt,p->dradius_dt,p->dspin_vec_dt[0],p->dspin_vec_dt[1],p->dspin_vec_dt[2],p->ospin_dot,p->spin_vec[0],p->spin_vec[1],p->spin_vec[2]);
            }
            #endif

        }
        else
        {
            if (p->integration_method==0)
            {
                //double factor_h_vec = p->child1_mass_dot/p->child1_mass + p->child2_mass_dot/p->child2_mass - (p->child1_mass_dot + p->child2_mass_dot)/p->child1_mass_plus_child2_mass;
                double factor_h_vec = 0.0; // now handled in mass_changes.cpp

                for (i=0; i<3; i++)
                {        
                    p->de_vec_dt[i] = 0.0;
                    p->dh_vec_dt[i] = p->h_vec[i]*factor_h_vec;
                }

               #ifdef VERBOSE
                if (verbose_flag > 2)
                {
                    printf("ODE_system.cpp -- reset_ODE_dots -- binary i %d de_vec_dt %g %g %g dh_vec_dt %g %g %g factor_h_vec %g\n",p->index,p->de_vec_dt[0],p->de_vec_dt[1],p->de_vec_dt[2],p->dh_vec_dt[0],p->dh_vec_dt[1],p->dh_vec_dt[2],factor_h_vec);
                }
                #endif

            }
            else if (p->integration_method==1)
            {
                for (i=0; i<4; i++)
                {        
                    p->KS_dalpha_vec_dt[i] = 0.0;
                    p->KS_dbeta_vec_dt[i] = 0.0;
                }

                for (i=0; i<3; i++)
                {        
                    p->a_vec[i] = 0.0; /* perturbing acceleration (without Keplerian part) */
                }

                p->KS_domega_dt = 0.0;
                p->KS_dE_dt = 0.0;

            }
            else if (p->integration_method==2)
            {

                double common_factor = -CONST_G*p->mass*p->r_pm3;
                for (i=0; i<3; i++)
                {        
                    p->a_vec[i] = common_factor*p->r_vec[i];
                }
            }
        }
    }
}


void write_ODE_variables_dots(ParticlesMap *particlesMap, N_Vector &y_dot)
{
    ParticlesMapIterator it_p;
    int k=1;
    int k_component;

    double spin_vec[3],e_vec[3],h_vec[3];
    for (it_p = particlesMap->begin(); it_p != particlesMap->end(); it_p++)
    {
        Particle *p = (*it_p).second;
        if (p->is_binary == false)
        {

            /* Limit to break-up rotation */
            double f_breakup = 1.0;
            double Omega_crit = compute_breakup_angular_frequency(p->mass,p->radius);

            if (p->spin_vec_norm > Omega_crit)
            {
                f_breakup = 0.0;
            }

            for (k_component=0; k_component<3; k_component++)
            {
                Ith(y_dot,k + k_component) = f_breakup * p->dspin_vec_dt[k_component];
                check_number(Ith(y_dot,k + k_component),                  "ODE_system.cpp -- write_ODE_variables_dots","spins", true);
            }
            
            Ith(y_dot,k + 2 + 1) = p->dmass_dt;
            Ith(y_dot,k + 2 + 2) = p->dradius_dt;
            check_number(Ith(y_dot,k + 2 + 1),                  "ODE_system.cpp -- write_ODE_variables_dots","dmass_dt", true);
            check_number(Ith(y_dot,k + 2 + 2),                  "ODE_system.cpp -- write_ODE_variables_dots","dradius_dt", true);
            
            k=k+5;

            #ifdef VERBOSE
            if (verbose_flag > 2)
            {
                printf("ODE_system.cpp -- write_ODE_variables_dots -- body i %d dspin_vec_dt %g %g %g dmass_dt %g dradius_dt %g\n",p->index,p->dspin_vec_dt[0],p->dspin_vec_dt[1],p->dspin_vec_dt[2],p->dmass_dt,p->dradius_dt);
            }
            #endif
        }
        if (p->is_binary == true) // particle is a binary
        {
            if (p->integration_method==0)
            {
                for (k_component=0; k_component<3; k_component++)
                {
                    Ith(y_dot,k + k_component)      = p->de_vec_dt[k_component];
                    Ith(y_dot,k + k_component + 3)  = p->dh_vec_dt[k_component];
                    check_number(Ith(y_dot,k + k_component),                  "ODE_system.cpp -- write_ODE_variables_dots","de_dt", true);
                    check_number(Ith(y_dot,k + k_component + 3),                  "ODE_system.cpp -- write_ODE_variables_dots","dh_dt", true);
                }

                #ifdef VERBOSE
                if (verbose_flag > 2)
                {
                    printf("ODE_system.cpp -- write_ODE_variables_dots -- binary i %d de_vec_dt %g %g %g dh_vec_dt %g %g %g\n",p->index,p->de_vec_dt[0],p->de_vec_dt[1],p->de_vec_dt[2],p->dh_vec_dt[0],p->dh_vec_dt[1],p->dh_vec_dt[2]);
                }
                #endif
                
                k=k+6;
            }
            else if (p->integration_method==1)
            {
                for (k_component=0; k_component<4; k_component++)
                {
                    Ith(y_dot,k + k_component)      = p->KS_dalpha_vec_dt[k_component];
                    Ith(y_dot,k + k_component + 4)  = p->KS_dbeta_vec_dt[k_component];

                }
                Ith(y_dot,k + 7 + 1)  = p->KS_domega_dt;
                Ith(y_dot,k + 7 + 2)  = p->KS_dE_dt;

                #ifdef VERBOSE
                if (verbose_flag > 2)
                {
                    printf("ODE_system.cpp -- write_ODE_variables_dots -- binary i %d KS_dalpha_vec_dt %g %g %g %g p->KS_dbeta_vec_dt %g %g %g %g\n",p->index,p->KS_dalpha_vec_dt[0],p->KS_dalpha_vec_dt[1],p->KS_dalpha_vec_dt[2],p->KS_dalpha_vec_dt[3], p->KS_dbeta_vec_dt[0], p->KS_dbeta_vec_dt[1], p->KS_dbeta_vec_dt[2], p->KS_dbeta_vec_dt[3]);
                    printf("ODE_system.cpp -- write_ODE_variables_dots -- binary i %d p->KS_domega_dt %g p->KS_dE_dt %g\n",p->index,p->KS_domega_dt,p->KS_dE_dt);
                }
                #endif
                
                k=k+10;
            }
            else if (p->integration_method==2)
            {
                for (k_component=0; k_component<3; k_component++)
                {
                    Ith(y_dot,k + k_component)      = p->v_vec[k_component];
                    Ith(y_dot,k + k_component + 3)  = p->a_vec[k_component];
                }

                #ifdef VERBOSE
                if (verbose_flag > 2)
                {
                    printf("ODE_system.cpp -- write_ODE_variables_dots -- binary i %d v_vec %g %g %g a_vec %g %g %g\n",p->index,p->v_vec[0],p->v_vec[1],p->v_vec[2],p->a_vec[0],p->a_vec[1],p->a_vec[2]);
                }
                #endif
                
                k=k+6;
            }
            
        }
    }
}


void set_initial_ODE_variables(ParticlesMap *particlesMap, N_Vector &y, N_Vector &y_abs_tol, double abs_tol_spin_vec, double abs_tol_e_vec, double abs_tol_h_vec)
{
    ParticlesMapIterator it_p;
    int k=1;
    int k_component;

    double spin_vec[3],e_vec[3],h_vec[3];
    double alpha_vec[4],beta_vec[4];
    
    for (it_p = particlesMap->begin(); it_p != particlesMap->end(); it_p++)
    {
        Particle *p = (*it_p).second;
        if (p->is_binary == false)
        {

            for (k_component=0; k_component<3; k_component++)
            {
                Ith(y,         k + k_component) = p->spin_vec[k_component];
                Ith(y_abs_tol, k + k_component) = abs_tol_spin_vec;
            }
            Ith(y,  k + 2 + 1) = p->mass;
            Ith(y,  k + 2 + 2) = p->radius;
            
            Ith(y_abs_tol, k + 2 + 1) = relative_tolerance*p->mass;
            Ith(y_abs_tol, k + 2 + 2) = relative_tolerance*p->radius;

            k=k+5;
            
            #ifdef VERBOSE
            if (verbose_flag > 2)
            {
                printf("ODE_system.cpp -- set_initial_ODE_variables -- body i %d spin_vec %g %g %g mass %g radius %g\n",p->index,p->spin_vec[0],p->spin_vec[1],p->spin_vec[2],p->mass,p->radius);
                printf("ODE_system.cpp -- set_initial_ODE_variables -- body i %d ABS TOL spin_vec %g %g %g mass %g radius %g\n",p->index,abs_tol_spin_vec,abs_tol_spin_vec,abs_tol_spin_vec,relative_tolerance*p->mass,relative_tolerance*p->radius);
            }
            #endif
        }
        if (p->is_binary == true) // particle is a binary
        {
            if (p->integration_method==0)
            {
                double h = norm3(p->h_vec);
        
                for (k_component=0; k_component<3; k_component++)
                {
                    Ith(y,k + k_component) = p->e_vec[k_component];
                    Ith(y,k + k_component + 3) = p->h_vec[k_component];
                    
                    Ith(y_abs_tol,k + k_component) = abs_tol_e_vec;
                    Ith(y_abs_tol,k + k_component + 3) = relative_tolerance*h;
                }
                
                k=k+6;
                
                #ifdef VERBOSE
                if (verbose_flag > 2)
                {
                    printf("ODE_system.cpp -- set_initial_ODE_variables -- binary i %d e_vec %g %g %g h_vec %g %g %g\n",p->index,p->e_vec[0],p->e_vec[1],p->e_vec[2],p->h_vec[0],p->h_vec[1],p->h_vec[2]);
                    printf("ODE_system.cpp -- set_initial_ODE_variables -- binary i %d ABS TOL e_vec %g %g %g h_vec %g %g %g\n",p->index,abs_tol_e_vec,abs_tol_e_vec,abs_tol_e_vec,relative_tolerance*h,relative_tolerance*h,relative_tolerance*h);
                }
                #endif
            }
            else if (p->integration_method==1)
            {

                for (k_component=0; k_component<4; k_component++)
                {
                    Ith(y,k + k_component) = p->KS_alpha_vec[k_component];
                    Ith(y,k + k_component + 4) = p->KS_beta_vec[k_component];
                    
                    Ith(y_abs_tol,k + k_component) = relative_tolerance;
                    Ith(y_abs_tol,k + k_component + 4) = relative_tolerance;
                }
                Ith(y,  k + 7 + 1) = p->KS_omega;
                Ith(y,  k + 7 + 2) = p->KS_E;

                Ith(y_abs_tol, k + 7 + 1) = relative_tolerance;
                Ith(y_abs_tol, k + 7 + 2) = relative_tolerance;
                
                k=k+10;

                #ifdef VERBOSE
                if (verbose_flag > 2)
                {
                    printf("ODE_system.cpp -- set_initial_ODE_variables -- binary i %d alpha_vec %g %g %g %g beta_vec %g %g %g %g\n",p->index,p->KS_alpha_vec[0],p->KS_alpha_vec[1],p->KS_alpha_vec[2],p->KS_alpha_vec[3],p->KS_beta_vec[0],p->KS_beta_vec[1],p->KS_beta_vec[2],p->KS_beta_vec[3]);
                }
                #endif
            }
            else if (p->integration_method==2)
            {
                for (k_component=0; k_component<3; k_component++)
                {
                    Ith(y,k + k_component) = p->r_vec[k_component];
                    Ith(y,k + k_component + 3) = p->v_vec[k_component];
                    
                    Ith(y_abs_tol,k + k_component) = relative_tolerance;
                    Ith(y_abs_tol,k + k_component + 3) = relative_tolerance;
                }
                k=k+6;

                #ifdef VERBOSE
                if (verbose_flag > 2)
                {
                    printf("ODE_system.cpp -- set_initial_ODE_variables -- binary i %d r %g %g %g v %g %g %g\n",p->index,p->r_vec[0],p->r_vec[1],p->r_vec[2],p->v_vec[0],p->v_vec[1],p->v_vec[2]);
                }
                #endif
            }
        }
    }
}


void extract_final_ODE_variables(ParticlesMap *particlesMap, N_Vector &y_out)
{
    ParticlesMapIterator it_p;
    int k=1;
    int k_component;

    double spin_vec[3],e_vec[3],h_vec[3];
    
    for (it_p = particlesMap->begin(); it_p != particlesMap->end(); it_p++)
    {
        Particle *p = (*it_p).second;
        if (p->is_binary == false)
        {
            for (k_component=0; k_component<3; k_component++)
            {
                p->spin_vec[k_component] = Ith(y_out,k + k_component);
            }

            p->mass = Ith(y_out,k + 2 + 1);
            p->radius = Ith(y_out,k + 2 + 2);

            k=k+5;
        }
        if (p->is_binary == true)
        {
            if (p->integration_method==0)
            {
                for (k_component=0; k_component<3; k_component++)
                {
                    p->e_vec[k_component] = Ith(y_out,k + k_component);
                    p->h_vec[k_component] = Ith(y_out,k + k_component + 3);
                }
                k=k+6;
                
                if (norm3(p->e_vec) > 1.0)
                {
                   
                    #ifdef VERBOSE
                    if (verbose_flag > 0)
                    {
                        printf("ODE_system.cpp -- extract_final_ODE_variables -- p %d e = %g > 1; original e_vec %g %g %g\n",p->index,norm3(p->e_vec),p->e_vec[0],p->e_vec[1],p->e_vec[2]); 
                    }
                    #endif

                    double e_norm = norm3(p->e_vec);
                    double e_overshoot = e_norm - 1.0;
                    for (k_component=0; k_component<3; k_component++)
                    {
                        p->e_vec[k_component] = fmod(1.0 - e_overshoot,1.0) * p->e_vec[k_component]/e_norm;
                    }

                    #ifdef VERBOSE
                    if (verbose_flag > 0)
                    {
                        printf("ODE_system.cpp -- extract_final_ODE_variables -- p %d; new adjusted e = %g; e_vec %g %g %g\n",p->index,norm3(p->e_vec),p->e_vec[0],p->e_vec[1],p->e_vec[2]); 
                    }
                    #endif


                }
                
            }
            else if (p->integration_method==1)
            {
                for (k_component=0; k_component<4; k_component++)
                {
                    p->KS_alpha_vec[k_component] = Ith(y_out,k + k_component);
                    p->KS_beta_vec[k_component] = Ith(y_out,k + k_component + 4);
                }
                
                p->KS_omega = Ith(y_out,k + 7 + 1);
                p->KS_E = Ith(y_out,k + 7 + 2);
                
                k=k+10;
            }
            else if (p->integration_method==2)
            {
                for (k_component=0; k_component<3; k_component++)
                {
                    p->r_vec[k_component] = Ith(y_out,k + k_component);
                    p->v_vec[k_component] = Ith(y_out,k + k_component + 3);
                }
               
                k=k+6;
            }
        }
    }
}

void check_for_integration_exclusion_orbits(ParticlesMap *particlesMap)
{
    set_up_derived_quantities(particlesMap);
    
    ParticlesMapIterator it_p;
    std::vector<int>::iterator it_parent_p,it_parent_q;
    int k=1;
    int k_component;

    for (it_p = particlesMap->begin(); it_p != particlesMap->end(); it_p++)
    {
        Particle *p = (*it_p).second;
        if (p->is_binary == true)
        {

            for (it_parent_p = p->parents.begin(); it_parent_p != p->parents.end(); it_parent_p++)
            {
               
                int i = std::distance(p->parents.begin(), it_parent_p);
                Particle *q = (*particlesMap)[(*it_parent_p)];
                int connecting_child_in_parent_q = p->connecting_child_in_parents[i];

                p->exclude_for_secular_integration = false;

                double t_sec = compute_order_of_magnitude_secular_timescale_for_pair(particlesMap,p->index,q->index,connecting_child_in_parent_q);

                /* 1PN apsidal motion */
                double t_1PN = compute_1PN_timescale(p->a,p->mass,p->e);
                if (t_1PN < t_sec * secular_integration_exclusion_safety_factor)
                {
                    p->exclude_for_secular_integration = true;
                    
                    #ifdef VERBOSE
                    if (verbose_flag > 1)
                    {
                        printf("ODE_system.cpp -- check_for_integration_exclusion_orbits -- 1PN excluding p %d parent %d t_sec %g t_1PN %g\n",p->index,q->index,t_sec,t_1PN);
                    }
                    #endif
                }
                
                /* Apsidal motion due to tides */
                double t_apsidal_motion1;
                double t_apsidal_motion2;
                Particle *child1 = (*particlesMap)[p->child1];
                Particle *child2 = (*particlesMap)[p->child2];
                compute_estimated_tidal_apsidal_motion_timescales(child1->mass,child2->mass,child1->radius,child1->apsidal_motion_constant,child1->spin_vec,p->a,p->e,p->e_vec,p->h_vec,&t_apsidal_motion1);
                compute_estimated_tidal_apsidal_motion_timescales(child2->mass,child1->mass,child2->radius,child2->apsidal_motion_constant,child2->spin_vec,p->a,p->e,p->e_vec,p->h_vec,&t_apsidal_motion2);
                double t_apsidal_motion_eff = CV_min(t_apsidal_motion1,t_apsidal_motion2);
                
                if (t_apsidal_motion_eff < t_sec * secular_integration_exclusion_safety_factor)
                {
                    p->exclude_for_secular_integration = true;
                    
                    #ifdef VERBOSE
                    if (verbose_flag > 1)
                    {
                        printf("ODE_system.cpp -- check_for_integration_exclusion_orbits -- tides excluding p %d parent %d t_sec %g t_apsidal_motion_eff %g\n",p->index,q->index,t_sec,t_apsidal_motion_eff);
                    }
                    #endif
                }
            }
        }
    }
}
            

/* function to check ODE solver-related function return values */
static int check_flag(void *flagvalue, char *funcname, int opt)
{
  int *errflag;

  /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
  if (opt == 0 && flagvalue == NULL) {
    fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n",
	    funcname);
    return(1); }

  /* Check if flag < 0 */
  else if (opt == 1) {
    errflag = (int *) flagvalue;
    if (*errflag < 0) {
      fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed with flag = %d\n\n",
	      funcname, *errflag);
      return(1); }}

  /* Check if function returned NULL pointer - no memory allocated */
  else if (opt == 2 && flagvalue == NULL) {
    fprintf(stderr, "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n",
	    funcname);
    return(1); }

  return 0;
}


}
