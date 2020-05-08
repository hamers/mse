#include "types.h"
#include "evolve.h"
#include "ODE_system.h"

extern "C"
{
int compute_y_dot(realtype time, N_Vector y, N_Vector y_dot, void *data_)
{
	UserData data;
	data = (UserData) data_;
    ParticlesMap *particlesMap = data->particlesMap;
    
    double start_time = data->start_time;
    double delta_time = time - start_time;

    #ifdef DEBUG
    printf("ODE_system.cpp -- compute_y_dot t=%g start_time=%g delta_time = %g\n",time,start_time,delta_time);
    #endif
    
    extract_ODE_variables(particlesMap, y, delta_time);
    reset_ODE_dots(particlesMap, y, delta_time);
    
    /****************************
     * compute right-hand sides *
     * **************************/

    double hamiltonian = 0.0;
    double KS_V = 0.0;
    ParticlesMapIterator it_p,it_f;
    std::vector<int>::iterator it_parent_p,it_parent_q;
    
    for (it_p = particlesMap->begin(); it_p != particlesMap->end(); it_p++)
    {
        Particle *p = (*it_p).second;
        if (p->is_binary == true)
        {

            /* Newtonian gravitational point mass dynamics (internal system) */
            compute_EOM_Newtonian_for_particle(particlesMap,p,&hamiltonian,&KS_V,false);
            
            /* Perturbations by flybys (external perturbations) */
            compute_EOM_Newtonian_external_perturbations(time,particlesMap,p,&hamiltonian,&KS_V,false); 

            /* VRR-related perturbations */

            if (p->VRR_model > 0) /* 0: no VRR perturbation */
            {
                hamiltonian += compute_VRR_perturbations(particlesMap,p->index,time);
            }


            /* Pairwise PN corrections */
            if (p->include_pairwise_1PN_terms == true)
            {
                hamiltonian += compute_EOM_pairwise_1PN(particlesMap,p->index,false);
            }
            if (p->include_pairwise_25PN_terms == true)
            {
                hamiltonian += compute_EOM_pairwise_25PN(particlesMap,p->index,false);
            }
            
            /* Tidal friction (ad hoc) */
            
            Particle *child1 = (*particlesMap)[p->child1];
            Particle *child2 = (*particlesMap)[p->child2];

            
            if (child1->include_tidal_friction_terms == true || child1->include_tidal_bulges_precession_terms == true || child1->include_rotation_precession_terms == true)
            {
                if (child1->tides_method == 0 || child1->tides_method == 1)
                {
                    compute_EOM_equilibrium_tide_BO_full(particlesMap,p->index,child1->index,child2->index,child1->include_tidal_friction_terms,child1->include_tidal_bulges_precession_terms,child1->include_rotation_precession_terms,child1->minimum_eccentricity_for_tidal_precession,child1->tides_method);
                }
                else if (child1->tides_method == 2)
                {
                    compute_EOM_equilibrium_tide(particlesMap,p->index,child1->index,child2->index,child1->include_tidal_friction_terms,child1->include_tidal_bulges_precession_terms,child1->include_rotation_precession_terms,child1->minimum_eccentricity_for_tidal_precession);
                }
            }
            if (child2->include_tidal_friction_terms == true || child2->include_tidal_bulges_precession_terms == true || child2->include_rotation_precession_terms == true)
            {
                if (child2->tides_method == 0 || child2->tides_method == 1)
                {
                    compute_EOM_equilibrium_tide_BO_full(particlesMap,p->index,child2->index,child1->index,child2->include_tidal_friction_terms,child2->include_tidal_bulges_precession_terms,child2->include_rotation_precession_terms,child2->minimum_eccentricity_for_tidal_precession,child2->tides_method);
                }
                else if (child2->tides_method == 2)
                {
                    compute_EOM_equilibrium_tide(particlesMap,p->index,child2->index,child1->index,child2->include_tidal_friction_terms,child2->include_tidal_bulges_precession_terms,child2->include_rotation_precession_terms,child2->minimum_eccentricity_for_tidal_precession);
                }

            }
            #ifdef IGNORE
            #endif
            

            /* Orbital changes due to mass changes */
            
            //ODE_handle_stellar_winds(p);
            //ODE_handle_RLOF(p,child1,child2);
            
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
    int highest_level = (*particlesMap)[0]->highest_level;
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
    set_up_derived_ODE_quantities(particlesMap);

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
                new_MA = remainder(new_MA,2.0*M_PI);

                p->true_anomaly = compute_true_anomaly_from_mean_anomaly(new_MA,norm3(p->e_vec));
                from_orbital_vectors_to_cartesian(p->child1_mass,p->child2_mass,p->e_vec,p->h_vec,p->true_anomaly,p->r_vec,p->v_vec);
                
                #ifdef DEBUG
                printf("ODE_system.cpp -- process_direct_integration_quantities -- n %g new_MA %g old_MA %g\n",n,new_MA,p->initial_mean_anomaly);
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
    set_up_derived_ODE_quantities(particlesMap);
}
    
void set_up_derived_ODE_quantities(ParticlesMap *particlesMap)
{
    /* These are derived quantities that are often used in the EOM, so they are calculated here once for speed up */
    int i;
    ParticlesMapIterator it_p;
    for (it_p = particlesMap->begin(); it_p != particlesMap->end(); it_p++)
    {
        Particle *p = (*it_p).second;
        if (p->is_binary == false)
        {
            p->spin_vec_norm = norm3(p->spin_vec);
        }
        else
        {
            p->e = norm3(p->e_vec);
            p->h = norm3(p->h_vec);

            for (i=0; i<3; i++)
            {        
                p->e_vec_unit[i] = p->e_vec[i]/p->e;
                p->h_vec_unit[i] = p->h_vec[i]/p->h;
            }
            
            cross3(p->h_vec_unit,p->e_vec_unit,p->q_vec_unit);
            
            p->e_p2 = p->e*p->e;
            p->j_p2 = 1.0 - p->e_p2;
            p->j = sqrt(p->j_p2);
            p->j_p3 = p->j*p->j_p2;
            p->j_p4 = p->j*p->j_p3;
            p->j_p5 = p->j*p->j_p4;

            p->a = p->h*p->h*p->child1_mass_plus_child2_mass/( CONST_G*p->child1_mass_times_child2_mass*p->child1_mass_times_child2_mass*p->j_p2 );

            p->r = norm3(p->r_vec);
            p->r_p2 = p->r*p->r;
            p->r_p3 = p->r*p->r_p2;
            p->r_pm1 = 1.0/p->r;
            p->r_pm2 = p->r_pm1*p->r_pm1;
            p->r_pm3 = p->r_pm1*p->r_pm2;
            
        }
    }
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
            p->dmass_dt = p->mass_dot_wind; //+ p->mass_dot_RLOF;
            p->dradius_dt = p->radius_dot + p->radius_ddot*delta_time;
            
            #ifdef DEBUG
            printf("ODE_system.cpp -- reset_ODE_dots -- body i %d dmass_dt %g dradius_dt %g dspin_vec_dt %g %g %g\n",p->index,p->dmass_dt,p->dradius_dt,p->dspin_vec_dt[0],p->dspin_vec_dt[1],p->dspin_vec_dt[2]);
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

                #ifdef DEBUG
                printf("ODE_system.cpp -- reset_ODE_dots -- binary i %d de_vec_dt %g %g %g dh_vec_dt %g %g %g factor_h_vec %g\n",p->index,p->de_vec_dt[0],p->de_vec_dt[1],p->de_vec_dt[2],p->dh_vec_dt[0],p->dh_vec_dt[1],p->dh_vec_dt[2],factor_h_vec);
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
            for (k_component=0; k_component<3; k_component++)
            {
                Ith(y_dot,k + k_component) = p->dspin_vec_dt[k_component];
            }
            
            Ith(y_dot,k + 2 + 1) = p->dmass_dt;
            Ith(y_dot,k + 2 + 2) = p->dradius_dt;
            
            k=k+5;
            
            #ifdef DEBUG
            printf("ODE_system.cpp -- write_ODE_variables_dots -- body i %d dspin_vec_dt %g %g %g dmass_dt %g dradius_dt %g\n",p->index,p->dspin_vec_dt[0],p->dspin_vec_dt[1],p->dspin_vec_dt[2],p->dmass_dt,p->dradius_dt);
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
                }

                #ifdef DEBUG
                printf("ODE_system.cpp -- write_ODE_variables_dots -- binary i %d de_vec_dt %g %g %g dh_vec_dt %g %g %g\n",p->index,p->de_vec_dt[0],p->de_vec_dt[1],p->de_vec_dt[2],p->dh_vec_dt[0],p->dh_vec_dt[1],p->dh_vec_dt[2]);
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

                #ifdef DEBUG
                printf("ODE_system.cpp -- write_ODE_variables_dots -- binary i %d KS_dalpha_vec_dt %g %g %g %g p->KS_dbeta_vec_dt %g %g %g %g\n",p->index,p->KS_dalpha_vec_dt[0],p->KS_dalpha_vec_dt[1],p->KS_dalpha_vec_dt[2],p->KS_dalpha_vec_dt[3], p->KS_dbeta_vec_dt[0], p->KS_dbeta_vec_dt[1], p->KS_dbeta_vec_dt[2], p->KS_dbeta_vec_dt[3]);
                printf("ODE_system.cpp -- write_ODE_variables_dots -- binary i %d p->KS_domega_dt %g p->KS_dE_dt %g\n",p->index,p->KS_domega_dt,p->KS_dE_dt);
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

                #ifdef DEBUG
                printf("ODE_system.cpp -- write_ODE_variables_dots -- binary i %d v_vec %g %g %g a_vec %g %g %g\n",p->index,p->v_vec[0],p->v_vec[1],p->v_vec[2],p->a_vec[0],p->a_vec[1],p->a_vec[2]);
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
            
            #ifdef DEBUG
            printf("ODE_system.cpp -- set_initial_ODE_variables -- body i %d spin_vec %g %g %g mass %g radius %g\n",p->index,p->spin_vec[0],p->spin_vec[1],p->spin_vec[2],p->mass,p->radius);
            printf("ODE_system.cpp -- set_initial_ODE_variables -- body i %d ABS TOL spin_vec %g %g %g mass %g radius %g\n",p->index,abs_tol_spin_vec,abs_tol_spin_vec,abs_tol_spin_vec,relative_tolerance*p->mass,relative_tolerance*p->radius);
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
                
                #ifdef DEBUG
                printf("ODE_system.cpp -- set_initial_ODE_variables -- binary i %d e_vec %g %g %g h_vec %g %g %g\n",p->index,p->e_vec[0],p->e_vec[1],p->e_vec[2],p->h_vec[0],p->h_vec[1],p->h_vec[2]);
                printf("ODE_system.cpp -- set_initial_ODE_variables -- binary i %d ABS TOL e_vec %g %g %g h_vec %g %g %g\n",p->index,abs_tol_e_vec,abs_tol_e_vec,abs_tol_e_vec,relative_tolerance*h,relative_tolerance*h,relative_tolerance*h);
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

                #ifdef DEBUG
                printf("ODE_system.cpp -- set_initial_ODE_variables -- binary i %d alpha_vec %g %g %g %g beta_vec %g %g %g %g\n",p->index,p->KS_alpha_vec[0],p->KS_alpha_vec[1],p->KS_alpha_vec[2],p->KS_alpha_vec[3],p->KS_beta_vec[0],p->KS_beta_vec[1],p->KS_beta_vec[2],p->KS_beta_vec[3]);
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

                #ifdef DEBUG
                printf("ODE_system.cpp -- set_initial_ODE_variables -- binary i %d r %g %g %g v %g %g %g\n",p->index,p->r_vec[0],p->r_vec[1],p->r_vec[2],p->v_vec[0],p->v_vec[1],p->v_vec[2]);
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

}
