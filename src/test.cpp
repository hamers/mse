/* MSE */

#include "evolve.h"
#include "test.h"
//#include "mstar/regularization.h"

extern "C"
{


/*********
 * Tools *
**********/

int test_tools()
{
    printf("test.cpp -- test_tools\n");
   
    int flag = 0;
    flag += test_a_P_orb_conversion();
    flag += test_a_h_conversion();
    flag += test_orbital_element_conversion();
    flag += test_kepler_equation_solver();
    flag += test_orbital_vectors_cartesian_conversion();
    
    if (flag == 0)
    {
        printf("test.cpp -- test_tools -- passed\n");
    }
    return flag;
}

int test_a_P_orb_conversion()
{
    printf("test.cpp -- test_a_P_orb_conversion\n");
    
    double a = 1.0;
    double m = 2.0;
    
    double P = compute_orbital_period_from_semimajor_axis(m,a);
    double an = compute_semimajor_axis_from_orbital_period(m,P);

    double tol = 1.0e-15;
    int flag = 0;
    if ( !equal_number(a,an,tol) )
    {
        printf("test.cpp -- error in test_a_P_orb_conversion!\n");
        flag = 1;
    }
    
    if (flag == 0)
    {
        printf("test.cpp -- test_a_P_orb_conversion -- passed\n");
    }
    return flag;
}

int test_a_h_conversion()
{
    printf("test.cpp -- test_a_h_conversion\n");
    
    double a = 1.0;
    double m1 = 2.0;
    double m2 = 3.0;
    double e = 0.5;
    
    double h = compute_h_from_a(m1,m2,a,e);
    double an = compute_a_from_h(m1,m2,h,e);
    
    double tol = 1.0e-15;
    int flag = 0;
    if ( !equal_number(a,an,tol) )
    {
        printf("test.cpp -- error in test_a_h_conversion!\n");
        flag = 1;
    }

    if (flag == 0)
    {
        printf("test.cpp -- test_a_h_conversion -- passed\n");
    }
    return flag;
}

int test_orbital_element_conversion()
{
    printf("test.cpp -- test_orbital_element_conversion\n");
   
    double a,m1,m2,e,INCL,AP,LAN;
    double an,en,INCLn,APn,LANn;
    double e_vec_x,e_vec_y,e_vec_z;
    double h_vec_x,h_vec_y,h_vec_z;
    double h_tot_vec[3];
    int N=1000;
    int flag = 0;
    double tol = 1.0e-8;
    for (int i=0; i<N; i++)
    {
        a = ((double) rand() / (RAND_MAX));
        m1 = ((double) rand() / (RAND_MAX));
        m2 = ((double) rand() / (RAND_MAX));
        e = ((double) rand() / (RAND_MAX));
        INCL = ((double) rand() / (RAND_MAX)) * M_PI;
        LAN = (2.0 * ((double) rand() / (RAND_MAX)) - 1.0) * M_PI;
        AP = (2.0 * ((double) rand() / (RAND_MAX)) - 1.0) * M_PI;

        compute_orbital_vectors_from_orbital_elements(m1,m2,a,e,INCL,AP,LAN,&e_vec_x,&e_vec_y,&e_vec_z,&h_vec_x,&h_vec_y,&h_vec_z);
        compute_orbital_elements_from_orbital_vectors(m1,m2,h_tot_vec,e_vec_x,e_vec_y,e_vec_z,h_vec_x,h_vec_y,h_vec_z,&an,&en,&INCLn,&APn,&LANn);
        
        if ( !equal_number(a,an,tol) or !equal_number(e,en,tol) or !equal_number(INCL,INCLn,tol) or !equal_number(AP,APn,tol) or !equal_number(LAN,LANn,tol) )
        {
            printf("test.cpp -- error in test_orbital_element_conversion! %g %g %g %g %g %g %g %g %g %g \n",a,an,e,en,INCL,INCLn,AP,APn,LAN,LANn);
            flag = 1;
        }
    }
    
    if (flag == 0)
    {
        printf("test.cpp -- test_orbital_element_conversion -- passed\n");
    }

    return flag;
}        

int test_orbital_vectors_cartesian_conversion()
{
    printf("test.cpp -- test_orbital_vectors_cartesian_conversion\n");
    
    double a,m1,m2,e,INCL,AP,LAN,TA,TAn;
    double e_vec[3],h_vec[3];
    double e_vecn[3],h_vecn[3];
    double r[3],v[3];
    double h_tot_vec[3];
    int N=1000;
    int flag = 0;
    double tol = 1.0e-10;
    for (int i=0; i<N; i++)
    {
        a = ((double) rand() / (RAND_MAX));
        m1 = ((double) rand() / (RAND_MAX));
        m2 = ((double) rand() / (RAND_MAX));
        e = ((double) rand() / (RAND_MAX));
        INCL = ((double) rand() / (RAND_MAX)) * M_PI;
        LAN = (2.0 * ((double) rand() / (RAND_MAX)) - 1.0) * M_PI;
        AP = (2.0 * ((double) rand() / (RAND_MAX)) - 1.0) * M_PI;
        TA = (2.0 * ((double) rand() / (RAND_MAX)) - 1.0) * M_PI;
        
        compute_orbital_vectors_from_orbital_elements(m1,m2,a,e,INCL,AP,LAN,&e_vec[0],&e_vec[1],&e_vec[2],&h_vec[0],&h_vec[1],&h_vec[2]);
        
        from_orbital_vectors_to_cartesian(m1,m2,e_vec,h_vec,TA,r,v);
        from_cartesian_to_orbital_vectors(m1,m2,r,v,e_vecn,h_vecn,&TAn);

        for (int j=0; j<3; j++)
        {
            if ( !equal_number(h_vec[j],h_vecn[j],tol) or !equal_number(e_vec[j],e_vecn[j],tol) or !equal_number(TA,TAn,tol))
            {
                printf("test.cpp -- error in test_orbital_vectors_cartesian_conversion!\n");
                flag = 1;
            }
        }
    }
    
    if (flag == 0)
    {
        printf("test.cpp -- test_orbital_vectors_cartesian_conversion -- passed\n");
    }

    return flag;
}        

int test_kepler_equation_solver()
{
    printf("test.cpp -- test_kepler_equation_solver\n");
    
    double e,MA,MAn,TA;
    int N=1000;
    int flag = 0;
    double tol = 1.0e-10;
    for (int i=0; i<N; i++)
    {
        e = ((double) rand() / (RAND_MAX));
        MA = (2.0 * ((double) rand() / (RAND_MAX)) - 1.0) * M_PI;

        TA = compute_true_anomaly_from_mean_anomaly(MA,e);
        MAn = compute_mean_anomaly_from_true_anomaly(TA,e);

        if ( !equal_number(MA,MAn,tol) )
        {
            printf("test.cpp -- error in test_kepler_equation_solver! %g %g\n",MA,MAn);
            flag = 1;
        }
    }
    
    if (flag == 0)
    {
        printf("test.cpp -- test_kepler_equation_solver -- passed\n");
    }

    return flag;
}


/**********
 * N-body *
***********/

int test_nbody()
{
    printf("test.cpp -- test_nbody\n");
    
    int flag = 0;
    flag += test_nbody_two_body_stopping_conditions();
    flag += test_nbody_two_body_kick();

    if (flag == 0)
    {
        printf("test.cpp -- test_nbody -- passed\n");
    }

    return flag;
    
}

int test_nbody_two_body_stopping_conditions()
{
    printf("test.cpp -- test_nbody_two_body_stopping_conditions\n");

    int flag = 0;    
    for (int stopping_condition_mode = 0; stopping_condition_mode < 2; stopping_condition_mode++)
    {
    
        initialize_mpi_or_serial(); 
        int N_bodies = 2;

        double m1 = 1.5;
        double m2 = 1.5;
        double R1 = 0.007;
        double R2 = 0.009;
        double a = 1.0;
        double e = 0.99;
        
        double a_init = a;
        double e_init = e;
        double gbs_tolerance = 1.0e-12;
        double stopping_condition_tolerance = 1.0e-12;
        struct RegularizedRegion *R = generate_binary_ICs(m1,m2,R1,R2,a,e,stopping_condition_mode,gbs_tolerance,stopping_condition_tolerance);


        double E_init = compute_nbody_total_energy(R);

        //print_state(R);

        int split_integration = 0;
        
        double tend = 1.0;
        
        int stopping_condition_occurred;
        double reached_dt;
        double t=0;

        if (split_integration==0)
        {
            run_integrator(R, tend, &reached_dt, &stopping_condition_occurred);
            t = reached_dt;
        }
        else
        {
            int Nsteps = 10;
            double dt = tend/( (double) Nsteps);
            while (t<=tend)
            {
                run_integrator(R, dt, &reached_dt, &stopping_condition_occurred);

                t+=reached_dt;
                
                //printf("t %g reached dt %g stopping_condition_occurred %d\n",t,reached_dt,stopping_condition_occurred);
                print_state(R);
                
                if (stopping_condition_occurred == 1)
                {
                    break;
                }
            }
        }
        double E_fin = compute_nbody_total_energy(R);

        //printf("Integration done reached t %g stopping_condition_occurred %d; t_end_an (collisions only) %g\n",t,stopping_condition_occurred,test_nbody_compute_time_of_collision(GCONST,m1+m2,R1+R2,a,e));
                
//        printf("Final state:\n");
//        print_state(R);
        
        /* Determine final orbital elements of binary for checking (independent of MSTAR) */
        double r[3],v[3];
        double M = m1+m2;
        for (int j=0; j<3; j++)
        {
            r[j] = R->Pos[3 * 0 + j] - R->Pos[3 * 1 + j];
            v[j] = R->Vel[3 * 0 + j] - R->Vel[3 * 1 + j];
        }
        
        //double r_RLOF1 = R1/fq_RLOF_Eggleton(m1,m2);
        //double r_RLOF2 = R2/fq_RLOF_Eggleton(m2,m1);
        double r_RLOF1 = R1/roche_radius_pericenter_eggleton(1.0,m1/m2);
        double r_RLOF2 = R2/roche_radius_pericenter_eggleton(1.0,m2/m1);
        
        double t_an = test_nbody_compute_time_of_collision(GCONST,m1+m2,R1+R2,a,e);
        
        test_nbody_compute_elements(GCONST,M,r,v,&a,&e);
        double r_norm = norm3(r);
        //printf("Final a %.15f e %.15f sep %.15f r_col %g r_RLOF1 %g r_RLOF2 %g\n",a,e,r_norm,R1+R2,r_RLOF1,r_RLOF2);

        free_data(R);

        double tol = 1.0e-10;
        if (!equal_number(E_init,E_fin,tol))
        {
            printf("test.cpp -- error in test_nbody_two_body_stopping_conditions -- E_init %g E_fin %g\n",E_init,E_fin);
            flag = 1;
        }
        
        if (!equal_number(a_init,a,tol))
        {
            printf("test.cpp -- error in test_nbody_two_body_stopping_conditions -- a %g a_init %g\n",a,a_init);
            flag = 1;
        }
        if (!equal_number(e_init,e,tol))
        {
            printf("test.cpp -- error in test_nbody_two_body_stopping_conditions -- e %g e_init %g\n",e,e_init);
            flag = 1;
        }

        if (stopping_condition_mode == 0)
        {
            if (!equal_number(r_norm,R1+R2,tol))
            {
                printf("test.cpp -- error in test_nbody_two_body_stopping_conditions -- r_norm %g R1+R2 %g\n",r_norm,R1+R2);
                flag = 1;
            }
            if (!equal_number(t,t_an,tol))
            {
                printf("test.cpp -- error in test_nbody_two_body_stopping_conditions -- t %g t_an %g\n",t,t_an);
                flag = 1;
            }
        }
        else if (stopping_condition_mode == 1)
        {
            if (!equal_number(r_norm,CV_max(r_RLOF1,r_RLOF2),tol))
            {
                printf("test.cpp -- error in test_nbody_two_body_stopping_conditions -- r_norm %g max(r_RLOF1,r_RLOF2) %g\n",r_norm,CV_max(r_RLOF1,r_RLOF2));
                flag = 1;
            }
        }
        
    }
    
    if (flag == 0)
    {
        printf("test.cpp -- test_nbody_two_body_stopping_conditions -- passed\n");
    }

    return flag;
}

double test_nbody_compute_time_of_collision(double CONST_G, double M, double r_col, double a, double e)
{
    double cos_theta = (1.0/e) * ( (a/r_col) * (1.0 - e*e) - 1.0 );
    double theta = acos(cos_theta);

    double mean_anomaly = M_PI - compute_mean_anomaly_from_true_anomaly(theta,e); // we started the orbit at apocenter!
    double n = sqrt(CONST_G*M / (a*a*a) );
    double Delta_t = mean_anomaly/n;

    return Delta_t;
}

int test_nbody_compute_elements(double CONST_G, double M, double *r, double *v, double *a, double *e)
{
    /* Not coded in the most efficient way, but suffices for testing */
    
    double E = 0.5*dot3(v,v) - CONST_G*M/norm3(r);
    *a = -CONST_G*M/(2.0*E);
    
    double e_vec[3],h[3],v_cross_h[3];
    cross3(r,v,h);
    cross3(v,h,v_cross_h);
    
    for (int i=0; i<3; i++)
    {
        e_vec[i] = v_cross_h[i]/(CONST_G*M) - r[i]/norm3(r);
    }
    *e = norm3(e_vec);
    
    return 0;
}

struct RegularizedRegion *generate_binary_ICs(double m1, double m2, double R1, double R2, double a, double e, int stopping_condition_mode, double gbs_tolerance, double stopping_condition_tolerance)
{
     
    struct RegularizedRegion *R;

    int i=0;
    int j;
   
    int N = 2;
    
    double R_cm_vec[3] = {0.1,0.0,0.0};
    double V_cm_vec[3] = {0.2,0.0,0.0};
    double r_vec[3],v_vec[3];
    
    double e_vec[3] = {1.0,0.0,0.0};
    double q_vec[3] = {0.0,1.0,0.0};
    
    double M = m1+m2;

    allocate_armst_structs(&R, N); // Initialize the data structure
    
    double R1_vec[3],R2_vec[3],V1_vec[3],V2_vec[3];
    
    double theta = M_PI;
    double r = a * (1.0 - e*e)/( 1.0 + e * cos(theta) );
    
    for (i=0; i<3; i++)
    {
        r_vec[i] = r * (cos(theta) * e_vec[i] + sin(theta) * q_vec[i]);
        v_vec[i] = sqrt(CONST_G * M/(a*(1.0-e*e))) * ( -sin(theta) * e_vec[i] + (e + cos(theta)) * q_vec[i] );

        R1_vec[i] = R_cm_vec[i] + (m2/M) * r_vec[i];
        V1_vec[i] = V_cm_vec[i] + (m2/M) * v_vec[i];

        R2_vec[i] = R_cm_vec[i] - (m1/M) * r_vec[i];
        V2_vec[i] = V_cm_vec[i] - (m1/M) * v_vec[i];

        R->Pos[3 * 0 + i] = R1_vec[i];
        R->Vel[3 * 0 + i] = V1_vec[i];

        R->Pos[3 * 1 + i] = R2_vec[i];
        R->Vel[3 * 1 + i] = V2_vec[i];
    }

    R->Mass[0] = m1;
    R->Mass[1] = m2;
    R->Radius[0] = R1;
    R->Radius[1] = R2;

    R->Stopping_Condition_Mode[0] = stopping_condition_mode;
    R->Stopping_Condition_Mode[1] = stopping_condition_mode;
    R->gbs_tolerance = gbs_tolerance;
    R->stopping_condition_tolerance = stopping_condition_tolerance;
    
   
    return R;
}

void compute_center_of_mass_position_and_velocity(struct RegularizedRegion *R, double R_cm[3], double V_cm[3])
{
    int i,j;
    
    for (j=0; j<3; j++)
    {
        R_cm[j] = 0.0;
        V_cm[j] = 0.0;
    }

    double m;
    double m_tot=0.0;
    for (i=0; i<R->NumVertex; i++)
    {
        m = R->Mass[i];
        m_tot += m;
        for (j=0; j<3; j++)
        {
            R_cm[j] += m * R->Pos[3 * i + j];
            V_cm[j] += m * R->Vel[3 * i + j];
        }
    }

    for (j=0; j<3; j++)
    {
        R_cm[j] = R_cm[j]/m_tot;
        V_cm[j] = V_cm[j]/m_tot;
    }
}

int test_nbody_two_body_kick()
{
    printf("test.cpp -- test_nbody_two_body_kick\n");

    int flag = 0;    
    int stopping_condition_mode = 0;
    
    initialize_mpi_or_serial(); 
    int N_bodies = 2;

    double m1 = 1.5;
    double m2 = 1.5;
    double R1 = 0.007;
    double R2 = 0.009;
    double a = 1.0;
    double e = 0.99;
    
    double gbs_tolerance = 1.0e-8;
    double stopping_condition_tolerance = 1.0e-6;
    struct RegularizedRegion *R = generate_binary_ICs(m1,m2,R1,R2,a,e,stopping_condition_mode,gbs_tolerance,stopping_condition_tolerance);

    
    /* Give kick to *2 */
    int i=2;
    int j=0;
    R->Vel[3 * j + i] += 1000.0*CONST_KM_PER_S;

    //double R_cm[3],V_cm[3];
    //compute_center_of_mass_position_and_velocity(R,R_cm,V_cm);
    //printf("pre int R_cm %g %g %g V_cm %g %g %g\n",R_cm[0],R_cm[1],R_cm[2],V_cm[0],V_cm[1],V_cm[2]);

    //compute_center_of_mass_position_and_velocity(R,R_cm,V_cm);
    //printf("post kick pre int %g %g %g %g %g %g\n",R_cm[0],R_cm[1],R_cm[2],V_cm[0],V_cm[1],V_cm[2]);

    //into_CoM_frame(R);
    double E_init = compute_nbody_total_energy(R);

    //print_state(R);

    double a_init,e_init;

    double r[3],v[3];
    double M = m1+m2;

    for (int j=0; j<3; j++)
    {
        r[j] = R->Pos[3 * 0 + j] - R->Pos[3 * 1 + j];
        v[j] = R->Vel[3 * 0 + j] - R->Vel[3 * 1 + j];
    }

    test_nbody_compute_elements(GCONST,M,r,v,&a_init,&e_init);
    
    int split_integration = 0;
    
    double tend = 10.0;
    
    int stopping_condition_occurred;
    double reached_dt;
    double t=0;

    if (split_integration==0)
    {
        run_integrator(R, tend, &reached_dt, &stopping_condition_occurred);
        t = reached_dt;
    }
    else
    {
        int Nsteps = 10;
        double dt = tend/( (double) Nsteps);
        while (t<=tend)
        {
            run_integrator(R, dt, &reached_dt, &stopping_condition_occurred);

            t+=reached_dt;
            
            //printf("t %g reached dt %g stopping_condition_occurred %d\n",t,reached_dt,stopping_condition_occurred);
            print_state(R);
            
            if (stopping_condition_occurred == 1)
            {
                break;
            }
        }
    }
    
    double E_fin = compute_nbody_total_energy(R);


    //compute_center_of_mass_position_and_velocity(R,R_cm,V_cm);
    //printf("post int R_cm %g %g %g V_cm %g %g %g\n",R_cm[0],R_cm[1],R_cm[2],V_cm[0],V_cm[1],V_cm[2]);
            
    //printf("Final state:\n");
    //print_state(R);
    
    /* Determine final orbital elements */
    for (int j=0; j<3; j++)
    {
        r[j] = R->Pos[3 * 0 + j] - R->Pos[3 * 1 + j];
        v[j] = R->Vel[3 * 0 + j] - R->Vel[3 * 1 + j];
    }
    
    test_nbody_compute_elements(GCONST,M,r,v,&a,&e);

    free_data(R);

    double tol = 1.0e-12;
    if (!equal_number(E_init,E_fin,tol))
    {
        printf("test.cpp -- error in test_nbody_two_body_kick -- E_init %g E_fin %g\n",E_init,E_fin);
        flag = 1;
    }

    if (!equal_number(a_init,a,tol))
    {
        printf("test.cpp -- error in test_nbody_two_body_kick -- a %g a_init %g\n",a,a_init);
        flag = 1;
    }
    if (!equal_number(e_init,e,tol))
    {
        printf("test.cpp -- error in test_nbody_two_body_kick -- e %g e_init %g\n",e,e_init);
        flag = 1;
    }

    if (flag == 0)
    {
        printf("test.cpp -- test_nbody_two_body_kick -- passed\n");
    }

    return flag;
}


/**********
 * Flybys *
***********/
    
int test_flybys()
{
    printf("test.cpp -- test_flybys\n");
    
    int flag=0;
    flag += test_flybys_integrals();
    flag += test_flybys_perturber_sampling();
    
    if (flag == 0)
    {
        printf("test.cpp -- test_flybys -- passed\n");
    }
    
    return flag;
}

int test_flybys_integrals()
{
    printf("test.cpp -- test_flybys_integrals\n");
    
    int flag = 0;
    flybys_mass_distribution = 0;
    flybys_mass_distribution_lower_value = 0.1;
    flybys_mass_distribution_upper_value = 100.0;
    flybys_mass_distribution = 0;
    flybys_encounter_sphere_radius = 1.0e5;
    flybys_stellar_density = 0.1*CONST_PER_PC3;
    flybys_stellar_relative_velocity_dispersion = 1.0*CONST_KM_PER_S;
    flybys_internal_mass = 20.0;
    
    double total_encounter_rate, stellar_density;
    compute_total_encounter_rate_and_density_at_R_enc(&total_encounter_rate,&stellar_density);
    
    double integral_rate = total_encounter_rate/(2.0*sqrt(2.0*M_PI)*flybys_encounter_sphere_radius*flybys_encounter_sphere_radius*flybys_stellar_density*flybys_stellar_relative_velocity_dispersion);
    double integral_density = stellar_density/flybys_stellar_density;
    
    if (!equal_number(integral_rate,1.181,1.0e-3))
    {
        printf("test.cpp -- test_flybys_integrals -- error in integral_rate %g\n",integral_rate);
        flag = 1;
    }
    if (!equal_number(integral_density,1.136,1.0e-3))
    {
        printf("test.cpp -- test_flybys_integrals -- error in integral_density %g\n",integral_density);
        flag = 1;
    }
    
    double W_max = compute_W_max();
    if (!equal_number(W_max,1.583,1.0e-3))
    {
        printf("test.cpp -- test_flybys_integrals -- error in W_max %g\n",W_max);
        flag = 1;
    }

    return flag;
}

int test_flybys_perturber_sampling()
{
    int flag = 0;
    flybys_mass_distribution = 0;
    flybys_mass_distribution_lower_value = 0.1;
    flybys_mass_distribution_upper_value = 100.0;
    flybys_mass_distribution = 0;
    flybys_encounter_sphere_radius = 1.0e5;
    flybys_stellar_density = 0.1*CONST_PER_PC3;
    flybys_stellar_relative_velocity_dispersion = 30.0*CONST_KM_PER_S;
    flybys_internal_mass = 20.0;
        
    ParticlesMap particlesMap;
    int N_bodies = 2;
    double masses[4] = {10.0,12.0};
    int stellar_types[4] = {1,1};
    double smas[3] = {10.0};
    double es[3] = {0.01};
    double TAs[3] = {0.01};
    double INCLs[3] = {0.01};
    double APs[3] = {0.01};
    double LANs[3] = {0.01};
    random_seed=0;
    create_nested_system(particlesMap,N_bodies,masses,stellar_types,smas,es,TAs,INCLs,APs,LANs);

    int N_points = 10000;
    bool apply_flyby;
    double t_next_encounter;
    int N_enc=0;
    int N_not_impulsive=0;
    double M_per;
    double b_vec[3],V_vec[3];

    double b_vec_mean[3];
    double M_per_mean = 0.0;

    for (int k=0; k<3; k++)
    {
        b_vec_mean[k] = 0.0;
    }

    for (int i=0; i<N_points; i++)
    {
        sample_next_flyby(&particlesMap, &apply_flyby, &t_next_encounter, &N_enc, &N_not_impulsive, &M_per, b_vec, V_vec);
        printf("T R %g %g %g V %g %g %g\n",b_vec[0],b_vec[1],b_vec[2],V_vec[0],V_vec[1],V_vec[2]);
        M_per_mean += M_per;
        for (int k=0; k<3; k++)
        {
            b_vec_mean[k] += b_vec[k];
        }
    }

    M_per_mean /= ((double) N_points);
    for (int k=0; k<3; k++)
    {
        b_vec_mean[k] /= ((double) N_points);
    }


    printf("D %g %g %g M %g\n",b_vec_mean[0],b_vec_mean[1],b_vec_mean[2],M_per_mean);
    return flag;
}


/*********************
 * Stellar evolution *
**********************/

int test_stellar_evolution()
{
    printf("test.cpp -- test_stellar_evolution\n");
    
    int flag;
    flag += test_apsidal_motion_constant();
    flag += test_sse();
    
    if (flag == 0)
    {
        printf("test.cpp -- test_stellar_evolution -- passed\n");
    }
    
    return flag;
}

int test_apsidal_motion_constant()
{
    /* Currently tests for NaNs in a few cases. */
    printf("test.cpp -- test_apsidal_motion_constant\n");

    int flag = 0;    
    int N_m=5;
    int N_st=15;
    double masses[N_m] = {0.08,0.5,1.0,10.0,100.0};
    int stellar_types[N_st] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14};
    
    double k_AM;
    for (int i=0; i<N_m; i++)
    {
        for (int j=0; j<N_st; j++)
        {
            ParticlesMap particlesMap;
            Particle *star = new Particle(0, false);
            particlesMap[0] = star;
            
            star->evolve_as_star = true;
            star->mass = masses[i];
            star->sse_initial_mass = masses[i];
            star->stellar_type = stellar_types[j];

            initialize_stars(&particlesMap);
        
            k_AM = compute_apsidal_motion_constant(star);
            if (k_AM != k_AM or k_AM < 0.0)
            {
                printf("test.cpp -- test_apsidal_motion_constant -- ERROR: k_AM %g is NaN or <0; m_init %g kw %d age %g\n",k_AM,star->mass,star->stellar_type,star->age);
                flag = 1;
            }
        }
    }

    if (flag == 0)
    {
        printf("test.cpp -- test_apsidal_motion_constant -- passed\n");
    }

    return flag;
}

int test_sse()
{
    printf("test.cpp -- test_sse\n");
    
    /* Test SSE implementation within MSE by evolving a few stars and comparing the output when evolving standalone in SSE. */

//Tev(Myr)    type      Mo        Mt      log10(L)  log10(R) log10(Teff)  Mc        Menv              epoch      spin     
//12000.0000    0.0000    0.1000    0.1000   -3.0181   -0.8710    3.4442    0.0000  5.0000E-02      0.0000  3.7150E+00
//12000.0000    3.0000    0.9995    0.9985    0.8274    0.5717    3.6843    0.1749  8.1269E-01    -21.9511  1.4691E+01
//12000.0000   13.0000    9.8891    1.3693   -9.7639   -4.8539    3.7492    1.3693  1.0000E-10     27.4622  2.0000E+08
//12000.0000   14.0000   15.3521   15.2432  -10.0000   -4.1896    3.3580   15.2432  1.0000E-10      3.9544  2.0000E+08

    int flag = 0;
    
    
    double z= 0.02;
    int N=4;
    double masses[N] = {0.1,1.0,10.0,100.0};
    
    int ref_kw_final[N] = {0,3,13,14};
    double ref_m_init_final[N] = {0.1,0.9995,9.8891,15.3521};
    double ref_m_final[N] = {0.1,0.9985,1.3693,15.2432};
    double ref_log10_L_final_LSun[N] = {-3.0181,0.8274,-9.7639,-10.0000};
    double ref_log10_R_final_RSun[N] = {-0.8710,0.5717,-4.8539,-4.1896};
    double ref_m_core_final[N] = {0.0,0.1749,1.3693,15.2432};
    double ref_m_env_final[N] = {5.0e-2,8.1269e-01,1.0e-10,1.0e-10};
    double ref_epoch_final[N] = {0.0,-21.9511,27.4622,3.9544};
    double ref_ospin_final[N] = {3.7150,1.4691e1,2.0e8,2.0e8};
            
    int kw_final;
    double m_init_final,m_final,R_final,ospin_final,L_final,m_core_final,m_env_final,epoch_final;
    
    for (int i=0; i<N; i++)
    {
        flag = test_sse_specific_model(masses[i],z,&kw_final,&m_init_final,&m_final,&R_final,&ospin_final,&L_final,&m_core_final,&m_env_final,&epoch_final);
        if (1==0)
        {
            printf("==========================\n");
            printf("kw_final %d %d\n",kw_final,ref_kw_final[i]);
            printf("m_init_final %g %g\n",m_init_final,ref_m_init_final[i]);
            printf("m_final %g %g\n",m_final,ref_m_final[i]);
            printf("log10_L_final_LSun %g %g\n",log10(L_final/CONST_L_SUN),ref_log10_L_final_LSun[i]);
            printf("log10_R_final_RSun %g %g\n",log10(R_final/CONST_R_SUN),ref_log10_R_final_RSun[i]);
            printf("m_core_final %g %g\n",m_core_final,ref_m_core_final[i]);
            printf("m_env_final %g %g\n",m_env_final,ref_m_env_final[i]);
            printf("epoch_final/Myr %g %g\n",epoch_final*1e-6,ref_epoch_final[i]);
            printf("ospin_final %g %g\n",ospin_final,ref_ospin_final[i]);
        }
        if (kw_final != ref_kw_final[i])
        {
            printf("test.cpp -- error in test_sse(); non-matching stellar types %d %d\n",kw_final,ref_kw_final[i]);
            flag = 1;
        }
        if ( !equal_number(m_init_final,ref_m_init_final[i],1.0e-2) )
        {
            printf("test.cpp -- error in test_sse(); non-matching m_init_final %g %g\n",m_init_final,ref_m_init_final[i]);
            flag = 1;
        }
        if ( !equal_number(m_final,ref_m_final[i],1.0e-2) )
        {
            printf("test.cpp -- error in test_sse(); non-matching m_final %g %g\n",m_final,ref_m_final[i]);
            flag = 1;
        }
        if ( !equal_number(log10(L_final/CONST_L_SUN),ref_log10_L_final_LSun[i],1.0e-2) )
        {
            printf("test.cpp -- error in test_sse(); non-matching log10_L_final_LSun %g %g\n",log10(L_final/CONST_L_SUN),ref_log10_L_final_LSun[i]);
            flag = 1;
        }
        if ( !equal_number(log10(R_final/CONST_R_SUN),ref_log10_R_final_RSun[i],1.0e-2) )
        {
            printf("test.cpp -- error in test_sse(); non-matching log10_R_final_RSun %g %g\n",log10(R_final/CONST_R_SUN),ref_log10_R_final_RSun[i]);
            flag = 1;
        }
        if ( !equal_number(m_core_final,ref_m_core_final[i],1.0e-2) and i!=0)
        {
            printf("test.cpp -- error in test_sse(); non-matching m_core_final %g %g\n",m_core_final,ref_m_core_final[i]);
            flag = 1;
        }
        if ( !equal_number(m_env_final,ref_m_env_final[i],1.0e-2) )
        {
            printf("test.cpp -- error in test_sse(); non-matching m_env_final %g %g\n",m_env_final,ref_m_env_final[i]);
            flag = 1;
        }
        if ( !equal_number(epoch_final*1e-6,ref_epoch_final[i],1.0e-2) and i!=0)
        {
            printf("test.cpp -- error in test_sse(); non-matching epoch_final/Myr %g %g\n",epoch_final*1e-6,ref_epoch_final[i]);
            flag = 1;
        }        
        if ( !equal_number(ospin_final,ref_ospin_final[i],1.0e-2) )
        {
            printf("test.cpp -- error in test_sse(); non-matching ospin_final %g %g\n",ospin_final,ref_ospin_final[i]);
            flag = 1;
        }
    }
    
    if (flag == 0)
    {
        printf("test.cpp -- test_sse -- passed\n");
    }

    return flag;
}

int test_sse_specific_model(double m, double z, int *kw_final, double *m_init_final, double *m_final, double *R_final, double *ospin_final, double *L_final, double *m_core_final, double *m_env_final, double *epoch_final)
{
    int flag = 0;
    
    ParticlesMap particlesMap;
    Particle *star = new Particle(0, false);
    particlesMap[0] = star;
    
    star->evolve_as_star = true;
    star->mass = m;
    star->sse_initial_mass = m;
    star->metallicity = z;
    star->stellar_type = 1;
    
    initialize_stars(&particlesMap);
    double t_old,t;
    t_old = t = 0.0;
    double dt = 1.0;
    
    double t_max = 1.2e10;
    int integration_flag=0;
    bool unbound_orbits;
    bool apply_SNe_effects;
    while (t < t_max)
    {
        t+=dt;

        evolve_stars(&particlesMap,t_old,t,&dt,false,&apply_SNe_effects);
        if (apply_SNe_effects == true)
        {
            handle_SNe_in_system(&particlesMap,&unbound_orbits,&integration_flag);
        }
        else
        {
            update_stellar_evolution_quantities_directly(&particlesMap,t-t_old);
        }

        t_old = t;

        if (t>t_max)
        {
            break;
        }
    }
    *kw_final = star->stellar_type;
    *m_final = star->mass;
    *m_init_final = star->sse_initial_mass;
    *R_final = star->radius;
    *ospin_final = norm3(star->spin_vec);
    *L_final = star->luminosity;
    *m_core_final = star->core_mass;
    *m_env_final = star->convective_envelope_mass;
    *epoch_final = star->epoch;
    
    return flag;
}

int test_kick_velocity(int kick_distribution, double m, int *kw, double *v_norm)
{
    ParticlesMap particlesMap;
    Particle *star = new Particle(0, false);
    particlesMap[0] = star;
    
    star->evolve_as_star = true;
    star->mass = m;
    star->sse_initial_mass = m;
    star->stellar_type = 1;
    star->kick_distribution = kick_distribution;
    
    initialize_stars(&particlesMap);
    double t_old,t;
    t_old = t = 0.0;
    double dt = 1.0;
    
    double vx,vy,vz;
    double t_max = 1.37e10;
    
    bool apply_SNe_effects = false;
    while (apply_SNe_effects == false)
    {
        t+=dt;
        //printf("t %g t_old %g dt %g\n",t,t_old,dt);
        evolve_stars(&particlesMap,t_old,t,&dt,false,&apply_SNe_effects);
        if (apply_SNe_effects == true)
        {
            sample_kick_velocity(star,&vx,&vy,&vz);
        }
        else
        {
            update_stellar_evolution_quantities_directly(&particlesMap,t-t_old);
        }

        t_old = t;
        if (dt<0)
        {
            printf("ERROR %g %\n",dt);
        }
        if (t>t_max)
        {
            break;
        }
        //printf("dt_new %g\n",dt_new);
    }
    *kw = star->stellar_type;
    double v_vec[3] = {vx,vy,vz};
    *v_norm = norm3(v_vec);
    
    if (apply_SNe_effects == false)
    {
        *v_norm = -1;
    }
    
    //printf("done m %g v %g\n",m,v_norm);
    return 0;
}
    
    
/*********************
 * Binary evolution *
**********************/
    
int test_binary_evolution()
{
    printf("test.cpp -- test_binary_evolution\n");
    
    int flag=0;
    flag += test_compute_Kelvin_Helmholtz_timescale();
    flag += test_compute_Eddington_accretion_rate();
    flag += test_handle_instantaneous_and_adiabatic_mass_changes_in_orbit();
    flag += test_wind_accretion();
    
    if (flag == 0)
    {
        printf("test.cpp -- test_binary_evolution -- passed\n");
    }
    
    exit(0);
    return flag;
}

int test_compute_Kelvin_Helmholtz_timescale()
{
    printf("test.cpp -- test_compute_Kelvin_Helmholtz_timescale\n");
    
    int kw = 1;
    double mass = 1.0;
    double core_mass = 0.0;
    double radius = 10.0 * CONST_R_SUN;
    double luminosity = 1.0 * CONST_L_SUN;
    double t_KH = compute_Kelvin_Helmholtz_timescale(kw,mass,core_mass,radius,luminosity);

    double num = 1.57e6;
    double tol = 1.0e-2;
    int flag = 0;
    if ( !equal_number(t_KH,num,tol) )
    {
        printf("test.cpp -- error in test_compute_Kelvin_Helmholtz_timescale!\n");
        flag = 1;
    }

    kw = 2;
    mass = 1.0;
    core_mass = 0.2;
    radius = 10.0 * CONST_R_SUN;
    luminosity = 1.0 * CONST_L_SUN;
    t_KH = compute_Kelvin_Helmholtz_timescale(kw,mass,core_mass,radius,luminosity);

    num = 1.25e6;
    tol = 1.0e-2;
    if ( !equal_number(t_KH,num,tol) )
    {
        printf("test.cpp -- error in test_compute_Kelvin_Helmholtz_timescale!\n");
        flag += 1;
    }

    if (flag == 0)
    {
        printf("test.cpp -- test_compute_Kelvin_Helmholtz_timescale -- passed\n");
    }

    return flag;
}

int test_compute_Eddington_accretion_rate()
{
    printf("test.cpp -- test_compute_Eddington_accretion_rate\n");
    
    double radius = 1.0 * CONST_R_SUN;
    double hydrogen_mass_fraction = 0.7;
    double m_dot = compute_Eddington_accretion_rate(radius,hydrogen_mass_fraction);
    
    int flag = 0;
    double num = 1.22e-2;
    double tol = 1e-2;
    if ( !equal_number(m_dot,num,tol) )
    {
        printf("test.cpp -- error in test_compute_Eddington_accretion_rate! %g %g\n",m_dot,num);
        flag = 1;
    }
    
    if (flag == 0)
    {
        printf("test.cpp -- test_compute_Eddington_accretion_rate -- passed\n");
    }

    return flag;
}

int test_handle_instantaneous_and_adiabatic_mass_changes_in_orbit()
{
    printf("test.cpp -- test_handle_instantaneous_and_adiabatic_mass_changes_in_orbit\n");
    
    ParticlesMap particlesMap;
    int N_bodies = 4;
    double m1 = 20.0;
    double m2 = 15.0;
    double m3 = 25.0;
    double m4 = 14.0;
    double a1 = 1.0;
    double a2 = 40.0;
    double a3 = 2000.0;
    double masses[4] = {m1,m2,m3,m4};
    int stellar_types[4] = {1,1,1,1};
    double smas[3] = {a1,a2,a3};
    double es[3] = {0.01,0.01,0.01};
    double TAs[3] = {0.01,0.01,0.01};
    double INCLs[3] = {0.01,0.01,0.01};
    double APs[3] = {0.01,0.01,0.01};
    double LANs[3] = {0.01,0.01,0.01};
    
    create_nested_system(particlesMap,N_bodies,masses,stellar_types,smas,es,TAs,INCLs,APs,LANs);// = create_nested_system();
//    printf("post s %d b %d %g r %g\n",particlesMap2.size(),particlesMap2[0]->is_binary,particlesMap2[0]->mass,particlesMap2[0]->radius);
    initialize_code(&particlesMap);
    
    Particle *star1 = particlesMap[0];
    Particle *star2 = particlesMap[1];
    double Delta_m1 = -5.0;
    double Delta_m2 = -2.0;
    double mass_loss_timescale = 1.0e2;
    int integration_flag = 0;
    
    handle_instantaneous_and_adiabatic_mass_changes_in_orbit(&particlesMap,star1,star2,Delta_m1,Delta_m2,mass_loss_timescale,&integration_flag);
    
    set_up_derived_quantities(&particlesMap);
    double m1_new = m1 + Delta_m1;
    double m2_new = m2 + Delta_m2;
    double a1_new = a1 * (m1 + m2)/(m1_new + m2_new);
    double a2_new = a2 * (m1 + m2 + m3)/(m1_new + m2_new + m3);
    //printf("Test a1n %g a2n %g %g %g\n",particlesMap[4]->a,particlesMap[5]->a,a1_new,a2_new);

    int flag=0;
    double tol = 1e-12;
    if ( !equal_number(particlesMap[4]->a,a1_new,tol) )
    {
        printf("test.cpp -- error in test_handle_instantaneous_and_adiabatic_mass_changes_in_orbit! %g %g\n",particlesMap[4]->a,a1_new);
        flag = 1;
    }
    if ( !equal_number(particlesMap[5]->a,a2_new,tol) )
    {
        printf("test.cpp -- error in test_handle_instantaneous_and_adiabatic_mass_changes_in_orbit! %g %g\n",particlesMap[5]->a,a2_new);
        flag = 1;
    }
    
    if (flag == 0)
    {
        printf("test.cpp -- test_handle_instantaneous_and_adiabatic_mass_changes_in_orbit -- passed\n");
    }

    return flag;
}

int test_wind_accretion()
{
    printf("test.cpp -- test_wind_accretion\n");
    
    ParticlesMap particlesMap;
    int N_bodies = 2;
    double m1 = 20.0;
    double m2 = 15.0;
    double a = 10.0;
    double e = 0.9;
    
    double masses[2] = {m1,m2};
    int stellar_types[4] = {1,1};
    double smas[1] = {a};
    double es[1] = {e};
    double TAs[1] = {0.01};
    double INCLs[1] = {0.01};
    double APs[1] = {0.01};
    double LANs[1] = {0.01};
    
    create_nested_system(particlesMap,N_bodies,masses,stellar_types,smas,es,TAs,INCLs,APs,LANs);

    //initialize_code(&particlesMap);
    initialize_stars(&particlesMap);
    
    Particle *star1 = particlesMap[0];
    Particle *star2 = particlesMap[1];
    Particle *orbit = particlesMap[2];
    int integration_flag = 0;
    double t_old = 0.0;
    double t = 0.0;
    double dt_binary_evolution;
    
    /* Let *1 have a wind which is partially accreted by *2 */
    double m1_dot = -1.0e-5;
    double m2_dot = 0.0;
    star1->mass_dot_wind = m1_dot;
    star2->mass_dot_wind = m2_dot;
    handle_wind_accretion(&particlesMap, t_old, t, &dt_binary_evolution, &integration_flag);
    
    double R1 = star1->radius;
    double mdot_test = star2->mass_dot_wind_accretion;
    
    double mdot_an;
    double v_orb_p2 = CONST_G * (m1+m2)/a;
    double v_wind_p2 = 2.0 * beta_wind_accretion * CONST_G * m1 / R1; /* squared wind speed from p */

    double factor = (1.0/sqrt(1.0 - e*e)) * pow( CONST_G * m2 / v_wind_p2, 2.0) * (alpha_wind_accretion / (2.0 * a*a)) * pow(1.0 + v_orb_p2/v_wind_p2, -1.5);

    if (factor >= 1.0)
    {
        mdot_an = -m1_dot;
    }
    else
    {
        mdot_an = -factor*m1_dot;
    }
            
    int flag = 0;
    double tol = 1e-12;
    if ( !equal_number(mdot_test,mdot_an,tol) )
    {
        printf("test.cpp -- test_wind_acretion -- error: mdot_test %g mdot_an %g\n",mdot_test,mdot_an);
        flag = 1;
    }

    if (flag == 0)
    {
        printf("test.cpp -- test_wind_accretion -- passed\n");
    }

    return flag;
}

int test_collisions()
{
    printf("test.cpp -- test_collisions\n");
    int flag;
    
    //flag = test_collision_MS_MS();
    //flag = test_collision_giant_MS();
    //flag = test_collision_star_MS(10.0,1,5.0);
    //flag = test_collision_star_MS(40.0,2,5.0);
    
    int i,j;
    for (i=1; i<=14; i++)
    {
        for (j=1; j<=14; j++)
        {
            printf("i %d j %d\n",i,j);
            //if (i!=1 or j!=4)
            //if (i>5 and i<12 or j>5 and j<12)
            if (j<6)
            {
                //continue;
            }
            flag += test_collision_stars(10.0,i,13,j,0);
            random_seed = 0;
            //flag = test_collision_stars(10.0,i,13,j,1);
            
        }
    }
    
    if (flag == 0)
    {
        printf("test.cpp -- test_collisions -- passed\n");
    }
    //flag = test_collision_stars(10.0,3,8,6);
    //flag = test_collision_stars(5.01,4,2.0);
    //flag = test_collision_stars(5.01,6,3.0);
    //flag = test_collision_stars(10.0,4,15);
    //flag = test_collision_stars(22.0,4,21.98);
    //flag = test_collision_stars(22.0,14,8.0);
    //flag = test_collision_stars(50.0,14,22.0);
    
//    double m=1.2;
//    double Omega = 1.8;
//    double chi = compute_spin_parameter_from_spin_frequency(m,Omega);
//    double Omega2 = compute_spin_frequency_from_spin_parameter(m,chi);
//    printf("test Omega %g %g\n",Omega,Omega2);
    
    return flag;
}

int test_collision_stars(double m1, int kw1, double m2, int kw2, int integration_flag)
{
    printf("*******************************************************************\n");
    printf("test_collision_star_MS m1 %g kw1 %d m2 %g kw2 %d integration_flag %d \n",m1,kw1,m2,kw2,integration_flag);
    printf("*******************************************************************\n");
    
    ParticlesMap particlesMap;
    int N_bodies = 4;
    double masses[4] = {m1,m2,95.0,40.0};
    int stellar_types[4] = {kw1,kw2,1,1};
    double smas[3] = {10.0,100.0,10000.0};
    double es[3] = {0.01,0.01,0.01};
    double TAs[3] = {0.01,0.01,0.01};
    double INCLs[3] = {0.01,0.01,0.01};
    double APs[3] = {0.01,0.01,0.01};
    double LANs[3] = {0.01,0.01,0.01};
    random_seed=6;
    
    create_nested_system(particlesMap,N_bodies,masses,stellar_types,smas,es,TAs,INCLs,APs,LANs);// = create_nested_system();
//    printf("post s %d b %d %g r %g\n",particlesMap2.size(),particlesMap2[0]->is_binary,particlesMap2[0]->mass,particlesMap2[0]->radius);
    initialize_code(&particlesMap);
    
    //int integration_flag_init = 0;
    printf("test_collision_stars -- pre merge\n");
    print_system(&particlesMap,integration_flag);

    double start_time = 0.0;
    double end_time = 1.0e1;
    double output_time,hamiltonian;
    int state,CVODE_flag,CVODE_error_code;
    //int integration_flag = 0;

    int kw = particlesMap[0]->stellar_type;
    double t_old=0.0;
    double t = 0.0;
    double dt;
    double dt_stev;
    bool apply_SNe_effects;
    bool unbound_orbits;


    int N_binaries,N_root_finding,N_ODE_equations;

    //printf("post merge\n");
    //set_up_derived_quantities(&particlesMap);
    //print_system(&particlesMap);
    if (integration_flag==0)
    {
        particlesMap[4]->merged = true;
    }
    else
    {
        particlesMap[0]->Stopping_Condition_Partner = 1;
        particlesMap[1]->Stopping_Condition_Partner = 0;
    }

    handle_collisions(&particlesMap,t,&integration_flag);
    //printf("pre s %d\n",particlesMap.size());
    //collision_product(&particlesMap, 4, 0, 1, &integration_flag);
    //printf("post s %d b %d %g r %g\n",particlesMap.size(),particlesMap[2]->is_binary,particlesMap[2]->mass,particlesMap[2]->radius);

    //printf("post1\n");
    //print_system(&particlesMap);



    printf("post merge integration_flag %d\n",integration_flag);
    
    print_system(&particlesMap,integration_flag);
    
    //#ifdef IGNORE
    printf("test_collision_stars -- pre evolve\n");
    evolve(&particlesMap,start_time,end_time,&output_time,&hamiltonian,&state,&CVODE_flag,&CVODE_error_code,&integration_flag);
    //initialize_code(&particlesMap);
    printf("test_collision_stars -- done evolve output_time %g\n",output_time);
    //evolve(&particlesMap,start_time,end_time,&output_time,&hamiltonian,&state,&CVODE_flag,&CVODE_error_code,&integration_flag);
    //printf("post evolve\n");
    //#endif
    
    
    print_system(&particlesMap,integration_flag);
    
    clear_particles(&particlesMap);
   
    int flag;
    
    
    return 0;
}





}
