/* MSE */

#include "evolve.h"
#include "test.h"

#include <iostream>
#include <fstream>

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
    flag += test_BH_spin_conversion();
    
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

    double m1_dot = -1.0e-10;
    double m2_dot = 1.0e-12;
    double a_dot = 1.0e-5;
    double e_dot = 1.0e-8;
    double h_dot_div_h = compute_h_dot_div_h(m1,m1_dot,m2,m2_dot,a,a_dot,e,e_dot);
    double h_dot = h*h_dot_div_h;
    double a_dot_new = a*compute_a_dot_div_a(m1,m1_dot,m2,m2_dot,h,h_dot,e,e_dot);

    if ( !equal_number(a_dot,a_dot_new,tol) )
    {
        printf("test.cpp -- error in test_a_h_conversion!\n");
        flag = 1;
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
    
    return flag;
}

int test_BH_spin_conversion()
{
    printf("test.cpp -- test_BH_spin_conversion\n");
    double m = 15.0;
    double Omega = 2.0;
    
    double chi = compute_spin_parameter_from_spin_frequency(m,Omega);
    double Omega2 = compute_spin_frequency_from_spin_parameter(m,chi);
    
    int flag = 0;
    double tol = 1.0e-10;

    if ( !equal_number(Omega,Omega2,tol) )
    {
        printf("test.cpp -- error in test_BH_spin_conversion! %g %g\n",Omega,Omega2);
        flag = 1;
    }
    
//    double m_dot = 1.0e-5;
//    double J = chi * CONST_G * m * m / CONST_C_LIGHT;
//    double J_dot = 1.0e-5 * J;
//    double Omega_dot = compute_spin_frequency_dot_BHs(m,Omega,J_dot,m_dot);
    
    return flag;
}


/**********
 * N-body *
***********/

int test_nbody(int mode)
{
    printf("test.cpp -- test_nbody\n");
    
    MSTAR_verbose = 0;
    
    int flag = 0;
    flag += test_nbody_two_body_stopping_conditions();
    flag += test_nbody_two_body_kick();
    flag += test_nbody_pythagorean();
    flag += test_nbody_inspiral(mode);
    flag += test_nbody_spin_orbit(mode);
    //flag += test_nbody_custom(mode);
    //flag += test_nbody_custom2(mode);
    
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
    for (int stopping_condition_mode = 0; stopping_condition_mode < 4; stopping_condition_mode++)
    {
        //printf("stopping_condition_mode %d\n",stopping_condition_mode);
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

        MSTAR_include_PN_acc_10 = false;
        MSTAR_include_PN_acc_20 = false;
        MSTAR_include_PN_acc_25 = false;
        MSTAR_include_PN_acc_30 = false;
        MSTAR_include_PN_acc_35 = false;

        MSTAR_include_PN_acc_SO = false;
        MSTAR_include_PN_acc_SS = false;
        MSTAR_include_PN_acc_Q = false;

        MSTAR_include_PN_spin_SO = false;
        MSTAR_include_PN_spin_SS = false;
        MSTAR_include_PN_spin_Q = false;
        
        double gbs_tolerance = 1.0e-12;
        double stopping_condition_tolerance = 1.0e-12;
        struct RegularizedRegion *R;
        if (stopping_condition_mode < 2)
        {
            
            R = generate_binary_ICs(m1,m2,R1,R2,a,e,stopping_condition_mode,stopping_condition_mode,gbs_tolerance,stopping_condition_tolerance);
        }
        else
        {
            if (stopping_condition_mode==2)
            {
                int stopping_condition_mode1 = 0;
                int stopping_condition_mode2 = 1;
                R = generate_binary_ICs(m1,m2,R1,R2,a,e,stopping_condition_mode1,stopping_condition_mode2,gbs_tolerance,stopping_condition_tolerance);
            }
            if (stopping_condition_mode==3)
            {
                int stopping_condition_mode1 = 1;
                int stopping_condition_mode2 = 0;
                R = generate_binary_ICs(m1,m2,R1,R2,a,e,stopping_condition_mode1,stopping_condition_mode2,gbs_tolerance,stopping_condition_tolerance);
            }
        }
        
        double E_init = compute_nbody_total_energy(R);

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
                
                mstar_print_state(R);
                
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
        
        double r_RLOF1 = R1/roche_radius_pericenter_eggleton(1.0,m1/m2);
        double r_RLOF2 = R2/roche_radius_pericenter_eggleton(1.0,m2/m1);
        
        double t_an = test_nbody_compute_time_of_collision(GCONST,m1+m2,R1+R2,a,e);
        
        test_nbody_compute_elements(GCONST,M,r,v,&a,&e);
        double r_norm = norm3(r);
        //printf("Final a %.15f e %.15f sep %.15f R1 %g R2 %g R1+R2 %g r_RLOF1 %g r_RLOF2 %g\n",a,e,r_norm,R1,R2,R1+R2,r_RLOF1,r_RLOF2);

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

        if (stopping_condition_mode == 0) // collision
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
        else if (stopping_condition_mode == 1) // RLOF
        {
            if (!equal_number(r_norm,CV_max(r_RLOF1,r_RLOF2),tol))
            {
                printf("test.cpp -- error in test_nbody_two_body_stopping_conditions -- r_norm %g max(r_RLOF1,r_RLOF2) %g\n",r_norm,CV_max(r_RLOF1,r_RLOF2));
                flag = 1;
            }
        }
        else if (stopping_condition_mode == 2) // col + RLOF
        {
            if (!equal_number(r_norm,CV_max(R1,r_RLOF2),tol))
            {
                printf("test.cpp -- error in test_nbody_two_body_stopping_conditions -- r_norm %g max(R1,r_RLOF2) %g\n",r_norm,CV_max(R1,r_RLOF2));
                flag = 1;
            }
        }
        else if (stopping_condition_mode == 3) // RLOF + col
        {
            if (!equal_number(r_norm,CV_max(r_RLOF1,R2),tol))
            {
                printf("test.cpp -- error in test_nbody_two_body_stopping_conditions -- r_norm %g max(R1,r_RLOF2) %g\n",r_norm,CV_max(r_RLOF1,R2));
                flag = 1;
            }
        }
        
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

struct RegularizedRegion *generate_binary_ICs(double m1, double m2, double R1, double R2, double a, double e, int stopping_condition_mode1, int stopping_condition_mode2, double gbs_tolerance, double stopping_condition_tolerance)
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
        R->Spin_S[3 * 0 + i] = 0.0;

        R->Pos[3 * 1 + i] = R2_vec[i];
        R->Vel[3 * 1 + i] = V2_vec[i];
        R->Spin_S[3 * 1 + i] = 0.0;
    }

    R->Mass[0] = m1;
    R->Mass[1] = m2;
    R->Radius[0] = R1;
    R->Radius[1] = R2;

    R->Stopping_Condition_Mode[0] = stopping_condition_mode1;
    R->Stopping_Condition_Mode[1] = stopping_condition_mode2;
    R->gbs_tolerance = gbs_tolerance;
    R->stopping_condition_tolerance = stopping_condition_tolerance;
    
   
    return R;
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
    
    MSTAR_include_PN_acc_10 = false;
    MSTAR_include_PN_acc_20 = false;
    MSTAR_include_PN_acc_25 = false;
    MSTAR_include_PN_acc_30 = false;
    MSTAR_include_PN_acc_35 = false;

    MSTAR_include_PN_acc_SO = false;
    MSTAR_include_PN_acc_SS = false;
    MSTAR_include_PN_acc_Q = false;

    MSTAR_include_PN_spin_SO = false;
    MSTAR_include_PN_spin_SS = false;
    MSTAR_include_PN_spin_Q = false;
    
    double gbs_tolerance = 1.0e-8;
    double stopping_condition_tolerance = 1.0e-6;
    struct RegularizedRegion *R = generate_binary_ICs(m1,m2,R1,R2,a,e,stopping_condition_mode,stopping_condition_mode,gbs_tolerance,stopping_condition_tolerance);

    
    /* Give kick to *2 */
    int i=2;
    int j=0;
    R->Vel[3 * j + i] += 1000.0*CONST_KM_PER_S;

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
            mstar_print_state(R);
            
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

    return flag;
}

struct RegularizedRegion *generate_pythagorean_ICs(double m1, double m2, double m3, int stopping_condition_mode, double gbs_tolerance, double stopping_condition_tolerance)
{
     
    struct RegularizedRegion *R;

    int i=0;
    int j;
   
    int N = 3;

    allocate_armst_structs(&R, N); // Initialize the data structure
    
    double R1_vec[3] = {1.0,3.0,0.0};
    double R2_vec[3] = {-2.0,-1.0,0.0};
    double R3_vec[3] = {1.0,-1.0,0.0};
    double V1_vec[3] = {0.0,0.0,0.0};
    double V2_vec[3] = {0.0,0.0,0.0};
    double V3_vec[3] = {0.0,0.0,0.0};
    
    for (i=0; i<3; i++)
    {
        R->Pos[3 * 0 + i] = R1_vec[i];
        R->Vel[3 * 0 + i] = V1_vec[i];
        R->Spin_S[3 * 0 + i] = 0.0;
        
        R->Pos[3 * 1 + i] = R2_vec[i];
        R->Vel[3 * 1 + i] = V2_vec[i];
        R->Spin_S[3 * 1 + i] = 0.0;

        R->Pos[3 * 2 + i] = R3_vec[i];
        R->Vel[3 * 2 + i] = V3_vec[i];
        R->Spin_S[3 * 2 + i] = 0.0;
    }

    R->Mass[0] = m1;
    R->Mass[1] = m2;
    R->Mass[2] = m3;

    R->Stopping_Condition_Mode[0] = stopping_condition_mode;
    R->Stopping_Condition_Mode[1] = stopping_condition_mode;
    R->gbs_tolerance = gbs_tolerance;
    R->stopping_condition_tolerance = stopping_condition_tolerance;
    
    return R;
}

int test_nbody_pythagorean()
{
    printf("test.cpp -- test_nbody_pythagorean\n");

    int flag = 0;    
    int stopping_condition_mode = 0;
    
    initialize_mpi_or_serial(); 
    int N_bodies = 2;

    double m1 = 3.0;
    double m2 = 4.0;
    double m3 = 5.0;

    MSTAR_include_PN_acc_10 = false;
    MSTAR_include_PN_acc_20 = false;
    MSTAR_include_PN_acc_25 = false;
    MSTAR_include_PN_acc_30 = false;
    MSTAR_include_PN_acc_35 = false;

    MSTAR_include_PN_acc_SO = false;
    MSTAR_include_PN_acc_SS = false;
    MSTAR_include_PN_acc_Q = false;

    MSTAR_include_PN_spin_SO = false;
    MSTAR_include_PN_spin_SS = false;
    MSTAR_include_PN_spin_Q = false;
    
    double gbs_tolerance = 1.0e-12;
    double stopping_condition_tolerance = 1.0e-6;
    struct RegularizedRegion *R = generate_pythagorean_ICs(m1,m2,m3,stopping_condition_mode,gbs_tolerance,stopping_condition_tolerance);

    double E_init = compute_nbody_total_energy(R);

    double tend = 20.0;
    double reached_dt;
    int stopping_condition_occurred;

    run_integrator(R, tend, &reached_dt, &stopping_condition_occurred);
    
    double E_fin = compute_nbody_total_energy(R);
    
    double tol = 1.0e-9;

    if (!equal_number(E_init,E_fin,tol))
    {
        printf("test.cpp -- error in test_nbody_pythagorean -- E_init %g E_fin %g err %g\n",E_init,E_fin,fabs( (E_init-E_fin)/E_init));
        flag = 1;
    }

    free_data(R);

    return flag;
}

int test_nbody_inspiral(int mode)
{
    printf("test.cpp -- test_nbody_inspiral\n");

    int flag = 0;    
    int stopping_condition_mode = 0;
    
    initialize_mpi_or_serial(); 
    int N_bodies = 2;

    //SPEEDOFLIGHT = 63239.72638679138;
    SPEEDOFLIGHT = CONST_C_LIGHT;

    MSTAR_include_PN_acc_10 = false;
    MSTAR_include_PN_acc_20 = false;
    MSTAR_include_PN_acc_25 = true;
    MSTAR_include_PN_acc_30 = false;
    MSTAR_include_PN_acc_35 = false;

    MSTAR_include_PN_acc_SO = false;
    MSTAR_include_PN_acc_SS = false;
    MSTAR_include_PN_acc_Q = false;

    MSTAR_include_PN_spin_SO = false;
    MSTAR_include_PN_spin_SS = false;
    MSTAR_include_PN_spin_Q = false;

    double m1 = 10.0;
    double m2 = 10.0;
    double a = 2.0e-5;
    double e = 1.0e-12;
    double R1 = 10.0*CONST_G*m1/(SPEEDOFLIGHT*SPEEDOFLIGHT);
    double R2 = 10.0*CONST_G*m2/(SPEEDOFLIGHT*SPEEDOFLIGHT);
    
    double gbs_tolerance = 1.0e-12;
    double stopping_condition_tolerance = 1.0e-12;
    struct RegularizedRegion *R = generate_binary_ICs(m1,m2,R1,R2,a,e,stopping_condition_mode,stopping_condition_mode,gbs_tolerance,stopping_condition_tolerance);
    
    double E_init = compute_nbody_total_energy(R);
    double P_orb = compute_orbital_period_from_semimajor_axis(m1+m2,a);
    
    double M = m1+m2;
    double mu = m1*m2/M;
    double L = mu * sqrt(GCONST* M * a);
    
    double beta = (64.0/5.0)*GCONST*GCONST*GCONST*m1*m2*(m1+m2)*pow(SPEEDOFLIGHT,-5.0);
    double t_GW = a*a*a*a/(4.0*beta); // circular orbits only
    
    double tend = 2*t_GW;
    //printf("t_GW/P %g\n",t_GW/P_orb);
    //printf("tend/P %g\n",tend/P_orb);
    
    double reached_dt;
    int stopping_condition_occurred;

    double t=0;
    double dt = tend/100.0;
    
    std::ofstream myfile;
    if (mode == 2)
    {
        myfile.open ("test_nbody_inspiral.txt");
    }
    
    double at,et;
    double r[3],v[3];
    while (t<tend)
    {
        run_integrator(R, dt, &reached_dt, &stopping_condition_occurred);
        t += reached_dt;
       
        for (int j=0; j<3; j++)
        {
            r[j] = R->Pos[3 * 0 + j] - R->Pos[3 * 1 + j];
            v[j] = R->Vel[3 * 0 + j] - R->Vel[3 * 1 + j];
        }
        
        test_nbody_compute_elements(GCONST,M,r,v,&at,&et);

        //printf("t/P_orb %g a %g e %g sep %g\n",t/P_orb,at,et,norm3(r));
        
        if (mode == 2)
        {
            myfile << t/P_orb << "," << norm3(r) << "\n";
        }
        
        if (stopping_condition_occurred==1)
        {
            break;
        }
    }

    if (mode == 2)
    {
        myfile.close();    
    }

    double tol = 1.0e-8;
//    if (!equal_number(E_init,E_fin,tol))
//    {
//        printf("test.cpp -- error in test_nbody_spin_orbit -- E_init %g E_fin %g err %g\n",E_init,E_fin,fabs( (E_init-E_fin)/E_init));
//        flag = 1;
//    }

    free_data(R);

    return flag;
}


int test_nbody_spin_orbit(int mode)
{
    printf("test.cpp -- test_nbody_spin_orbit\n");

    int flag = 0;    
    int stopping_condition_mode = 0;
    
    initialize_mpi_or_serial(); 
    int N_bodies = 2;

    double m1 = 10.0;
    double m2 = 10.0;
    double R1 = 1.0e-10;
    double R2 = 1.0e-10;
    double a = 1.0e-2;
    double e = 0.99;
    
    //SPEEDOFLIGHT = 63239.72638679138;
    SPEEDOFLIGHT = CONST_C_LIGHT;
    
    MSTAR_include_PN_acc_10 = false;
    MSTAR_include_PN_acc_20 = false;
    MSTAR_include_PN_acc_25 = false;
    MSTAR_include_PN_acc_30 = false;
    MSTAR_include_PN_acc_35 = false;

    MSTAR_include_PN_acc_SO = false;
    MSTAR_include_PN_acc_SS = false;
    MSTAR_include_PN_acc_Q = false;

    MSTAR_include_PN_spin_SO = true;
    MSTAR_include_PN_spin_SS = false;
    MSTAR_include_PN_spin_Q = false;
    
    double gbs_tolerance = 1.0e-12;
    double stopping_condition_tolerance = 1.0e-12;
    struct RegularizedRegion *R = generate_binary_ICs(m1,m2,R1,R2,a,e,stopping_condition_mode,stopping_condition_mode,gbs_tolerance,stopping_condition_tolerance);
    
    double eps = 1.0e-3;
    double spin1[3] = {eps,0.0,eps};
    double spin2[3] = {eps,0.0,eps};
    for (int i=0; i<3; i++)
    {
        R->Spin_S[3 * 0 + i] = spin1[i];
        R->Spin_S[3 * 1 + i] = spin2[i];
    }
    
    double E_init = compute_nbody_total_energy(R);
    double P_orb = compute_orbital_period_from_semimajor_axis(m1+m2,a);
    
    double M = m1+m2;
    double mu = m1*m2/M;
    double L = mu * sqrt(GCONST* M * a);
    double P_SO = 2.0*M_PI * SPEEDOFLIGHT*SPEEDOFLIGHT * a*a*a * (1.0 - e*e) * (1.0/(2.0*GCONST)) * (1.0/(1.0 + 0.75*m2/m1)) * (1.0/L);
    double tend = P_SO;
    
    double reached_dt;
    int stopping_condition_occurred;

    double t=0;
    double dt = tend/100.0;
    double precession_angle_old;
    
    std::ofstream myfile;
    if (mode == 2)
    {
        myfile.open ("test_nbody_spin_orbit.txt");
    }
    
    double at,et;
    double r[3],v[3];
    while (t<tend)
    {
        run_integrator(R, dt, &reached_dt, &stopping_condition_occurred);
        t += reached_dt;

        double precession_angle = atan2(R->Spin_S[3 * 0 + 1],R->Spin_S[3 * 0 + 0]);
        if (precession_angle < 0.0)
        {
            precession_angle += 2.0*M_PI;
        }

        for (int j=0; j<3; j++)
        {
            r[j] = R->Pos[3 * 0 + j] - R->Pos[3 * 1 + j];
            v[j] = R->Vel[3 * 0 + j] - R->Vel[3 * 1 + j];
        }
        
        test_nbody_compute_elements(GCONST,M,r,v,&at,&et);

        //printf("t/P_orb %g precession angle/deg %g delta/deg %g a %g e %g sep %g\n",t/P_orb,precession_angle*180.0/M_PI,(precession_angle-precession_angle_old)*180.0/M_PI,at,et,norm3(r));
        //printf("spin 1 %g %g %g\n",R->Spin_S[3 * 0 + 0],R->Spin_S[3 * 0 + 1],R->Spin_S[3 * 0 + 2]);
        //printf("spin 2 %g %g %g\n",R->Spin_S[3 * 1 + 0],R->Spin_S[3 * 1 + 1],R->Spin_S[3 * 1 + 2]);
        
        if (mode == 2)
        {
            myfile << t/P_orb << "," << precession_angle*(180.0/M_PI) << "\n";
        }

        precession_angle_old = precession_angle;
        
        if (stopping_condition_occurred==1)
        {
            break;
        }

    }
    
    if (mode == 2)
    {
        myfile.close();    
    }
   
    double tol = 1.0e-8;
//    printf("test.cpp -- error in test_nbody_spin_orbit -- E_init %g E_fin %g err %g\n",E_init,E_fin,fabs( (E_init-E_fin)/E_init));
//    if (!equal_number(E_init,E_fin,tol))
//    {
//        printf("test.cpp -- error in test_nbody_spin_orbit -- E_init %g E_fin %g err %g\n",E_init,E_fin,fabs( (E_init-E_fin)/E_init));
//        flag = 1;
//    }

    free_data(R);
    
    return flag;
}


int test_nbody_custom(int mode)
{
    printf("test.cpp -- test_nbody_custom\n");

    int flag = 0;    
    int stopping_condition_mode = 0;
    
    initialize_mpi_or_serial(); 
   
    MSTAR_verbose = 1;
    SPEEDOFLIGHT = CONST_C_LIGHT;

    bool include = false;
    MSTAR_include_PN_acc_10 = include;
    MSTAR_include_PN_acc_20 = include;
    MSTAR_include_PN_acc_25 = include;
    MSTAR_include_PN_acc_30 = include;
    MSTAR_include_PN_acc_35 = include;

    MSTAR_include_PN_acc_SO = include;
    MSTAR_include_PN_acc_SS = include;
    MSTAR_include_PN_acc_Q = include;

    MSTAR_include_PN_spin_SO = include;
    MSTAR_include_PN_spin_SS = include;
    MSTAR_include_PN_spin_Q = include;
    
    
    double gbs_tolerance = 1.0e-12;
    double stopping_condition_tolerance = 1.0e-12;

     
    struct RegularizedRegion *R;

    int i=0;
    int j;
   
    int N = 3;

    allocate_armst_structs(&R, N); // Initialize the data structure
    
    // Issue #16
    R->Mass[0] = 1.33919;
    R->Mass[1] = 1.03431;
    R->Mass[2] = 0.18581;

    double R1_vec[3] = {334.275, 78.5216 ,270.365};
    double R2_vec[3] = {-2352.3, -528.001, -1809.79};
    double R3_vec[3] = {-2352.27, -528.187 ,-1809.92};
    double V1_vec[3] = {103.938 ,58.4821 ,-35.9006};
    double V2_vec[3] = {-1.26751, -1.49463, -0.332985};
    double V3_vec[3] = { 8.35301, 7.92899, 0.320513};
    
    #ifdef IGNORE
    R->Mass[0] = 4.24276;
    R->Mass[1] = 5.05712;
    R->Mass[2] = 16.5278;

    double R1_vec[3] = {-471.064, 832.573 ,6284.28};
    double R2_vec[3] = {-478.381, 817.815 ,6261.61};
    double R3_vec[3] = { 249.561, -255.505 ,-2196.88};
    double V1_vec[3] = {0.577908, -0.0394511, 1.41379};
    double V2_vec[3] = { -0.779966, -0.152874 ,-1.60409};
    double V3_vec[3] = { 14.1263 ,-1.08172, 7.26681};
    #endif    
    
    for (i=0; i<3; i++)
    {
        R->Pos[3 * 0 + i] = R1_vec[i];
        R->Vel[3 * 0 + i] = V1_vec[i];

        R->Pos[3 * 1 + i] = R2_vec[i];
        R->Vel[3 * 1 + i] = V2_vec[i];

        R->Pos[3 * 2 + i] = R3_vec[i];
        R->Vel[3 * 2 + i] = V3_vec[i];
    }


    R->Stopping_Condition_Mode[0] = stopping_condition_mode;
    R->Stopping_Condition_Mode[1] = stopping_condition_mode;
    R->Stopping_Condition_Mode[2] = stopping_condition_mode;
    R->gbs_tolerance = gbs_tolerance;
    R->stopping_condition_tolerance = stopping_condition_tolerance;

    for (int i=0; i<3; i++)
    {
        R->Spin_S[3 * 0 + i] = 0.0;
        R->Spin_S[3 * 1 + i] = 0.0;
    }
        
    double reached_dt;
    int stopping_condition_occurred;

    double t=0;
    double dt = 1756.4;
    //double dt = 20000;
    
    printf("pre\n");
    mstar_print_state(R);
    double E_init = compute_nbody_total_energy(R);
    
    printf("running\n");
    run_integrator(R, dt, &reached_dt, &stopping_condition_occurred);
    
    double E_fin = compute_nbody_total_energy(R);
    
    printf("post\n");
    mstar_print_state(R);
    
    printf("E_init %g E_fin %g E_err %g\n",E_init,E_fin,fabs((E_init-E_fin)/E_init));
    
    free_data(R);
    
    return flag;
}

int test_nbody_custom2(int mode)
{
    printf("test.cpp -- test_nbody_custom\n");

    int flag = 0;    
    int stopping_condition_mode = 0;
    
    initialize_mpi_or_serial(); 
   
    MSTAR_verbose = 1;
    SPEEDOFLIGHT = CONST_C_LIGHT;

    bool include = false;
    MSTAR_include_PN_acc_10 = include;
    MSTAR_include_PN_acc_20 = include;
    MSTAR_include_PN_acc_25 = include;
    MSTAR_include_PN_acc_30 = include;
    MSTAR_include_PN_acc_35 = include;

    MSTAR_include_PN_acc_SO = include;
    MSTAR_include_PN_acc_SS = include;
    MSTAR_include_PN_acc_Q = include;

    MSTAR_include_PN_spin_SO = include;
    MSTAR_include_PN_spin_SS = include;
    MSTAR_include_PN_spin_Q = include;
    
    
    double gbs_tolerance = 1.0e-12;
    double stopping_condition_tolerance = 1.0e-12;

     
    struct RegularizedRegion *R;

    int i=0;
    int j;
   
    int N = 4;

    allocate_armst_structs(&R, N); // Initialize the data structure
    
    // Issue #34
    R->Mass[0] = 4.72225;
    R->Mass[1] = 3.75456;
    R->Mass[2] = 1.56984;
    R->Mass[3] = 2.459;

    double R1_vec[3] = {381.896, 581.169, 572.307};
    double R2_vec[3] = {381.506, 581.167 ,572.323};
    double R3_vec[3] = {-3397.39, -5128.91, -4872.17};
    double R4_vec[3] = {-3394.78, -5128.14 ,-4872.58};
    double V1_vec[3] = {-4.98801, -14.9323 ,13.7696};
    double V2_vec[3] = {0.0712215 ,3.93403, -27.1107};
    double V3_vec[3] = {1.36039, -3.57182 ,3.06344};
    double V4_vec[3] = {-0.533719, 2.65386, -2.18125};
    
    for (i=0; i<3; i++)
    {
        R->Pos[3 * 0 + i] = R1_vec[i];
        R->Vel[3 * 0 + i] = V1_vec[i];

        R->Pos[3 * 1 + i] = R2_vec[i];
        R->Vel[3 * 1 + i] = V2_vec[i];

        R->Pos[3 * 2 + i] = R3_vec[i];
        R->Vel[3 * 2 + i] = V3_vec[i];

        R->Pos[3 * 3 + i] = R4_vec[i];
        R->Vel[3 * 3 + i] = V4_vec[i];

    }


    R->Stopping_Condition_Mode[0] = stopping_condition_mode;
    R->Stopping_Condition_Mode[1] = stopping_condition_mode;
    R->Stopping_Condition_Mode[2] = stopping_condition_mode;
    R->gbs_tolerance = gbs_tolerance;
    R->stopping_condition_tolerance = stopping_condition_tolerance;

    for (int i=0; i<3; i++)
    {
        R->Spin_S[3 * 0 + i] = 0.0;
        R->Spin_S[3 * 1 + i] = 0.0;
    }
        
    double reached_dt;
    int stopping_condition_occurred;

    double t=0;
    double dt = 999986; // gives error almost immediately
    //double dt = 200000; // takes a while, but runs without errors
    
    printf("pre\n");
    mstar_print_state(R);
    double E_init = compute_nbody_total_energy(R);
    
    printf("running\n");
    run_integrator(R, dt, &reached_dt, &stopping_condition_occurred);
    
    double E_fin = compute_nbody_total_energy(R);
    
    printf("post\n");
    mstar_print_state(R);
    
    printf("E_init %g E_fin %g E_err %g\n",E_init,E_fin,fabs((E_init-E_fin)/E_init));
    
    free_data(R);
    
    return flag;
}

/**********
 * Flybys *
***********/
    
int test_flybys()
{
    printf("test.cpp -- test_flybys\n");
    
    int flag = 0;
    flag += test_flybys_integrals();
    /* Note: test_flybys_perturber_sampling is carried out within Python (test_mse.py) */
    flag += test_flybys_compute_effects_of_flyby_on_system();
    
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

int test_flybys_perturber_sampling(double R_enc, double n_star, double sigma_rel, double M_int, double *M_per, double *b_vec_x, double *b_vec_y, double *b_vec_z, double *V_vec_x, double *V_vec_y, double *V_vec_z) 
{
    int flag = 0;
    flybys_mass_distribution = 0;
    flybys_mass_distribution_lower_value = 0.1;
    flybys_mass_distribution_upper_value = 100.0;
    flybys_mass_distribution = 0;
    flybys_encounter_sphere_radius = R_enc;
    flybys_stellar_density = n_star;
    flybys_stellar_relative_velocity_dispersion = sigma_rel;
    flybys_internal_mass = M_int;
        
    ParticlesMap particlesMap;
    int N_bodies = 2;
    double masses[2] = {10.0,12.0};
    int stellar_types[2] = {1,1};
    int object_types[2] = {1,1};
    double smas[1] = {10.0};
    double es[1] = {0.01};
    double TAs[1] = {0.01};
    double INCLs[1] = {0.01};
    double APs[1] = {0.01};
    double LANs[1] = {0.01};

    create_nested_system(particlesMap,N_bodies,masses,stellar_types,object_types,smas,es,TAs,INCLs,APs,LANs);

    bool apply_flyby;
    double t_next_encounter;
    double b_vec[3],V_vec[3],e_vec_unit_per[3],h_vec_unit_per[3];
    int N_enc=0;
    int N_not_impulsive=0;
    
    //sample_next_flyby(&particlesMap, &apply_flyby, &t_next_encounter, &N_enc, &N_not_impulsive, M_per, b_vec, V_vec);
    int flyby_type;
    double e_per,Q_per;
    sample_next_flyby(&particlesMap, &apply_flyby, &flyby_type, &flybys_t_next_encounter, &N_enc, &N_not_impulsive, M_per, b_vec, V_vec, &e_per, &Q_per, e_vec_unit_per, h_vec_unit_per);
    *b_vec_x = b_vec[0];
    *b_vec_y = b_vec[1];
    *b_vec_z = b_vec[2];
    *V_vec_x = V_vec[0];
    *V_vec_y = V_vec[1];
    *V_vec_z = V_vec[2];
    
    return flag;
}

int test_flybys_compute_effects_of_flyby_on_system()
{
    printf("test.cpp -- test_flybys_compute_effects_of_flyby_on_system\n");
    int flag = 0;
    flybys_mass_distribution = 0;
    flybys_mass_distribution_lower_value = 0.1;
    flybys_mass_distribution_upper_value = 100.0;
    flybys_mass_distribution = 0;
    flybys_encounter_sphere_radius = 1.0e5;
    flybys_stellar_density = 0.1*CONST_PER_PC3;
    flybys_stellar_relative_velocity_dispersion = 30.0*CONST_KM_PER_S;
    flybys_reference_binary = -1;
    
    ParticlesMap particlesMap;
    int N_bodies = 2;
    double m1 = 10.0;
    double m2 = 8.0;
    flybys_internal_mass = m1+m2;
    
    double masses[2] = {m1,m2};
    int stellar_types[2] = {1,1};
    int object_types[2] = {1,1};
    double smas[1] = {10.0};
    double es[1] = {0.01};
    double TAs[1] = {0.01};
    double INCLs[1] = {0.01};
    double APs[1] = {0.01};
    double LANs[1] = {0.01};

    create_nested_system(particlesMap,N_bodies,masses,stellar_types,object_types,smas,es,TAs,INCLs,APs,LANs);

    Particle *s1 = particlesMap[0];
    Particle *s2 = particlesMap[1];
    
    s1->R_vec[0] = 1.0e2;
    s1->R_vec[1] = 1.0e1;
    s1->R_vec[2] = 2.0e1;

    s2->R_vec[0] = -1.0e2;
    s2->R_vec[1] = -1.0e1;
    s2->R_vec[2] = -4.0e1;
        
    double M_per = 100.0;
    double b_per_vec[3] = {1.0e3,1.0e3,1.0e2};
    double V_per_vec[3] = {10.0*CONST_KM_PER_S,5.0*CONST_KM_PER_S,15.0*CONST_KM_PER_S};
    bool unbound_orbits;
    int integration_flag;
    int flyby_type = 1; /* impulsive */
    double e_per,Q_per;
    double e_vec_unit_per[3],h_vec_unit_per[3];
    
    compute_effects_of_flyby_on_system(&particlesMap, flyby_type, M_per,b_per_vec,V_per_vec,e_per,Q_per,e_vec_unit_per,h_vec_unit_per,&unbound_orbits,false,&integration_flag);
    
    double tol = 1.0e-2; /* Note: smaller tolerance than 1e-2 yields discrepancies; most likely due to slightly different values of constants used */
    if (!equal_number(s1->instantaneous_perturbation_delta_VX,0.818813,tol))
    {
        printf("test.cpp -- test_compute_effects_of_flyby_on_system -- error in s1 Delta VX %g\n",s1->instantaneous_perturbation_delta_VX);
        flag = 1;
    }
    if (!equal_number(s1->instantaneous_perturbation_delta_VY,1.35592,tol))
    {
        printf("test.cpp -- test_compute_effects_of_flyby_on_system -- error in s1 Delta VY %g\n",s1->instantaneous_perturbation_delta_VY);
        flag = 1;
    }
    if (!equal_number(s1->instantaneous_perturbation_delta_VZ,-0.99785,tol))
    {
        printf("test.cpp -- test_compute_effects_of_flyby_on_system -- error in s1 Delta VZ %g\n",s1->instantaneous_perturbation_delta_VZ);
        flag = 1;
    }
    if (!equal_number(s2->instantaneous_perturbation_delta_VX,0.888464,tol))
    {
        printf("test.cpp -- test_compute_effects_of_flyby_on_system -- error in s2 Delta VX %g\n",s2->instantaneous_perturbation_delta_VX);
        flag = 1;
    }
    if (!equal_number(s2->instantaneous_perturbation_delta_VY,1.14714,tol))
    {
        printf("test.cpp -- test_compute_effects_of_flyby_on_system -- error in s2 Delta VY %g\n",s2->instantaneous_perturbation_delta_VY);
        flag = 1;
    }
    if (!equal_number(s2->instantaneous_perturbation_delta_VZ,-0.97469,tol))
    {
        printf("test.cpp -- test_compute_effects_of_flyby_on_system -- error in s2 Delta VZ %g\n",s2->instantaneous_perturbation_delta_VZ);
        flag = 1;
    }
    
//    printf("S1 %g %g %g\n",s1->instantaneous_perturbation_delta_VX,s1->instantaneous_perturbation_delta_VY,s1->instantaneous_perturbation_delta_VZ);
//    printf("S2 %g %g %g\n",s2->instantaneous_perturbation_delta_VX,s2->instantaneous_perturbation_delta_VY,s2->instantaneous_perturbation_delta_VZ);
    return flag;
}


/*********************
 * Stellar evolution *
**********************/

int test_stellar_evolution(int mode)
{
    printf("test.cpp -- test_stellar_evolution\n");
    
    int flag = 0;
    flag += test_spin_conversion();
    flag += test_apsidal_motion_constant();
    flag += test_sse();
    //flag += test_sse_custom();
    flag += test_remove_massless_remnants_from_system();
    flag += test_NS_models(mode);
    flag += test_determine_sse_compact_object_radius(mode);
    
    if (flag == 0)
    {
        printf("test.cpp -- test_stellar_evolution -- passed\n");
    }
    
    return flag;
}

int test_spin_conversion()
{
    printf("test.cpp -- test_spin_conversion\n");
    int flag = 0;

    double stellar_types[2] = {1,10};
    double Omega = 1.0;
    double mass = 1.0;
    double core_mass = 0.5;
    double radius = 1.0;
    double core_radius = 0.5;
    double k2 = 0.02;
    double k3 = 0.02;
    
    double stellar_type;
    double S;
    double Omega2;
    double tol = 1.0e-14;
    for (int i=0; i<2; i++)
    {
        stellar_type = stellar_types[i];
    
        S = compute_spin_angular_momentum_from_spin_frequency(Omega, stellar_type, 1, mass, core_mass, radius, core_radius, k2, k3);
        Omega2 = compute_spin_frequency_from_spin_angular_momentum(S, stellar_type, 1, mass, core_mass, radius, core_radius, k2, k3);

        if (!equal_number(Omega,Omega2,tol))
        {
            printf("test.cpp -- test_spin_conversion -- error %g %g\n",Omega,Omega2);
            flag = 1;
        }
    }

    return flag;
}

int test_apsidal_motion_constant()
{
    /* Currently tests for NaNs in a few cases. */
    printf("test.cpp -- test_apsidal_motion_constant\n");

    int flag = 0;    
    int N_m=6;
    int N_st=15;
    double masses[N_m] = {0.08,0.5,1.0,10.0,20.5,100.0};
    int stellar_types[N_st] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14};
    
    double k_AM;
    for (int i=0; i<N_m; i++)
    {
        for (int j=0; j<N_st; j++)
        {
            ParticlesMap particlesMap;
            Particle *star = new Particle(0, false);
            particlesMap[0] = star;
            
            star->object_type = 1;
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
    
    int NS_model_old = NS_model;
    NS_model = 0; /* With other NS models, final spins of NSs will not be consistent with vanilla SSE */
    
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
        test_sse_specific_model(masses[i],z,&kw_final,&m_init_final,&m_final,&R_final,&ospin_final,&L_final,&m_core_final,&m_env_final,&epoch_final);

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
        
        test_sse_specific_model_evolve_function(masses[i],z,&kw_final,&m_init_final,&m_final,&R_final,&ospin_final,&L_final,&m_core_final,&m_env_final,&epoch_final);
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

        test_sse_specific_model_stopping_conditions(masses[i],z,&kw_final,&m_init_final,&m_final,&R_final,&ospin_final,&L_final,&m_core_final,&m_env_final,&epoch_final);
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
    
    NS_model = NS_model_old;
    return flag;
}

int test_sse_specific_model(double m, double z, int *kw_final, double *m_init_final, double *m_final, double *R_final, double *ospin_final, double *L_final, double *m_core_final, double *m_env_final, double *epoch_final)
{

    int flag = 0;
    
    ParticlesMap particlesMap;
    Particle *star = new Particle(0, false);
    particlesMap[0] = star;
    
    star->object_type = 1;
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

        evolve_stars(&particlesMap,t_old,t,&dt,false,&apply_SNe_effects,&integration_flag);
        if (apply_SNe_effects == true)
        {
            handle_SNe_in_system(&particlesMap,&unbound_orbits,&integration_flag);
        }
        else
        {
            update_stellar_evolution_quantities_directly_nbody(&particlesMap,t,t-t_old);
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

int test_sse_specific_model_evolve_function(double m, double z, int *kw_final, double *m_init_final, double *m_final, double *R_final, double *ospin_final, double *L_final, double *m_core_final, double *m_env_final, double *epoch_final)
{

    int flag = 0;
    
    ParticlesMap particlesMap;
    Particle *star = new Particle(0, false);
    particlesMap[0] = star;
    
    star->object_type = 1;
    star->mass = m;
    star->sse_initial_mass = m;
    star->metallicity = z;
    star->stellar_type = 1;
    
    initialize_stars(&particlesMap);
    double t_old,t;
    t_old = t = 0.0;
    double dt = 1.0e9;
    
    double t_max = 1.2e10;

    bool unbound_orbits;
    bool apply_SNe_effects;
    
    double output_time,hamiltonian;
    int state, CVODE_flag, CVODE_error_code;
    int integration_flag = 0;
    initialize_code(&particlesMap);
    while (t < t_max)
    {
        evolve(&particlesMap,t,t+dt, &output_time, &hamiltonian, &state, &CVODE_flag, &CVODE_error_code, &integration_flag);
        
        t = output_time;

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

int test_sse_specific_model_stopping_conditions(double m, double z, int *kw_final, double *m_init_final, double *m_final, double *R_final, double *ospin_final, double *L_final, double *m_core_final, double *m_env_final, double *epoch_final)
{

    int flag = 0;
    
    ParticlesMap particlesMap;
    Particle *star = new Particle(0, false);
    particlesMap[0] = star;
    
    star->object_type = 1;
    star->mass = m;
    star->sse_initial_mass = m;
    star->metallicity = z;
    star->stellar_type = 1;
    
    initialize_stars(&particlesMap);
    double t_old,t;
    t_old = t = 0.0;
    double dt = 1.0;
    
    int sec_flag;
    double t_max = 1.2e10;
    int integration_flag=0;
    bool unbound_orbits;
    bool apply_SNe_effects;
    double t_out;
    double hamiltonian;
    int CVODE_flag,CVODE_error_code;
    double spin_vec_norm;
    double dt_stev,dt_stev_new,dt_true;
    double t_stev = 1.0;
    int kw_old,kw_new;
    
    int k=0;
    int k_max=4;
    while (t_old < t_max)
    {
        kw_old = star->stellar_type;
        evolve_stars(&particlesMap,t_old,t_stev,&dt_stev_new,false,&apply_SNe_effects,&integration_flag);
        kw_new = star->stellar_type;

        dt_stev = t_stev-t_old;

        dt = dt_stev;
        
        if (kw_old == kw_new)
        {
            k += 1;
        }
        
        if (kw_old == kw_new and kw_new < 10 and k == k_max)
        {
            dt = dt_stev * 0.1;
            //printf("adjust dt %g k %d\n",dt,k);
        }
        if (k>k_max)
        {
            k=0;
        }

        if (apply_SNe_effects == true)
        {
            handle_SNe_in_system(&particlesMap,&unbound_orbits,&integration_flag);
        }
        else
        {

            star->mass += star->mass_dot_wind * dt;
            star->radius += star->radius_dot * dt;

            spin_vec_norm = norm3(star->spin_vec);
            for (int i=0; i<3; i++)
            {
                star->spin_vec[i] += dt * star->ospin_dot * (star->spin_vec[i]/spin_vec_norm);
            }
            
            update_stellar_evolution_quantities_directly_secular(&particlesMap, dt_stev, t_old, t_old + dt);
        
        }

        t_old += dt;

        dt = dt_stev_new;
        t_stev+=dt;
        
        if (t_old>t_max)
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

int test_sse_custom()
{
    printf("test.cpp -- test_sse_custom\n");
    struct value1__ value1_;
    struct value2__ value2_;
    struct value3__ value3_;
    struct value4__ value4_;
    struct value5__ value5_;
    struct flags__ flags_;
    struct points__ points_;
    struct sse_error_output__ sse_error_output_;

    
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
    

    int i;
    int kw,kw_desired;
    double tm,tn;
    double age;
    double *GB,*tscls,*lums;
    double r,lum,mc,rc,menv,renv,k2;
            
    double z = 0.02;
            
    double *zpars;
    zpars = new double[20];
    //double zpars[20];
    zcnsts_(&z,zpars);
    
    double mt = 12.385436420676200;
    //double sse_initial_mass = 6.125196978052001;
    double sse_initial_mass = 6.125196978052001;
    double epoch = 22.550149349880343;
    //double epoch = 0;
    //age = 107.105369463623362;
    age = 0;
    
    GB = new double[10];
    tscls = new double[20];
    lums = new double[10];    
    
//            star_(&kw, &mass, &mt, &tm, &tn, tscls, lums, GB, zpars);
//            hrdiag_(&mass,&age,&mt,&tm,&tn,tscls,lums,GB,zpars,
//                &r,&lum,&kw,&mc,&rc,&menv,&renv,&k2);
    
    double dtp=0.0;
    double ospin=1490.002987504306020;
    double tms;
    double tphys = 129.655518813503704;
    double tphysf = 129.657412089438111;
    //double tphysf = 1500;
    double dt;
    
    kw = 2;

    printf("pre evolv1_() call --  kw %d mt %.15f sse_initial_mass %.15f r %.15f epoch %.15f age %.15f tphys %.15f tphysf %.15f ospin %.15f tms %.15f epoch %.15f k2 %.15f rc %.15f mc %.15f menv %.15f renv %.15f \n",kw,mt,sse_initial_mass,r,epoch,age,tphys,tphysf,ospin,tms,epoch,k2,rc,mc,menv,renv);

    evolv1_(&kw,&sse_initial_mass,&mt,&r,&lum,&mc,&rc,&menv,&renv,&ospin,&epoch,&tms,&tphys,&tphysf,&dtp,&z,zpars,&k2);
    
    printf("post evolv1_() call --  kw %d mt %.15f sse_initial_mass %.15f r %.15f epoch %.15f age %.15f tphys %.15f tphysf %.15f ospin %.15f tms %.15f epoch %.15f k2 %.15f rc %.15f mc %.15f menv %.15f renv %.15f \n",kw,mt,sse_initial_mass,r,epoch,age,tphys,tphysf,ospin,tms,epoch,k2,rc,mc,menv,renv);

    
    return 0;
}

int test_kick_velocity(int kick_distribution, double m, int *kw, double *v_norm)
{
    ParticlesMap particlesMap;
    Particle *star = new Particle(0, false);
    particlesMap[0] = star;
    
    star->object_type = 1;
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
    
    int integration_flag = 0;
    bool apply_SNe_effects = false;
    while (apply_SNe_effects == false)
    {
        t+=dt;
        //printf("t %g t_old %g dt %g\n",t,t_old,dt);
        evolve_stars(&particlesMap,t_old,t,&dt,false,&apply_SNe_effects,&integration_flag);
        if (apply_SNe_effects == true)
        {
            sample_kick_velocity(star,&vx,&vy,&vz);
        }
        else
        {
            update_stellar_evolution_quantities_directly_nbody(&particlesMap,t,t-t_old);
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
    
int test_remove_massless_remnants_from_system()
{
    printf("test.cpp -- test_remove_massless_remnants_from_system\n");
    
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
    int object_types[4] = {1,1,1,1};
    double smas[3] = {a1,a2,a3};
    double es[3] = {0.01,0.01,0.01};
    double TAs[3] = {0.01,0.01,0.01};
    double INCLs[3] = {0.01,0.01,0.01};
    double APs[3] = {0.01,0.01,0.01};
    double LANs[3] = {0.01,0.01,0.01};
    
    create_nested_system(particlesMap,N_bodies,masses,stellar_types,object_types,smas,es,TAs,INCLs,APs,LANs);// = create_nested_system();
//    printf("post s %d b %d %g r %g\n",particlesMap2.size(),particlesMap2[0]->is_binary,particlesMap2[0]->mass,particlesMap2[0]->radius);
    initialize_code(&particlesMap);
    
    particlesMap[0]->stellar_type = 15;
    
    int integration_flag = 1;
    remove_massless_remnants_from_system(&particlesMap,&integration_flag);
    
    return 0;
}

int test_NS_models(int mode)
{

    printf("test.cpp -- test_NS_models\n");

    int NS_model_old = NS_model;
    NS_model = 1;

    ParticlesMap particlesMap;
    Particle *star = new Particle(0, false);
    particlesMap[0] = star;
    double m = 15.0;
    
    star->object_type = 1;
    star->mass = m;
    star->sse_initial_mass = m;
    star->stellar_type = 13;
    
    //random_seed = 3;
    initialize_code(&particlesMap);
    //particlesMap[0]->magnetic_field_strength_gauss = 1.0e10;
    double P_s0 = compute_spin_period_from_spin_angular_frequency(norm3(particlesMap[0]->spin_vec)) * yr_to_s;
    double B_G0 = particlesMap[0]->magnetic_field_strength_gauss;
    //printf("P_s0 %g B_G0 %g\n",P_s0,B_G0);
    
    double t_old,t;
    t_old = t = 0.0;
    double dt = 1.0e8;
    
    double tend = 5e9;

    std::ofstream myfile;
    if (mode == 2)
    {
        myfile.open ("test_NS_models.txt");
    }

    int integration_flag = 0;

    double output_time;
    double hamiltonian;
    int state,CVODE_flag,CVODE_error_code;
    
    while (t < tend)
    {
        evolve(&particlesMap,t,t+dt,&output_time,&hamiltonian,&state,&CVODE_flag,&CVODE_error_code,&integration_flag);
        double P_s = compute_spin_period_from_spin_angular_frequency(norm3(particlesMap[0]->spin_vec)) * yr_to_s;
        double B_G = particlesMap[0]->magnetic_field_strength_gauss;

        t = output_time;

        if (mode == 2)
        {
            myfile << t << "," << P_s << "," << B_G << "," << P_s0 << "," << B_G0 << "\n";
        }

    }

    if (mode == 2)
    {
        myfile.close();    
    }

    NS_model = NS_model_old;
    
    return 0;
}

int test_determine_sse_compact_object_radius(int mode)
{
    printf("test.cpp -- test_determine_sse_compact_object_radius\n");
    int flag = 0;
    
    std::ofstream myfile;
    if (mode == 2)
    {
        myfile.open ("test_determine_sse_compact_object_radius.txt");
    }
    
    double m = 0.01;
    double dm = 0.001;
    double r;
    
    while (m < chandrasekhar_mass)
    {
        r = determine_sse_compact_object_radius_RSun(10,m);
        
        if (mode == 2)
        {
            myfile << m << "," << r << "\n";
        }
        
        m += dm;
    }

    if (mode == 2)
    {
        myfile.close();    
    }

    r = determine_sse_compact_object_radius_RSun(14,m);
    double r_an = (1.0/CONST_R_SUN) * 2.0*CONST_G*m/(CONST_C_LIGHT_P2);
    if ( !equal_number(r,r_an,1.0e-2) )
    {
        printf("test.cpp -- error in test_determine_sse_compact_object_radius! %g %g\n",r,r_an);
        flag += 1;
    }

    return flag;
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
    flag += test_mass_accretion_events_with_degenerate_objects();
    flag += test_compute_bse_mass_transfer_amount_averaged();
    flag += test_binary_common_envelope_evolution();
    
    if (flag == 0)
    {
        printf("test.cpp -- test_binary_evolution -- passed\n");
    }
    
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
    int object_types[4] = {1,1,1,1};
    double smas[3] = {a1,a2,a3};
    double es[3] = {0.01,0.01,0.01};
    double TAs[3] = {0.01,0.01,0.01};
    double INCLs[3] = {0.01,0.01,0.01};
    double APs[3] = {0.01,0.01,0.01};
    double LANs[3] = {0.01,0.01,0.01};
    
    create_nested_system(particlesMap,N_bodies,masses,stellar_types,object_types,smas,es,TAs,INCLs,APs,LANs);// = create_nested_system();
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
    int stellar_types[2] = {1,1};
    int object_types[2] = {1,1};
    double smas[1] = {a};
    double es[1] = {e};
    double TAs[1] = {0.01};
    double INCLs[1] = {0.01};
    double APs[1] = {0.01};
    double LANs[1] = {0.01};
    
    create_nested_system(particlesMap,N_bodies,masses,stellar_types,object_types,smas,es,TAs,INCLs,APs,LANs);

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

    return flag;
}

int test_mass_accretion_events_with_degenerate_objects()
{
    printf("test.cpp -- test_mass_accretion_events_with_degenerate_objects\n");

    int flag = 0;

    //verbose_flag = 1;

    /* CO WD exceeding M_Ch */
    ParticlesMap particlesMap;
    int N_bodies = 2;
    double m1 = 20.0;
    double m2 = 6.0;
    double a = 1.0;
    double e = 0.95;
    
    double masses[2] = {m1,m2};
    int stellar_types[2] = {4,11};
    int object_types[2] = {1,1};
    double smas[1] = {a};
    double es[1] = {e};
    double TAs[1] = {0.01};
    double INCLs[1] = {0.01};
    double APs[1] = {0.01};
    double LANs[1] = {0.01};
    
    create_nested_system(particlesMap,N_bodies,masses,stellar_types,object_types,smas,es,TAs,INCLs,APs,LANs);

    initialize_code(&particlesMap);
    
    int donor_index = 0;
    int accretor_index = 1;
    int parent_index = 2;
        
    Particle *star1 = particlesMap[donor_index];
    Particle *star2 = particlesMap[accretor_index];
    Particle *orbit = particlesMap[parent_index];

    int integration_flag = 0;
    double t_old = 0.0;
    double t = 1.0e1;
    double dt_binary_evolution;

    star2->mass = 1.4;

    star2->mass_dot_RLOF = 0.1;
    handle_mass_accretion_events_with_degenerate_objects(&particlesMap, t_old, t, &integration_flag, &dt_binary_evolution);
    
    //Log_type &last_entry = logData.back();
    //printf("log SNe type %d info %d\n",last_entry.log_info.SNe_type,last_entry.log_info.SNe_info);
    
    double start_time = 0.0;
    double end_time = 1.0e1;
    double output_time,hamiltonian;
    int state,CVODE_flag,CVODE_error_code;
    evolve(&particlesMap,start_time,end_time,&output_time,&hamiltonian,&state,&CVODE_flag,&CVODE_error_code,&integration_flag);

    if (particlesMap.size() != 1)
    {
        printf("test.cpp -- test_mass_accretion_events_with_degenerate_objects -- error: incorrect number (%d) of particles after SNe \n",particlesMap.size());
        flag = 1;
    }

    reset_interface();

    /* ONe WD electron capture SNe */
    int ECSNe_model_old = ECSNe_model;
    int NS_model_old = NS_model;
    ECSNe_model = 1;
    NS_model = 1;
    stellar_types[1] = 12;
    masses[1] = 7.0;
    
    create_nested_system(particlesMap,N_bodies,masses,stellar_types,object_types,smas,es,TAs,INCLs,APs,LANs);
    
    initialize_code(&particlesMap);

    integration_flag = 0;
    t_old = 0.0;
    t = 1.0e1;

    star1 = particlesMap[donor_index];
    star2 = particlesMap[accretor_index];
    orbit = particlesMap[parent_index];
    star2->mass = 1.4;
    star2->mass_dot_RLOF = 0.1;

    handle_mass_accretion_events_with_degenerate_objects(&particlesMap, t_old, t, &integration_flag, &dt_binary_evolution);
    
    if (particlesMap.size() != 3 or star2->stellar_type != 13) /* Check that NS was formed */
    {
        printf("test.cpp -- test_mass_accretion_events_with_degenerate_objects -- error: incorrect number (%d) of particles after ECSN or incorrect stellar type (%d) \n",particlesMap.size(),star2->stellar_type);
        flag = 1;
    }

    reset_interface();

    ECSNe_model = ECSNe_model_old;
    NS_model = NS_model_old;

    
    /* Multiple events in a system (2+2) -- skip by default; uncomment the next line to include */
    //return flag;
    
    /* CO WD M>M_Ch */
    int N_bodies2 = 4;
    double masses2[4] = {1.0,6.0,1.5,6.1}; // CO WDs
    int stellar_types2[4] = {1,11,1,11}; // CO WDs
    
    int object_types2[4] = {1,1,1,1};
    double smas2[3] = {10.0,20.0,10000.0};
    double es2[3] = {0.01,0.01,0.01};
    double TAs2[3] = {0.01,0.01,0.01};
    double INCLs2[3] = {0.01,0.01,0.01};
    double APs2[3] = {0.01,0.01,0.01};
    double LANs2[3] = {0.01,0.01,0.01};
    
    create_2p2_quadruple_system(particlesMap,masses2,stellar_types2,object_types2,smas2,es2,TAs2,INCLs2,APs2,LANs2);
    
    initialize_code(&particlesMap);

    integration_flag = 0;
    t_old = 0.0;
    t = 1.0e1;

    star1 = particlesMap[0];
    star2 = particlesMap[1];
    Particle *star3 = particlesMap[2];
    Particle *star4 = particlesMap[3];

    star2->mass = 1.4;
    star2->mass_dot_RLOF = 0.1;
    star4->mass = 1.4;
    star4->mass_dot_RLOF = 0.1;

    handle_mass_accretion_events_with_degenerate_objects(&particlesMap, t_old, t, &integration_flag, &dt_binary_evolution);
    //print_system(&particlesMap,integration_flag);
    
    if (particlesMap.size() != 2) /* Check that two WDs were destroyed */
    {
        printf("test.cpp -- test_mass_accretion_events_with_degenerate_objects -- error: incorrect number (%d) of particles after two SNe\n",particlesMap.size());
        flag = 1;
    }

    evolve(&particlesMap,start_time,end_time,&output_time,&hamiltonian,&state,&CVODE_flag,&CVODE_error_code,&integration_flag);
    reset_interface();
    

    /* ONeWD ECSN */
    masses2[1] = 7.0;
    masses2[3] = 7.1;
    stellar_types2[1] = 12;
    stellar_types2[3] = 12;

    create_2p2_quadruple_system(particlesMap,masses2,stellar_types2,object_types2,smas2,es2,TAs2,INCLs2,APs2,LANs2);
    
    ECSNe_model = 1;
    NS_model = 1;

    initialize_code(&particlesMap);

    integration_flag = 0;
    t_old = 0.0;
    t = 1.0e1;

    star1 = particlesMap[0];
    star2 = particlesMap[1];
    star3 = particlesMap[2];
    star4 = particlesMap[3];

    star2->mass = 1.4;
    star2->mass_dot_RLOF = 0.1;
    star4->mass = 1.4;
    star4->mass_dot_RLOF = 0.1;

    handle_mass_accretion_events_with_degenerate_objects(&particlesMap, t_old, t, &integration_flag, &dt_binary_evolution);
    //print_system(&particlesMap,integration_flag);
    
    if (particlesMap.size() != 4 or star2->stellar_type != 13 or star4->stellar_type != 13) /* Check that two NSs were formed */
    {
        printf("test.cpp -- test_mass_accretion_events_with_degenerate_objects -- error: incorrect number (%d) of particles after two ECSN or incorrect stellar types (%d; %d) \n",particlesMap.size(),star2->stellar_type,star4->stellar_type);
        flag = 1;
    }

    evolve(&particlesMap,start_time,end_time,&output_time,&hamiltonian,&state,&CVODE_flag,&CVODE_error_code,&integration_flag);
    reset_interface();


    ECSNe_model = ECSNe_model_old;
    NS_model = NS_model_old;

    return flag;
}

int test_compute_bse_mass_transfer_amount_averaged()
{
    printf("test.cpp -- test_compute_bse_mass_transfer_amount_averaged\n");
    
    double m_donor = 20.0;
    double m_accretor = 10.0;
    double core_mass_donor = 1.0;
    double R_donor = 1.0;
    double a = 1.0;
    double e = 1.0e-10;
    double dt = 1.0;
    double t_dyn_donor = 1.0;
    double t_KH_donor = 1.0;
    int kw1 = 1;
    
    int binary_evolution_numerical_mass_transfer_rate_number_of_points_old = binary_evolution_numerical_mass_transfer_rate_number_of_points;
    
    /* Test circular case with default value of binary_evolution_numerical_mass_transfer_rate_number_of_points */
    double dm_av = compute_bse_mass_transfer_amount_averaged(kw1,m_donor,core_mass_donor,R_donor,m_accretor,a,e,dt,t_dyn_donor,t_KH_donor);
    double q = m_donor/m_accretor;
    double dm_av_circ = compute_bse_mass_transfer_amount(kw1, m_donor, core_mass_donor, R_donor, roche_radius_pericenter_eggleton(a*(1.0-e),q), dt, t_dyn_donor, t_KH_donor);

    int flag = 0;
    double tol = 1e-8;
    if ( !equal_number(dm_av,dm_av_circ,tol) )
    {
        printf("test.cpp -- test_compute_bse_mass_transfer_amount_averaged -- error: dm_av %g dm_av_circ %g\n",dm_av,dm_av_circ);
        flag = 1;
    }

    /* Test convergence (eccentric case, with different tolerances up to and including the default setting) */
    e = 0.999;
    int N_s = 7;
    int N;
    double fs[N_s] = {0.1,0.2,0.4,0.6,0.8,0.9,1.0};
    double eps;
    double dm_av_old;
    for (int i=0; i<N_s; i++)
    {
        N = int(binary_evolution_numerical_mass_transfer_rate_number_of_points_old * fs[i] );
        binary_evolution_numerical_mass_transfer_rate_number_of_points = N;
        dm_av = compute_bse_mass_transfer_amount_averaged(kw1,m_donor,core_mass_donor,R_donor,m_accretor,a,e,dt,t_dyn_donor,t_KH_donor);

        eps = fabs((dm_av_old-dm_av)/dm_av_old);
        
        dm_av_old = dm_av;
       
        //printf("N %d dm_av %g eps %g\n",N,dm_av,eps);
    }
    tol = 1e-4;
    if ( eps > tol )
    {
        printf("test.cpp -- test_compute_bse_mass_transfer_amount_averaged -- error (convergence in eccentric case): eps %g\n",eps);
        flag = 1;
    }


    binary_evolution_numerical_mass_transfer_rate_number_of_points = binary_evolution_numerical_mass_transfer_rate_number_of_points_old;

    return flag;
}



int test_mass_transfer_special_cases_old()
{
    printf("test.cpp -- test_mass_transfer_special_cases\n");

    int flag = 0;

    /* CO WD exceeding M_Ch */
    ParticlesMap particlesMap;
    int N_bodies = 2;
    double m1 = 20.0;
    double m2 = 6.0;
    double a = 1.0;
    double e = 0.95;
    
    double masses[2] = {m1,m2};
    int stellar_types[2] = {4,11};
    int object_types[2] = {1,1};
    double smas[1] = {a};
    double es[1] = {e};
    double TAs[1] = {0.01};
    double INCLs[1] = {0.01};
    double APs[1] = {0.01};
    double LANs[1] = {0.01};
    
    create_nested_system(particlesMap,N_bodies,masses,stellar_types,object_types,smas,es,TAs,INCLs,APs,LANs);

    initialize_code(&particlesMap);
    
    int donor_index = 0;
    int accretor_index = 1;
    int parent_index = 2;
        
    Particle *star1 = particlesMap[donor_index];
    Particle *star2 = particlesMap[accretor_index];
    Particle *orbit = particlesMap[parent_index];

    int integration_flag = 0;
    double t_old = 0.0;
    double t = 1.0e1;
    double dt_binary_evolution;

    star2->mass = 1.4;
    //verbose_flag = 0;
    binary_stable_mass_transfer_evolution(&particlesMap, parent_index, donor_index, accretor_index, t_old, t, &integration_flag, &dt_binary_evolution);
    
    double start_time = 0.0;
    double end_time = 1.0e1;
    double output_time,hamiltonian;
    int state,CVODE_flag,CVODE_error_code;
    evolve(&particlesMap,start_time,end_time,&output_time,&hamiltonian,&state,&CVODE_flag,&CVODE_error_code,&integration_flag);

    if (particlesMap.size() != 1)
    {
        printf("test.cpp -- test_mass_transfer_special_cases -- error: incorrect number (%d) of particles after SNe \n",particlesMap.size());
        flag = 1;
    }
    
    reset_interface();
    
    /* ONe WD electron capture SNe */
    int ECSNe_model_old = ECSNe_model;
    int NS_model_old = NS_model;
    ECSNe_model = 1;
    NS_model = 1;
    stellar_types[1] = 12;
    masses[1] = 7.0;
    
    create_nested_system(particlesMap,N_bodies,masses,stellar_types,object_types,smas,es,TAs,INCLs,APs,LANs);
    
    initialize_code(&particlesMap);

    integration_flag = 0;
    t_old = 0.0;
    t = 1.0e1;

    star1 = particlesMap[donor_index];
    star2 = particlesMap[accretor_index];
    orbit = particlesMap[parent_index];
    star2->mass = 1.4;
    //verbose_flag = 0;

    binary_stable_mass_transfer_evolution(&particlesMap, parent_index, donor_index, accretor_index, t_old, t, &integration_flag, &dt_binary_evolution);
    
    if (particlesMap.size() != 3 or star2->stellar_type != 13) /* Check that NS was formed */
    {
        printf("test.cpp -- test_mass_transfer_special_cases -- error: incorrect number (%d) of particles after ECSN or incorrect stellar type (%d) \n",particlesMap.size(),star2->stellar_type);
        flag = 1;
    }

    ECSNe_model = ECSNe_model_old;
    NS_model = NS_model_old;

    return flag;
}




int test_collisions()
{
    printf("test.cpp -- test_collisions\n");
    int flag;
    
    int i,j;
    for (i=1; i<=14; i++)
    {
        for (j=1; j<=14; j++)
        {
            if (j<i)
            {
                continue;
            }

            printf("i %d j %d\n",i,j);
            //if (i!=1 or j!=4)
            //if (i>5 and i<12 or j>5 and j<12)
            if (i!=2 or j!= 5)
            {
                //continue;
            }
            flag += test_collision_stars(10.0,i,13,j,0);
            random_seed = 0;
            
        }
    }
    
    if (flag == 0)
    {
        printf("test.cpp -- test_collisions -- passed\n");
    }

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
    int object_types[4] = {1,1,1,1};
    double smas[3] = {10.0,100.0,10000.0};
    double es[3] = {0.01,0.01,0.01};
    double TAs[3] = {0.01,0.01,0.01};
    double INCLs[3] = {0.01,0.01,0.01};
    double APs[3] = {0.01,0.01,0.01};
    double LANs[3] = {0.01,0.01,0.01};
    random_seed=6;
    
    create_nested_system(particlesMap,N_bodies,masses,stellar_types,object_types,smas,es,TAs,INCLs,APs,LANs);// = create_nested_system();
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
    
    print_system(&particlesMap,integration_flag);
    
        ParticlesMapIterator it_p,it_q;
    for (it_p = particlesMap.begin(); it_p != particlesMap.end(); it_p++)
    {
        Particle *P_p = (*it_p).second;

        P_p->parent = -1;
    }
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


int test_binary_common_envelope_evolution()
{
    printf("test.cpp -- test_binary_common_envelope_evolution\n");
    
    bool verbose_testing = false;
    if (verbose_testing == true)
    {
        verbose_flag = 1;
    }
    
    int flag = 0;
    
    for (int i=0; i<2; i++)
    {
        //printf("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n");
        if (i==0)
        {
            //effective_radius_multiplication_factor_for_collisions_stars = 3.0;
            
            ParticlesMap particlesMap;
            int N_bodies = 3;
            double masses[3] = {5.0,0.8,1.0};
            int stellar_types[3] = {5,1,1};
            int object_types[3] = {1,1,1};
            double smas[2] = {0.1,5.0};
            double es[2] = {0.01,0.01};
            double TAs[2] = {0.01,0.01};
            double INCLs[2] = {0.01,0.01};
            double APs[2] = {0.01,0.01};
            double LANs[2] = {0.01,0.01};
            random_seed=6;
            
            create_nested_system(particlesMap,N_bodies,masses,stellar_types,object_types,smas,es,TAs,INCLs,APs,LANs);// = create_nested_system();
            initialize_code(&particlesMap);

            int integration_flag = 0;
            
            if (verbose_testing == true)
            {
                printf("test_binary_common_envelope_evolution -- pre CE\n");
                print_system(&particlesMap,integration_flag);
            }
            
            Particle *star1 = particlesMap[0];
            double M1 = star1->mass;
            double M2 = particlesMap[1]->mass;

            int KW1 = star1->stellar_type;
            double M01 = star1->sse_initial_mass;
            double MC1 = star1->core_mass;
            double R1 = star1->radius;
            double R1_RSUN = star1->radius/CONST_R_SUN;

            double TM1,TN1;
            double *GB1,*TSCLS1,*LUMS1;
            GB1 = new double[10];
            TSCLS1 = new double[20];
            LUMS1 = new double[10];    
            double *ZPARS1 = star1->zpars;   
            
            star_(&KW1, &M01, &M1, &TM1, &TN1, TSCLS1, LUMS1, GB1, ZPARS1);

            double MENV1 = star1->convective_envelope_mass;
            
            double MENVD1 = MENV1 / (M1 - MC1);
            double RZAMS1 = rzamsf_(&M01);
            double fac = binary_evolution_CE_recombination_fraction;
            double L1 = star1->luminosity/CONST_L_SUN;
            double LAMB1 = celamf_(&KW1,&M01,&L1,&R1_RSUN,&RZAMS1,&MENVD1,&fac);
            
            if (verbose_testing == true)
            {
                printf("test_binary_common_envelope_evolution -- pre CE\n");
                print_system(&particlesMap,integration_flag);
            }
            
            double start_time = 0.0;
            double end_time = 1.0e3;
            double output_time,hamiltonian;
            int state,CVODE_flag,CVODE_error_code;

            double alpha=2.0;
            particlesMap[0]->common_envelope_alpha = alpha;
            particlesMap[0]->common_envelope_timescale = 10.0;
            int kw = particlesMap[0]->stellar_type;
            double t_old=0.0;
            double t = 0.0;
            double dt;
            double dt_stev;
            bool apply_SNe_effects;
            bool unbound_orbits;


            int N_binaries,N_root_finding,N_ODE_equations;

            int binary_index = 3;
            int index1 = 0;
            int index2 = 1;
            binary_common_envelope_evolution(&particlesMap, binary_index, index1, index2, t, &integration_flag);

            if (verbose_testing == true)
            {
                printf("test_binary_common_envelope_evolution -- post CE; pre evolve\n");
                print_system(&particlesMap,integration_flag);
            }
            
            evolve(&particlesMap,start_time,end_time,&output_time,&hamiltonian,&state,&CVODE_flag,&CVODE_error_code,&integration_flag);

            if (verbose_testing == true)
            {
                printf("test_binary_common_envelope_evolution --post evolve\n");
                print_system(&particlesMap,integration_flag);
            }
            
            double num = particlesMap[3]->a;
            
            double EBINDI = M1 * (M1 - MC1)/(LAMB1 * R1);
            double EORBI = MC1 * M2/(2.0 * (smas[0] ));
            double EORBF = EORBI + EBINDI/alpha;
            double a_f =  MC1 * M2/(2.0 * EORBF); // secondary does not have a core; assume it did not lose mass
            
            double tol = 1e-4;
            if ( !equal_number(a_f,num,tol) )
            {
                printf("test.cpp -- error in test_binary_common_envelope_evolution! %g %g\n",a_f,num);
                flag = 1;
            }

            
            clear_particles(&particlesMap);
        }
            
    }
        
    return 0;


}


/*********************
 * Triple interactions
**********************/

int test_triple_interactions()
{
    printf("test.cpp -- test_triple_interactions\n");
    
    int flag=0;
    flag += test_triple_common_envelope_evolution();
    
    if (flag == 0)
    {
        printf("test.cpp -- test_triple_interactions -- passed\n");
    }
    
    return flag;
}


int test_triple_common_envelope_evolution()
{
    printf("test.cpp -- test_triple_common_envelope_evolution\n");
    
    int flag = 0;
    
    bool verbose_testing = false;
    if (verbose_testing == true)
    {
        verbose_flag = 1;
    }
    
    for (int i=0; i<4; i++)
    {
        if (verbose_testing == true)
        {
            printf("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n");
        }
        if (i==0) // new outer orbit stable
        {
            ParticlesMap particlesMap;
            int N_bodies = 3;
            double masses[3] = {1.0,0.8,4.0};
            int stellar_types[3] = {1,1,5};
            int object_types[3] = {1,1,1};
            double smas[2] = {0.1,5.0};
            double es[2] = {0.01,0.01};
            double TAs[2] = {0.01,0.01};
            double INCLs[2] = {0.01,0.01};
            double APs[2] = {0.01,0.01};
            double LANs[2] = {0.01,0.01};
            random_seed=6;
            
            create_nested_system(particlesMap,N_bodies,masses,stellar_types,object_types,smas,es,TAs,INCLs,APs,LANs);// = create_nested_system();
            initialize_code(&particlesMap);

            
            int integration_flag = 0;
            if (verbose_testing == true)
            {
                printf("test_triple_common_envelope_evolution -- pre triple CE\n");
                print_system(&particlesMap,integration_flag);
            }

            double start_time = 0.0;
            double end_time = 1.0e-10;
            double output_time,hamiltonian;
            int state,CVODE_flag,CVODE_error_code;

            particlesMap[2]->triple_common_envelope_alpha = 100.0;
            int kw = particlesMap[0]->stellar_type;
            double t_old=0.0;
            double t = 0.0;
            double dt;
            double dt_stev;
            bool apply_SNe_effects;
            bool unbound_orbits;

            int N_binaries,N_root_finding,N_ODE_equations;

            int binary_index = 4;
            int index1 = 2;
            int index2 = 3;
            triple_common_envelope_evolution(&particlesMap, binary_index, index1, index2, t, &integration_flag);
            
            if (verbose_testing == true)
            {
                printf("test_triple_common_envelope_evolution -- post triple CE\n");
                print_system(&particlesMap,integration_flag);
                printf("test_triple_common_envelope_evolution -- pre evolve\n");
            }
            

            evolve(&particlesMap,start_time,end_time,&output_time,&hamiltonian,&state,&CVODE_flag,&CVODE_error_code,&integration_flag);

            if (verbose_testing == true)
            {

                printf("test_triple_common_envelope_evolution -- done evolve output_time %g\n",output_time);
                print_system(&particlesMap,integration_flag);
            }
            
            double a1_f = particlesMap[3]->a;
            double a2_f = particlesMap[4]->a;
            
            double tol = 1e-3;
            double num = smas[0];

            if ( !equal_number(a1_f,num,tol) )
            {
                printf("test.cpp -- error in test_triple_common_envelope_evolution! %g %g\n",a1_f,num);
                flag = 1;
            }

            num = 0.367243;

            if ( !equal_number(a2_f,num,tol) )
            {
                printf("test.cpp -- error in test_triple_common_envelope_evolution! %g %g\n",a2_f,num);
                flag = 1;
            }
            
            clear_particles(&particlesMap);
        }
        else if (i==1) // new outer orbit unstable
        {
            ParticlesMap particlesMap;
            int N_bodies = 3;
            double masses[3] = {1.0,0.8,4.0};
            int stellar_types[3] = {1,1,5};
            int object_types[3] = {1,1,1};
            double smas[2] = {0.1,5.0};
            double es[2] = {0.01,0.01};
            double TAs[2] = {0.01,0.01};
            double INCLs[2] = {0.01,0.01};
            double APs[2] = {0.01,0.01};
            double LANs[2] = {0.01,0.01};
            random_seed=6;
            
            create_nested_system(particlesMap,N_bodies,masses,stellar_types,object_types,smas,es,TAs,INCLs,APs,LANs);// = create_nested_system();
            initialize_code(&particlesMap);

            
            int integration_flag = 0;
            if (verbose_testing == true)
            {
                printf("test_triple_common_envelope_evolution -- pre triple CE\n");
                print_system(&particlesMap,integration_flag);
            }

            double start_time = 0.0;
            double end_time = 1.0e-10;
            double output_time,hamiltonian;
            int state,CVODE_flag,CVODE_error_code;

            particlesMap[2]->triple_common_envelope_alpha = 1.0;
            int kw = particlesMap[0]->stellar_type;
            double t_old=0.0;
            double t = 0.0;
            double dt;
            double dt_stev;
            bool apply_SNe_effects;
            bool unbound_orbits;

            int N_binaries,N_root_finding,N_ODE_equations;

            int binary_index = 4;
            int index1 = 2;
            int index2 = 3;
            triple_common_envelope_evolution(&particlesMap, binary_index, index1, index2, t, &integration_flag);
            
            if (verbose_testing == true)
            {
                printf("test_triple_common_envelope_evolution -- post triple CE\n");
                print_system(&particlesMap,integration_flag);
                printf("test_triple_common_envelope_evolution -- pre evolve\n");
            }
            

            evolve(&particlesMap,start_time,end_time,&output_time,&hamiltonian,&state,&CVODE_flag,&CVODE_error_code,&integration_flag);

            if (verbose_testing == true)
            {

                printf("test_triple_common_envelope_evolution -- done evolve output_time %g\n",output_time);
                print_system(&particlesMap,integration_flag);
            }
            
            double a1_f = particlesMap[3]->a;
            double a2_f = particlesMap[4]->a;
            
            double tol = 1e-4;
            double num = smas[0];

            if ( !equal_number(a1_f,num,tol) )
            {
                printf("test.cpp -- error in test_triple_common_envelope_evolution! %g %g\n",a1_f,num);
                flag = 1;
            }

            num = 0.329554;

            if ( !equal_number(a2_f,num,tol) )
            {
                printf("test.cpp -- error in test_triple_common_envelope_evolution! %g %g\n",a2_f,num);
                flag = 1;
            }
            
            clear_particles(&particlesMap);
        }

        else if (i==2) // with external (fourth) body; stable
        {
            ParticlesMap particlesMap;
            int N_bodies = 4;
            double masses[4] = {1.0,0.8,4.0,2.0};
            int stellar_types[4] = {1,1,5,1};
            int object_types[4] = {1,1,1,1};
            double smas[3] = {0.1,5.0,100.0};
            double es[3] = {0.01,0.01,0.1};
            double TAs[3] = {0.01,0.01,0.1};
            double INCLs[3] = {0.01,0.01,0.1};
            double APs[3] = {0.01,0.01,0.1};
            double LANs[3] = {0.01,0.01,0.1};
            random_seed=6;
            
            create_nested_system(particlesMap,N_bodies,masses,stellar_types,object_types,smas,es,TAs,INCLs,APs,LANs);// = create_nested_system();
            initialize_code(&particlesMap);

            
            int integration_flag = 0;
            if (verbose_testing == true)
            {
                printf("test_triple_common_envelope_evolution -- pre triple CE\n");
                print_system(&particlesMap,integration_flag);
            }

            double start_time = 0.0;
            double end_time = 1.0e-10;
            double output_time,hamiltonian;
            int state,CVODE_flag,CVODE_error_code;

            particlesMap[2]->triple_common_envelope_alpha = 100.0;
            int kw = particlesMap[0]->stellar_type;
            double t_old=0.0;
            double t = 0.0;
            double dt;
            double dt_stev;
            bool apply_SNe_effects;
            bool unbound_orbits;

            int N_binaries,N_root_finding,N_ODE_equations;

            int binary_index = 5;
            int index1 = 2;
            int index2 = 4;
            triple_common_envelope_evolution(&particlesMap, binary_index, index1, index2, t, &integration_flag);
            
            if (verbose_testing == true)
            {
                printf("test_triple_common_envelope_evolution -- post triple CE\n");
                print_system(&particlesMap,integration_flag);
                printf("test_triple_common_envelope_evolution -- pre evolve\n");
            }
            

            evolve(&particlesMap,start_time,end_time,&output_time,&hamiltonian,&state,&CVODE_flag,&CVODE_error_code,&integration_flag);

            if (verbose_testing == true)
            {

                printf("test_triple_common_envelope_evolution -- done evolve output_time %g\n",output_time);
                print_system(&particlesMap,integration_flag);
            }
            
            double a1_f = particlesMap[4]->a;
            double a2_f = particlesMap[5]->a;
            
            double tol = 1e-4;
            double num = smas[0];

            if ( !equal_number(a1_f,num,tol) )
            {
                printf("test.cpp -- error in test_triple_common_envelope_evolution! %g %g\n",a1_f,num);
                flag = 1;
            }

            num = 0.367243;

            if ( !equal_number(a2_f,num,tol) )
            {
                printf("test.cpp -- error in test_triple_common_envelope_evolution! %g %g\n",a2_f,num);
                flag = 1;
            }
            
            clear_particles(&particlesMap);
        }

        else if (i==3) // with external (fourth) body; unstable
        {
            ParticlesMap particlesMap;
            int N_bodies = 4;
            double masses[4] = {1.0,0.8,4.0,2.0};
            int stellar_types[4] = {1,1,5,1};
            int object_types[4] = {1,1,1,1};
            double smas[3] = {0.1,5.0,100.0};
            double es[3] = {0.01,0.01,0.1};
            double TAs[3] = {0.01,0.01,0.1};
            double INCLs[3] = {0.01,0.01,0.1};
            double APs[3] = {0.01,0.01,0.1};
            double LANs[3] = {0.01,0.01,0.1};
            random_seed=6;
            
            create_nested_system(particlesMap,N_bodies,masses,stellar_types,object_types,smas,es,TAs,INCLs,APs,LANs);// = create_nested_system();
            initialize_code(&particlesMap);

            
            int integration_flag = 0;
            if (verbose_testing == true)
            {
                printf("test_triple_common_envelope_evolution -- pre triple CE\n");
                print_system(&particlesMap,integration_flag);
            }

            double start_time = 0.0;
            double end_time = 1.0e-4;
            double output_time,hamiltonian;
            int state,CVODE_flag,CVODE_error_code;

            particlesMap[2]->triple_common_envelope_alpha = 1.0;
            int kw = particlesMap[0]->stellar_type;
            double t_old=0.0;
            double t = 0.0;
            double dt;
            double dt_stev;
            bool apply_SNe_effects;
            bool unbound_orbits;

            int N_binaries,N_root_finding,N_ODE_equations;

            int binary_index = 5;
            int index1 = 2;
            int index2 = 4;
            triple_common_envelope_evolution(&particlesMap, binary_index, index1, index2, t, &integration_flag);
            
            if (verbose_testing == true)
            {
                printf("test_triple_common_envelope_evolution -- post triple CE\n");
                print_system(&particlesMap,integration_flag);
                printf("test_triple_common_envelope_evolution -- pre evolve\n");
            }
            

            evolve(&particlesMap,start_time,end_time,&output_time,&hamiltonian,&state,&CVODE_flag,&CVODE_error_code,&integration_flag);

            if (verbose_testing == true)
            {

                printf("test_triple_common_envelope_evolution -- done evolve output_time %g\n",output_time);
                print_system(&particlesMap,integration_flag);
            }
            
            double a1_f = particlesMap[4]->a;
            double a2_f = particlesMap[5]->a;
            
            double tol = 1e-4;
            double num = smas[0];

            if ( !equal_number(a1_f,num,tol) )
            {
                printf("test.cpp -- error in test_triple_common_envelope_evolution! %g %g\n",a1_f,num);
                flag = 1;
            }

            num = 0.329554;

            if ( !equal_number(a2_f,num,tol) )
            {
                printf("test.cpp -- error in test_triple_common_envelope_evolution! %g %g\n",a2_f,num);
                flag = 1;
            }
            
            clear_particles(&particlesMap);
        }


        else if (i==4) // custom
        {
            ParticlesMap particlesMap;
            int N_bodies = 4;
            double masses[4] = {1.0,1.8,4.0,2.0};
            int stellar_types[4] = {1,1,5,1};
            int object_types[4] = {1,1,1,1};
            double smas[3] = {0.2,20.0,1000.0};
            double es[3] = {0.01,0.01,0.1};
            double TAs[3] = {0.01,0.01,0.1};
            double INCLs[3] = {0.01,0.01,0.1};
            double APs[3] = {0.01,0.01,0.1};
            double LANs[3] = {0.01,0.01,0.1};
            random_seed=6;
            
            create_nested_system(particlesMap,N_bodies,masses,stellar_types,object_types,smas,es,TAs,INCLs,APs,LANs);// = create_nested_system();
            initialize_code(&particlesMap);

            
            int integration_flag = 0;
            if (verbose_testing == true)
            {
                printf("test_triple_common_envelope_evolution -- pre triple CE\n");
                print_system(&particlesMap,integration_flag);
            }

            double start_time = 0.0;
            double end_time = 1.0e0;
            double output_time,hamiltonian;
            int state,CVODE_flag,CVODE_error_code;

            particlesMap[2]->triple_common_envelope_alpha = 10000.0;
            particlesMap[2]->common_envelope_timescale = 1000.0;
            int kw = particlesMap[0]->stellar_type;
            double t_old=0.0;
            double t = 0.0;
            double dt;
            double dt_stev;
            bool apply_SNe_effects;
            bool unbound_orbits;

            int N_binaries,N_root_finding,N_ODE_equations;

            int binary_index = 5;
            int index1 = 2;
            int index2 = 4;
            triple_common_envelope_evolution(&particlesMap, binary_index, index1, index2, t, &integration_flag);
            
            if (verbose_testing == true)
            {
                printf("test_triple_common_envelope_evolution -- post triple CE\n");
                print_system(&particlesMap,integration_flag);
                printf("test_triple_common_envelope_evolution -- pre evolve\n");
            }
            

            evolve(&particlesMap,start_time,end_time,&output_time,&hamiltonian,&state,&CVODE_flag,&CVODE_error_code,&integration_flag);

            if (verbose_testing == true)
            {

                printf("test_triple_common_envelope_evolution -- done evolve output_time %g\n",output_time);
                print_system(&particlesMap,integration_flag);
            }
            
            double a1_f = particlesMap[4]->a;
            double a2_f = particlesMap[5]->a;
            
            double tol = 1e-4;
            double num = smas[0];

            if ( !equal_number(a1_f,num,tol) )
            {
                printf("test.cpp -- error in test_triple_common_envelope_evolution! %g %g\n",a1_f,num);
                flag = 1;
            }

            num = 0.367243;

            if ( !equal_number(a2_f,num,tol) )
            {
                printf("test.cpp -- error in test_triple_common_envelope_evolution! %g %g\n",a2_f,num);
                flag = 1;
            }
            
            clear_particles(&particlesMap);
        }

            
    }
       
    
    return flag;


}

}
