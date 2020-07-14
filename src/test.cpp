/* MSE */

#include "evolve.h"
#include "test.h"

extern "C"
{

int test_tools()
{
   
    int flag = 0;
    flag += test_a_P_orb_conversion();
    flag += test_a_h_conversion();
    flag += test_orbital_element_conversion();
    flag += test_kepler_equation_solver();
    flag += test_orbital_vectors_cartesian_conversion();
    //flag += test_kroupa_imf_sampling();
    
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

#ifdef IGNORE
int test_kroupa_imf_sampling()
{
    printf("test.cpp -- test_kroupa_imf_sampling\n");
    
    int flag = 0;
    
    //double m = sample_from_Kroupa_93_imf();
    printf("m %g\n",m);
    //printf("K C1 %g C2 %g C3 %g x1 %g x2 %g x3 %g\n",kroupa_C1,kroupa_C2,kroupa_C3,kroupa_x1,kroupa_x2,kroupa_x3);
    
    
    return flag;
}
#endif

int test_binary_evolution()
{
    printf("test.cpp -- test_binary_evolution\n");
    
    int flag=0;
    flag += test_compute_Kelvin_Helmholtz_timescale();
    flag += test_compute_Eddington_accretion_rate();
    flag += test_handle_instantaneous_and_adiabatic_mass_changes_in_orbit();
    
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
    
    //print_system(&particlesMap);
    
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
            //flag = test_collision_stars(10.0,i,13,j,0);
            random_seed = 0;
            flag = test_collision_stars(10.0,i,13,j,1);
            
        }
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
        particlesMap[0]->Collision_Partner = 1;
        particlesMap[1]->Collision_Partner = 0;
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


#ifdef IGNORE
int test_collision_MS_MS()
{
    ParticlesMap particlesMap;
    int N_bodies = 4;
    double masses[4] = {10.0,5.0,1.0,1.0};
    int stellar_types[4] = {1,1,1,1};
    double smas[3] = {1.0,100.0,10000.0};
    double es[3] = {0.01,0.01,0.01};
    double TAs[3] = {0.01,0.01,0.01};
    double INCLs[3] = {0.01,0.01,0.01};
    double APs[3] = {0.01,0.01,0.01};
    double LANs[3] = {0.01,0.01,0.01};
    
    create_nested_system(particlesMap,N_bodies,masses,stellar_types,smas,es,TAs,INCLs,APs,LANs);// = create_nested_system();
//    printf("post s %d b %d %g r %g\n",particlesMap2.size(),particlesMap2[0]->is_binary,particlesMap2[0]->mass,particlesMap2[0]->radius);

    particlesMap[4]->merged = true;

    printf("pre\n");
    print_system(&particlesMap);
    
    int integration_flag = 0;
    handle_collisions(&particlesMap,&integration_flag);
    //printf("pre s %d\n",particlesMap.size());
    //collision_product(&particlesMap, 4, 0, 1, &integration_flag);
    //printf("post s %d b %d %g r %g\n",particlesMap.size(),particlesMap[2]->is_binary,particlesMap[2]->mass,particlesMap[2]->radius);

    printf("post1\n");
    print_system(&particlesMap);

    int N_binaries,N_root_finding,N_ODE_equations;

    determine_binary_parents_and_levels(&particlesMap, &N_bodies, &N_binaries, &N_root_finding,&N_ODE_equations);

    printf("post2\n");
    print_system(&particlesMap);

    
    int flag;
    
    
    return 0;
}
#endif


    
}
