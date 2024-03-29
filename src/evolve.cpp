/* MSE */

#include "evolve.h"

extern "C"
{

int initialize_code(ParticlesMap *particlesMap)
{

    #ifdef VERBOSE
    if (verbose_flag > 0)
    {
        printf("evolve.cpp -- initialize_code; set_up_flybys and initialize stars\n");
        print_system(particlesMap,0);
    }
    #endif

    update_structure(particlesMap, 0);

    random_number_generator.seed(random_seed);
    flybys_t_next_encounter = 0.0;
    
    error_code = 0;

    int integration_flag = 0;

    initialize_stars(particlesMap);
    set_positions_and_velocities(particlesMap);

    if (include_flybys == true)
    {
        handle_next_flyby(particlesMap,true,&integration_flag);
    }
 
    #ifdef LOGGING
    Log_info_type log_info;
    update_log_data(particlesMap, 0.0, 0, LOG_INIT, log_info);
    #endif

    signal(SIGSEGV, custom_segfault_handler); // set up custom handler for segfaults
    
    time(&wall_time_start);

    return 0;
}

void custom_segfault_handler(int s)
{
    switch(s)
    {
        case SIGSEGV:
        error_code = -1;
        printf("tools.cpp -- Segmentation fault signal caught! System index %d \n",system_index);
        
        FILE *fp;
        fp = fopen("segfault.txt", "w+");
        fprintf(fp,"Segmentation fault signal caught! System index %d \n",system_index);
        fclose(fp);
        
        longjmp(jump_buf,1);
        break;
    }
}

int evolve(ParticlesMap *particlesMap, double start_time, double end_time, double *output_time, double *hamiltonian, int *state, int *CVODE_flag, int *CVODE_error_code, int *integration_flag)
{

    /* Error handling -- begin */
    int setjmp_return_val = setjmp(jump_buf);
    if (setjmp_return_val == 1) // error occurred somewhere in the code; will return control to Python
    {
        *output_time = start_time;
        *hamiltonian = 0.0;
        *state = -1;
        *CVODE_flag = -1;
        *CVODE_error_code = -1;
        *integration_flag = 0;
        return error_code;
    }

    /* Some lines to test segfault handling (should normally NEVER be enabled) */
    //if (start_time > 1e4)
    //{
        //int *p = NULL;
        //*p=0xdead; // DELIBERATE segfault
    //}

    /* Error handling -- end */
    
    int N_bodies,N_binaries,N_root_finding,N_ODE_equations;

    #ifdef VERBOSE
    if (verbose_flag > 1)
    {
        printf("evolve.cpp -- evolve -- start \n");
        print_system(particlesMap,*integration_flag);
    }
    #endif

    if (start_time == end_time)
    {
        *output_time = start_time;
        *hamiltonian = 0.0;
        *state = 0.0;
        *CVODE_flag = 0;
        *CVODE_error_code = 0;

        return 0;
    }
    
    /* Set up integration variables */
    int flag,binary_flag;
    double t = start_time;
    double t_old = t;
    double t_end = end_time;
    double dt,dt_stev,dt_binary_evolution,dt_nbody;
    double t_out;

    bool apply_SNe_effects = false;
    bool unbound_orbits = false;

    /* Initial timestep */
    double min_dt = ODE_min_dt;
    double max_dt = end_time - start_time;

    dt = CV_min(min_dt,max_dt);
    dt_binary_evolution = max_dt;

    /* Time loop */

    int i=0;
    *state = 0;
    bool last_iteration = false;
    while (t <= t_end)
    {
        t += dt;

        /* Stellar evolution */
        flag = evolve_stars(particlesMap,t_old,t,&dt_stev,false,&apply_SNe_effects,integration_flag);

        /* SNe */
        if (apply_SNe_effects == true)
        {
            flag = handle_SNe_in_system(particlesMap,&unbound_orbits,integration_flag);
        }
        
        /* Binary evolution */
        binary_flag = handle_binary_evolution(particlesMap,t_old,t,&dt_binary_evolution,integration_flag);

        /* Dynamical evolution -- will also update masses, radii, and spins */
        if (*integration_flag == 0) // Secular
        {
            flag = integrate_ODE_system(particlesMap,t_old,t,&t_out,hamiltonian,CVODE_flag,CVODE_error_code);
        }
        else // Direct N-body
        {
            integrate_nbody_system(particlesMap, integration_flag, t_old, t, &t_out, &dt_nbody);
        }

        t = t_out;

        #ifdef LOGGING
        Log_type &last_entry = logData.back();
        if (last_entry.event_flag == LOG_CE_START)
        {
            Log_info_type &last_log_info = last_entry.log_info;
            Log_info_type log_info;
            
            log_info.binary_index = last_log_info.binary_index;
            log_info.index1 = last_log_info.index1;
            log_info.index2 = last_log_info.index2;
            update_log_data(particlesMap, t, *integration_flag, LOG_CE_END, log_info);
        }
        if (last_entry.event_flag == LOG_COL_START)
        {
            Log_info_type &last_log_info = last_entry.log_info;
            Log_info_type log_info;

            log_info.binary_index = last_log_info.binary_index;
            log_info.index1 = last_log_info.index1;
            log_info.index2 = last_log_info.index2;
            update_log_data(particlesMap, t, *integration_flag, LOG_COL_END, log_info);
        }
        #endif

        #ifdef VERBOSE
        if (verbose_flag > 1)
        {
            printf("evolve.cpp -- evolve -- t %g ODE dt %g t_old - t %g t - t_out %g\n",t,dt,t_old-t,t-t_out);
        }
        #endif

        /* Handle roots */
        if (*CVODE_flag==2)
        {
            
            if (stop_after_root_found == true) // Used for root finding implementation testing only
            {
                *CVODE_flag = 2;
                break;
            }

            #ifdef VERBOSE
            if (verbose_flag > 0)
            {
                printf("evolve.cpp -- evolve -- ROOT occurred; setting t = t_out; t %g t_out %g\n",t,t_out);
                print_system(particlesMap,*integration_flag);
            }
            #endif
            
            flag = investigate_roots_in_system(particlesMap, t, *integration_flag);
            handle_roots(particlesMap, flag, integration_flag, CVODE_flag, t, &dt_stev, &dt_binary_evolution);
            
        }
       
        /* Time step (phase 1) */
        dt = CV_min(dt_stev,dt_binary_evolution);

        #ifdef VERBOSE
        if (verbose_flag > 1)
        {
            printf("evolve.cpp -- evolve -- t %g dt_stev %g dt_binary_evolution %g\n",t,dt_stev,dt_binary_evolution);
            print_system(particlesMap,*integration_flag);
        }
        #endif

        if (*integration_flag > 0)
        {
            dt = CV_min(dt,dt_nbody);
        }

        if (*integration_flag == 2) // semisecular; use stellar evolution timesteps
        {
            dt = dt_stev;
        }

        /* Flybys */
        if (include_flybys == true)
        {
            bool apply_flyby = flyby_criterion(particlesMap, integration_flag);
            if (t + dt >= flybys_t_next_encounter and last_iteration == false and apply_flyby == true)
            {
                dt = flybys_t_next_encounter - t;
                
                handle_next_flyby(particlesMap, false, integration_flag);
            }
        }

        /* Time step (phase 2) */
        if (t + dt >= end_time)
        {
            //printf("AD t %g dt %g t+dt %g end_time %g\n",t,dt,t+dt,end_time);
            dt = end_time - t;
            last_iteration = true;
            
            #ifdef VERBOSE
            if (verbose_flag > 1)
            {
                printf("evolve.cpp -- evolve -- adjust dt to reach end time; t %g end_time %g new dt %g\n",t,end_time,dt);
                print_system(particlesMap,*integration_flag);
            }
            #endif
            
        }

        dt = CV_min(dt,max_dt);
        dt = CV_max(dt,min_dt);

        t_old = t;
//        printf("dt %g dt_stev %g dt_binary_evolution %g\n",dt,dt_stev,dt_binary_evolution);
        if (t >= end_time)
        {
            break;
        }
        i+=1;
        
        
        #ifdef LOGGING
        if (logData.size() > 0)
        {
            Log_type last_entry = logData.back();
            Log_info_type last_entry_info = last_entry.log_info;
            if (last_entry.event_flag == LOG_SNE_START and *integration_flag == 0) // end of SNe phase
            {
                Log_info_type log_info;
                log_info.index1 = last_entry_info.index1;
                update_log_data(particlesMap, t, *integration_flag, LOG_SNE_END, log_info);
            }
            if (last_entry.event_flag == LOG_WD_KICK_START and *integration_flag == 0) // end of WD kick phase
            {
                Log_info_type log_info;
                log_info.index1 = last_entry_info.index1;
                update_log_data(particlesMap, t, *integration_flag, LOG_WD_KICK_END, log_info);
            }
        }
        #endif
        
        if (error_code != 0)
        {
            break;
        }
        
    }

    *output_time = t;
    *hamiltonian = 0.0;

    #ifdef VERBOSE
    if (verbose_flag > 1)
    {
        printf("evolve.cpp -- evolve -- end \n");
        print_system(particlesMap,*integration_flag);
    }
    #endif


    return error_code;
}

}
