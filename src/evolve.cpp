/* MSE */

#include "evolve.h"

extern "C"
{

int initialize_code(ParticlesMap *particlesMap)
{
    #ifdef DEBUG
    printf("evolve.cpp -- initialize_code; set_up_flybys and initialize stars\n");
    #endif

    update_structure(particlesMap, 0);

    srand(random_seed);
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
    
    return 0;
}

int evolve(ParticlesMap *particlesMap, double start_time, double end_time, double *output_time, double *hamiltonian, int *state, int *CVODE_flag, int *CVODE_error_code, int *integration_flag)
{
    int N_bodies,N_binaries,N_root_finding,N_ODE_equations;
    //if (*integration_flag == 0)
    {
        //determine_binary_parents_and_levels(particlesMap,&N_bodies,&N_binaries,&N_root_finding,&N_ODE_equations);
        //printf("beginning of evolve start_time %g  end_time  %g %d %d %d %d\n",start_time,end_time,N_bodies,N_binaries,N_root_finding,N_ODE_equations);
        //print_system(particlesMap,*integration_flag);
    }

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
    int flag;
    double t = start_time;
    double t_old = t;
    double t_end = end_time;
    double dt,dt_stev,dt_binary_evolution,dt_nbody;
    double t_out;
    //double t_next_encounter = 0.0;

    bool apply_SNe_effects = false;
    bool unbound_orbits = false;

    /* Initial timestep */
    double min_dt = 1.0e0;
    double max_dt = end_time - start_time;

    dt = CV_min(min_dt,max_dt);
    dt_binary_evolution = max_dt;
    
    /* Time loop */
    //printf("evolve.cpp -- start integration_flag %d\n",*integration_flag);

    int i=0;
    *state = 0;
    bool last_iteration = false;
    while (t <= t_end)
    {
        t += dt;

        /* Stellar evolution */
        //printf("pre ev\n");
        //print_system(particlesMap,*integration_flag);
        flag = evolve_stars(particlesMap,t_old,t,&dt_stev,false,&apply_SNe_effects);
        //printf("post ev\n");
        //print_system(particlesMap,*integration_flag);

        /* SNe */
        if (apply_SNe_effects == true)
        {
            flag = handle_SNe_in_system(particlesMap,&unbound_orbits,integration_flag);
        }
        //printf("post SNe\n");
        
        /* Binary evolution */
        flag = handle_binary_evolution(particlesMap,t_old,t,&dt_binary_evolution,integration_flag);
        //printf("post bin\n");
        //print_system(particlesMap,*integration_flag);

        /* Dynamical evolution -- will update masses/radii */
        if (*integration_flag == 0) // Secular
        {
            flag = integrate_ODE_system(particlesMap,t_old,t,&t_out,hamiltonian,CVODE_flag,CVODE_error_code);
        }
        else // Direct N-body
        {
            integrate_nbody_system(particlesMap, integration_flag, t_old, t, &t_out, &dt_nbody);
        }

        //printf("post dyn\n");
        //print_system(particlesMap,*integration_flag);
        
        t = t_out;
        
        #ifdef DEBUG
        printf("ODE dt %g t_old - t %g t - t_out %g\n",dt,t_old-t,t-t_out);
        #endif
        
        /* Handle roots */
        if (*CVODE_flag==2)
        {
            printf("ROOT occurred; setting t = t_out; t %g t_out %g\n",t,t_out);
            print_system(particlesMap,*integration_flag);
            
            flag = investigate_roots_in_system(particlesMap, t, *integration_flag);
            handle_roots(particlesMap, flag, integration_flag, CVODE_flag, t, &dt_stev, &dt_binary_evolution);
        }
       
        /* Time step (phase 1) */
        dt = CV_min(dt_stev,dt_binary_evolution);
        
        if (*integration_flag > 0)
        {
            dt = CV_min(dt,dt_nbody);
        }
        
        /* Flybys */
        if (include_flybys == true)
        {
            bool apply_flyby = flyby_criterion(particlesMap, integration_flag);
            if (t + dt >= flybys_t_next_encounter and last_iteration == false and apply_flyby == true)// and *integration_flag==0)
            {
                //printf("NE t+dt %g flybys_t_next_encounter %g \n",t+dt,flybys_t_next_encounter);
                //printf("flybys dt %g %g\n",dt,flybys_t_next_encounter - t);
                dt = flybys_t_next_encounter - t;
                
                handle_next_flyby(particlesMap, false, integration_flag);
            }
        }

        /* Time step (phase 2) */
        if (t + dt >= end_time)
        {
            dt = end_time - t;
            last_iteration = true;
            //printf("adjust dt to reach end time \n");
        }

        dt = CV_min(dt,max_dt);
        dt = CV_max(dt,min_dt);

        t_old = t;
        
        if (t >= t_end)
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
            if (last_entry.event_flag == LOG_CE_START and *integration_flag == 0) // end of CE phase
            {
                Log_info_type log_info;
                log_info.binary_index = last_entry_info.binary_index;
                log_info.index1 = last_entry_info.index1;
                log_info.index2 = last_entry_info.index2;
                update_log_data(particlesMap, t, *integration_flag, LOG_CE_END, log_info);
            }
            if (last_entry.event_flag == LOG_COL_START and *integration_flag == 0) // end of collision phase
            {
                Log_info_type log_info;
                log_info.binary_index = last_entry_info.binary_index;
                log_info.index1 = last_entry_info.index1;
                log_info.index2 = last_entry_info.index2;
                update_log_data(particlesMap, t, *integration_flag, LOG_COL_END, log_info);
            }
        }
        #endif
        
        //*integration_flag=1;

    }

  //if (*integration_flag == 0)
    {
        //int N_bodies_new,N_binaries_new,N_root_finding_new,N_ODE_equations_new;
        //determine_binary_parents_and_levels(particlesMap,&N_bodies_new,&N_binaries_new,&N_root_finding_new,&N_ODE_equations_new);
        //if (N_binaries != N_binaries_new and *integration_flag == 0)
        {
        //    printf("Restructuring of system! %d %d\n",N_binaries,N_binaries_new);
            //*state = 1;
        }
    }

    //printf("end of evolve\n");
    //print_system(particlesMap,*integration_flag);
    
    //update_log_data(particlesMap, t, *integration_flag, 0);
    
    *output_time = t;
    *hamiltonian = 0.0;
    //*output_flag = 0;
    //*error_code = error_code;
    
    return 0;
}

}
