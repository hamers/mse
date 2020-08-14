/* MSE */

#include "evolve.h"
#include "logging.h"

extern "C"
{

void update_log_data(ParticlesMap *particlesMap, double time, int integration_flag, int event_flag, Log_info_type log_info)
/* Event flags *
 * 0: initial state
 * 1: stellar type change
 * 2: SNe
 * 3: onset of mass transfer
 * 4: stop mass transfer
 * 5: CE
 * 6: collision
 * 7: dynamical instability
 * 8: secular breakdown
 */
{
    Log_type new_entry;
    new_entry.index = logData.size();
    new_entry.time = time;
    new_entry.event_flag = event_flag;
    new_entry.integration_flag = integration_flag;
    
    ParticlesMap log_particlesMap;
    copy_particlesMap(particlesMap,&log_particlesMap);
    new_entry.particlesMap = log_particlesMap;

    new_entry.log_info = log_info;
    //if (event_flag == 1)
    //{
     //   printf("ST CHANGE %d\n",log_info.index1);
    //}
    //new_entry.N_particles = particlesMap->size();
    //if (integration_flag == 0)
    //{
        //int N_bodies,N_binaries,N_root_finding,N_ODE_equations;
        //determine_binary_parents_and_levels(particlesMap,&N_bodies,&N_binaries,&N_root_finding,&N_ODE_equations);
        //new_entry.N_bodies = N_bodies;
        //new_entry.N_binaries = N_binaries;
        //printf("T %d %d\n",N_bodies,N_binaries);
    //}
    
    logData.push_back(new_entry);
    
    //printf("logging.cpp -- event_flag %d\n",event_flag,logData.size());
}

}
