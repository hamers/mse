/* MSE */

#include "evolve.h"
#include "logging.h"

extern "C"
{

void update_log_data(ParticlesMap *particlesMap, double time, int integration_flag, int event_flag)
{
    Log_type new_entry;
    new_entry.index = logData.size();
    new_entry.time = time;
    new_entry.event_flag = event_flag;

    ParticlesMap log_particlesMap;
    copy_particlesMap(particlesMap,&log_particlesMap);
    new_entry.particlesMap = log_particlesMap;

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
    
    printf("logging.cpp -- event_flag %d N_log %d %g\n",event_flag,logData.size(),logData[0].particlesMap[0]->mass);
}

}
