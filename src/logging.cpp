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
    
    //ParticlesMap log_particlesMap;
    //copy_particlesMap(particlesMap,&log_particlesMap);
    
    //log_particlesMap.insert(particlesMap->begin(), particlesMap->end());
    
    ParticlesMap log_particlesMap = copy_particlesMap_for_logging(particlesMap);
    new_entry.particlesMap = log_particlesMap;

    new_entry.log_info = log_info;
    
    //printf("TEST %g\n",log_particlesMap[0]->mass);
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


ParticlesMap copy_particlesMap_for_logging(ParticlesMap *particlesMap)
{
    ParticlesMap log_particlesMap;

    int index;
    int j=0;
    bool is_binary;
    ParticlesMapIterator it_p;
    for (it_p = particlesMap->begin(); it_p != particlesMap->end(); it_p++)
    {
        Particle *p = (*it_p).second;
        //Particle *q;
        Particle *q = new Particle(p->index, p->is_binary);
        //printf("? %d\n",p->index);
        //q->index = p->index;
        //q->is_binary = p->is_binary;
        q->mass = p->mass;
        q->parent = p->parent;
        
        if (q->is_binary == false)
        {
            q->radius = p->radius;
            q->stellar_type = p->stellar_type;
            q->core_mass = p->core_mass;
            q->sse_initial_mass = p->sse_initial_mass;
            q->convective_envelope_mass = p->convective_envelope_mass;
            q->epoch = p->epoch;
            q->age = p->age;
            q->core_radius = p->core_radius;
            q->convective_envelope_radius = p->convective_envelope_radius;
            q->luminosity = p->luminosity;
            
            for (int k=0; k<3; k++)
            {
                q->spin_vec[k] = p->spin_vec[k];
            }
        }
        else
        {
            q->child1 = p->child1;
            q->child2 = p->child2;
            q->true_anomaly = p->true_anomaly;
        
            for (int k=0; k<3; k++)
            {
                q->h_vec[k] = p->h_vec[k];
                q->e_vec[k] = p->e_vec[k];
            }
        }
        
        log_particlesMap[q->index] = q;

    }

    return log_particlesMap;
}

}
