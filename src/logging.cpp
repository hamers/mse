/* MSE */

#include "evolve.h"
#include "logging.h"

extern "C"
{

void update_log_data(ParticlesMap *particlesMap, double time, int integration_flag, int event_flag, Log_info_type log_info)
/* Event flags *
 * 0: initial state
 * 1: stellar type change
 * 2: SNe start
 * 3: SNe end
 * 4: onset of mass transfer
 * 5: stop mass transfer
 * 6: CE start
 * 7: CE end
 * 8: collision start
 * 9: collision end
 * 10: dynamical instability
 * 11: secular breakdown
 * 12: WD kick start
 * 13: WD kick end
 * 14: Triple CE start
 * 15: Triple CE end
 * 16: MSP formation
 * 17: final state (written when write_final_log_entry is called within Python)
 */
{
    Log_type new_entry;
    new_entry.index = logData.size();
    new_entry.time = time;
    new_entry.event_flag = event_flag;
    new_entry.integration_flag = integration_flag;
    
    ParticlesMap log_particlesMap = copy_particlesMap_for_logging(particlesMap);
    new_entry.particlesMap = log_particlesMap;

    new_entry.log_info = log_info;
    
    logData.push_back(new_entry);
    
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
        Particle *q = new Particle(p->index, p->is_binary);

        q->mass = p->mass;
        q->parent = p->parent;
        
        if (q->is_binary == false)
        {
            q->radius = p->radius;
            q->stellar_type = p->stellar_type;
            q->object_type = p->object_type;

            q->metallicity = p->metallicity;            
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
                q->R_vec[k] = p->R_vec[k];
                q->V_vec[k] = p->V_vec[k];
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
                q->R_vec[k] = p->R_vec[k];
                q->V_vec[k] = p->V_vec[k];
            }
        }
        
        log_particlesMap[q->index] = q;

    }

    return log_particlesMap;
}

}
