#include "types.h"
extern "C"
{
void update_log_data(ParticlesMap *particlesMap, double time, int integration_flag, int event_flag, Log_info_type log_info);
ParticlesMap copy_particlesMap_for_logging(ParticlesMap *particlesMap);
}
