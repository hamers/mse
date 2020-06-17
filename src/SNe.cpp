/* MSE */

#include "evolve.h"
#include "SNe.h"

extern "C"
{

int handle_SNe_in_system(ParticlesMap *particlesMap, bool *unbound_orbits, int *integration_flag)
{
    int flag;
    double VX,VY,VZ;
    ParticlesMapIterator it_p;
    //std::vector<int>::iterator it_parent_p,it_parent_q;

//    int seed = orbital_phases_random_seed;
    int index=0;
    
    for (it_p = particlesMap->begin(); it_p != particlesMap->end(); it_p++)
    {
        Particle *p = (*it_p).second;
        if (p->is_binary == false and p->evolve_as_star == true)
        {
            /* p's instantaneous_perturbation_delta_mass is assumed to be set before calling handle_SNe_in_system() */
            p->instantaneous_perturbation_delta_X = 0.0;
            p->instantaneous_perturbation_delta_Y = 0.0;
            p->instantaneous_perturbation_delta_Z = 0.0;
            p->instantaneous_perturbation_delta_VX = 0.0;
            p->instantaneous_perturbation_delta_VY = 0.0;
            p->instantaneous_perturbation_delta_VZ = 0.0;

            //if (fabs(p->instantaneous_perturbation_delta_mass) > 0.0)
            if (p->apply_kick == true)
            {
                flag = sample_kick_velocity(p,&VX,&VY,&VZ);
                p->instantaneous_perturbation_delta_VX = VX;
                p->instantaneous_perturbation_delta_VY = VY;
                p->instantaneous_perturbation_delta_VZ = VZ;
                index+=1;
            }
            printf("SNe.cpp -- delta_m %g vk %g % g %g\n",p->instantaneous_perturbation_delta_mass,VX,VY,VZ);
            
        }
    }
    if (*integration_flag == 0) /* secular case */
    {
        apply_user_specified_instantaneous_perturbation(particlesMap);
    }
    else
    {
        apply_user_specified_instantaneous_perturbation_nbody(particlesMap);
    }
    
    *unbound_orbits = check_for_unbound_orbits(particlesMap);

    if (*unbound_orbits == true)
    {
        printf("SNe.cpp -- handle_SNe_in_system -- Unbound orbits in system due to supernova!\n");
        *integration_flag = 3;
    }

    reset_instantaneous_perturbation_quantities(particlesMap);
            
    return 0;
}

int sample_kick_velocity(Particle *p, double *vx, double *vy, double *vz)
{
    //srand(seed);
    double x;
    x = ((double) rand() / (RAND_MAX));
    double theta = 2.0*M_PI*x - M_PI;
    x = ((double) rand() / (RAND_MAX));
    double phi = 2.0*M_PI*x;
    *vx = sin(theta)*cos(phi);
    *vy = sin(theta)*sin(phi);
    *vz = cos(theta);

    double vnorm;
    
    if (p->kick_distribution == 0) // no kicks
    {
        vnorm = 0.0;
    }
    else if (p->kick_distribution == 1) // Maxwellian for NS; 0 otherwise
    {
        if (p->stellar_type == 13)
        {
            double sigma = p->kick_distribution_sigma;
            double v[3];
            sample_from_3d_maxwellian_distribution(sigma, v);
            vnorm = norm3(v);
        }
    }
    
  //    std::default_random_engine generator (seed);

//  std::normal_distribution<double> distribution (0.0,1.0);

//    vnorm = 0.0;
    printf("SNe.cpp -- i %d apply_kick %d distr %d sigma %g vnorm %g\n",p->index,p->apply_kick,p->kick_distribution,p->kick_distribution_sigma,vnorm);
    *vx *= vnorm;
    *vy *= vnorm;
    *vz *= vnorm;
    
    return 0;
}

bool check_for_unbound_orbits(ParticlesMap *particlesMap)
{
    ParticlesMapIterator it_p;
//    double h_vec[3],e_vec[3];
    double e;
    
    bool unbound_orbits = false;
    
    for (it_p = particlesMap->begin(); it_p != particlesMap->end(); it_p++)
    {
        Particle *p = (*it_p).second;
        if (p->is_binary == 1)
        {
            //get_e_and_h_vectors_from_particle(p,e_vec,h_vec);
            e = norm3(p->e_vec);

            if (e<0 or e >= 1.0)
            {
                unbound_orbits = true;
            }
            //printf("test e %.15f unbound_orbits %d\n",e,unbound_orbits);
        }
    }
    
    return unbound_orbits;
}

}
