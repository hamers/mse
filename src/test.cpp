/* MSE */
/* Adrian Hamers November 2019 */

#include "evolve.h"
#include "merger.h"

extern "C"
{

int test_collisions()
{
    
    int flag;
    
    //flag = test_collision_MS_MS();
    //flag = test_collision_giant_MS();
    //flag = test_collision_star_MS(10.0,1,5.0);
    //flag = test_collision_star_MS(40.0,2,5.0);
    
    flag = test_collision_stars(5.01,1,5);
    //flag = test_collision_stars(5.01,4,2.0);
    //flag = test_collision_stars(5.01,6,3.0);
    //flag = test_collision_stars(10.0,4,15);
    //flag = test_collision_stars(22.0,4,21.98);
    //flag = test_collision_stars(22.0,14,8.0);
    //flag = test_collision_stars(50.0,14,22.0);
    
    return flag;
}

int test_collision_stars(double m1, int kw1, double m2)
{
    printf("*************************************\n");
    printf("test_collision_star_MS m1 %g kw1 %d\n",m1,kw1);
    printf("*************************************\n");
    
    ParticlesMap particlesMap;
    int N_bodies = 4;
    double masses[4] = {m1,m2,1.0,1.0};
    double smas[3] = {1.0,100.0,10000.0};
    double es[3] = {0.01,0.01,0.01};
    double TAs[3] = {0.01,0.01,0.01};
    double INCLs[3] = {0.01,0.01,0.01};
    double APs[3] = {0.01,0.01,0.01};
    double LANs[3] = {0.01,0.01,0.01};
    
    create_nested_system(particlesMap,N_bodies,masses,smas,es,TAs,INCLs,APs,LANs);// = create_nested_system();
//    printf("post s %d b %d %g r %g\n",particlesMap2.size(),particlesMap2[0]->is_binary,particlesMap2[0]->mass,particlesMap2[0]->radius);

    printf("pre evolve\n");
    print_system(&particlesMap);
    
    int kw = particlesMap[0]->stellar_type;
    double t_old=0.0;
    double t = 0.0;
    double dt;
    double dt_stev;
    bool apply_SNe_effects;
    bool unbound_orbits;

    int integration_flag = 0;

    int N_binaries,N_root_finding,N_ODE_equations;
            
    while (kw < kw1)
    {
        dt = t-t_old;
        dt_stev=0.0;
        evolve_stars(&particlesMap,t_old,t,&dt_stev,false,&apply_SNe_effects);
        kw = particlesMap[0]->stellar_type;
        if (dt != 0.0)
        {
            update_stellar_evolution_quantities_directly(&particlesMap,dt);
        }
        if (apply_SNe_effects == true)
        {
            
            handle_SNe_in_system(&particlesMap,&unbound_orbits,&integration_flag);
            
            /* Or: do not actually take into account effect of mass loss on orbit; just update the masses */
            ParticlesMapIterator it;
            for (it = particlesMap.begin(); it != particlesMap.end(); it++)
            {
                Particle *body = (*it).second;
                if (body->is_binary == false)
                {
                    //body->mass += body->instantaneous_perturbation_delta_mass;
                    //body->instantaneous_perturbation_delta_mass = 0.0;
                }
            }
        }
        t_old = t;
        
        //dt_stev *= 0.4;
        t += dt_stev;
        //t += 1000.0;
        
       // printf("t %g kw %d\n",t,kw);
    }
    update_structure(&particlesMap);
        
    printf("post evolve\n");
    print_system(&particlesMap);
    particlesMap[4]->merged = true;


    
    handle_collisions(&particlesMap,&integration_flag);
    //printf("pre s %d\n",particlesMap.size());
    //collision_product(&particlesMap, 4, 0, 1, &integration_flag);
    //printf("post s %d b %d %g r %g\n",particlesMap.size(),particlesMap[2]->is_binary,particlesMap[2]->mass,particlesMap[2]->radius);

    //printf("post1\n");
    //print_system(&particlesMap);



    printf("post merge integration_flag %d\n",integration_flag);
    print_system(&particlesMap);

    
    int flag;
    
    
    return 0;
}



int test_collision_MS_MS()
{
    ParticlesMap particlesMap;
    int N_bodies = 4;
    double masses[4] = {10.0,5.0,1.0,1.0};
    double smas[3] = {1.0,100.0,10000.0};
    double es[3] = {0.01,0.01,0.01};
    double TAs[3] = {0.01,0.01,0.01};
    double INCLs[3] = {0.01,0.01,0.01};
    double APs[3] = {0.01,0.01,0.01};
    double LANs[3] = {0.01,0.01,0.01};
    
    create_nested_system(particlesMap,N_bodies,masses,smas,es,TAs,INCLs,APs,LANs);// = create_nested_system();
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


}
