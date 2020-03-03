#include "types.h"
#include "evolve.h"
#include "ODE_system.h"

extern "C"
{
int compute_y_dot(realtype time, N_Vector y, N_Vector y_dot, void *data_)
{
	UserData data;
	data = (UserData) data_;
    ParticlesMap *particlesMap = data->particlesMap;
    
    double start_time = data->start_time;
    double delta_time = time - start_time;

    
    extract_ODE_variables(particlesMap, y, delta_time);
    reset_ODE_dots(particlesMap, y, delta_time);
    
    /****************************
     * compute right-hand sides *
     * **************************/

    double hamiltonian = 0.0;
    ParticlesMapIterator it_p,it_f;
    std::vector<int>::iterator it_parent_p,it_parent_q;
    
    

    for (it_p = particlesMap->begin(); it_p != particlesMap->end(); it_p++)
    {
        Particle *p = (*it_p).second;
        if (p->is_binary == 1)
        {
            Particle *child1 = (*particlesMap)[p->child1];
            Particle *child2 = (*particlesMap)[p->child2];

            /* Newtonian gravitational point mass dynamics */
            
            /* binary pairs */
            for (it_parent_p = p->parents.begin(); it_parent_p != p->parents.end(); it_parent_p++)
            {
                int i = std::distance(p->parents.begin(), it_parent_p);
                Particle *P_q = (*particlesMap)[(*it_parent_p)];
                int connecting_child_in_parent_q = p->connecting_child_in_parents[i];
                hamiltonian += compute_EOM_binary_pairs(particlesMap,p->index,P_q->index,connecting_child_in_parent_q,false);
                
                /* binary triplets */
                for (it_parent_q = P_q->parents.begin(); it_parent_q != P_q->parents.end(); it_parent_q++)
                {
                    int j = std::distance(P_q->parents.begin(), it_parent_q);
                    Particle *P_u = (*particlesMap)[(*it_parent_q)];
                    int connecting_child_in_parent_u = P_q->connecting_child_in_parents[j];
                    hamiltonian += compute_EOM_binary_triplets(particlesMap,p->index,P_q->index,P_u->index,connecting_child_in_parent_q,connecting_child_in_parent_u,false);
                    //printf("cross applied %d %d %d %d %d\n",P_p->index,P_q->index,P_u->index,connecting_child_in_parent_q,connecting_child_in_parent_u);

                }
            }
            
            
            /* Perturbations by flybys */
            for (it_f = particlesMap->begin(); it_f != particlesMap->end(); it_f++)
            {
                Particle *f = (*it_f).second;
                //if (f->mode == 0)
                if (f->is_external == true)
                {
                    hamiltonian += compute_EOM_binary_pairs_external_perturbation(particlesMap,p->index,f->index,time,false);
                }
            }

            /* VRR-related perturbations */

            if (p->VRR_model > 0) /* 0: no VRR perturbation */
            {
                hamiltonian += compute_VRR_perturbations(particlesMap,p->index,time);
            }


            /* Pairwise PN corrections */
            if (p->include_pairwise_1PN_terms == 1)
            {
                hamiltonian += compute_EOM_pairwise_1PN(particlesMap,p->index,false);
            }
            if (p->include_pairwise_25PN_terms == 1)
            {
                hamiltonian += compute_EOM_pairwise_25PN(particlesMap,p->index,false);
            }
            
            /* Tidal friction (ad hoc) */
            if (child1->include_tidal_friction_terms == 1 || child1->include_tidal_bulges_precession_terms == 1 || child1->include_rotation_precession_terms == 1)
            {
                if (child1->tides_method == 0 || child1->tides_method == 1)
                {
                    compute_EOM_equilibrium_tide_BO_full(particlesMap,p->index,child1->index,child2->index,child1->include_tidal_friction_terms,child1->include_tidal_bulges_precession_terms,child1->include_rotation_precession_terms,child1->minimum_eccentricity_for_tidal_precession,child1->tides_method);
                }
                else if (child1->tides_method == 2)
                {
                    compute_EOM_equilibrium_tide(particlesMap,p->index,child1->index,child2->index,child1->include_tidal_friction_terms,child1->include_tidal_bulges_precession_terms,child1->include_rotation_precession_terms,child1->minimum_eccentricity_for_tidal_precession);
                }
            }
            if (child2->include_tidal_friction_terms == 1 || child2->include_tidal_bulges_precession_terms == 1 || child2->include_rotation_precession_terms == 1)
            {
                if (child2->tides_method == 0 || child2->tides_method == 1)
                {
                    compute_EOM_equilibrium_tide_BO_full(particlesMap,p->index,child2->index,child1->index,child2->include_tidal_friction_terms,child2->include_tidal_bulges_precession_terms,child2->include_rotation_precession_terms,child2->minimum_eccentricity_for_tidal_precession,child2->tides_method);
                }
                else if (child2->tides_method == 2)
                {
                    compute_EOM_equilibrium_tide(particlesMap,p->index,child2->index,child1->index,child2->include_tidal_friction_terms,child2->include_tidal_bulges_precession_terms,child2->include_rotation_precession_terms,child2->minimum_eccentricity_for_tidal_precession);
                }

            }
            

            /* Orbital changes due to mass changes */
            
            ODE_handle_stellar_winds(p);
            ODE_handle_RLOF(p,child1,child2);
            
        }
            
    }
    
    write_ODE_variables_dots(particlesMap,y_dot);

    data->hamiltonian = hamiltonian;

    return 0;
}


void extract_ODE_variables(ParticlesMap *particlesMap, N_Vector &y, double delta_time)
{
    ParticlesMapIterator it_p;
    int k=1;
    int k_component;

//    double spin_vec[3],e_vec[3],h_vec[3];
//    double mass,radius;
    
    for (it_p = particlesMap->begin(); it_p != particlesMap->end(); it_p++)
    {
        Particle *p = (*it_p).second;
        if (p->is_binary == 0) // particle is a body
        {
            for (k_component=0; k_component<3; k_component++)
            {
                p->spin_vec[k_component] = Ith(y,k + k_component);
            }
            p->mass = Ith(y,k + 2 + 1);
            p->radius = Ith(y,k + 2 + 2);
            //printf("m %g index %d dt %g\n",p->mass,p->index,delta_time);
            k=k+5;
        }
        if (p->is_binary == 1) // particle is a binary
        {
            for (k_component=0; k_component<3; k_component++)
            {
                p->e_vec[k_component] = Ith(y,k + k_component);
                p->h_vec[k_component] = Ith(y,k + k_component + 3);
            }
            //printf("testttt %g %g %g\n",P_p->h_vec_x_dot_external,P_p->h_vec_y_dot_external,P_p->h_vec_z_dot_external);
            k=k+6;
        }
//        P_p->set_ODE_quantities(delta_time);
        
//        if (reset_ODE_quantities == true)
        //{
          //  P_p->reset_ODE_quantities();
        //}        
//        printf("p %d a %g\n",P_p->index,P_p->a);
//        printf("particle %d mass %g e_vec_x %g h_vec_x %g\n",P_p->index,P_p->mass,P_p->e_vec_x,P_p->h_vec_x);
    }
    
    set_binary_masses_from_body_masses(particlesMap);

    for (it_p = particlesMap->begin(); it_p != particlesMap->end(); it_p++)
    {
        Particle *p = (*it_p).second;
        if (p->is_binary == 0) // particle is a body
        {
            p->spin_vec_norm = norm3(p->spin_vec);
        }
        else
        {
            p->e = norm3(p->e_vec);
            p->h = norm3(p->h_vec);

            for (int i=0; i<3; i++)
            {        
                p->e_vec_unit[i] = p->e_vec[i]/p->e;
                p->h_vec_unit[i] = p->h_vec[i]/p->h;
            }
            cross3(p->h_vec_unit,p->e_vec_unit,p->q_vec_unit);
            
            p->e_p2 = p->e*p->e;
            p->j_p2 = 1.0 - p->e_p2;
            p->j = sqrt(p->j_p2);
            p->j_p3 = p->j*p->j_p2;
            p->j_p4 = p->j*p->j_p3;
            p->j_p5 = p->j*p->j_p4;

            p->a = p->h*p->h*p->child1_mass_plus_child2_mass/( CONST_G*p->child1_mass_times_child2_mass*p->child1_mass_times_child2_mass*p->j_p2 );
        }
    }
}

void reset_ODE_dots(ParticlesMap *particlesMap, N_Vector &y, double delta_time)
{
    ParticlesMapIterator it_p;
    int i;
    double spin_vec_norm;
//    double spin_vec[3],e_vec[3],h_vec[3];
//    double mass,radius;
    
    for (it_p = particlesMap->begin(); it_p != particlesMap->end(); it_p++)
    {
        Particle *p = (*it_p).second;
        if (p->is_binary == 0) // particle is a body
        {

            /* Externally imposed */
            spin_vec_norm = norm3(p->spin_vec);
            //printf("spin_vec_norm %g\n",spin_vec_norm);
            for (i=0; i<3; i++)
            {
                p->dspin_vec_dt[i] = p->ospin_dot * (p->spin_vec[i]/spin_vec_norm);
                //p->dspin_vec_dt[i] = 0.0;
            }

            p->dmass_dt = p->mass_dot_wind; //+ p->mass_dot_RLOF;
    //printf("dmass_dt %g %g\n",dmass_dt,radius_dot);
            p->dradius_dt = p->radius_dot + p->radius_ddot*delta_time;
        }
        else
        {
            for (i=0; i<3; i++)
            {
                p->de_vec_dt[i] = 0.0;
                p->dh_vec_dt[i] = 0.0;
            }
        }
    }
}

void write_ODE_variables_dots(ParticlesMap *particlesMap, N_Vector &y_dot)
{
    ParticlesMapIterator it_p;
    int k=1;
    int k_component;

    double spin_vec[3],e_vec[3],h_vec[3];
    for (it_p = particlesMap->begin(); it_p != particlesMap->end(); it_p++)
    {
        Particle *P_p = (*it_p).second;
        if (P_p->is_binary == 0) // particle is a body
        {
            for (k_component=0; k_component<3; k_component++)
            {
                Ith(y_dot,k + k_component) = P_p->dspin_vec_dt[k_component];
            }
            
            Ith(y_dot,k + 2 + 1) = P_p->dmass_dt;
            Ith(y_dot,k + 2 + 2) = P_p->dradius_dt;
            //printf("ODE dR/dt %g\n",P_p->dradius_dt);
            
            k=k+5;
        }
        if (P_p->is_binary == 1) // particle is a binary
        {
            for (k_component=0; k_component<3; k_component++)
            {
                Ith(y_dot,k + k_component)      = P_p->de_vec_dt[k_component];
                Ith(y_dot,k + k_component + 3)  = P_p->dh_vec_dt[k_component];
                //printf("k %d out %g\n",k,Ith(y_dot,k + k_component));
                //printf("out %g\n",Ith(y_dot,k + k_component+3));
                
            }

            //printf("P_p->child1_mass_dot_external %g P_p->child2_mass_dot_external %g\n",P_p->child1_mass_dot_external,P_p->child2_mass_dot_external);
            //printf("out %g\n",2.0*dot3(P_p->h_vec_unit,P_p->dh_vec_dt)/norm3(P_p->h_vec) + (P_p->child1_mass_dot_external + P_p->child2_mass_dot_external)/P_p->child1_mass_plus_child2_mass \
                    - 2.0*(P_p->child1_mass_dot_external/P_p->child1_mass + P_p->child2_mass_dot_external/P_p->child2_mass) );
            
                        
            k=k+6;
        }
    }
}

void set_initial_ODE_variables(ParticlesMap *particlesMap, N_Vector &y, N_Vector &y_abs_tol, double abs_tol_spin_vec, double abs_tol_e_vec, double abs_tol_h_vec)
{
    ParticlesMapIterator it_p;
    int k=1;
    int k_component;

    double spin_vec[3],e_vec[3],h_vec[3];
    for (it_p = particlesMap->begin(); it_p != particlesMap->end(); it_p++)
    {
        Particle *P_p = (*it_p).second;
        if (P_p->is_binary == 0) // particle is a body
        {
            spin_vec[0] = P_p->spin_vec_x;
            spin_vec[1] = P_p->spin_vec_y;
            spin_vec[2] = P_p->spin_vec_z;

            for (k_component=0; k_component<3; k_component++)
            {
                Ith(y,         k + k_component) = spin_vec[k_component];
                Ith(y_abs_tol, k + k_component) = abs_tol_spin_vec;
            }
            Ith(y,  k + 2 + 1) = P_p->mass;
            Ith(y,  k + 2 + 2) = P_p->radius;
            
            Ith(y_abs_tol, k + 2 + 1) = relative_tolerance*P_p->mass;
            Ith(y_abs_tol, k + 2 + 2) = relative_tolerance*P_p->radius;
            //Ith(y_abs_tol, k + 2 + 2) = 1.0e-10;

            k=k+5;
        }
        if (P_p->is_binary == 1) // particle is a binary
        {
            e_vec[0] = P_p->e_vec_x;
            e_vec[1] = P_p->e_vec_y;
            e_vec[2] = P_p->e_vec_z;
            h_vec[0] = P_p->h_vec_x;
            h_vec[1] = P_p->h_vec_y;
            h_vec[2] = P_p->h_vec_z;
            double h = norm3(h_vec);
    
            for (k_component=0; k_component<3; k_component++)
            {
                Ith(y,k + k_component) = e_vec[k_component];
                Ith(y,k + k_component + 3) = h_vec[k_component];
                
                Ith(y_abs_tol,k + k_component) = abs_tol_e_vec;
//                Ith(y_abs_tol,k + k_component + 3) = abs_tol_h_vec;
                Ith(y_abs_tol,k + k_component + 3) = relative_tolerance*h;
            }
            k=k+6;
        }
    }
}

void extract_final_ODE_variables(ParticlesMap *particlesMap, N_Vector &y_out)
{
    ParticlesMapIterator it_p;
    int k=1;
    int k_component;

    double spin_vec[3],e_vec[3],h_vec[3];
    for (it_p = particlesMap->begin(); it_p != particlesMap->end(); it_p++)
    {
        Particle *P_p = (*it_p).second;
        if (P_p->is_binary == 0) // particle is a body
        {
            for (k_component=0; k_component<3; k_component++)
            {
                spin_vec[k_component] = Ith(y_out,k + k_component);
            }
            P_p->spin_vec_x = spin_vec[0];
            P_p->spin_vec_y = spin_vec[1];
            P_p->spin_vec_z = spin_vec[2];

            P_p->mass = Ith(y_out,k + 2 + 1);
            P_p->radius = Ith(y_out,k + 2 + 2);
            //printf("actual R %g m %g\n",P_p->radius,P_p->mass);

            k=k+5;
        }
        if (P_p->is_binary == 1) // particle is a binary
        {
            for (k_component=0; k_component<3; k_component++)
            {
                e_vec[k_component] = Ith(y_out,k + k_component);
                h_vec[k_component] = Ith(y_out,k + k_component + 3);
            }

            P_p->e_vec_x = e_vec[0];
            P_p->e_vec_y = e_vec[1];
            P_p->e_vec_z = e_vec[2];
            P_p->h_vec_x = h_vec[0];
            P_p->h_vec_y = h_vec[1];
            P_p->h_vec_z = h_vec[2];

            k=k+6;
        }
    }
}
}
