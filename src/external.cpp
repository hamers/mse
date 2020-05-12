/*
*/

#include "types.h"
#include "external.h"
#include "evolve.h"
//#include "../interface.h" /* for parameters */
#include "structure.h" /* for determine_binary_parents_and_levels */
#include <stdio.h>

extern "C"
{

void apply_instantaneous_mass_changes_and_kicks(ParticlesMap *particlesMap, int *integration_flag)
{

    apply_user_specified_instantaneous_perturbation(particlesMap);
    reset_instantaneous_perturbation_quantities(particlesMap);
    
    bool unbound_orbits = check_for_unbound_orbits(particlesMap);
    if (unbound_orbits == true)
    {
        *integration_flag = 5;
        printf("external.cpp -- apply_instantaneous_mass_changes_and_kicks -- Unbound orbits in system; switching to integration_flag %d\n",*integration_flag);
    }
}

void apply_user_specified_instantaneous_perturbation(ParticlesMap *particlesMap)
{
    set_positions_and_velocities(particlesMap);
    update_masses_positions_and_velocities_of_all_bodies_instantaneous_perturbation(particlesMap);
    update_masses_positions_and_velocities_of_all_binaries(particlesMap);
    update_orbital_vectors_in_binaries_from_positions_and_velocities(particlesMap);
}

void apply_user_specified_instantaneous_perturbation_nbody(ParticlesMap *particlesMap)
{
    update_masses_positions_and_velocities_of_all_bodies_instantaneous_perturbation(particlesMap);
    set_binary_masses_from_body_masses(particlesMap);
}


void reset_instantaneous_perturbation_quantities(ParticlesMap *particlesMap)
{
    
    ParticlesMapIterator it;
    for (it = particlesMap->begin(); it != particlesMap->end(); it++)
    {
        Particle *p = (*it).second;
        {
            p->instantaneous_perturbation_delta_mass = 0.0;
            p->instantaneous_perturbation_delta_X = 0.0;
            p->instantaneous_perturbation_delta_Y = 0.0;
            p->instantaneous_perturbation_delta_Z = 0.0;

            p->instantaneous_perturbation_delta_VX = 0.0;
            p->instantaneous_perturbation_delta_VY = 0.0;
            p->instantaneous_perturbation_delta_VZ = 0.0;
        }
    }
}

void update_masses_positions_and_velocities_of_all_bodies_instantaneous_perturbation(ParticlesMap *particlesMap)
{
    /* Update the masses, positions and velocities of the bodies */
    
    ParticlesMapIterator it;
    for (it = particlesMap->begin(); it != particlesMap->end(); it++)
    {
        Particle *body = (*it).second;
        if (body->is_binary == false)
        {
            body->mass += body->instantaneous_perturbation_delta_mass;
            body->R_vec[0] += body->instantaneous_perturbation_delta_X;
            body->R_vec[1] += body->instantaneous_perturbation_delta_Y;
            body->R_vec[2] += body->instantaneous_perturbation_delta_Z;

            body->V_vec[0] += body->instantaneous_perturbation_delta_VX;
            body->V_vec[1] += body->instantaneous_perturbation_delta_VY;
            body->V_vec[2] += body->instantaneous_perturbation_delta_VZ;
            
            //printf("test instantaneous_perturbation_delta_mass body %d delta m %g\n",body->index,body->instantaneous_perturbation_delta_mass);
          //printf("x %g y %g z %g \n",body->position_X,body->position_Y,body->position_Z);
            //printf("vx %g vy %g vz %g \n",body->instantaneous_perturbation_delta_VX,body->instantaneous_perturbation_delta_VY,body->instantaneous_perturbation_delta_VZ);
        }
    }
}

void compute_EOM_Newtonian_external_perturbations(double time, ParticlesMap *particlesMap, Particle *p, double *hamiltonian, double *KS_V, bool compute_hamiltonian_only)
{
    ParticlesMapIterator it_f;
    for (it_f = particlesMap->begin(); it_f != particlesMap->end(); it_f++)
    {
        Particle *f = (*it_f).second;
        if (f->is_external == true)
        {
            compute_EOM_binary_pairs_external_perturbation(particlesMap,p->index,f->index,time,hamiltonian,KS_V,compute_hamiltonian_only);
        }
    }
}

void compute_EOM_binary_pairs_external_perturbation(ParticlesMap *particlesMap, int binary_index, int perturber_index, double time, double *hamiltonian, double *KS_V, bool compute_hamiltonian_only)
{
    /* last checked 20-03-20 */

    /* stop if no pairwise terms are to be computed */
    if ((include_quadrupole_order_terms == false) && (include_octupole_order_binary_pair_terms == false) && (include_hexadecupole_order_binary_pair_terms == false) && (include_dotriacontupole_order_binary_pair_terms == false) )
    {
        return;
    }

    Particle *binary = (*particlesMap)[binary_index];
    Particle *perturber = (*particlesMap)[perturber_index];

    #ifdef DEBUG
    printf("external.cpp -- compute_EOM_binary_pairs_external_perturbation -- binary_index %d perturber_index %d time %g\n",binary_index,perturber_index,time);
    #endif
    
    double e = binary->e;
    double e_P2 = binary->e_p2;
    
    double *e_vec = binary->e_vec;
    double *h_vec = binary->h_vec;
    
    double *e_vec_unit = binary->e_vec_unit;
    double *h_vec_unit = binary->h_vec_unit;
    
    double h = binary->h;
    double j = binary->j;
        
    double j_vec[3];
    for (int i=0; i<3; i++)
    {
        j_vec[i] = j*h_vec_unit[i];
    }
    
    double a = binary->a;
    
    double t_per_ref = perturber->external_t_ref;
    double dt = time - t_per_ref;
    double M_per = perturber->mass;
    //printf("test t %.15g %.15g %.15g\n",t_per_ref,time,dt);
    //printf("test %g %g %g %g %g %g %g %g\n",perturber->r0_vec_x,perturber->r0_vec_y,perturber->r0_vec_z,perturber->rdot_vec_x,perturber->rdot_vec_y,perturber->rdot_vec_z,perturber->t_ref,perturber->mass);
    
    double r_per_vec[3];
    double r_per;
    
    compute_position_vectors_external_particle(particlesMap,perturber,time,&r_per,r_per_vec);
    
    double r_per_p2 = r_per*r_per;
    double r_per_pm1 = 1.0/r_per;
    
    //printf("test r %.15g %.15g %.15g %.15g\n",r_per,r_per_vec[0],r_per_vec[1],r_per_vec[2]);
    double e_vec_dot_r_per_vec = dot3(e_vec,r_per_vec);
    double j_vec_dot_r_per_vec = dot3(j_vec,r_per_vec);
    
    double m1 = binary->child1_mass;
    double m2 = binary->child2_mass;

    //printf("m1 %g m2 %g m3 %g\n",m1,m2,M_per);
//    printf("a_in %g a_out %g\n",a_in,a_out);
        
    double m1_plus_m2 = binary->child1_mass_plus_child2_mass;
    double m1_minus_m2 = binary->child1_mass_minus_child2_mass;
    double m1_times_m2 = binary->child1_mass_times_child2_mass;

    int n,m,i1,i2; /* n: order of expansion, starts at n = 2; 0 <= m <= n; i1 + i2 <= m */
    double constant_Hamiltonian_factor = (m1_times_m2/m1_plus_m2)*(CONST_G*M_per/r_per); /* constant factor in Hamiltonian*/
    
    double M_bin_pnm1; /* M_bin to power n - 1 */
    double M_bin_child1_pnm1; /* M_bin_child1 to power n - 1 */
    double M_bin_child2_pnm1; /* M_bin_child2 to power n - 1 */
    double minusone_pnp1; /* minus one to power n + 1 */
    double r_per_pmn; /* r_per to power 1/n */
    double a_pn; /* binary a to power n */
    
    double mass_factor_children = 0.0;
    double binary_pair_hamiltonian = 0.0;
    double hamiltonian_factor = 0.0;
    double A_n_m = 0.0;
    double B_n_m_i1_i2 = 0.0;
    double dB_n_m_i1_i2_de = 0.0; /* derivative of B-function w.r.t. e */
    

    int index_B_eccp2;
    double e_Peven,e_Podd; // e to a power which is even or odd
    double e_P3 = e * e_P2;
    double e_P4 = e * e_P3;
    double e_P5 = e * e_P4;
    double e_P6 = e * e_P5;
    
    double grad_e_vec_H[3],grad_j_vec_H[3];
    for (int i=0; i<3; i++)
    {
        grad_e_vec_H[i] = grad_j_vec_H[i] = 0.0;
    }
    
    int index_A,index_B;
    int n_lookup,m_lookup;
    int n_old = 2;
    
    bool continue_to_B_table;
    double B_lookup = 0.0;
    double r_per_pow_mi1mi2,e_vec_dot_r_per_vec_pi1,j_vec_dot_r_per_vec_pi2,e_vec_dot_r_per_vec_pi1m1,j_vec_dot_r_per_vec_pi2m1;
    for (index_B=0; index_B<TABLELENGTH_B; index_B++)
    {
        B_n_m_i1_i2 = 0.0;
        dB_n_m_i1_i2_de = 0.0;
        
        n = int(B_TABLE[index_B][0]);
        m = int(B_TABLE[index_B][1]);
        i1 = int(B_TABLE[index_B][2]);
        i2 = int(B_TABLE[index_B][3]);
        
        /* include only user-specific orders */
        if (n==2 && include_quadrupole_order_terms==false) continue;
        if (n==3 && include_octupole_order_binary_pair_terms==false) continue;
        if (n==4 && include_hexadecupole_order_binary_pair_terms==false) continue;
        if (n==5 && include_dotriacontupole_order_binary_pair_terms==false) continue;
        
        /* retrieve A from table */
        continue_to_B_table = false;
        for (index_A=0; index_A<TABLELENGTH_A; index_A++)
        {
            n_lookup = int(A_TABLE[index_A][0]);
            m_lookup = int(A_TABLE[index_A][1]);
            
            if ((n_lookup == n) && (m_lookup == m))
            {
                A_n_m = A_TABLE[index_A][2];
                continue_to_B_table = true;
            }
        }     
        //printf("look %d %d %d %d\n",n,m,i1,i2);
        
        /* continue in B table if there is no entry for given n and m (i.e. A_n_m = 0) */   
        if (continue_to_B_table == false) continue;
        
        /* construct B-function from table */
        B_n_m_i1_i2 = 0.0;
        dB_n_m_i1_i2_de = 0.0;
        for (index_B_eccp2=0; index_B_eccp2<= HIGHEST_POWER_ECCP2_IN_B_TABLE; index_B_eccp2++)
        {
            B_lookup = B_TABLE[index_B][4+index_B_eccp2]; /* take into account offset for n,m,i1,i2 */

            if (index_B_eccp2 == 0)
            {
                e_Peven = 1.0;
                e_Podd = 0.0;
            }
            else if (index_B_eccp2 == 1)
            {
                e_Peven = e_P2;
                e_Podd = e;
            }
            else if (index_B_eccp2 == 2)
            {
                e_Peven = e_P4;
                e_Podd = e_P3;
            }
            else if (index_B_eccp2 == 3)
            {
                e_Peven = e_P6;
                e_Podd = e_P5;
            }
            else
            {
                printf("external.cpp -- FATAL ERROR in constructing B-function \n");
                exit(-1);
            }
            
            B_n_m_i1_i2 += B_lookup * e_Peven;
            dB_n_m_i1_i2_de += 2.0* ( (double) index_B_eccp2) * B_lookup * e_Podd;
        }
        //printf("test... n %d m %d i1 %d i2 %d e %g B %g \n",n,m,i1,i2,e,B_n_m_i1_i2);
        //printf("test D... n %d m %d i1 %d i2 %d e %g B %g \n",n,m,i1,i2,e,dB_n_m_i1_i2_de);
        
        M_bin_pnm1 = pow(m1_plus_m2,n-1.0);
        M_bin_child1_pnm1 = pow(m1,n-1.0);
        M_bin_child2_pnm1 = pow(m2,n-1.0);
        a_pn = pow(a,n);
        r_per_pmn = pow(r_per_pm1,n);    
        minusone_pnp1 = pow(-1.0,n+1.0);
    
        //printf("test2 %d %d %d %d %g %g %g\n",n,m,i1,i2,A_n_m,B_n_m_i1_i2,dB_n_m_i1_i2_de);
        
        /* compute the Hamiltonian */
        r_per_pow_mi1mi2 = pow(r_per,-i1-i2);
        e_vec_dot_r_per_vec_pi1 = pow(e_vec_dot_r_per_vec,i1);
        j_vec_dot_r_per_vec_pi2 = pow(j_vec_dot_r_per_vec,i2);

        e_vec_dot_r_per_vec_pi1m1 = pow(e_vec_dot_r_per_vec,i1-1.0);
        j_vec_dot_r_per_vec_pi2m1 = pow(j_vec_dot_r_per_vec,i2-1.0);

        mass_factor_children = M_bin_child2_pnm1 - minusone_pnp1 * M_bin_child1_pnm1;
        hamiltonian_factor = minusone_pnp1*constant_Hamiltonian_factor*(mass_factor_children/M_bin_pnm1)*a_pn*r_per_pmn*A_n_m*r_per_pow_mi1mi2;
        
        binary_pair_hamiltonian += hamiltonian_factor*e_vec_dot_r_per_vec_pi1*j_vec_dot_r_per_vec_pi2*B_n_m_i1_i2;

        /* compute EOM */
        if (compute_hamiltonian_only == false)
        {
            for (int i=0; i<3; i++)
            {
                grad_e_vec_H[i] += hamiltonian_factor*j_vec_dot_r_per_vec_pi2*( double(i1)*B_n_m_i1_i2*e_vec_dot_r_per_vec_pi1m1*r_per_vec[i] + e_vec_dot_r_per_vec_pi1*dB_n_m_i1_i2_de*e_vec_unit[i] );
                grad_j_vec_H[i] += hamiltonian_factor*e_vec_dot_r_per_vec_pi1*B_n_m_i1_i2*double(i2)*j_vec_dot_r_per_vec_pi2m1*r_per_vec[i];
            }
        }
    }
    
    
    double j_vec_cross_grad_j_vec_H[3],j_vec_cross_grad_e_vec_H[3];
    double e_vec_cross_grad_e_vec_H[3],e_vec_cross_grad_j_vec_H[3];
    
    cross3(j_vec,        grad_j_vec_H,              j_vec_cross_grad_j_vec_H);
    cross3(j_vec,        grad_e_vec_H,              j_vec_cross_grad_e_vec_H);
    cross3(e_vec,        grad_e_vec_H,              e_vec_cross_grad_e_vec_H);
    cross3(e_vec,        grad_j_vec_H,              e_vec_cross_grad_j_vec_H);

    double Lambda = h/j;
    for (int i=0; i<3; i++)
    {
        binary->de_vec_dt[i] += (-1.0/(Lambda))*( e_vec_cross_grad_j_vec_H[i] \
            + j_vec_cross_grad_e_vec_H[i] );
        binary->dh_vec_dt[i] += -1.0*( j_vec_cross_grad_j_vec_H[i] \
            + e_vec_cross_grad_e_vec_H[i] );
    }

    *hamiltonian = binary_pair_hamiltonian;
}

void compute_position_vectors_external_particle(ParticlesMap *particlesMap, Particle *perturber, double time, double *r_per, double r_per_vec[3])
{
    /* last checked 20-03-20 */
    
    double dt = time - perturber->external_t_ref;

    double e_f = perturber->external_e;
    double rp_f = perturber->external_r_p;
    double M_f = perturber->mass;
    double abs_a_f = rp_f/(e_f-1.0);
    double total_internal_system_mass = (*particlesMap)[0]->total_system_mass;
    double n_f = sqrt(CONST_G*(total_internal_system_mass + M_f)/(abs_a_f*abs_a_f*abs_a_f));
    double mean_anomaly = n_f*dt;
    double cos_true_anomaly,sin_true_anomaly;
    compute_true_anomaly_from_mean_anomaly_hyperbolic_orbit(mean_anomaly,e_f,&cos_true_anomaly,&sin_true_anomaly);

    *r_per = rp_f*(1.0 + e_f)/(1.0 + e_f*cos_true_anomaly);            
        
    double e_f_hat_vec[3] = {perturber->e_vec[0],perturber->e_vec[1],perturber->e_vec[2]}; // e_vec for external particles is assumed to be a unit vector!
    double h_f_hat_vec[3] = {perturber->h_vec[0],perturber->h_vec[1],perturber->h_vec[2]}; // h_vec for external particles is assumed to be a unit vector!
    double q_f_hat_vec[3];
    cross3(h_f_hat_vec,e_f_hat_vec,q_f_hat_vec);

    #ifdef DEBUG
    printf("external.cpp -- compute_position_vectors_external_particle -- r_per %g mean_anomaly %g total_internal_system_mass %g\n",r_per,mean_anomaly,total_internal_system_mass);
    #endif

    for (int i=0; i<3; i++)
    {        
        r_per_vec[i] = *r_per*(cos_true_anomaly*e_f_hat_vec[i] + sin_true_anomaly*q_f_hat_vec[i]);
        perturber->r_vec[i] = r_per_vec[i];
        
        #ifdef DEBUG
        printf("external.cpp -- compute_position_vectors_external_particle -- r_vec %g\n",r_per_vec[i]);
        #endif
    }  
}


int apply_external_perturbation_assuming_integrated_orbits(ParticlesMap *particlesMap)
{
    ParticlesMapIterator it_p,it_f;
    
    /* make sure that total_system_mass is updated */
    int N_particles = particlesMap->size();
    int N_bodies, N_binaries;
    int N_root_finding;
    int N_ODE_equations;
    
    determine_binary_parents_and_levels(particlesMap,&N_bodies,&N_binaries,&N_root_finding,&N_ODE_equations);
    set_binary_masses_from_body_masses(particlesMap);
    set_up_derived_ODE_quantities(particlesMap);
    
    /* compute and apply perturbations */
    for (it_p = particlesMap->begin(); it_p != particlesMap->end(); it_p++)
    {
        Particle *p = (*it_p).second;
        
        if (p->is_binary == true)
        {
            for (it_f = particlesMap->begin(); it_f != particlesMap->end(); it_f++)
            {
                Particle *f = (*it_f).second;
                if (f->is_external == true)
                {
                    {
                        apply_external_perturbation_assuming_integrated_orbits_binary_pair(particlesMap,p->index,f->index);
                    }
                }
            }
        }
    }
    
    return 0;
}

void apply_external_perturbation_assuming_integrated_orbits_binary_pair(ParticlesMap *particlesMap, int binary_index, int perturber_index)
{
    /* last checked 20-03-20 */
    
    /* stop if no pairwise terms are to be computed */
    if ((include_quadrupole_order_terms == false) && (include_octupole_order_binary_pair_terms == false) && (include_hexadecupole_order_binary_pair_terms == false) && (include_dotriacontupole_order_binary_pair_terms == false) )
    {
        return;
    }

    Particle *binary = (*particlesMap)[binary_index];
    Particle *perturber = (*particlesMap)[perturber_index];

    /* unsubscripted elements/vectors refer to the binary */
    double e = binary->e;
    double e_p2 = binary->e_p2;
    
    double *e_vec = binary->e_vec;
    double *h_vec = binary->h_vec;
    
    double *e_vec_unit = binary->e_vec_unit;
    double *h_vec_unit = binary->h_vec_unit;
    
    double h = binary->h;
    double j = binary->j;
    double a = binary->a;
    
    double j_vec[3];
    for (int i=0; i<3; i++)
    {
        j_vec[i] = j*h_vec_unit[i];
    }

    double M_per = perturber->mass;
    
    double grad_e_vec_H[3],grad_j_vec_H[3];
    for (int i=0; i<3; i++)
    {
        grad_e_vec_H[i] = grad_j_vec_H[i] = 0.0;
    }

    double e_f = perturber->external_e;
    double q_f = perturber->external_r_p;
    double abs_a_f = q_f/(e_f-1.0);
    double total_internal_system_mass = binary->total_system_mass; /* make sure this is up to date */
    double n_f = sqrt(CONST_G*(total_internal_system_mass + M_per)/(abs_a_f*abs_a_f*abs_a_f));

    /* note: e & h vectors of external particles are assumed to be unit vectors */
    double e_f_hat_vec[3] = {perturber->e_vec[0],perturber->e_vec[1],perturber->e_vec[2]};
    double j_f_hat_vec[3] = {perturber->h_vec[0],perturber->h_vec[1],perturber->h_vec[2]}; 

    double e_f_p2 = e_f*e_f;
    double e_f_p4 = e_f_p2*e_f_p2;
    double one_div_e_f_p1 = 1.0/e_f;
    double one_div_e_f_p2 = one_div_e_f_p1*one_div_e_f_p1;
    double one_div_e_f_p3 = one_div_e_f_p1*one_div_e_f_p2;
    double e_f_p2_minus_one = e_f_p2 - 1.0;
    double sqrt_e_f_p2_minus_one = sqrt(e_f_p2_minus_one);
    double e_f_p2_minus_one_p3div2 = e_f_p2_minus_one*sqrt_e_f_p2_minus_one;
    double asec_minus_e_f = acos(-1.0/e_f); /* asec(x) = acos(1/x) */

    double one_plus_e_f = 1.0 + e_f;
    double one_plus_e_f_pm1 = 1.0/one_plus_e_f;
    double one_plus_e_f_pmn;

    
    double e_vec_dot_e_f_hat_vec = dot3(e_vec,e_f_hat_vec);
    double j_vec_dot_j_f_hat_vec = dot3(j_vec,j_f_hat_vec);
    double e_vec_dot_j_f_hat_vec = dot3(e_vec,j_f_hat_vec);
    double j_vec_dot_e_f_hat_vec = dot3(j_vec,e_f_hat_vec);
    
    double e_vec_dot_e_f_hat_vec_pl1,e_vec_dot_e_f_hat_vec_pl1m1;
    double j_vec_dot_j_f_hat_vec_pl2,j_vec_dot_j_f_hat_vec_pl2m1;
    double e_vec_dot_j_f_hat_vec_pl3,e_vec_dot_j_f_hat_vec_pl3m1;
    double j_vec_dot_e_f_hat_vec_pl4,j_vec_dot_e_f_hat_vec_pl4m1;
    
    double m1 = binary->child1_mass;
    double m2 = binary->child2_mass;

    double m1_plus_m2 = binary->child1_mass_plus_child2_mass;
    double m1_minus_m2 = binary->child1_mass_minus_child2_mass;
    double m1_times_m2 = binary->child1_mass_times_child2_mass;

    int n,m,i1,i2,l1,l2,l3,l4; /* n: order of expansion, starts at n = 2; 0 <= m <= n; i1 + i2 <= m */
    double constant_integrated_hamiltonian_factor = (m1_times_m2/m1_plus_m2)*(CONST_G*M_per/(q_f*one_plus_e_f))*(1.0/n_f)*e_f_p2_minus_one_p3div2; /* constant factor in Hamiltonian*/

    
    double M_bin_pnm1; /* M_bin to power n - 1 */
    double M_bin_child1_pnm1; /* M_bin_child1 to power n - 1 */
    double M_bin_child2_pnm1; /* M_bin_child2 to power n - 1 */
    double minusone_pnp1; /* minus one to power n + 1 */
    double a_pn; /* binary a to power n */
    double q_f_pm1 = 1.0/q_f;
    double q_f_pmn;
    
    double mass_factor_children = 0.0;
    double binary_pair_integrated_hamiltonian = 0.0;
    double integrated_hamiltonian_factor = 0.0;
    double A_n_m = 0.0;
    double B_n_m_i1_i2 = 0.0;
    double dB_n_m_i1_i2_de = 0.0; /* derivative of B function w.r.t. e */
    int D_n_i1_i2_l1_l2_l3_l4_function_index;
    double D_n_i1_i2_l1_l2_l3_l4 = 0.0;
    double dD_n_i1_i2_l1_l2_l3_l4_de = 0.0; /* derivative of D function w.r.t. e */
    
    double e_p_array_even[HIGHEST_POWER_ECCP2_IN_B_TABLE+1];
    double e_p_array_odd[HIGHEST_POWER_ECCP2_IN_B_TABLE+1];
    e_p_array_even[0] = 1.0;
    e_p_array_odd[1] = e;
    int index_B_eccp2;
    for (index_B_eccp2=1; index_B_eccp2 <= HIGHEST_POWER_ECCP2_IN_B_TABLE; index_B_eccp2++)
    {

        e_p_array_even[index_B_eccp2] = e_p_array_even[index_B_eccp2-1]*e_p2;
    
        if (index_B_eccp2>1)
        {
            e_p_array_odd[index_B_eccp2] = e_p_array_odd[index_B_eccp2-1]*e_p2;
        }
    }
    
    int index_A,index_B,index_D;
    int n_lookup,m_lookup,i1_lookup,i2_lookup;
    
    bool continue_after_A_table;
    double B_lookup = 0.0;
    for (index_B=0; index_B<TABLELENGTH_B; index_B++)
    {
        B_n_m_i1_i2 = 0.0;
        dB_n_m_i1_i2_de = 0.0;
        
        n = int(B_TABLE[index_B][0]);
        m = int(B_TABLE[index_B][1]);
        i1 = int(B_TABLE[index_B][2]);
        i2 = int(B_TABLE[index_B][3]);
        
        /* include only user-specific orders */
        if (n==2 && include_quadrupole_order_terms==false) continue;
        if (n==3 && include_octupole_order_binary_pair_terms==false) continue;
        if (n==4 && include_hexadecupole_order_binary_pair_terms==false) continue;
        if (n==5 && include_dotriacontupole_order_binary_pair_terms==false) continue;
        
        /* retrieve A from table */
        continue_after_A_table = false;
        for (index_A=0; index_A<TABLELENGTH_A; index_A++)
        {
            n_lookup = int(A_TABLE[index_A][0]);
            m_lookup = int(A_TABLE[index_A][1]);
            
            if ((n_lookup == n) && (m_lookup == m))
            {
                A_n_m = A_TABLE[index_A][2];
                continue_after_A_table = true;
            }
        }     
        
        /* continue in B table if there is no entry for given n and m (i.e. A_n_m = 0) */   
        if (continue_after_A_table == false) continue;
        
        /* construct B function from table */
        B_n_m_i1_i2 = 0.0;
        dB_n_m_i1_i2_de = 0.0;
        for (index_B_eccp2=0; index_B_eccp2 <= HIGHEST_POWER_ECCP2_IN_B_TABLE; index_B_eccp2++)
        {
            B_lookup = B_TABLE[index_B][4+index_B_eccp2]; /* take into account offset for n,m,i1,i2 */
            
            B_n_m_i1_i2 += B_lookup*e_p_array_even[index_B_eccp2];
            if (index_B_eccp2>0)
            {
                dB_n_m_i1_i2_de += 2.0*index_B_eccp2*B_lookup*e_p_array_odd[index_B_eccp2];
            }
        }

        /* compute quantities depending on n */
        minusone_pnp1 = pow(-1.0,n+1.0);
        M_bin_child1_pnm1 = pow(m1,n-1.0);
        M_bin_child2_pnm1 = pow(m2,n-1.0);
        M_bin_pnm1 = pow(m1_plus_m2,n-1.0);            
        a_pn = pow(a,n);
        q_f_pmn = pow(q_f_pm1,n);
        one_plus_e_f_pmn = pow(one_plus_e_f_pm1,n);
        
        //mass_factor_children = fabs(M_bin_child1_pnm1 - minusone_pnp1*M_bin_child2_pnm1);
        mass_factor_children = M_bin_child2_pnm1 - minusone_pnp1 * M_bin_child1_pnm1;
        integrated_hamiltonian_factor = minusone_pnp1*constant_integrated_hamiltonian_factor*(mass_factor_children/M_bin_pnm1)*a_pn*q_f_pmn*one_plus_e_f_pmn*A_n_m;
        
        /* construct D function from table */
        for (index_D=0; index_D<TABLELENGTH_D; index_D++)
        {
            n_lookup = int(D_TABLE[index_D][0]);
            i1_lookup = int(D_TABLE[index_D][1]);
            i2_lookup = int(D_TABLE[index_D][2]);
            if ((n_lookup == n) && (i1_lookup == i1) && (i2_lookup == i2))
            {
                l1 = int(D_TABLE[index_D][3]);
                l2 = int(D_TABLE[index_D][4]);
                l3 = int(D_TABLE[index_D][5]);
                l4 = int(D_TABLE[index_D][6]);
                D_n_i1_i2_l1_l2_l3_l4_function_index = int(D_TABLE[index_D][7]);

                D_n_i1_i2_l1_l2_l3_l4 = retrieve_D_function(D_n_i1_i2_l1_l2_l3_l4_function_index,e,e_p2,e_f,e_f_p2,e_f_p4,one_div_e_f_p1,one_div_e_f_p2,one_div_e_f_p3,e_f_p2_minus_one,sqrt_e_f_p2_minus_one,asec_minus_e_f);
                dD_n_i1_i2_l1_l2_l3_l4_de = retrieve_D_function_e_derivative(D_n_i1_i2_l1_l2_l3_l4_function_index,e,e_p2,e_f,e_f_p2,e_f_p4,one_div_e_f_p1,one_div_e_f_p2,one_div_e_f_p3,e_f_p2_minus_one,sqrt_e_f_p2_minus_one,asec_minus_e_f);
                    
                //printf("test... n %d m %d i1 %d i2 %d e %g B %g \n",n,m,i1,i2,e,B_n_m_i1_i2);
                //printf("============================\n");
                //printf("n %d m %d i1 %d i2 %d l1 %d l2 %d l3 %d l4 %d \n",n,m,i1,i2,l1,l2,l3,l4);
                
                //printf("D_n_i1_i2_l1_l2_l3_l4_function_index %d\n",D_n_i1_i2_l1_l2_l3_l4_function_index);
                //printf("D %g dDde %g\n",D_n_i1_i2_l1_l2_l3_l4,dD_n_i1_i2_l1_l2_l3_l4_de);

                /* compute quantities depending on l1,l2,l3,l4 */
                e_vec_dot_e_f_hat_vec_pl1 = pow(e_vec_dot_e_f_hat_vec,l1);
                j_vec_dot_j_f_hat_vec_pl2 = pow(j_vec_dot_j_f_hat_vec,l2);
                e_vec_dot_j_f_hat_vec_pl3 = pow(e_vec_dot_j_f_hat_vec,l3);
                j_vec_dot_e_f_hat_vec_pl4 = pow(j_vec_dot_e_f_hat_vec,l4);

                e_vec_dot_e_f_hat_vec_pl1m1 = pow(e_vec_dot_e_f_hat_vec,l1-1.0);
                j_vec_dot_j_f_hat_vec_pl2m1 = pow(j_vec_dot_j_f_hat_vec,l2-1.0);
                e_vec_dot_j_f_hat_vec_pl3m1 = pow(e_vec_dot_j_f_hat_vec,l3-1.0);
                j_vec_dot_e_f_hat_vec_pl4m1 = pow(j_vec_dot_e_f_hat_vec,l4-1.0);
                
                /* compute integrated Hamiltonian */
                binary_pair_integrated_hamiltonian += integrated_hamiltonian_factor*D_n_i1_i2_l1_l2_l3_l4*e_vec_dot_e_f_hat_vec_pl1*j_vec_dot_j_f_hat_vec_pl2*e_vec_dot_j_f_hat_vec_pl3*j_vec_dot_e_f_hat_vec_pl4;
                
                //printf("integrated_hamiltonian_factor %g\n",integrated_hamiltonian_factor);
                //printf("binary_pair_integrated_hamiltonian %g\n",binary_pair_integrated_hamiltonian);

                /* compute gradients of integrated Hamiltonian */
                for (int i=0; i<3; i++)
                {
                    grad_e_vec_H[i] += integrated_hamiltonian_factor*j_vec_dot_j_f_hat_vec_pl2*j_vec_dot_e_f_hat_vec_pl4*( \
                        (dB_n_m_i1_i2_de*D_n_i1_i2_l1_l2_l3_l4 + B_n_m_i1_i2*dD_n_i1_i2_l1_l2_l3_l4_de)*e_vec_dot_e_f_hat_vec_pl1*e_vec_dot_j_f_hat_vec_pl3*e_vec_unit[i] \
                        + B_n_m_i1_i2*D_n_i1_i2_l1_l2_l3_l4*( double(l1)*e_vec_dot_e_f_hat_vec_pl1m1*e_vec_dot_j_f_hat_vec_pl3*e_f_hat_vec[i] \
                        + e_vec_dot_e_f_hat_vec_pl1*double(l3)*e_vec_dot_j_f_hat_vec_pl3m1*j_f_hat_vec[i]) \
                    );
                    
                    grad_j_vec_H[i] += integrated_hamiltonian_factor*B_n_m_i1_i2*D_n_i1_i2_l1_l2_l3_l4*e_vec_dot_e_f_hat_vec_pl1*e_vec_dot_j_f_hat_vec_pl3*( \
                        double(l2)*j_vec_dot_j_f_hat_vec_pl2m1*j_vec_dot_e_f_hat_vec_pl4*j_f_hat_vec[i] + j_vec_dot_j_f_hat_vec_pl2*double(l4)*j_vec_dot_e_f_hat_vec_pl4m1*e_f_hat_vec[i] );
                    
                }
            }
        }
    }

    double j_vec_cross_grad_j_vec_H[3],j_vec_cross_grad_e_vec_H[3];
    double e_vec_cross_grad_e_vec_H[3],e_vec_cross_grad_j_vec_H[3];
    
    cross3(j_vec,        grad_j_vec_H,              j_vec_cross_grad_j_vec_H);
    cross3(j_vec,        grad_e_vec_H,              j_vec_cross_grad_e_vec_H);
    cross3(e_vec,        grad_e_vec_H,              e_vec_cross_grad_e_vec_H);
    cross3(e_vec,        grad_j_vec_H,              e_vec_cross_grad_j_vec_H);

    double Lambda = h/j;
    /* note: `H' is the time-integrated Hamiltonian */
    for (int i=0; i<3; i++)
    {
        binary->e_vec[i] += (-1.0/(Lambda))*( e_vec_cross_grad_j_vec_H[i] \
            + j_vec_cross_grad_e_vec_H[i] );
        binary->h_vec[i] += -1.0*( j_vec_cross_grad_j_vec_H[i] \
            + e_vec_cross_grad_e_vec_H[i] );
    }
}

double retrieve_D_function(int function_index, double ep, double ep_p2, double ef, double ef_p2, double ef_p4, double one_div_ef_p1, double one_div_ef_p2, double one_div_ef_p3, double ef_p2_minus_one, double sqrt_ef_p2_minus_one, double asec_minus_ef)
{
    double result = 0.0;
    switch (function_index)
    {
        case 1: result = D_TABLE_FUNC1(ep,ep_p2,ef,ef_p2,ef_p4,one_div_ef_p1,one_div_ef_p2,one_div_ef_p3,ef_p2_minus_one,sqrt_ef_p2_minus_one,asec_minus_ef); break;
        case 2: result = D_TABLE_FUNC2(ep,ep_p2,ef,ef_p2,ef_p4,one_div_ef_p1,one_div_ef_p2,one_div_ef_p3,ef_p2_minus_one,sqrt_ef_p2_minus_one,asec_minus_ef); break;
        case 3: result = D_TABLE_FUNC3(ep,ep_p2,ef,ef_p2,ef_p4,one_div_ef_p1,one_div_ef_p2,one_div_ef_p3,ef_p2_minus_one,sqrt_ef_p2_minus_one,asec_minus_ef); break;
        case 4: result = D_TABLE_FUNC4(ep,ep_p2,ef,ef_p2,ef_p4,one_div_ef_p1,one_div_ef_p2,one_div_ef_p3,ef_p2_minus_one,sqrt_ef_p2_minus_one,asec_minus_ef); break;
        case 5: result = D_TABLE_FUNC5(ep,ep_p2,ef,ef_p2,ef_p4,one_div_ef_p1,one_div_ef_p2,one_div_ef_p3,ef_p2_minus_one,sqrt_ef_p2_minus_one,asec_minus_ef); break;
        case 6: result = D_TABLE_FUNC6(ep,ep_p2,ef,ef_p2,ef_p4,one_div_ef_p1,one_div_ef_p2,one_div_ef_p3,ef_p2_minus_one,sqrt_ef_p2_minus_one,asec_minus_ef); break;
        case 7: result = D_TABLE_FUNC7(ep,ep_p2,ef,ef_p2,ef_p4,one_div_ef_p1,one_div_ef_p2,one_div_ef_p3,ef_p2_minus_one,sqrt_ef_p2_minus_one,asec_minus_ef); break;
        case 8: result = D_TABLE_FUNC8(ep,ep_p2,ef,ef_p2,ef_p4,one_div_ef_p1,one_div_ef_p2,one_div_ef_p3,ef_p2_minus_one,sqrt_ef_p2_minus_one,asec_minus_ef); break;
        case 9: result = D_TABLE_FUNC9(ep,ep_p2,ef,ef_p2,ef_p4,one_div_ef_p1,one_div_ef_p2,one_div_ef_p3,ef_p2_minus_one,sqrt_ef_p2_minus_one,asec_minus_ef); break;
        case 10: result = D_TABLE_FUNC10(ep,ep_p2,ef,ef_p2,ef_p4,one_div_ef_p1,one_div_ef_p2,one_div_ef_p3,ef_p2_minus_one,sqrt_ef_p2_minus_one,asec_minus_ef); break;
        case 11: result = D_TABLE_FUNC11(ep,ep_p2,ef,ef_p2,ef_p4,one_div_ef_p1,one_div_ef_p2,one_div_ef_p3,ef_p2_minus_one,sqrt_ef_p2_minus_one,asec_minus_ef); break;
        case 12: result = D_TABLE_FUNC12(ep,ep_p2,ef,ef_p2,ef_p4,one_div_ef_p1,one_div_ef_p2,one_div_ef_p3,ef_p2_minus_one,sqrt_ef_p2_minus_one,asec_minus_ef); break;
        case 13: result = D_TABLE_FUNC13(ep,ep_p2,ef,ef_p2,ef_p4,one_div_ef_p1,one_div_ef_p2,one_div_ef_p3,ef_p2_minus_one,sqrt_ef_p2_minus_one,asec_minus_ef); break;
        case 14: result = D_TABLE_FUNC14(ep,ep_p2,ef,ef_p2,ef_p4,one_div_ef_p1,one_div_ef_p2,one_div_ef_p3,ef_p2_minus_one,sqrt_ef_p2_minus_one,asec_minus_ef); break;
        case 15: result = D_TABLE_FUNC15(ep,ep_p2,ef,ef_p2,ef_p4,one_div_ef_p1,one_div_ef_p2,one_div_ef_p3,ef_p2_minus_one,sqrt_ef_p2_minus_one,asec_minus_ef); break;
    }
    //printf("function_index %d %g %g %g result %g \n",function_index,ep,ep_p2,ef,result);    
    return result;
}

double retrieve_D_function_e_derivative(int function_index, double ep, double ep_p2, double ef, double ef_p2, double ef_p4, double one_div_ef_p1, double one_div_ef_p2, double one_div_ef_p3, double ef_p2_minus_one, double sqrt_ef_p2_minus_one, double asec_minus_ef)
{
    double result = 0.0;
    switch (function_index)
    {
        case 1: result = D_TABLE_FUNC1_DER(ep,ep_p2,ef,ef_p2,ef_p4,one_div_ef_p1,one_div_ef_p2,one_div_ef_p3,ef_p2_minus_one,sqrt_ef_p2_minus_one,asec_minus_ef); break;
        case 2: result = D_TABLE_FUNC2_DER(ep,ep_p2,ef,ef_p2,ef_p4,one_div_ef_p1,one_div_ef_p2,one_div_ef_p3,ef_p2_minus_one,sqrt_ef_p2_minus_one,asec_minus_ef); break;
        case 3: result = D_TABLE_FUNC3_DER(ep,ep_p2,ef,ef_p2,ef_p4,one_div_ef_p1,one_div_ef_p2,one_div_ef_p3,ef_p2_minus_one,sqrt_ef_p2_minus_one,asec_minus_ef); break;
        case 4: result = D_TABLE_FUNC4_DER(ep,ep_p2,ef,ef_p2,ef_p4,one_div_ef_p1,one_div_ef_p2,one_div_ef_p3,ef_p2_minus_one,sqrt_ef_p2_minus_one,asec_minus_ef); break;
        case 5: result = D_TABLE_FUNC5_DER(ep,ep_p2,ef,ef_p2,ef_p4,one_div_ef_p1,one_div_ef_p2,one_div_ef_p3,ef_p2_minus_one,sqrt_ef_p2_minus_one,asec_minus_ef); break;
        case 6: result = D_TABLE_FUNC6_DER(ep,ep_p2,ef,ef_p2,ef_p4,one_div_ef_p1,one_div_ef_p2,one_div_ef_p3,ef_p2_minus_one,sqrt_ef_p2_minus_one,asec_minus_ef); break;
        case 7: result = D_TABLE_FUNC7_DER(ep,ep_p2,ef,ef_p2,ef_p4,one_div_ef_p1,one_div_ef_p2,one_div_ef_p3,ef_p2_minus_one,sqrt_ef_p2_minus_one,asec_minus_ef); break;
        case 8: result = D_TABLE_FUNC8_DER(ep,ep_p2,ef,ef_p2,ef_p4,one_div_ef_p1,one_div_ef_p2,one_div_ef_p3,ef_p2_minus_one,sqrt_ef_p2_minus_one,asec_minus_ef); break;
        case 9: result = D_TABLE_FUNC9_DER(ep,ep_p2,ef,ef_p2,ef_p4,one_div_ef_p1,one_div_ef_p2,one_div_ef_p3,ef_p2_minus_one,sqrt_ef_p2_minus_one,asec_minus_ef); break;
        case 10: result = D_TABLE_FUNC10_DER(ep,ep_p2,ef,ef_p2,ef_p4,one_div_ef_p1,one_div_ef_p2,one_div_ef_p3,ef_p2_minus_one,sqrt_ef_p2_minus_one,asec_minus_ef); break;
        case 11: result = D_TABLE_FUNC11_DER(ep,ep_p2,ef,ef_p2,ef_p4,one_div_ef_p1,one_div_ef_p2,one_div_ef_p3,ef_p2_minus_one,sqrt_ef_p2_minus_one,asec_minus_ef); break;
        case 12: result = D_TABLE_FUNC12_DER(ep,ep_p2,ef,ef_p2,ef_p4,one_div_ef_p1,one_div_ef_p2,one_div_ef_p3,ef_p2_minus_one,sqrt_ef_p2_minus_one,asec_minus_ef); break;
        case 13: result = D_TABLE_FUNC13_DER(ep,ep_p2,ef,ef_p2,ef_p4,one_div_ef_p1,one_div_ef_p2,one_div_ef_p3,ef_p2_minus_one,sqrt_ef_p2_minus_one,asec_minus_ef); break;
        case 14: result = D_TABLE_FUNC14_DER(ep,ep_p2,ef,ef_p2,ef_p4,one_div_ef_p1,one_div_ef_p2,one_div_ef_p3,ef_p2_minus_one,sqrt_ef_p2_minus_one,asec_minus_ef); break;
        case 15: result = D_TABLE_FUNC15_DER(ep,ep_p2,ef,ef_p2,ef_p4,one_div_ef_p1,one_div_ef_p2,one_div_ef_p3,ef_p2_minus_one,sqrt_ef_p2_minus_one,asec_minus_ef); break;
    }
    return result;
}
}
