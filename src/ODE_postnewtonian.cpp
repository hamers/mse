/* MSE */

#include "types.h"
#include "ODE_postnewtonian.h"
#include "evolve.h"
//#include <stdio.h>
extern "C"
{
    
void compute_EOM_Post_Newtonian_for_particle(ParticlesMap *particlesMap, Particle *p, double *hamiltonian, bool compute_hamiltonian_only)
{
    if (p->include_pairwise_1PN_terms == true)
    {
        *hamiltonian += compute_EOM_pairwise_1PN(particlesMap,p->index,compute_hamiltonian_only);
    }
    if (p->include_pairwise_25PN_terms == true)
    {
        *hamiltonian += compute_EOM_pairwise_25PN(particlesMap,p->index,compute_hamiltonian_only);
    }
    
    Particle *child1 = (*particlesMap)[p->child1];
    Particle *child2 = (*particlesMap)[p->child2];
    if (child1->is_binary == false && child1->include_spin_orbit_1PN_terms == true)
    {
        *hamiltonian += compute_EOM_spin_orbit_coupling_1PN(particlesMap,p->index,child1->index,child2->index,false);
    }
    if (child2->is_binary == false && child2->include_spin_orbit_1PN_terms == true)
    {
        *hamiltonian += compute_EOM_spin_orbit_coupling_1PN(particlesMap,p->index,child2->index,child1->index,false);
    }

}

double compute_EOM_pairwise_1PN(ParticlesMap *particlesMap, int binary_index, bool compute_hamiltonian_only)
{
    #ifdef VERBOSE
    if (verbose_flag > 1)
    {
        printf("ODE_postnewtonian.cpp -- compute_EOM_pairwise_1PN\n");
    }
    #endif

    Particle *binary = (*particlesMap)[binary_index];
    double e = binary->e;
    double a = binary->a;
    double *e_vec_unit = binary->e_vec_unit;
    double *h_vec_unit = binary->h_vec_unit;
    double j = binary->j;
    double j_p2 = binary->j_p2;
    double m1 = binary->child1_mass;
    double m2 = binary->child2_mass;
    double mt = m1+m2;

    if (binary->exclude_for_secular_integration == true)
    {
        #ifdef VERBOSE
        if (verbose_flag > 1)
        {
            printf("ODE_postnewtonian.cpp -- compute_EOM_pairwise_1PN -- not applying 1PN terms for particle %d\n",binary->index);
        }
        #endif

        return 0;
    }


    double hamiltonian_1PN = -3.0*CONST_G_P2*m1*m2*mt/(a*a*CONST_C_LIGHT_P2*j);
    if (compute_hamiltonian_only == true)
    {
        return hamiltonian_1PN;
    }
    
    double GMdiva = CONST_G*mt/a;
    double Z_1PN = 3.0*sqrt(GMdiva)*GMdiva/(a*CONST_C_LIGHT_P2*j_p2);
    
    if (binary->parent == -1 and binary->exclude_1PN_precession_in_case_of_isolated_binary == true)
    {
        Z_1PN = 0.0;
    }
    
    for (int i=0; i<3; i++)
    {
        binary->de_vec_dt[i] += e*Z_1PN*binary->q_vec_unit[i];
    }
    
    return hamiltonian_1PN;
}

double compute_EOM_pairwise_25PN(ParticlesMap *particlesMap, int binary_index, bool compute_hamiltonian_only)
{
    #ifdef VERBOSE
    if (verbose_flag > 1)
    {
        printf("ODE_postnewtonian.cpp -- compute_EOM_pairwise_25PN\n");
    }
    #endif
    
    Particle *binary = (*particlesMap)[binary_index];
    double e = binary->e;
    double e_p2 = e*e;
    double a = binary->a;
    double *e_vec_unit = binary->e_vec_unit;
    double *h_vec_unit = binary->h_vec_unit;
    double j = binary->j;
    double j_p2 = binary->j_p2;
    double j_p4 = binary->j_p4;
    double m1 = binary->child1_mass;
    double m2 = binary->child2_mass;
    double mt = m1+m2;

    double a_p3 = a*a*a;
    double GMdiva = CONST_G*mt/a;
    double c_common = CONST_G_P3*m1*m2/(CONST_C_LIGHT_P5*a_p3*j_p4);
    double f_e = 1.0 + c_121div304*e_p2;
    double f_h = 1.0 + c_7div8*e_p2;

    double de_dt = -c_304div15*c_common*mt*e*f_e/(a*j);
    double dh_dt = -c_32div5*c_common*m1*m2*sqrt(GMdiva)*f_h;

    for (int i=0; i<3; i++)
    {
        binary->de_vec_dt[i] += de_dt*e_vec_unit[i];
        binary->dh_vec_dt[i] += dh_dt*h_vec_unit[i];
    }
    
    return 0.0; // N/A
}

double compute_EOM_spin_orbit_coupling_1PN(ParticlesMap *particlesMap, int binary_index, int body_index, int companion_index, bool compute_hamiltonian_only)
{
    /* 1PN spin-orbit coupling */
    /* See, e.g., https://arxiv.org/abs/1711.07142 Eqs. (2) and (3) or https://journals.aps.org/prd/abstract/10.1103/PhysRevD.12.329
     * NOTE: the equations of motion strictly apply to the spin angular momentum vector S. Here,
     * we apply the equations of motion for the angular frequency vector \Omega. This should be correct
     * as long as S and \Omega are parallel. */
     
    Particle *binary = (*particlesMap)[binary_index];
    Particle *body = (*particlesMap)[body_index];
    Particle *companion = (*particlesMap)[companion_index];
    
    double mp = body->mass;
    double mc = companion->mass;
    
    double h = binary->h;
    double *h_vec = binary->h_vec;
    double a = binary->a;
    double j_p3 = binary->j_p3;
    
    double constant_factor = (2.0*CONST_G/(CONST_C_LIGHT_P2*a*a*a*j_p3))*(1.0 + c_3div4*mc/mp);
    
    double *spin_vec = body->spin_vec;
    double h_vec_cross_spin_vec[3];
    double z_vec[3];
    cross3(h_vec,spin_vec,h_vec_cross_spin_vec);

    for (int i=0; i<3; i++)
    {
        body->dspin_vec_dt[i] += constant_factor*h_vec_cross_spin_vec[i];
        
        /* include effect on orbit? */
        
    }

    #ifdef VERBOSE
    if (verbose_flag > 1)
    {
        printf("ODE_postnewtonian -- compute_EOM_spin_orbit_coupling_1PN -- bin %d body %d comp %d constant_factor %g body->dspin_vec_dt %g %g %g\n",binary_index,body_index,companion_index,constant_factor,body->dspin_vec_dt[0],body->dspin_vec_dt[1],body->dspin_vec_dt[2]);
    }
    #endif
    
    return 0.0;

}

double compute_spin_parameter_from_spin_frequency(double m, double Omega)
{
    double A = 2.0 * CONST_G * m * Omega;
    double chi = 2.0 * CONST_C_LIGHT_P3 * A / (CONST_C_LIGHT_P6 + A*A);
    
    return chi;
}

double compute_spin_frequency_from_spin_parameter(double m, double chi)
{
    double Omega = CONST_C_LIGHT_P3 * chi / ( 2.0 * CONST_G * m * (1.0 + sqrt(1.0 - chi*chi)) );
    
    return Omega;
}

double compute_spin_frequency_dot_BHs(double m, double Omega, double J_dot, double m_dot)
{
    double chi = compute_spin_parameter_from_spin_frequency(m, Omega);
    double J = chi * CONST_G * m * m / CONST_C_LIGHT;
    double J_dot_div_J = J_dot/J;
    double m_dot_div_m = m_dot/m;
    
    double chi_sq = chi*chi;
    double sqrt_one_minus_chi_sq = sqrt(1.0 - chi_sq);
    
    double Omega_dot = Omega*( J_dot_div_J - 3.0*m_dot_div_m + (chi_sq/(sqrt_one_minus_chi_sq * (1.0 + sqrt_one_minus_chi_sq))) * (J_dot_div_J - 2.0*m_dot_div_m) );
    
    return Omega_dot;
}

double compute_1PN_timescale(double a, double M, double e)
{
    double j_p2 = 1.0 - e*e;
    
    double P_orb = compute_orbital_period_from_semimajor_axis(M, a);
    double rg = CONST_G*M/(CONST_C_LIGHT_P2);
    double t_1PN = c_1div3 * P_orb * j_p2 * (a/rg);
    return t_1PN;
}

}
