/* MSE */

#include "types.h"
#include "ODE_newtonian.h"
#include "../interface.h" /* for parameters */
#include "evolve.h"
//#include <stdio.h>

extern "C"
{
    
void compute_EOM_Newtonian_for_particle(ParticlesMap *particlesMap, Particle *p, double *hamiltonian, double *KS_V, bool compute_hamiltonian_only)
{
    
    if (p->exclude_for_secular_integration == true)
    {
        #ifdef VERBOSE
        if (verbose_flag > 1)
        {
            printf("ODE_newtonian.cpp -- compute_EOM_Newtonian_for_particle -- not applying terms for particle %d\n",p->index);
        }
        #endif

        return;
    }
    
    std::vector<int>::iterator it_parent_p,it_parent_q;

    /* binary pairs */
    for (it_parent_p = p->parents.begin(); it_parent_p != p->parents.end(); it_parent_p++)
    {
        
        int i = std::distance(p->parents.begin(), it_parent_p);
        Particle *P_q = (*particlesMap)[(*it_parent_p)];
        int connecting_child_in_parent_q = p->connecting_child_in_parents[i];
        
        compute_EOM_binary_pairs(particlesMap,p->index,P_q->index,connecting_child_in_parent_q,hamiltonian,KS_V,compute_hamiltonian_only);
        
        /* binary triplets */
        for (it_parent_q = P_q->parents.begin(); it_parent_q != P_q->parents.end(); it_parent_q++)
        {
            int j = std::distance(P_q->parents.begin(), it_parent_q);
            Particle *P_u = (*particlesMap)[(*it_parent_q)];
            int connecting_child_in_parent_u = P_q->connecting_child_in_parents[j];
            compute_EOM_binary_triplets(particlesMap,p->index,P_q->index,P_u->index,connecting_child_in_parent_q,connecting_child_in_parent_u,hamiltonian,KS_V,compute_hamiltonian_only);

        }
    }
}            

void compute_EOM_binary_pairs(ParticlesMap *particlesMap, int inner_binary_index, int outer_binary_index, int connecting_child_in_outer_binary, double *hamiltonian, double *KS_V, bool compute_hamiltonian_only)
{


    Particle *inner_binary = (*particlesMap)[inner_binary_index];
    Particle *outer_binary = (*particlesMap)[outer_binary_index];
    
    if (inner_binary->integration_method>0)
    {
        if (outer_binary->integration_method>0)
        {
            compute_EOM_binary_pairs_non_averaged(particlesMap,inner_binary_index,outer_binary_index,connecting_child_in_outer_binary,hamiltonian,KS_V,compute_hamiltonian_only);
        }
        else if (outer_binary->integration_method==0)
        {
            printf("FATAL ERROR in newtonian.cpp: an outer orbit cannot be averaged if the inner one is not! inner_binary_index=%d; outer_binary_index=%d\n",inner_binary_index,outer_binary_index);
            exit(-1);
        }   
    }
    else if (inner_binary->integration_method==0)
    {
        if (outer_binary->integration_method>0)
        {
            compute_EOM_binary_pairs_single_averaged(particlesMap,inner_binary_index,outer_binary_index,connecting_child_in_outer_binary,hamiltonian,KS_V,compute_hamiltonian_only);
        }
        else if (outer_binary->integration_method==0)
        {
            compute_EOM_binary_pairs_double_averaged(particlesMap,inner_binary_index,outer_binary_index,connecting_child_in_outer_binary,hamiltonian,KS_V,compute_hamiltonian_only);
        }
    }
    
}

void compute_EOM_binary_pairs_non_averaged(ParticlesMap *particlesMap, int inner_binary_index, int outer_binary_index, int child_connected_to_p_in_parent_k, double *hamiltonian, double *KS_V, bool compute_hamiltonian_only)
{
    /* last checked/updated 05/03/20 */
    
    /* stop if no pairwise terms are to be computed */
    if ((include_quadrupole_order_terms == false) && (include_octupole_order_binary_pair_terms == false) && (include_hexadecupole_order_binary_pair_terms == false) && (include_dotriacontupole_order_binary_pair_terms == false) )
    {
        return;
    }
    
    Particle *p = (*particlesMap)[inner_binary_index];
    Particle *k = (*particlesMap)[outer_binary_index];
    
    Particle *sibling_of_child_connected_to_p_in_k;
    if (child_connected_to_p_in_parent_k==1)
    {
        sibling_of_child_connected_to_p_in_k = (*particlesMap)[k->child2];
    }
    else if (child_connected_to_p_in_parent_k==2)
    {
        sibling_of_child_connected_to_p_in_k = (*particlesMap)[k->child1];
    }

    double alpha = 1.0;
    
    double M_p = p->mass;
    double M_pC1 = p->child1_mass;
    double M_pC2 = p->child2_mass;
    double M_kCSp = sibling_of_child_connected_to_p_in_k->mass;

    double *r_p_vec = p->r_vec;
    double *r_k_vec = k->r_vec;
    //double *v_p_vec = p->v_vec;
    //double *v_k_vec = k->v_vec;
    
    double r_p_vec_dot_r_k_vec = dot3(r_p_vec,r_k_vec);
    //double r_p_vec_dot_v_p_vec = dot3(r_p_vec,v_p_vec);
    //double r_p_vec_dot_v_k_vec = dot3(r_p_vec,v_k_vec);
    //double v_p_vec_dot_r_k_vec = dot3(v_p_vec,r_k_vec);
    //double r_k_vec_dot_v_k_vec = dot3(r_k_vec,v_k_vec);
    
    double r_p_P1 = norm3(r_p_vec);
    double r_k_P1 = norm3(r_k_vec);
    double r_p_P2 = r_p_P1*r_p_P1;
    double r_k_P2 = r_k_P1*r_k_P1;

    double r_p_div_r_k = r_p_P1/r_k_P1;
    double r_p_div_r_k_P2 = r_p_div_r_k*r_p_div_r_k;

    double mu_p = M_pC1*M_pC2/M_p;
    double mu_k = k->child1_mass*k->child2_mass/k->mass;
    
    double C_1 = CONST_G*M_kCSp;
    double C_n,r_p_div_r_k_Pn;
    double n, n_minus_1;
    double M_p_n;
    
    for (int int_n=2; int_n<=5; int_n++)
    {
        /* include only user-specific orders */
        if (int_n==2 && include_quadrupole_order_terms==false) continue;
        if (int_n==3 && include_octupole_order_binary_pair_terms==false) continue;
        if (int_n==4 && include_hexadecupole_order_binary_pair_terms==false) continue;
        if (int_n==5 && include_dotriacontupole_order_binary_pair_terms==false) continue;
        
        n = (double) int_n;
        n_minus_1 = n - 1.0;
        
        r_p_div_r_k_Pn = pow(r_p_div_r_k,n);
        M_p_n = ( pow(M_pC2,n_minus_1) - pow(-M_pC1,n_minus_1) )/pow(M_p,n_minus_1);
        
        C_n = C_1 * pow(-1.0, n) * pow(alpha,n) * M_p_n * r_p_div_r_k_Pn;
        
        get_gravity_binary_pairs_order_n(p,k,int_n,n,C_n, \
            hamiltonian, KS_V, \
            mu_p, mu_k, \
            r_p_vec,r_k_vec, \
            r_p_P1,r_k_P1,r_p_P2,r_k_P2, \
            r_p_div_r_k_P2,r_p_div_r_k_Pn,r_p_vec_dot_r_k_vec);
    }
}

void get_gravity_binary_pairs_order_n(
    Particle *p, Particle *k, \
    int int_n, double n, double C_n, \
    double *hamiltonian, double *KS_V, \
    double mu_p, double mu_k, \
    double *r_p_vec, double *r_k_vec, \
    double r_p_P1, double r_k_P1, double r_p_P2, double r_k_P2, \
    double r_p_div_r_k_P2, double r_p_div_r_k_Pn, double r_p_vec_dot_r_k_vec)
{
    int c;
    double m,m_minus_1,n_minus_m,n_plus_m_plus_1,n_plus_m_plus_3;
    double A_mn;
    double r_p_vec_dot_r_k_vec_Pm_minus_1;
    double r_p_Pminus_m_minus_2;
    double r_k_Pminus_m_minus_1;
    double c_a_p,c_a_k;
    double H;
    
    for (int int_m=0; int_m<=int_n; int_m++)
    {
        A_mn = A_mn_table[int_n-2][int_m];
        if (A_mn==0.0)
        {
            continue;
        }
        
        m = (double) int_m;
        m_minus_1 = m - 1.0;
        n_minus_m = n - m;
        n_plus_m_plus_1 = n + m + 1.0;
        n_plus_m_plus_3 = n_plus_m_plus_1 + 2.0;
        
        r_p_vec_dot_r_k_vec_Pm_minus_1 = pow(r_p_vec_dot_r_k_vec,m_minus_1);

        r_p_Pminus_m_minus_2 = pow(r_p_P1, -m - 2.0);
        r_k_Pminus_m_minus_1 = pow(r_k_P1, -m - 1.0);
        
        c_a_p = C_n * A_mn * r_p_vec_dot_r_k_vec_Pm_minus_1 * r_p_Pminus_m_minus_2 * r_k_Pminus_m_minus_1;
        c_a_k = c_a_p * (mu_p/mu_k) * r_p_div_r_k_P2;

        H = -c_a_p * mu_p * r_p_P2 * r_p_vec_dot_r_k_vec;
        *hamiltonian += H;
        *KS_V += H;
        
        for (c=0; c<3; c++)
        {
            p->a_vec[c] += c_a_p * (m*r_p_P2*r_k_vec[c] + n_minus_m*r_p_vec_dot_r_k_vec*r_p_vec[c]);
            k->a_vec[c] += c_a_k * (m*r_k_P2*r_p_vec[c] - n_plus_m_plus_1*r_p_vec_dot_r_k_vec*r_k_vec[c]);
        
            //p->j_ext[c] += c_j_ext_p*( m*r_p_norm_pow_2*r_k_norm_pow_2*(v_p_dot_r_k+r_p_dot_v_k)*(m_minus_1*r_p_norm_pow_2*r_k[c] + n_minus_m*r_p_dot_r_k*r_p[c]) \
                + r_k_norm_pow_2*r_p_dot_r_k*r_p_dot_v_p*(3.0*m*r_p_norm_pow_2*r_k[c] + n_minus_m*r_p_dot_r_k*r_p[c]) \
                - n_plus_m_plus_1*r_p_norm_pow_2*r_p_dot_r_k*r_k_dot_v_k*(m*r_p_norm_pow_2*r_k[c] + n_minus_m*r_p_dot_r_k*r_p[c]) \
                + m*r_p_norm_pow_2*r_p_norm_pow_2*r_k_norm_pow_2*r_p_dot_r_k*v_k[c] \
                + n_minus_m*r_p_norm_pow_2*r_k_norm_pow_2*r_p_dot_r_k*v_p[c] );

            //k->j_ext[c] += c_j_ext_k*( m*r_p_norm_pow_2*r_k_norm_pow_2*(v_p_dot_r_k+r_p_dot_v_k)*(m_minus_1*r_k_norm_pow_2*r_p[c] - n_plus_m_plus_1*r_p_dot_r_k*r_k[c]) \
                + n_minus_m*r_k_norm_pow_2*r_p_dot_r_k*r_p_dot_v_p*(m*r_k_norm_pow_2*r_p[c] - n_plus_m_plus_1*r_p_dot_r_k*r_k[c]) \
                - n_plus_m_plus_1*n_plus_m_plus_3*r_p_norm_pow_2*r_p_dot_r_k*r_k_dot_v_k*(r_k_norm_pow_2*r_p[c] - r_p_dot_r_k*r_k[c]) \
                + m*r_p_norm_pow_2*r_k_norm_pow_2*r_k_norm_pow_2*r_p_dot_r_k*v_p[c] \
                + n_plus_m_plus_1*r_p_norm_pow_2*r_k_norm_pow_2*r_p_dot_r_k*v_k[c] );
        }
    }
}


void compute_EOM_binary_pairs_single_averaged(ParticlesMap *particlesMap, int inner_binary_index, int outer_binary_index, int child_connected_to_p_in_parent_k, double *hamiltonian, double *KS_V, bool compute_hamiltonian_only)
{
    /* last checked/updated 05/03/20 */

    /* stop if no pairwise terms are to be computed */
    if ((include_quadrupole_order_terms == false) && (include_octupole_order_binary_pair_terms == false) && (include_hexadecupole_order_binary_pair_terms == false) && (include_dotriacontupole_order_binary_pair_terms == false) )
    {
        return;
    }
    
    Particle *p = (*particlesMap)[inner_binary_index];
    Particle *k = (*particlesMap)[outer_binary_index];
    
    Particle *sibling_of_child_connected_to_p_in_k;
    if (child_connected_to_p_in_parent_k==1)
    {
        sibling_of_child_connected_to_p_in_k = (*particlesMap)[k->child2];
    }
    else if (child_connected_to_p_in_parent_k==2)
    {
        sibling_of_child_connected_to_p_in_k = (*particlesMap)[k->child1];
    }

    double e = p->e;
    double e_P2 = p->e_p2;
    
    double *e_vec = p->e_vec;
    double *h_vec = p->h_vec;
    
    double *e_vec_unit = p->e_vec_unit;
    double *h_vec_unit = p->h_vec_unit;

    double h = p->h;
    double j = p->j;
        
    double j_vec[3];
    for (int i=0; i<3; i++)
    {
        j_vec[i] = j*h_vec_unit[i];
    }
    
    double a = p->a;

    double M_p = p->mass;
    double M_pC1 = p->child1_mass;
    double M_pC2 = p->child2_mass;    
    double M_kCSp = sibling_of_child_connected_to_p_in_k->mass;

    double mu_p = M_pC1*M_pC2/M_p;
    double mu_k = k->child1_mass*k->child2_mass/k->mass;
    
    double *r_k_vec = k->r_vec;
    double r_k_P1 = norm3(r_k_vec);
    double r_k_P2 = r_k_P1*r_k_P1;
    double r_k_P3 = r_k_P1*r_k_P2;
    double r_k_Pminus_1 = 1.0/r_k_P1;
    double r_k_Pminus_3 = 1.0/r_k_P3;

    double r_k_Pminus_i1_minus_i2;
    
    double a_div_r_k = a*r_k_Pminus_1;
    
    double e_vec_dot_r_k_vec = dot3(e_vec,r_k_vec);
    double j_vec_dot_r_k_vec = dot3(j_vec,r_k_vec);
    double e_vec_dot_r_k_vec_Pi1,j_vec_dot_r_k_vec_Pi2;
    double e_vec_dot_r_k_vec_Pi1_minus_1,j_vec_dot_r_k_vec_Pi2_minus_1;
    
    if (e_vec_dot_r_k_vec==0.0)
    {
        e_vec_dot_r_k_vec = epsilon;
    }
    if (j_vec_dot_r_k_vec==0.0)
    {
        j_vec_dot_r_k_vec = epsilon;
    }
    
    double n,m,i1,i2; /* n: order of expansion, starts at n = 2; 0 <= m <= n; i1 + i2 <= m */
    int int_n,int_m,int_i1,int_i2;
    
    double M_p_n;
    double minus_1_Pn;
    double a_div_r_k_Pn;
    
    double C_1 = CONST_G*M_kCSp;
    
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

    double C_r_k_n,C_r_k_n2,C_H_n;
    
    double grad_e_vec_H[3],grad_j_vec_H[3];
    for (int i=0; i<3; i++)
    {
        grad_e_vec_H[i] = grad_j_vec_H[i] = 0.0;
    }
    
    int index_A,index_B;
    int n_lookup,m_lookup;
    double n_minus_1,n_plus_1_plus_i1_plus_i2;
    
    int i;
    bool continue_to_B_table;
    double B_lookup = 0.0;
    double r_per_pow_mi1mi2,e_vec_dot_r_per_vec_pi1,j_vec_dot_r_per_vec_pi2,e_vec_dot_r_per_vec_pi1m1,j_vec_dot_r_per_vec_pi2m1;
    for (index_B=0; index_B<TABLELENGTH_B; index_B++)
    {
        B_n_m_i1_i2 = 0.0;
        dB_n_m_i1_i2_de = 0.0;
        
        int_n = int(B_TABLE[index_B][0]);
        int_m = int(B_TABLE[index_B][1]);
        int_i1 = int(B_TABLE[index_B][2]);
        int_i2 = int(B_TABLE[index_B][3]);
        
        n = (double) int_n;
        m = (double) int_m;
        i1 = (double) int_i1;
        i2 = (double) int_i2;
        
        /* include only user-specific orders */
        if (int_n==2 && include_quadrupole_order_terms==false) continue;
        if (int_n==3 && include_octupole_order_binary_pair_terms==false) continue;
        if (int_n==4 && include_hexadecupole_order_binary_pair_terms==false) continue;
        if (int_n==5 && include_dotriacontupole_order_binary_pair_terms==false) continue;
        
        /* retrieve A from table */
        continue_to_B_table = false;
        for (index_A=0; index_A<TABLELENGTH_A; index_A++)
        {
            n_lookup = int(A_TABLE[index_A][0]);
            m_lookup = int(A_TABLE[index_A][1]);
            
            if ((n_lookup == int_n) && (m_lookup == int_m))
            {
                A_n_m = A_TABLE[index_A][2];
                continue_to_B_table = true;
            }
        }     
        //printf("look %d %d %d %d %g\n",int_n,int_m,int_i1,int_i2,A_n_m);
        
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
                printf("newtonian.cpp -- FATAL ERROR in constructing B-function \n");
                exit(-1);
            }
            
            B_n_m_i1_i2 += B_lookup * e_Peven;
            dB_n_m_i1_i2_de += 2.0* ( (double) index_B_eccp2) * B_lookup * e_Podd;
        }

        #ifdef VERBOSE
        if (verbose_flag > 2)
        {
            printf("ODE_newtonian.cpp -- compute_EOM_binary_pairs_single_averaged -- n %g m %g i1 %g i2 %g A_n_m %g e %g B %g dB/de %g\n",n,m,i1,i2,A_n_m,e,B_n_m_i1_i2,dB_n_m_i1_i2_de);
        }
        #endif
        
        n_minus_1 = n - 1.0;
        n_plus_1_plus_i1_plus_i2 = n + 1.0 + i1 + i2;
        
        a_div_r_k_Pn = pow(a_div_r_k,n);
        minus_1_Pn = pow(-1.0,n);
        
        M_p_n = ( pow(M_pC2,n_minus_1) - pow(-M_pC1,n_minus_1) )/pow(M_p,n_minus_1);

        r_k_Pminus_i1_minus_i2 = pow(r_k_P1, -i1 - i2);
        e_vec_dot_r_k_vec_Pi1 = pow(e_vec_dot_r_k_vec,i1);
        j_vec_dot_r_k_vec_Pi2 = pow(j_vec_dot_r_k_vec,i2);

        e_vec_dot_r_k_vec_Pi1_minus_1 = pow(e_vec_dot_r_k_vec, i1 - 1.0);
        j_vec_dot_r_k_vec_Pi2_minus_1 = pow(j_vec_dot_r_k_vec, i2 - 1.0);

        C_r_k_n = minus_1_Pn * C_1 * mu_p * M_p_n * a_div_r_k_Pn * A_n_m;
        C_H_n = -C_r_k_n * r_k_Pminus_i1_minus_i2 * r_k_Pminus_1;
        
        #ifdef VERBOSE
        if (verbose_flag > 2)
        {
            printf("ODE_newtonian.cpp -- compute_EOM_binary_pairs_single_averaged -- C_r_k_n %g C_H_n %g r_k_vec %g %g %g e_vec %g %g %g j_vec %g %g %g\n",C_r_k_n,C_H_n,r_k_vec[0],r_k_vec[1],r_k_vec[2],e_vec[0],e_vec[1],e_vec[2],j_vec[0],j_vec[1],j_vec[2]);
        }
        #endif

        /* compute the Hamiltonian */
        binary_pair_hamiltonian += C_H_n * e_vec_dot_r_k_vec_Pi1 * j_vec_dot_r_k_vec_Pi2 * B_n_m_i1_i2;

        /* compute EOM */
        if (compute_hamiltonian_only == false)
        {
            for (i=0; i<3; i++)
            {
                grad_e_vec_H[i] += C_H_n * j_vec_dot_r_k_vec_Pi2 * ( i1 * B_n_m_i1_i2 * e_vec_dot_r_k_vec_Pi1_minus_1 * r_k_vec[i] + e_vec_dot_r_k_vec_Pi1 * dB_n_m_i1_i2_de * e_vec_unit[i] );
                grad_j_vec_H[i] += C_H_n * e_vec_dot_r_k_vec_Pi1 * B_n_m_i1_i2 * i2 * j_vec_dot_r_k_vec_Pi2_minus_1 * r_k_vec[i];
            }
            
            C_r_k_n2 = C_r_k_n * (1.0/mu_k) * B_n_m_i1_i2 * e_vec_dot_r_k_vec_Pi1_minus_1 * j_vec_dot_r_k_vec_Pi2_minus_1 * r_k_Pminus_3 * r_k_Pminus_i1_minus_i2;
            for (i=0; i<3; i++)
            {
                k->a_vec[i] += C_r_k_n2 * ( i1 * j_vec_dot_r_k_vec * r_k_P2 * e_vec[i] + i2 * e_vec_dot_r_k_vec * r_k_P2 * j_vec[i] - n_plus_1_plus_i1_plus_i2 * e_vec_dot_r_k_vec * j_vec_dot_r_k_vec * r_k_vec[i] );

            #ifdef VERBOSE
            if (verbose_flag > 2)
            {
                printf("ODE_newtonian.cpp -- compute_EOM_binary_pairs_single_averaged -- a_vec[i] += %g C_r_k_n2 %g\n",C_r_k_n2 * ( i1 * j_vec_dot_r_k_vec * r_k_P2 * e_vec[i] + i2 * e_vec_dot_r_k_vec * r_k_P2 * j_vec[i] - n_plus_1_plus_i1_plus_i2 * e_vec_dot_r_k_vec * j_vec_dot_r_k_vec * r_k_vec[i]),C_r_k_n2);
            }
            #endif

            }
        }
    }
    *hamiltonian += binary_pair_hamiltonian;
    *KS_V += binary_pair_hamiltonian;
    //binary_pair_hamiltonian += mu_k*( c_1div2*dot3(k->v_vec,k->v_vec) - CONST_G*k->mass/norm3(k->r_vec) );
    
    double j_vec_cross_grad_j_vec_H[3],j_vec_cross_grad_e_vec_H[3];
    double e_vec_cross_grad_e_vec_H[3],e_vec_cross_grad_j_vec_H[3];
    
    cross3(j_vec,        grad_j_vec_H,              j_vec_cross_grad_j_vec_H);
    cross3(j_vec,        grad_e_vec_H,              j_vec_cross_grad_e_vec_H);
    cross3(e_vec,        grad_e_vec_H,              e_vec_cross_grad_e_vec_H);
    cross3(e_vec,        grad_j_vec_H,              e_vec_cross_grad_j_vec_H);

    double Lambda = h/j;

    for (int i=0; i<3; i++)
    {
        p->de_vec_dt[i] += (-1.0/(Lambda))*( e_vec_cross_grad_j_vec_H[i] \
            + j_vec_cross_grad_e_vec_H[i] );
        p->dh_vec_dt[i] += -1.0*( j_vec_cross_grad_j_vec_H[i] \
            + e_vec_cross_grad_e_vec_H[i] );
    }
}


void compute_EOM_binary_pairs_double_averaged(ParticlesMap *particlesMap, int inner_binary_index, int outer_binary_index, int connecting_child_in_outer_binary, double *hamiltonian, double *KS_V, bool compute_hamiltonian_only)
{
    /* last checked/updated 06-03-20 */
    
    /* stop if no triple terms are to be computed for the binary pair */
    if ((include_quadrupole_order_terms == false) && (include_octupole_order_binary_pair_terms == false) && (include_hexadecupole_order_binary_pair_terms == false) && (include_dotriacontupole_order_binary_pair_terms == false) )
    {
        return;
    }


    /*********************
     * preamble          *
     ********************/
    Particle *inner_binary = (*particlesMap)[inner_binary_index];
    Particle *outer_binary = (*particlesMap)[outer_binary_index];
    
    Particle *P_child1 = (*particlesMap)[inner_binary->child1];
    Particle *P_child2 = (*particlesMap)[inner_binary->child2];
    Particle *P_sibling;
    if (connecting_child_in_outer_binary==1)
    {
        P_sibling = (*particlesMap)[outer_binary->child2];
    }
    else if (connecting_child_in_outer_binary==2)
    {
        P_sibling = (*particlesMap)[outer_binary->child1];
    }
    
    #ifdef VERBOSE
    if (verbose_flag > 2)
    {
        printf("ODE_newtonian.cpp -- compute_EOM_binary_pairs_double_averaged -- compute_EOM_binary_pairs inner_binary_index %d outer_binary_index %d connecting_child_in_outer_binary %d P_sibling %d sibling_mass %g\n",inner_binary_index,outer_binary_index,connecting_child_in_outer_binary,P_sibling->index,P_sibling->mass);
    }
    #endif

    
    double e_in = inner_binary->e;
    double e_in_p2 = inner_binary->e_p2;
    double e_in_p4 = e_in_p2*e_in_p2;
    double e_out = outer_binary->e;
    double e_out_p2 = outer_binary->e_p2;
    
    double *e_in_vec = inner_binary->e_vec;
    double *e_out_vec = outer_binary->e_vec;
    double *h_in_vec = inner_binary->h_vec;
    double *h_out_vec = outer_binary->h_vec;
    
    double *e_in_vec_unit = inner_binary->e_vec_unit;
    double *e_out_vec_unit = outer_binary->e_vec_unit;
    double *h_in_vec_unit = inner_binary->h_vec_unit;
    double *h_out_vec_unit = outer_binary->h_vec_unit;
    
    double *q_out_vec_unit = outer_binary->q_vec_unit;
    
    double h_in = inner_binary->h;
    double h_out = outer_binary->h;
    
    double j_in = inner_binary->j;
    double j_in_p2 = inner_binary->j_p2;
    double j_in_p3 = inner_binary->j_p3;
    double j_out = outer_binary->j;
    double j_out_p2 = outer_binary->j_p2;
    double j_out_p3 = outer_binary->j_p3;
    double j_out_p4 = outer_binary->j_p4;
    double j_out_p5 = outer_binary->j_p5;
    double j_out_p6 = j_out*j_out_p5;
    double j_out_p7 = j_out*j_out_p6;
    double j_out_p8 = j_out*j_out_p7;
    double j_out_p9 = j_out*j_out_p8;
    double j_out_p10 = j_out*j_out_p9;
    double j_out_p11 = j_out*j_out_p10;
    double j_out_p13 = j_out_p2*j_out_p11;
    
    double j_out_p2_inv = 1.0/j_out_p2;
    double j_out_p3_inv = 1.0/j_out_p3;
    double j_out_p5_inv = 1.0/j_out_p5;
    double j_out_p7_inv = 1.0/j_out_p7;
    double j_out_p9_inv = 1.0/j_out_p9;
    double j_out_p11_inv = 1.0/j_out_p11;
    double j_out_p13_inv = 1.0/j_out_p13;
    
    double j_in_vec[3],j_out_vec[3];
    for (int i=0; i<3; i++)
    {
        j_in_vec[i] = j_in*h_in_vec_unit[i];
        j_out_vec[i] = j_out*h_out_vec_unit[i];
    }
    
    double a_in = inner_binary->a;
    double a_out = outer_binary->a;
    
    /* set alpha = +1 */
    double m1 = P_child1->mass;
    double m2 = P_child2->mass;
    double m3 = P_sibling->mass;
        
    double m1_plus_m2 = inner_binary->child1_mass_plus_child2_mass;
    double m1_minus_m2 = inner_binary->child1_mass_minus_child2_mass;
    double m1_times_m2 = inner_binary->child1_mass_times_child2_mass;

    double A_quad = c_1div8*CONST_G*(a_in*a_in/(a_out*a_out*a_out))*m1_times_m2*m3/m1_plus_m2;
    double A_oct = A_quad*c_15div8*(a_in/a_out)*(m1_minus_m2)/m1_plus_m2;
    double A_hd = 0.0;
    double A_tc = 0.0;

    if (include_quadrupole_order_terms == false)
    {
        A_quad = 0.0;
    }
    if (include_octupole_order_binary_pair_terms == false)
    {
        A_oct = 0.0;
    }
    if (include_hexadecupole_order_binary_pair_terms == true)
    {
        A_hd = c_3div1024*CONST_G*(a_in*a_in*a_in*a_in/(a_out*a_out*a_out*a_out*a_out))*(m1_times_m2*m3*(m1*m1 - m1_times_m2 + m2*m2)/(m1_plus_m2*m1_plus_m2*m1_plus_m2));
//        A_hd = c_3div1024*CONST_G*(a_in*a_in*a_in*a_in/(a_out*a_out*a_out*a_out*a_out))*(m1_times_m2*m3*(m1*m1*m1 + m2*m2*m2)/(m1_plus_m2*m1_plus_m2*m1_plus_m2*m1_plus_m2));
        /* the above two expressions are mathematically identical */
    }
    if (include_dotriacontupole_order_binary_pair_terms == true)
    {
        A_tc = -c_105div4096*CONST_G*(a_in*a_in*a_in*a_in*a_in/(a_out*a_out*a_out*a_out*a_out*a_out))*(m1_times_m2*m3*(m1_minus_m2)*(m1*m1 + m2*m2)/(m1_plus_m2*m1_plus_m2*m1_plus_m2*m1_plus_m2));  /* changed 06-05-20 */
    }

    double Lambda_in = h_in/j_in;
    double Lambda_out = h_out/j_out;

    double e_in_vec_dot_e_out_vec = dot3(e_in_vec,e_out_vec);
    double j_in_vec_dot_j_out_vec = dot3(j_in_vec,j_out_vec);
    double e_in_vec_dot_j_out_vec = dot3(e_in_vec,j_out_vec);
    double j_in_vec_dot_e_out_vec = dot3(j_in_vec,e_out_vec);
    
    double e_in_vec_dot_e_out_vec_p2 = e_in_vec_dot_e_out_vec*e_in_vec_dot_e_out_vec;
    double j_in_vec_dot_j_out_vec_p2 = j_in_vec_dot_j_out_vec*j_in_vec_dot_j_out_vec;
    double e_in_vec_dot_j_out_vec_p2 = e_in_vec_dot_j_out_vec*e_in_vec_dot_j_out_vec;
    double j_in_vec_dot_e_out_vec_p2 = j_in_vec_dot_e_out_vec*j_in_vec_dot_e_out_vec;

    double j_in_vec_dot_j_out_vec_p4 = j_in_vec_dot_j_out_vec_p2*j_in_vec_dot_j_out_vec_p2;
    double e_in_vec_dot_j_out_vec_p4 = e_in_vec_dot_j_out_vec_p2*e_in_vec_dot_j_out_vec_p2;

    /* dotriacontupole */
    double e_in_vec_dot_e_out_vec_p3 = e_in_vec_dot_e_out_vec*e_in_vec_dot_e_out_vec_p2;
    double j_in_vec_dot_j_out_vec_p3 = j_in_vec_dot_j_out_vec*j_in_vec_dot_j_out_vec_p2;
    double e_in_vec_dot_j_out_vec_p3 = e_in_vec_dot_j_out_vec*e_in_vec_dot_j_out_vec_p2;
    

    /***************************
     * compute the Hamiltonian *
     **************************/

    double f1 = (1.0-6.0*e_in_p2)*j_out_p2 + 15.0*e_in_vec_dot_j_out_vec_p2 - 3.0*j_in_vec_dot_j_out_vec_p2;
	double f2 = (1.0-8.0*e_in_p2)*j_out_p2 + 35.0*e_in_vec_dot_j_out_vec_p2 - 5.0*j_in_vec_dot_j_out_vec_p2;
	double f3 = -10.0*e_in_vec_dot_j_out_vec*j_in_vec_dot_e_out_vec*j_in_vec_dot_j_out_vec;
    double f4,f5,f6,f7,f8,f9,f10,f11,f12; /* hexadecupole */
    double g,g1,g2,g3,h1,h2,h3,h4,h5; /* dotriacontupole */
    if (include_hexadecupole_order_binary_pair_terms == true)
    {
        f4 = -6.0 + e_out_p2 + 40.0*e_in_p2*(1.0 + 8.0*e_out_p2) - 20.0*e_in_p4*(8.0+15.0*e_out_p2);
        f5 = -2.0*j_out_p2 - e_in_p2*j_out_p2 + 21.0*e_in_vec_dot_j_out_vec_p2;
        f6 = (1.0 - 10.0*e_in_p2)*(4.0 + 3.0*e_out_p2);
        f7 = 8.0 + 6.0*e_out_p2 + e_in_p2*(6.0 + 29.0*e_out_p2);
        f8 = j_out_p2 + 13.0*e_in_p2*j_out_p2 - 7.0*j_in_vec_dot_j_out_vec_p2;
        f9 = -2.0 - 3.0*e_out_p2 + 4.0*e_in_p2*(5.0 + 3.0*e_out_p2);
        f10 = j_out_p2 - e_in_p2*j_out_p2 + 7.0*e_in_vec_dot_j_out_vec_p2;
        f11 = 2.0 + e_out_p2;
        f12 = 3.0*f4*j_out_p4 + 420.0*e_in_vec_dot_e_out_vec_p2*j_out_p2*f5 \
            - 5880.0*j_out_p2*e_in_vec_dot_e_out_vec*j_in_vec_dot_e_out_vec*e_in_vec_dot_j_out_vec*j_in_vec_dot_j_out_vec \
            + 5.0*( 28.0*j_out_p2*f6*e_in_vec_dot_j_out_vec_p2 - 6.0*j_out_p2*f7*j_in_vec_dot_j_out_vec_p2 \
                - 12.0*j_out_p2*f8*j_in_vec_dot_e_out_vec_p2 + 98.0*j_out_p2*f9*e_in_vec_dot_j_out_vec_p2 \
                - 441.0*f11*e_in_vec_dot_j_out_vec_p4 + 42.0*f11*f10*j_in_vec_dot_j_out_vec_p2 \
                - 21.0*f11*j_in_vec_dot_j_out_vec_p4);
    }
    if (include_dotriacontupole_order_binary_pair_terms == true)
    {
        h1 = (1.0 - 4.0*e_in_p2)*(8.0 + e_out_p2);
        h2 = 8.0 + 3.0*e_out_p2;
        h3 = -8.0 + e_out_p2 - 4.0*e_in_p4*(80.0 + 179.0*e_out_p2) + e_in_p2*(64.0 + 748.0*e_out_p2);
        h4 = -8.0 - 19.0*e_out_p2 + 6.0*e_in_p2*(16.0 + 5.0*e_out_p2);
        h5 = 8.0 + e_out_p2 - 2.0*e_in_p2*(16.0 + 29.0*e_out_p2); 
        
        g1 = (-26.0 + 15.0*e_in_p2)*j_out_p2 + 18.0*j_in_vec_dot_j_out_vec_p2 + 99.0*e_in_vec_dot_j_out_vec_p2;
        g2 = h1*j_out_p2 + 9.0*h2*e_in_vec_dot_j_out_vec_p2 + 6.0*j_out_p2*j_in_vec_dot_e_out_vec_p2 - 3.0*h2*j_in_vec_dot_j_out_vec_p2;
        g3 = h3*j_out_p4 - 693.0*h2*e_in_vec_dot_j_out_vec_p4 + 42.0*e_in_vec_dot_j_out_vec_p2*(h4*j_out_p2 + 9.0*h2*j_in_vec_dot_j_out_vec_p2) \
            + 14.0*h5*j_in_vec_dot_j_out_vec_p2*j_out_p2 - 21.0*h2*j_in_vec_dot_j_out_vec_p4 \
            - 28.0*j_in_vec_dot_e_out_vec_p2*j_out_p2*( (1.0 + 23.0*e_in_p2)*j_out_p2 - 9.0*j_in_vec_dot_j_out_vec_p2 );
        
        g = -3024.0*e_in_vec_dot_e_out_vec_p2*e_in_vec_dot_j_out_vec*j_in_vec_dot_e_out_vec*j_in_vec_dot_j_out_vec*j_out_p2 \
            + 28.0*j_out_p2*e_in_vec_dot_e_out_vec_p3*g1 + 28.0*e_in_vec_dot_j_out_vec*j_in_vec_dot_e_out_vec*j_in_vec_dot_j_out_vec*g2 \
            + e_in_vec_dot_e_out_vec*g3;

    }
    
    double binary_pair_hamiltonian = A_quad*j_out_p5_inv*f1 - A_oct*j_out_p7_inv*( e_in_vec_dot_e_out_vec*f2 + f3 ) \
        + A_hd*j_out_p11_inv*f12 + A_tc*j_out_p13_inv*g;


    *hamiltonian = binary_pair_hamiltonian;
    *KS_V = binary_pair_hamiltonian;

    if (compute_hamiltonian_only == true)
    {
        return;
    }

    
    /****************************************
     * compute gradients of the Hamiltonian *
     ***************************************/    
    double grad_e_in_vec_phi[3],    grad_j_in_vec_phi[3];
    double grad_e_out_vec_phi[3],   grad_j_out_vec_phi[3];
    
    double grad_j_in_vec_f1[3],     grad_j_in_vec_f2[3],        grad_j_in_vec_f3[3];
    double grad_j_out_vec_f1[3],    grad_j_out_vec_f2[3],       grad_j_out_vec_f3[3];    
    double grad_e_in_vec_f1[3],     grad_e_in_vec_f2[3],        grad_e_in_vec_f3[3];
    double grad_e_out_vec_f3[3];

    /* dotriacontadupole */
    double grad_j_in_vec_g1[3],     grad_j_in_vec_g2[3],        grad_j_in_vec_g3[3];
    double grad_j_out_vec_g1[3],    grad_j_out_vec_g2[3],       grad_j_out_vec_g3[3];    
    double grad_e_in_vec_g1[3],     grad_e_in_vec_g2[3],        grad_e_in_vec_g3[3];
    double grad_e_out_vec_g2[3],    grad_e_out_vec_g3[3];


    for (int i=0; i<3; i++)
    {
        /* separate terms occurring in the gradients */
        if (include_quadrupole_order_terms == true)
        {
            grad_j_in_vec_f1[i] = -6.0*j_in_vec_dot_j_out_vec*j_out_vec[i];
            grad_j_out_vec_f1[i] = -6.0*j_in_vec_dot_j_out_vec*j_in_vec[i] + 2.0*(1.0-6.0*e_in_p2)*j_out_vec[i] \
                + 30.0*e_in_vec_dot_j_out_vec*e_in_vec[i];
            grad_e_in_vec_f1[i] = -12.0*j_out_p2*e_in_vec[i] + 30.0*e_in_vec_dot_j_out_vec*j_out_vec[i];
        }
        if (include_octupole_order_binary_pair_terms == true)
        {
            grad_j_in_vec_f2[i] = -10.0*j_in_vec_dot_j_out_vec*j_out_vec[i];
            grad_j_in_vec_f3[i] = -10.0*e_in_vec_dot_j_out_vec*( j_in_vec_dot_j_out_vec*e_out_vec[i] + j_in_vec_dot_e_out_vec*j_out_vec[i] );
            grad_j_out_vec_f2[i] = -10.0*j_in_vec_dot_j_out_vec*j_in_vec[i] + 2.0*(1.0-8.0*e_in_p2)*j_out_vec[i] \
                + 70.0*e_in_vec_dot_j_out_vec*e_in_vec[i];
            grad_j_out_vec_f3[i] = -10.0*j_in_vec_dot_e_out_vec*( j_in_vec_dot_j_out_vec*e_in_vec[i] + e_in_vec_dot_j_out_vec*j_in_vec[i] );
            grad_e_in_vec_f2[i] = -16.0*j_out_p2*e_in_vec[i] + 70.0*e_in_vec_dot_j_out_vec*j_out_vec[i];
            grad_e_in_vec_f3[i] = -10.0*j_in_vec_dot_e_out_vec*j_in_vec_dot_j_out_vec*j_out_vec[i];
            grad_e_out_vec_f3[i] = -10.0*e_in_vec_dot_j_out_vec*j_in_vec_dot_j_out_vec*j_in_vec[i];
        }
        if (include_dotriacontupole_order_binary_pair_terms == true)
        {
            grad_j_in_vec_g1[i] = 36.0*j_in_vec_dot_j_out_vec*j_out_vec[i];
            grad_j_in_vec_g2[i] = 12.0*j_out_p2*j_in_vec_dot_e_out_vec*e_out_vec[i] - 6.0*h2*j_in_vec_dot_j_out_vec*j_out_vec[i];
            grad_j_in_vec_g3[i] = 756.0*e_in_vec_dot_j_out_vec_p2*h2*j_in_vec_dot_j_out_vec*j_out_vec[i] + 28.0*h5*j_in_vec_dot_j_out_vec*j_out_p2*j_out_vec[i] \
                - 84.0*h2*j_in_vec_dot_j_out_vec_p3*j_out_vec[i] - 56.0*j_in_vec_dot_e_out_vec*j_out_p2*((1.0 + 23.0*e_in_p2)*j_out_p2 - 9.0*j_in_vec_dot_j_out_vec_p2)*e_out_vec[i] \
                + 504.0*j_in_vec_dot_e_out_vec_p2*j_out_p2*j_in_vec_dot_j_out_vec*j_out_vec[i];

            grad_j_out_vec_g1[i] = 2.0*(-26.0 + 15.0*e_in_p2)*j_out_vec[i] + 36.0*j_in_vec_dot_j_out_vec*j_in_vec[i] + 198.0*e_in_vec_dot_j_out_vec*e_in_vec[i];
            grad_j_out_vec_g2[i] = 2.0*h1*j_out_vec[i] + 18.0*h2*e_in_vec_dot_j_out_vec*e_in_vec[i] + 12.0*j_in_vec_dot_e_out_vec_p2*j_out_vec[i] \
                - 6.0*h2*j_in_vec_dot_j_out_vec*j_in_vec[i];
            grad_j_out_vec_g3[i] = 4.0*h3*j_out_p2*j_out_vec[i] - 2772.0*h2*e_in_vec_dot_j_out_vec_p3*e_in_vec[i] \
                + 84.0*e_in_vec_dot_j_out_vec*( h4*j_out_p2 + 9.0*h2*j_in_vec_dot_j_out_vec_p2 )*e_in_vec[i] \
                + 42.0*e_in_vec_dot_j_out_vec_p2*( 2.0*h4*j_out_vec[i] + 18.0*h2*j_in_vec_dot_j_out_vec*j_in_vec[i] ) \
                + 28.0*h5*( j_out_p2*j_in_vec_dot_j_out_vec*j_in_vec[i] + j_in_vec_dot_j_out_vec_p2*j_out_vec[i] ) \
                - 84.0*h2*j_in_vec_dot_j_out_vec_p3*j_in_vec[i] \
                - 28.0*j_in_vec_dot_e_out_vec_p2*( 2.0*j_out_vec[i]*( (1.0 + 23.0*e_in_p2)*j_out_p2 - 9.0*j_in_vec_dot_j_out_vec_p2 ) \
                    + j_out_p2*( 2.0*(1.0 + 23.0*e_in_p2)*j_out_vec[i] - 18.0*j_in_vec_dot_j_out_vec*j_in_vec[i]) );

            grad_e_in_vec_g1[i] = 30.0*j_out_p2*e_in_vec[i] + 198.0*e_in_vec_dot_j_out_vec*j_out_vec[i];
            grad_e_in_vec_g2[i] = -8.0*(8.0 + e_out_p2)*j_out_p2*e_in_vec[i] + 18.0*e_in_vec_dot_j_out_vec*h2*j_out_vec[i];
            grad_e_in_vec_g3[i] = j_out_p4*( -16.0*e_in_p2*(80.0 + 179.0*e_out_p2) + 2.0*(64.0 + 748.0*e_out_p2) )*e_in_vec[i] \
                - 2772.0*h2*e_in_vec_dot_j_out_vec_p3*j_out_vec[i] + 84.0*e_in_vec_dot_j_out_vec*( h4*j_out_p2 + 9.0*h2*j_in_vec_dot_j_out_vec_p2 )*j_out_vec[i] \
                + 504.0*e_in_vec_dot_j_out_vec_p2*j_out_p2*(16.0 + 5.0*e_out_p2)*e_in_vec[i] - 56.0*j_in_vec_dot_j_out_vec_p2*j_out_p2*(16.0 + 29.0*e_out_p2)*e_in_vec[i] \
                - 1288.0*j_in_vec_dot_e_out_vec_p2*j_out_p4*e_in_vec[i];
        
            grad_e_out_vec_g2[i] = 2.0*(1.0 - 4.0*e_in_p2)*j_out_p2*e_out_vec[i] + 54.0*e_in_vec_dot_j_out_vec_p2*e_out_vec[i] \
                + 12.0*j_out_p2*j_in_vec_dot_e_out_vec*j_in_vec[i] - 18.0*j_in_vec_dot_j_out_vec_p2*e_out_vec[i];
            grad_e_out_vec_g3[i] = j_out_p4*( 2.0 + 1496.0*e_in_p2 -1432.0*e_in_p4 )*e_out_vec[i] - 4158.0*e_in_vec_dot_j_out_vec_p4*e_out_vec[i] \
                + 42.0*e_in_vec_dot_j_out_vec_p2*( (-38.0 + 60.0*e_in_p2)*j_out_p2 + 54.0*j_in_vec_dot_j_out_vec_p2 )*e_out_vec[i] \
                + 14.0*j_in_vec_dot_j_out_vec_p2*j_out_p2*(2.0 - 116.0*e_in_p2)*e_out_vec[i] - 126.0*j_in_vec_dot_j_out_vec_p4*e_out_vec[i] \
                - 56.0*j_in_vec_dot_e_out_vec*j_out_p2*( (1.0 + 23.0*e_in_p2)*j_out_p2 - 9.0*j_in_vec_dot_j_out_vec_p2 )*j_in_vec[i];
        }
            
        /* complete gradients */
        grad_j_in_vec_phi[i] = 0.0;
        grad_j_out_vec_phi[i] = 0.0;
        grad_e_in_vec_phi[i] = 0.0;
        grad_e_out_vec_phi[i] = 0.0;
        
        if (include_quadrupole_order_terms == true)
        {
            grad_j_in_vec_phi[i] += A_quad*j_out_p5_inv*grad_j_in_vec_f1[i];
            grad_j_out_vec_phi[i] += -5.0*A_quad*j_out_p7_inv*j_out_vec[i]*f1 + A_quad*j_out_p5_inv*grad_j_out_vec_f1[i];
            grad_e_in_vec_phi[i] += A_quad*j_out_p5_inv*grad_e_in_vec_f1[i];
        }    
        if (include_octupole_order_binary_pair_terms == true)
        {
            grad_j_in_vec_phi[i] += -A_oct*j_out_p7_inv*( e_in_vec_dot_e_out_vec*grad_j_in_vec_f2[i] + grad_j_in_vec_f3[i] );
            grad_j_out_vec_phi[i] += 7.0*A_oct*j_out_p9_inv*j_out_vec[i]*( e_in_vec_dot_e_out_vec*f2 + f3 ) \
                - A_oct*j_out_p7_inv*( e_in_vec_dot_e_out_vec*grad_j_out_vec_f2[i] + grad_j_out_vec_f3[i] );
            grad_e_in_vec_phi[i] += -A_oct*j_out_p7_inv*( e_out_vec[i]*f2 + e_in_vec_dot_e_out_vec*grad_e_in_vec_f2[i] \
                + grad_e_in_vec_f3[i] );
            grad_e_out_vec_phi[i] += -A_oct*j_out_p7_inv*( e_in_vec[i]*f2 + grad_e_out_vec_f3[i] );
        }
        if (include_hexadecupole_order_binary_pair_terms == true)
        {
            grad_j_in_vec_phi[i] += A_hd*j_out_p11_inv*( \
                - 5880.0*j_out_p2*e_in_vec_dot_e_out_vec*e_in_vec_dot_j_out_vec*(j_in_vec_dot_j_out_vec*e_out_vec[i] + j_in_vec_dot_e_out_vec*j_out_vec[i]) \
                + 5.0*( -12.0*j_out_p2*f7*j_in_vec_dot_j_out_vec*j_out_vec[i] - 12.0*j_out_p2*(2.0*f8*j_in_vec_dot_e_out_vec*e_out_vec[i] \
                        - 14.0*j_in_vec_dot_e_out_vec_p2*j_in_vec_dot_j_out_vec*j_out_vec[i]) \
                    + 84.0*f11*f10*j_in_vec_dot_j_out_vec*j_out_vec[i] \
                    - 84.0*f11*j_in_vec_dot_j_out_vec_p2*j_in_vec_dot_j_out_vec*j_out_vec[i] ) \
                );
            grad_j_out_vec_phi[i] += -11.0*A_hd*j_out_p11_inv*j_out_p2_inv*f12*j_out_vec[i] \
                + A_hd*j_out_p11_inv*(12.0*f4*j_out_p2*j_out_vec[i] \
                    + 420.0*e_in_vec_dot_e_out_vec_p2*(2.0*f5*j_out_vec[i] + j_out_p2*(-4.0*j_out_vec[i] - 2.0*e_in_p2*j_out_vec[i] + 42.0*e_in_vec_dot_j_out_vec*e_in_vec[i])) \
                    - 5880.0*e_in_vec_dot_e_out_vec*j_in_vec_dot_e_out_vec*(2.0*e_in_vec_dot_j_out_vec*j_in_vec_dot_j_out_vec*j_out_vec[i] \
                        + j_out_p2*j_in_vec_dot_j_out_vec*e_in_vec[i] + j_out_p2*e_in_vec_dot_j_out_vec*j_in_vec[i]) \
                    + 5.0*( \
                        + 56.0*f6*(e_in_vec_dot_j_out_vec_p2*j_out_vec[i] + j_out_p2*e_in_vec_dot_j_out_vec*e_in_vec[i]) \
                        - 12.0*f7*(j_in_vec_dot_j_out_vec_p2*j_out_vec[i] + j_out_p2*j_in_vec_dot_j_out_vec*j_in_vec[i]) \
                        - 12.0*j_in_vec_dot_e_out_vec_p2*(2.0*f8*j_out_vec[i] + j_out_p2*(2.0*j_out_vec[i] + 26.0*e_in_p2*j_out_vec[i] - 14.0*j_in_vec_dot_j_out_vec*j_in_vec[i]) ) \
                        + 196.0*f9*(e_in_vec_dot_j_out_vec_p2*j_out_vec[i] + j_out_p2*e_in_vec_dot_j_out_vec*e_in_vec[i]) \
                        - 1764.0*f11*e_in_vec_dot_j_out_vec_p2*e_in_vec_dot_j_out_vec*e_in_vec[i] \
                        + 42.0*f11*( j_in_vec_dot_j_out_vec_p2*(2.0*j_out_vec[i] - 2.0*e_in_p2*j_out_vec[i] + 14.0*e_in_vec_dot_j_out_vec*e_in_vec[i]) \
                            + 2.0*f10*j_in_vec_dot_j_out_vec*j_in_vec[i] ) - 84.0*f11*j_in_vec_dot_j_out_vec_p2*j_in_vec_dot_j_out_vec*j_in_vec[i] ) \
                );
            grad_e_in_vec_phi[i] += A_hd*j_out_p11_inv*( \
                + 240.0*j_out_p4*(1.0 + 8.0*e_out_p2 - e_in_p2*(8.0 + 15.0*e_out_p2))*e_in_vec[i] + 840.0*e_in_vec_dot_e_out_vec*j_out_p2*f5*e_out_vec[i] \
                + 420.0*e_in_vec_dot_e_out_vec_p2*j_out_p2*(-2.0*j_out_p2*e_in_vec[i] + 42.0*e_in_vec_dot_j_out_vec*j_out_vec[i]) \
                - 5880.0*j_out_p2*j_in_vec_dot_e_out_vec*j_in_vec_dot_j_out_vec*( e_in_vec_dot_j_out_vec*e_out_vec[i] + e_in_vec_dot_e_out_vec*j_out_vec[i] ) \
                + 5.0*( \
                    + 28.0*j_out_p2*(4.0 + 3.0*e_out_p2)*( -20.0*e_in_vec_dot_j_out_vec_p2*e_in_vec[i] + 2.0*(1.0 - 10.0*e_in_p2)*e_in_vec_dot_j_out_vec*j_out_vec[i] ) \
                    - 12.0*j_out_p2*(6.0 + 29.0*e_out_p2)*j_in_vec_dot_j_out_vec_p2*e_in_vec[i] - 312.0*j_out_p4*j_in_vec_dot_e_out_vec_p2*e_in_vec[i] \
                    + 98.0*j_out_p2*(8.0*e_in_vec_dot_j_out_vec_p2*(5.0 + 3.0*e_out_p2)*e_in_vec[i] + 2.0*e_in_vec_dot_j_out_vec*f9*j_out_vec[i]) \
                    - 1764.0*f11*e_in_vec_dot_j_out_vec_p2*e_in_vec_dot_j_out_vec*j_out_vec[i] \
                    + 42.0*f11*j_in_vec_dot_j_out_vec_p2*(-2.0*j_out_p2*e_in_vec[i] + 14.0*e_in_vec_dot_j_out_vec*j_out_vec[i]) ) \
                );
            grad_e_out_vec_phi[i] += A_hd*j_out_p11_inv*( \
                + 6.0*j_out_p4*(1.0 + 320.0*e_in_p2 - 300.0*e_in_p4)*e_out_vec[i] + 840.0*e_in_vec_dot_e_out_vec*j_out_p2*f5*e_in_vec[i] \
                - 5880.0*j_out_p2*e_in_vec_dot_j_out_vec*j_in_vec_dot_j_out_vec*(j_in_vec_dot_e_out_vec*e_in_vec[i] + e_in_vec_dot_e_out_vec*j_in_vec[i]) \
                + 5.0*( \
                    + 168.0*j_out_p2*(1.0 - 10.0*e_in_p2)*e_in_vec_dot_j_out_vec_p2*e_out_vec[i] - 6.0*j_out_p2*j_in_vec_dot_j_out_vec_p2*(12.0 + 58.0*e_in_p2)*e_out_vec[i] \
                    - 24.0*j_out_p2*f8*j_in_vec_dot_e_out_vec*j_in_vec[i] + 98.0*j_out_p2*e_in_vec_dot_j_out_vec_p2*(-6.0 + 24.0*e_in_p2)*e_out_vec[i] \
                    - 882.0*e_in_vec_dot_j_out_vec_p4*e_out_vec[i] + 84.0*f10*j_in_vec_dot_j_out_vec_p2*e_out_vec[i] - 42.0*j_in_vec_dot_j_out_vec_p4*e_out_vec[i]) \
                );
        }
        if (include_dotriacontupole_order_binary_pair_terms == true)
        {
            grad_j_in_vec_phi[i] += A_tc*j_out_p13_inv*( \
                - 3024.0*e_in_vec_dot_e_out_vec_p2*e_in_vec_dot_j_out_vec*j_out_p2*( j_in_vec_dot_j_out_vec*e_out_vec[i] + j_in_vec_dot_e_out_vec*j_out_vec[i] ) \
                + 28.0*j_out_p2*e_in_vec_dot_e_out_vec_p3*grad_j_in_vec_g1[i] + 28.0*e_in_vec_dot_j_out_vec*( j_in_vec_dot_j_out_vec*g2*e_out_vec[i] \
                    + j_in_vec_dot_e_out_vec*g2*j_out_vec[i] + j_in_vec_dot_e_out_vec*j_in_vec_dot_j_out_vec*grad_j_in_vec_g2[i] ) \
                    + e_in_vec_dot_e_out_vec*grad_j_in_vec_g3[i] );
            grad_j_out_vec_phi[i] += -13.0*A_tc*j_out_p13_inv*j_out_p2_inv*g*j_out_vec[i] + A_tc*j_out_p13_inv*( \
                - 3024.0*e_in_vec_dot_e_out_vec_p2*j_in_vec_dot_e_out_vec*( j_in_vec_dot_j_out_vec*j_out_p2*e_in_vec[i] + e_in_vec_dot_j_out_vec*j_out_p2*j_in_vec[i] \
                    + 2.0*e_in_vec_dot_j_out_vec*j_in_vec_dot_j_out_vec*j_out_vec[i] ) + 28.0*e_in_vec_dot_e_out_vec_p3*( 2.0*g1*j_out_vec[i] + j_out_p2*grad_j_out_vec_g1[i] ) \
                    + 28.0*j_in_vec_dot_e_out_vec*( j_in_vec_dot_j_out_vec*g2*e_in_vec[i] + e_in_vec_dot_j_out_vec*g2*j_in_vec[i] \
                        + e_in_vec_dot_j_out_vec*j_in_vec_dot_j_out_vec*grad_j_out_vec_g2[i] ) + e_in_vec_dot_e_out_vec*grad_j_out_vec_g3[i] );
            grad_e_in_vec_phi[i] += A_tc*j_out_p13_inv*( \
                - 3024.0*j_in_vec_dot_e_out_vec*j_in_vec_dot_j_out_vec*j_out_p2*( 2.0*e_in_vec_dot_e_out_vec*e_in_vec_dot_j_out_vec*e_out_vec[i] \
                    + e_in_vec_dot_e_out_vec_p2*j_out_vec[i] ) + 84.0*j_out_p2*e_in_vec_dot_e_out_vec_p2*g1*e_out_vec[i] \
                + 28.0*j_out_p2*e_in_vec_dot_e_out_vec_p3*grad_e_in_vec_g1[i] + 28.0*j_in_vec_dot_e_out_vec*j_in_vec_dot_j_out_vec*( g2*j_out_vec[i] \
                    + e_in_vec_dot_j_out_vec*grad_e_in_vec_g2[i] ) + g3*e_out_vec[i] + e_in_vec_dot_e_out_vec*grad_e_in_vec_g3[i] );
            grad_e_out_vec_phi[i] += A_tc*j_out_p13_inv*( \
                - 3024.0*e_in_vec_dot_j_out_vec*j_in_vec_dot_j_out_vec*j_out_p2*( 2.0*j_in_vec_dot_e_out_vec*e_in_vec_dot_e_out_vec*e_in_vec[i] \
                    + e_in_vec_dot_e_out_vec_p2*j_in_vec[i] ) + 84.0*j_out_p2*e_in_vec_dot_e_out_vec_p2*g1*e_in_vec[i] \
                + 28.0*e_in_vec_dot_j_out_vec*j_in_vec_dot_j_out_vec*( g2*j_in_vec[i] + j_in_vec_dot_e_out_vec*grad_e_out_vec_g2[i] ) \
                + g3*e_in_vec[i] + e_in_vec_dot_e_out_vec*grad_e_out_vec_g3[i] );
        }
    }
    
    double j_in_vec_cross_grad_j_in_vec_phi[3],             j_in_vec_cross_grad_e_in_vec_phi[3];
    double j_out_vec_cross_grad_j_out_vec_phi[3],           j_out_vec_cross_grad_e_out_vec_phi[3];
    
    double e_in_vec_cross_grad_e_in_vec_phi[3],             e_in_vec_cross_grad_j_in_vec_phi[3];
    double e_out_vec_cross_grad_e_out_vec_phi[3],           e_out_vec_cross_grad_j_out_vec_phi[3];
    

    cross3(j_in_vec,        grad_j_in_vec_phi,              j_in_vec_cross_grad_j_in_vec_phi);
    cross3(j_in_vec,        grad_e_in_vec_phi,              j_in_vec_cross_grad_e_in_vec_phi);
    cross3(e_in_vec,        grad_e_in_vec_phi,              e_in_vec_cross_grad_e_in_vec_phi);
    cross3(e_in_vec,        grad_j_in_vec_phi,              e_in_vec_cross_grad_j_in_vec_phi);
    
    cross3(j_out_vec,       grad_j_out_vec_phi,             j_out_vec_cross_grad_j_out_vec_phi);
    cross3(j_out_vec,       grad_e_out_vec_phi,             j_out_vec_cross_grad_e_out_vec_phi);
    cross3(e_out_vec,       grad_e_out_vec_phi,             e_out_vec_cross_grad_e_out_vec_phi);
    cross3(e_out_vec,       grad_j_out_vec_phi,             e_out_vec_cross_grad_j_out_vec_phi);
    


    for (int i=0; i<3; i++)
    {
        inner_binary->de_vec_dt[i] += (-1.0/(Lambda_in))*( e_in_vec_cross_grad_j_in_vec_phi[i] \
            + j_in_vec_cross_grad_e_in_vec_phi[i] );
        inner_binary->dh_vec_dt[i] += -1.0*( j_in_vec_cross_grad_j_in_vec_phi[i] \
            + e_in_vec_cross_grad_e_in_vec_phi[i] );

        outer_binary->de_vec_dt[i] += (-1.0/(Lambda_out))*( e_out_vec_cross_grad_j_out_vec_phi[i] \
            + j_out_vec_cross_grad_e_out_vec_phi[i] );
        outer_binary->dh_vec_dt[i] += -1.0*( j_out_vec_cross_grad_j_out_vec_phi[i] \
            + e_out_vec_cross_grad_e_out_vec_phi[i] );  
        
    }

    /* DA corrections */
    if (include_double_averaging_corrections==true)
    {
        /* Based on eqs of https://ui.adsabs.harvard.edu/abs/2016MNRAS.458.3060L/abstract (test particle limit). 
         * The latter paper assumes a particular (fixed) orientation of the perturber: j along z, e along x. This is not
         * the case here, so e & j of the inner orbit need to be projected onto e_out and j_out. After computing the EOM for
         * the projected e_in and j_in, these quantities need to be transformed back into the original frame. */
        
        double b_per = a_out * j_out;
        double t_sec = sqrt(m1_plus_m2) * b_per * b_per * b_per / ( sqrt(CONST_G * a_in * a_in * a_in) * m3 );
        double eps_SA = pow(a_in/a_out,3.0/2.0) * m3 * (1.0/sqrt( (m1_plus_m2 + m3) * m1_plus_m2 ) ) * j_out_p3_inv;
        
        /* Project e & j vectors to outer orbital plane (new frame; primed) */
        double *x_hat_prime = e_out_vec_unit;
        double *y_hat_prime = q_out_vec_unit;
        double *z_hat_prime = h_out_vec_unit;
        
        double ex = dot3(e_in_vec,x_hat_prime);
        double ey = dot3(e_in_vec,y_hat_prime);
        double ez = dot3(e_in_vec,z_hat_prime);
        double jx = dot3(j_in_vec,x_hat_prime);
        double jy = dot3(j_in_vec,y_hat_prime);
        double jz = dot3(j_in_vec,z_hat_prime);

        double ex_p2 = ex * ex;
        double ey_p2 = ey * ey;
        double ez_p2 = ez * ez;
        double jx_p2 = jx * jx;
        double jy_p2 = jy * jy;
        double jz_p2 = jz * jz;
        
        /* Compute CDA EOM */
        double c0 = (c_1div3 + 8.0*ex_p2 + 8.0*ey_p2 + 3.0*ez_p2);
        double cj = c0 - jz_p2;
        double ce = c0 - 17.0*jz_p2;
        
        double dj_dtau_CDA[3];
        double de_dtau_CDA[3];
        
        dj_dtau_CDA[0] =  c_27div64*( -10.0*ey*ez*jz + jy*cj ) + e_out_p2*( -c_9div64*(10.0*ey*ez*jz + jy*(-c_2div3 - 21.0*ex_p2 + 9.0*ey_p2 - 16.0*ez_p2 - jx_p2 + jy_p2)) );
        dj_dtau_CDA[1] = -c_27div64*( -10.0*ex*ez*jz + jx*cj ) + e_out_p2*( c_9div64 * (30.0*ex*ez*jz + jx*(c_10div3 - 45.0*ex_p2 - 15.0*ey_p2 - 5.0*jx_p2 - 3.0*jy_p2)) ); 
        dj_dtau_CDA[2] = -e_out_p2 * c_9div16*( 5.0*ey*ez*jx + 5.0*ex*ez*jy + 5.0*ex*ey*jz + jx*jy*jz );
        
        de_dtau_CDA[0] =  c_27div64*( 6.0*ez*jy*jz + ey*ce) + e_out_p2*(  c_9div64*( 14.0*ez*jy*jz + ey*( c_35div3 + 10.0*ex_p2 + 5.0*ez_p2 - 10.0*jx_p2 - 32.0*jy_p2 - 35.0*jz_p2) ) );
        de_dtau_CDA[1] = -c_27div64*( 6.0*ez*jx*jz + ex*ce) + e_out_p2*( -c_9div64*( 10.0*ez*jx*jz + ex*( c_64div3 - 10.0*ey_p2 - 25.0*ez_p2 - 22.0*jy_p2 - 65.0*jz_p2) ) );
        de_dtau_CDA[2] =  c_27div4*(ey*jx - ex*jy)*jz + e_out_p2*( -c_9div16*( 5.0*ex*ey*ez + 5.0*ez*jx*jy - 5.0*ey*jx*jz + 11.0*ex*jy*jz) );
        
        double de_dt_prime[3];
        double dh_dt_prime[3];
        
        for (int i=0; i<3; i++)
        {
            de_dt_prime[i] = eps_SA * de_dtau_CDA[i] / t_sec;
            dh_dt_prime[i] = Lambda_in * eps_SA * dj_dtau_CDA[i] / t_sec;
        }

        /* Determine time derivatives of x_hat_prime, etc. */
        double dx_hat_prime_dt[3];
        double dy_hat_prime_dt[3];
        double dz_hat_prime_dt[3];
        
        double cx = dot3(x_hat_prime,outer_binary->de_vec_dt);
        double cz = dot3(z_hat_prime,outer_binary->dh_vec_dt);
        
        for (int i=0; i<3; i++)
        {
            dx_hat_prime_dt[i] = (1.0/e_out)*( outer_binary->de_vec_dt[i] - cx * x_hat_prime[i] );
            dz_hat_prime_dt[i] = (1.0/h_out)*( outer_binary->dh_vec_dt[i] - cz * z_hat_prime[i] );
        }
        
        double cy1[3],cy2[3];
        cross3(dz_hat_prime_dt,x_hat_prime,cy1);
        cross3(z_hat_prime,dx_hat_prime_dt,cy2);

        for (int i=0; i<3; i++)
        {
            dy_hat_prime_dt[i] = cy1[i] + cy2[i];
        }
        
        /* Convert EOM of primed e & h vectors to unprimed ones (in original frame) */
        double cedotx = de_dt_prime[0] - dot3(e_in_vec,dx_hat_prime_dt);
        double cedoty = de_dt_prime[1] - dot3(e_in_vec,dy_hat_prime_dt);
        double cedotz = de_dt_prime[2] - dot3(e_in_vec,dz_hat_prime_dt);
        double chdotx = dh_dt_prime[0] - dot3(h_in_vec,dx_hat_prime_dt);
        double chdoty = dh_dt_prime[1] - dot3(h_in_vec,dy_hat_prime_dt);
        double chdotz = dh_dt_prime[2] - dot3(h_in_vec,dz_hat_prime_dt);

        for (int i=0; i<3; i++)
        {
            inner_binary->de_vec_dt[i] += cedotx * x_hat_prime[i] + cedoty * y_hat_prime[i] + cedotz * z_hat_prime[i];
            inner_binary->dh_vec_dt[i] += chdotx * x_hat_prime[i] + chdoty * y_hat_prime[i] + chdotz * z_hat_prime[i];
            //inner_binary->dh_vec_dt[i] += Lambda_in * eps_SA * dj_dtau_CDA[i] / t_sec;
            //printf("CDA %g %g\n",eps_SA * de_dtau_CDA[i] / t_sec,Lambda_in * eps_SA * dj_dtau_CDA[i] / t_sec);
        }
        
        //printf("test %g %g\n",e_in_vec[0],dot3(e_in_vec,x_hat_prime)*x_hat_prime[0] + dot3(e_in_vec,y_hat_prime)*y_hat_prime[0] + dot3(e_in_vec,z_hat_prime)*z_hat_prime[0]);
    }
    

    #ifdef VERBOSE
    if (verbose_flag > 2)
    {
        printf("ODE_newtonian.cpp -- compute_EOM_binary_pairs_double_averaged -- begin print\n");
        printf("e_in %g %g %g\n",e_in_vec[0],e_in_vec[1],e_in_vec[2]);
        printf("e_out %g %g %g\n",e_out_vec[0],e_out_vec[1],e_out_vec[2]);    
        printf("h_in %g %g %g\n",h_in_vec[0],h_in_vec[1],h_in_vec[2]);
        printf("h_out %g %g %g\n",h_out_vec[0],h_out_vec[1],h_out_vec[2]);    

        printf("grad1 %g %g %g\n",grad_e_in_vec_f1[0],grad_e_in_vec_f1[1],grad_e_in_vec_f1[2]);
        printf("grad2 %g %g %g\n",grad_e_in_vec_f2[0],grad_e_in_vec_f2[1],grad_e_in_vec_f2[2]);    
        printf("grad3 %g %g %g\n",grad_e_in_vec_f3[0],grad_e_in_vec_f3[1],grad_e_in_vec_f3[2]);    

        printf("de_in_dt %g %g %g\n",inner_binary->de_vec_dt[0],inner_binary->de_vec_dt[1],inner_binary->de_vec_dt[2]);
        printf("de_out_dt %g %g %g\n",outer_binary->de_vec_dt[0],outer_binary->de_vec_dt[1],outer_binary->de_vec_dt[2]);    
        printf("dh_in_dt %g %g %g\n",inner_binary->dh_vec_dt[0],inner_binary->dh_vec_dt[1],inner_binary->dh_vec_dt[2]);
        printf("dh_out_dt %g %g %g\n",outer_binary->dh_vec_dt[0],outer_binary->dh_vec_dt[1],outer_binary->dh_vec_dt[2]); 
        printf("ODE_newtonian.cpp -- compute_EOM_binary_pairs_double_averaged -- end print\n");
    }
    #endif
}



void compute_EOM_binary_triplets(ParticlesMap *particlesMap, int binary_A_index, int binary_B_index, int binary_C_index, int connecting_child_in_binary_B_to_binary_A, int connecting_child_in_binary_C_to_binary_B, double *hamiltonian, double *KS_V, bool compute_hamiltonian_only)
{
    /* last checked/updated 06-03-20 */
    
    if (include_octupole_order_binary_triplet_terms == false)
    {
        return;
    }

    /*********************
     * preamble          *
     ********************/

    Particle *binary_A = (*particlesMap)[binary_A_index];
    Particle *binary_B = (*particlesMap)[binary_B_index];
    Particle *binary_C = (*particlesMap)[binary_C_index];

    Particle *binary_A_child1 = (*particlesMap)[binary_A->child1];
    Particle *binary_A_child2 = (*particlesMap)[binary_A->child2];

    Particle *binary_B_child1 = (*particlesMap)[binary_B->child1];
    Particle *binary_B_child2 = (*particlesMap)[binary_B->child2];

    Particle *binary_C_child1 = (*particlesMap)[binary_C->child1];
    Particle *binary_C_child2 = (*particlesMap)[binary_C->child2];
   
   /* set alpha = +1 */
    double B_ijB = 0.0;
    
    if (connecting_child_in_binary_B_to_binary_A==1)
    {
        B_ijB = binary_B_child2->mass/binary_B->mass;
    }
    else if (connecting_child_in_binary_B_to_binary_A==2)
    {
        B_ijB = -binary_B_child1->mass/binary_B->mass;
    }

    double M_C_CS_B = 0.0;
    
    if (connecting_child_in_binary_C_to_binary_B==1)
    {
        M_C_CS_B = binary_C_child2->mass;
    }
    else if (connecting_child_in_binary_C_to_binary_B==2)
    {
        M_C_CS_B = binary_C_child1->mass;
    }

    double e_A = binary_A->e;
    double e_B = binary_B->e;
    double e_C = binary_C->e;
    double e_A_P2 = binary_A->e_p2;
    
    double *e_A_vec = binary_A->e_vec;
    double *e_B_vec = binary_B->e_vec;
    double *e_C_vec = binary_C->e_vec;
    
    double *h_A_vec = binary_A->h_vec;
    double *h_B_vec = binary_B->h_vec;
    double *h_C_vec = binary_C->h_vec;

    double *e_A_vec_unit = binary_A->e_vec_unit;
    double *e_B_vec_unit = binary_B->e_vec_unit;
    double *e_C_vec_unit = binary_C->e_vec_unit;
    
    double *h_A_vec_unit = binary_A->h_vec_unit;
    double *h_B_vec_unit = binary_B->h_vec_unit;
    double *h_C_vec_unit = binary_C->h_vec_unit;
    
    double *j_A_vec_unit = h_A_vec_unit;
    double *j_B_vec_unit = h_B_vec_unit;    
    double *j_C_vec_unit = h_C_vec_unit;    

    double h_A = binary_A->h;
    double h_B = binary_B->h;
    double h_C = binary_C->h;
    
    double j_A = binary_A->j;
    double j_A_p2 = binary_A->j_p2;
    double j_B = binary_B->j;

    double j_C = binary_C->j;
    double j_C_p2 = binary_C->j_p2;
    double j_C_p4 = binary_C->j_p4;
    double j_C_p7 = j_C*j_C_p2*j_C_p4;
    double j_C_p9 = j_C_p7*j_C_p2;
        
    double j_C_p7_inv = 1.0/j_C_p7;
    double j_C_p9_inv = 1.0/j_C_p9;
    
    double j_A_vec[3],j_B_vec[3],j_C_vec[3];
    for (int i=0; i<3; i++)
    {
        j_A_vec[i] = j_A*h_A_vec_unit[i];
        j_B_vec[i] = j_B*h_B_vec_unit[i];
        j_C_vec[i] = j_C*h_C_vec_unit[i];
    }
    
    double a_A = binary_A->a;
    double a_B = binary_B->a;
    double a_C = binary_C->a;    
    
    double mu_A = binary_A_child1->mass * binary_A_child2->mass / binary_A->mass;
    double mu_B = binary_B_child1->mass * binary_B_child2->mass / binary_B->mass;
    double mu_C = binary_C_child1->mass * binary_C_child2->mass / binary_C->mass;

    double Lambda_A = h_A/j_A;
    double Lambda_B = h_B/j_B;
    double Lambda_C = h_C/j_C;
   
    double C_tr = mu_A * B_ijB * M_C_CS_B;
    double H;
    int i;

    if (binary_A->integration_method==0 && binary_B->integration_method==0 && binary_C->integration_method==0)
    {
        /* A, B, C averaged */
       
        #ifdef VERBOSE
        if (verbose_flag > 2)
        {
            printf("ODE_newtonian.cpp -- compute_EOM_binary_triplets -- A %d, B %d, C %d averaged\n",binary_A->index,binary_B->index,binary_C->index);
        }
        #endif
        
        double A_cross = -(c_9div32) * C_tr *(a_A*a_A*a_B/(a_C*a_C*a_C*a_C));

        double e_A_vec_dot_e_B_vec = dot3(e_A_vec,e_B_vec);
        double e_B_vec_dot_e_C_vec = dot3(e_B_vec,e_C_vec);
        double e_A_vec_dot_e_C_vec = dot3(e_A_vec,e_C_vec);

        double e_A_vec_dot_j_C_vec = dot3(e_A_vec,j_C_vec);
        double e_B_vec_dot_j_C_vec = dot3(e_B_vec,j_C_vec);
        double e_B_vec_dot_j_A_vec = dot3(e_B_vec,j_A_vec);
        double e_C_vec_dot_j_A_vec = dot3(e_C_vec,j_A_vec);
        double j_A_vec_dot_j_C_vec = dot3(j_A_vec,j_C_vec);
        
        double e_A_vec_dot_j_C_vec_p2 = e_A_vec_dot_j_C_vec*e_A_vec_dot_j_C_vec;
        double j_A_vec_dot_j_C_vec_p2 = j_A_vec_dot_j_C_vec*j_A_vec_dot_j_C_vec;
        
        /***************************
         * compute the Hamiltonian *
         **************************/
        
        double f1 = j_C_p2*(1.0 - 6.0*e_A_P2) + 25.0*e_A_vec_dot_j_C_vec_p2 - 5.0*j_A_vec_dot_j_C_vec_p2;
        double f0 = -10.0*e_A_vec_dot_e_B_vec*e_A_vec_dot_e_C_vec*j_C_p2 + 50.0*e_A_vec_dot_e_C_vec*e_A_vec_dot_j_C_vec*e_B_vec_dot_j_C_vec \
            + 2.0*e_C_vec_dot_j_A_vec*e_B_vec_dot_j_A_vec*j_C_p2 - 10.0*e_B_vec_dot_j_C_vec*e_C_vec_dot_j_A_vec*j_A_vec_dot_j_C_vec \
            + e_B_vec_dot_e_C_vec*f1;
                
        double binary_triplet_hamiltonian = A_cross*j_C_p7_inv*f0;
        
        *hamiltonian += binary_triplet_hamiltonian;
        *KS_V += binary_triplet_hamiltonian;
        
        if (compute_hamiltonian_only == true)
        {
            return;
        }

        /****************************************
         * compute gradients of the Hamiltonian *
         ***************************************/
        double grad_e_A_vec_H[3],     grad_j_A_vec_H[3];
        double grad_e_B_vec_H[3],     grad_j_B_vec_H[3];
        double grad_e_C_vec_H[3],     grad_j_C_vec_H[3];    
        
        for (i=0; i<3; i++)
        {
            
            /* gradient w.r.t. e_A */
            grad_e_A_vec_H[i] = A_cross*j_C_p7_inv*( \
                - 10.0*j_C_p2*(e_A_vec_dot_e_C_vec*e_B_vec[i] + e_A_vec_dot_e_B_vec*e_C_vec[i]) \
                + 50.0*e_B_vec_dot_j_C_vec*(e_A_vec_dot_j_C_vec*e_C_vec[i] + e_A_vec_dot_e_C_vec*j_C_vec[i]) \
                + e_B_vec_dot_e_C_vec*(50.0*e_A_vec_dot_j_C_vec*j_C_vec[i] - 12.0*j_C_p2*e_A_vec[i]) );
            
            /* gradient w.r.t. j_A */
            grad_j_A_vec_H[i] = A_cross*j_C_p7_inv*( \
                + 2.0*e_B_vec_dot_j_A_vec*j_C_p2*e_C_vec[i] + 2.0*e_C_vec_dot_j_A_vec*j_C_p2*e_B_vec[i] \
                - 10.0*e_B_vec_dot_j_C_vec*(j_A_vec_dot_j_C_vec*e_C_vec[i] + e_C_vec_dot_j_A_vec*j_C_vec[i]) \
                - 10.0*e_B_vec_dot_e_C_vec*j_A_vec_dot_j_C_vec*j_C_vec[i] );
                
            /* gradient w.r.t. e_B */
            grad_e_B_vec_H[i] = A_cross*j_C_p7_inv*( \
                - 10.0*e_A_vec_dot_e_C_vec*j_C_p2*e_A_vec[i] + 50.0*e_A_vec_dot_e_C_vec*e_A_vec_dot_j_C_vec*j_C_vec[i] \
                + 2.0*e_C_vec_dot_j_A_vec*j_C_p2*j_A_vec[i] - 10.0*e_C_vec_dot_j_A_vec*j_A_vec_dot_j_C_vec*j_C_vec[i] \
                + f1*e_C_vec[i] );
                
            /* gradient w.r.t. j_B */
            grad_j_B_vec_H[i] = 0.0;

            /* gradient w.r.t. e_C */
            grad_e_C_vec_H[i] = A_cross*j_C_p7_inv*( \
                - 10.0*e_A_vec_dot_e_B_vec*j_C_p2*e_A_vec[i] + 50.0*e_A_vec_dot_j_C_vec*e_B_vec_dot_j_C_vec*e_A_vec[i] \
                + 2.0*e_B_vec_dot_j_A_vec*j_C_p2*j_A_vec[i] - 10.0*e_B_vec_dot_j_C_vec*j_A_vec_dot_j_C_vec*j_A_vec[i] \
                + f1*e_B_vec[i] );

            /* gradient w.r.t. j_C */
            grad_j_C_vec_H[i] = -7.0*A_cross*j_C_p9_inv*f0*j_C_vec[i] + A_cross*j_C_p7_inv*( \
                - 20.0*e_A_vec_dot_e_B_vec*e_A_vec_dot_e_C_vec*j_C_vec[i] \
                + 50.0*e_A_vec_dot_e_C_vec*(e_B_vec_dot_j_C_vec*e_A_vec[i] + e_A_vec_dot_j_C_vec*e_B_vec[i]) \
                + 4.0*e_C_vec_dot_j_A_vec*e_B_vec_dot_j_A_vec*j_C_vec[i] \
                - 10.0*e_C_vec_dot_j_A_vec*(j_A_vec_dot_j_C_vec*e_B_vec[i] + e_B_vec_dot_j_C_vec*j_A_vec[i]) \
                + e_B_vec_dot_e_C_vec*(2.0*(1.0 - 6.0*e_A_P2)*j_C_vec[i] + 50.0*e_A_vec_dot_j_C_vec*e_A_vec[i] \
                    - 10.0*j_A_vec_dot_j_C_vec*j_A_vec[i]) );
                
        }
        
        double j_A_vec_cross_grad_j_A_vec_H[3],                   j_A_vec_cross_grad_e_A_vec_H[3];
        double j_B_vec_cross_grad_j_B_vec_H[3],                   j_B_vec_cross_grad_e_B_vec_H[3];    
        double j_C_vec_cross_grad_j_C_vec_H[3],                   j_C_vec_cross_grad_e_C_vec_H[3];        
        
        double e_A_vec_cross_grad_e_A_vec_H[3],                   e_A_vec_cross_grad_j_A_vec_H[3];
        double e_B_vec_cross_grad_e_B_vec_H[3],                   e_B_vec_cross_grad_j_B_vec_H[3];
        double e_C_vec_cross_grad_e_C_vec_H[3],                   e_C_vec_cross_grad_j_C_vec_H[3];
        
        cross3(j_A_vec,             grad_j_A_vec_H,               j_A_vec_cross_grad_j_A_vec_H);
        cross3(j_A_vec,             grad_e_A_vec_H,               j_A_vec_cross_grad_e_A_vec_H);
        cross3(j_B_vec,             grad_j_B_vec_H,               j_B_vec_cross_grad_j_B_vec_H);
        cross3(j_B_vec,             grad_e_B_vec_H,               j_B_vec_cross_grad_e_B_vec_H);
        cross3(j_C_vec,             grad_j_C_vec_H,               j_C_vec_cross_grad_j_C_vec_H);
        cross3(j_C_vec,             grad_e_C_vec_H,               j_C_vec_cross_grad_e_C_vec_H);
        
        cross3(e_A_vec,             grad_e_A_vec_H,               e_A_vec_cross_grad_e_A_vec_H);
        cross3(e_A_vec,             grad_j_A_vec_H,               e_A_vec_cross_grad_j_A_vec_H);    
        cross3(e_B_vec,             grad_e_B_vec_H,               e_B_vec_cross_grad_e_B_vec_H);
        cross3(e_B_vec,             grad_j_B_vec_H,               e_B_vec_cross_grad_j_B_vec_H);    
        cross3(e_C_vec,             grad_e_C_vec_H,               e_C_vec_cross_grad_e_C_vec_H);
        cross3(e_C_vec,             grad_j_C_vec_H,               e_C_vec_cross_grad_j_C_vec_H);
        
        for (i=0; i<3; i++)
        {
            binary_A->de_vec_dt[i] += (-1.0/(Lambda_A))*( e_A_vec_cross_grad_j_A_vec_H[i] \
                + j_A_vec_cross_grad_e_A_vec_H[i] );
            binary_A->dh_vec_dt[i] += -1.0*( j_A_vec_cross_grad_j_A_vec_H[i] \
                + e_A_vec_cross_grad_e_A_vec_H[i] );

            binary_B->de_vec_dt[i] += (-1.0/(Lambda_B))*( e_B_vec_cross_grad_j_B_vec_H[i] \
                + j_B_vec_cross_grad_e_B_vec_H[i] );
            binary_B->dh_vec_dt[i] += -1.0*( j_B_vec_cross_grad_j_B_vec_H[i] \
                + e_B_vec_cross_grad_e_B_vec_H[i] );

            binary_C->de_vec_dt[i] += (-1.0/(Lambda_C))*( e_C_vec_cross_grad_j_C_vec_H[i] \
                + j_C_vec_cross_grad_e_C_vec_H[i] );
            binary_C->dh_vec_dt[i] += -1.0*( j_C_vec_cross_grad_j_C_vec_H[i] \
                + e_C_vec_cross_grad_e_C_vec_H[i] );
        }
    }
    else if (binary_A->integration_method==0 && binary_B->integration_method==0 && binary_C->integration_method>0)
    {
        /* A, B averaged; C direct */
        #ifdef VERBOSE
        if (verbose_flag > 2)
        {
            printf("ODE_newtonian.cpp -- compute_EOM_binary_triplets -- A %d B %d averaged, C %d direct\n",binary_A->index,binary_B->index,binary_C->index);
        }
        #endif

        double *r_C_vec = binary_C->r_vec;
        
        double r_C_P1 = norm3(r_C_vec);
        double r_C_P2 = r_C_P1 * r_C_P1;
        double r_C_Pminus_1 = 1.0/r_C_P1;
        double r_C_Pminus_2 = r_C_Pminus_1 * r_C_Pminus_1;
        double r_C_Pminus_4 = r_C_Pminus_2 * r_C_Pminus_2;
        double r_C_Pminus_5 = r_C_Pminus_1 * r_C_Pminus_4;
        double r_C_Pminus_7 = r_C_Pminus_2 * r_C_Pminus_5;
        double r_C_Pminus_9 = r_C_Pminus_2 * r_C_Pminus_7;
        
        double e_A_vec_dot_r_C_vec = dot3(e_A_vec,r_C_vec);
        double e_A_vec_dot_r_C_vec_P2 = e_A_vec_dot_r_C_vec * e_A_vec_dot_r_C_vec;        
        double j_A_vec_dot_r_C_vec = dot3(j_A_vec,r_C_vec);
        double j_A_vec_dot_r_C_vec_P2 = j_A_vec_dot_r_C_vec * j_A_vec_dot_r_C_vec;

        double e_A_vec_dot_e_B_vec = dot3(e_A_vec,e_B_vec);
        double e_B_vec_dot_r_C_vec = dot3(e_B_vec,r_C_vec);
        double j_A_vec_dot_e_B_vec = dot3(j_A_vec,e_B_vec);
        
        double C_V = -c_9div8 * C_tr * a_A * a_A * a_B;
        double f_V11 = ( (1.0 - e_A_P2) * r_C_P2 + 5.0 * e_A_vec_dot_r_C_vec_P2 - j_A_vec_dot_r_C_vec_P2 );
        double f_V1 = f_V11 * e_B_vec_dot_r_C_vec;
        double f_V2 = (1.0 - e_A_P2) * e_B_vec_dot_r_C_vec + 5.0 * e_A_vec_dot_e_B_vec * e_A_vec_dot_r_C_vec - j_A_vec_dot_e_B_vec * j_A_vec_dot_r_C_vec;
        double f_V31 = (2.0 + 3.0 * e_A_P2);
        double f_V3 = f_V31 * e_B_vec_dot_r_C_vec;
        
        H = C_V * ( f_V1 * r_C_Pminus_7 - 2.0 * f_V2 * r_C_Pminus_5 - f_V3 * r_C_Pminus_5 );
        *hamiltonian += H;
        *KS_V += H;

        if (compute_hamiltonian_only == true)
        {
            return;
        }

        double C_a_C = -C_V / mu_C;
        
        double grad_e_A_vec_H[3],grad_j_A_vec_H[3];
        double grad_e_B_vec_H[3],grad_j_B_vec_H[3];
        
        for (i=0; i<3; i++)
        {
            grad_e_A_vec_H[i] = C_V * ( e_B_vec_dot_r_C_vec * (-2.0 * r_C_P2 * e_A_vec[i] + 10.0 * e_A_vec_dot_r_C_vec * r_C_vec[i] ) * r_C_Pminus_7 - 2.0 * ( -2.0 * e_B_vec_dot_r_C_vec * e_A_vec[i] + 5.0 * e_A_vec_dot_r_C_vec * e_B_vec[i] + 5.0 * e_A_vec_dot_e_B_vec * r_C_vec[i] ) * r_C_Pminus_5 - 6.0 * e_B_vec_dot_r_C_vec * r_C_Pminus_5 * e_A_vec[i] );
            grad_j_A_vec_H[i] = C_V * ( -2.0 * e_B_vec_dot_r_C_vec * j_A_vec_dot_r_C_vec * r_C_Pminus_7 * r_C_vec[i] + 2.0 * ( j_A_vec_dot_r_C_vec * e_B_vec[i] + j_A_vec_dot_e_B_vec * r_C_vec[i] ) * r_C_Pminus_5 );
            
            grad_e_B_vec_H[i] = C_V * ( f_V11 * r_C_Pminus_7 * r_C_vec[i] - 2.0 * ( (1.0 - e_A_P2) * r_C_vec[i] + 5.0 * e_A_vec_dot_r_C_vec * e_A_vec[i] - j_A_vec_dot_r_C_vec * j_A_vec[i] ) * r_C_Pminus_5 - f_V31 * r_C_Pminus_5 * r_C_vec[i] );
            grad_j_B_vec_H[i] = 0.0;
            
            binary_C->a_vec[i] += C_a_C * ( ( f_V11 * e_B_vec[i] + e_B_vec_dot_r_C_vec * ( 2.0 * (1.0 - e_A_P2) * r_C_vec[i] + 10.0 * e_A_vec_dot_r_C_vec * e_A_vec[i] - 2.0 * j_A_vec_dot_r_C_vec * j_A_vec[i] ) ) * r_C_Pminus_7 - 7.0 * f_V1 * r_C_Pminus_9 * r_C_vec[i] - 2.0 * ( (1.0 - e_A_P2) * e_B_vec[i] + 5.0 * e_A_vec_dot_e_B_vec * e_A_vec[i] - j_A_vec_dot_e_B_vec * j_A_vec[i] ) * r_C_Pminus_5  + 10.0 * f_V2 * r_C_Pminus_7 * r_C_vec[i] - f_V31 * r_C_Pminus_5 * e_B_vec[i] + 5.0 * f_V3 * r_C_Pminus_7 * r_C_vec[i] );
        }

        double j_A_vec_cross_grad_j_A_vec_H[3],j_A_vec_cross_grad_e_A_vec_H[3];
        double e_A_vec_cross_grad_e_A_vec_H[3],e_A_vec_cross_grad_j_A_vec_H[3];

        double j_B_vec_cross_grad_j_B_vec_H[3],j_B_vec_cross_grad_e_B_vec_H[3];
        double e_B_vec_cross_grad_e_B_vec_H[3],e_B_vec_cross_grad_j_B_vec_H[3];
        
        cross3(j_A_vec,             grad_j_A_vec_H,               j_A_vec_cross_grad_j_A_vec_H);
        cross3(j_A_vec,             grad_e_A_vec_H,               j_A_vec_cross_grad_e_A_vec_H);
        cross3(e_A_vec,             grad_e_A_vec_H,               e_A_vec_cross_grad_e_A_vec_H);
        cross3(e_A_vec,             grad_j_A_vec_H,               e_A_vec_cross_grad_j_A_vec_H);    

        cross3(j_B_vec,             grad_j_B_vec_H,               j_B_vec_cross_grad_j_B_vec_H);
        cross3(j_B_vec,             grad_e_B_vec_H,               j_B_vec_cross_grad_e_B_vec_H);
        cross3(e_B_vec,             grad_e_B_vec_H,               e_B_vec_cross_grad_e_B_vec_H);
        cross3(e_B_vec,             grad_j_B_vec_H,               e_B_vec_cross_grad_j_B_vec_H);    

        for (i=0; i<3; i++)
        {
            binary_A->de_vec_dt[i] += (-1.0/(Lambda_A))*( e_A_vec_cross_grad_j_A_vec_H[i] \
                + j_A_vec_cross_grad_e_A_vec_H[i] );
            binary_A->dh_vec_dt[i] += -1.0*( j_A_vec_cross_grad_j_A_vec_H[i] \
                + e_A_vec_cross_grad_e_A_vec_H[i] );

            binary_B->de_vec_dt[i] += (-1.0/(Lambda_B))*( e_B_vec_cross_grad_j_B_vec_H[i] \
                + j_B_vec_cross_grad_e_B_vec_H[i] );
            binary_B->dh_vec_dt[i] += -1.0*( j_B_vec_cross_grad_j_B_vec_H[i] \
                + e_B_vec_cross_grad_e_B_vec_H[i] );
        }
    }
    else if (binary_A->integration_method==0 && binary_B->integration_method>0 && binary_C->integration_method>0)
    {
        /* A averaged; B & C direct */
        #ifdef VERBOSE
        if (verbose_flag > 2)
        {
            printf("ODE_newtonian.cpp -- compute_EOM_binary_triplets -- A %d averaged, B %d C %d direct\n",binary_A->index,binary_B->index,binary_C->index);
        }
        #endif
        
        double *r_B_vec = binary_B->r_vec;
        double *r_C_vec = binary_C->r_vec;
        
        double r_C_P1 = norm3(r_C_vec);
        double r_C_P2 = r_C_P1 * r_C_P1;
        double r_C_Pminus_1 = 1.0/r_C_P1;
        double r_C_Pminus_2 = r_C_Pminus_1 * r_C_Pminus_1;
        double r_C_Pminus_4 = r_C_Pminus_2 * r_C_Pminus_2;
        double r_C_Pminus_5 = r_C_Pminus_1 * r_C_Pminus_4;
        double r_C_Pminus_7 = r_C_Pminus_2 * r_C_Pminus_5;
        double r_C_Pminus_9 = r_C_Pminus_2 * r_C_Pminus_7;
        
        double e_A_vec_dot_r_B_vec = dot3(e_A_vec,r_B_vec);
        double e_A_vec_dot_r_C_vec = dot3(e_A_vec,r_C_vec);
        double e_A_vec_dot_r_C_vec_P2 = e_A_vec_dot_r_C_vec * e_A_vec_dot_r_C_vec;        
        double j_A_vec_dot_r_B_vec = dot3(j_A_vec,r_B_vec);
        double j_A_vec_dot_r_C_vec = dot3(j_A_vec,r_C_vec);
        double j_A_vec_dot_r_C_vec_P2 = j_A_vec_dot_r_C_vec * j_A_vec_dot_r_C_vec;
        double r_B_vec_dot_r_C_vec = dot3(r_B_vec,r_C_vec);
        
        double C_V = c_3div4 * C_tr * a_A * a_A;
        double f_V11 = ( (1.0 - e_A_P2) * r_C_P2 + 5.0 * e_A_vec_dot_r_C_vec_P2 - j_A_vec_dot_r_C_vec_P2 );
        double f_V1 = f_V11 * r_B_vec_dot_r_C_vec;
        double f_V2 = (1.0 - e_A_P2) * r_B_vec_dot_r_C_vec + 5.0 * e_A_vec_dot_r_B_vec * e_A_vec_dot_r_C_vec - j_A_vec_dot_r_B_vec * j_A_vec_dot_r_C_vec;
        double f_V31 = (2.0 + 3.0 * e_A_P2);
        double f_V3 = f_V31 * r_B_vec_dot_r_C_vec;
        
        H = C_V * ( f_V1 * r_C_Pminus_7 - 2.0 * f_V2 * r_C_Pminus_5 - f_V3 * r_C_Pminus_5 );
        *hamiltonian += H;
        *KS_V += H;
        
        if (compute_hamiltonian_only == true)
        {
            return;
        }

        double C_a_B = -C_V / mu_B;
        double C_a_C = -C_V / mu_C;
        
        double grad_e_A_vec_H[3],     grad_j_A_vec_H[3];
        
        for (i=0; i<3; i++)
        {
            grad_e_A_vec_H[i] = C_V * ( r_B_vec_dot_r_C_vec * (-2.0 * r_C_P2 * e_A_vec[i] + 10.0 * e_A_vec_dot_r_C_vec * r_C_vec[i] ) * r_C_Pminus_7 - 2.0 * ( -2.0 * r_B_vec_dot_r_C_vec * e_A_vec[i] + 5.0 * e_A_vec_dot_r_C_vec * r_B_vec[i] + 5.0 * e_A_vec_dot_r_B_vec * r_C_vec[i] ) * r_C_Pminus_5 - 6.0 * r_B_vec_dot_r_C_vec * r_C_Pminus_5 * e_A_vec[i] );
            grad_j_A_vec_H[i] = C_V * ( -2.0 * r_B_vec_dot_r_C_vec * j_A_vec_dot_r_C_vec * r_C_Pminus_7 * r_C_vec[i] + 2.0 * ( j_A_vec_dot_r_C_vec * r_B_vec[i] + j_A_vec_dot_r_B_vec * r_C_vec[i] ) * r_C_Pminus_5 );
            
            binary_B->a_vec[i] += C_a_B * ( f_V11 * r_C_Pminus_7 * r_C_vec[i] - 2.0 * ( (1.0 - e_A_P2) * r_C_vec[i] + 5.0 * e_A_vec_dot_r_C_vec * e_A_vec[i] - j_A_vec_dot_r_C_vec * j_A_vec[i] ) * r_C_Pminus_5 - f_V31 * r_C_Pminus_5 * r_C_vec[i] );
            binary_C->a_vec[i] += C_a_C * ( ( f_V11 * r_B_vec[i] + r_B_vec_dot_r_C_vec * ( 2.0 * (1.0 - e_A_P2) * r_C_vec[i] + 10.0 * e_A_vec_dot_r_C_vec * e_A_vec[i] - 2.0 * j_A_vec_dot_r_C_vec * j_A_vec[i] ) ) * r_C_Pminus_7 - 7.0 * f_V1 * r_C_Pminus_9 * r_C_vec[i] - 2.0 * ( (1.0 - e_A_P2) * r_B_vec[i] + 5.0 * e_A_vec_dot_r_B_vec * e_A_vec[i] - j_A_vec_dot_r_B_vec * j_A_vec[i] ) * r_C_Pminus_5  + 10.0 * f_V2 * r_C_Pminus_7 * r_C_vec[i] - f_V31 * r_C_Pminus_5 * r_B_vec[i] + 5.0 * f_V3 * r_C_Pminus_7 * r_C_vec[i] );
        }

        double j_A_vec_cross_grad_j_A_vec_H[3],j_A_vec_cross_grad_e_A_vec_H[3];
        double e_A_vec_cross_grad_e_A_vec_H[3],e_A_vec_cross_grad_j_A_vec_H[3];
        
        cross3(j_A_vec,             grad_j_A_vec_H,               j_A_vec_cross_grad_j_A_vec_H);
        cross3(j_A_vec,             grad_e_A_vec_H,               j_A_vec_cross_grad_e_A_vec_H);
        cross3(e_A_vec,             grad_e_A_vec_H,               e_A_vec_cross_grad_e_A_vec_H);
        cross3(e_A_vec,             grad_j_A_vec_H,               e_A_vec_cross_grad_j_A_vec_H);    

        for (i=0; i<3; i++)
        {
            binary_A->de_vec_dt[i] += (-1.0/(Lambda_A))*( e_A_vec_cross_grad_j_A_vec_H[i] \
                + j_A_vec_cross_grad_e_A_vec_H[i] );
            binary_A->dh_vec_dt[i] += -1.0*( j_A_vec_cross_grad_j_A_vec_H[i] \
                + e_A_vec_cross_grad_e_A_vec_H[i] );
        }
        
    }
    else if (binary_A->integration_method>0 && binary_B->integration_method>0 && binary_C->integration_method>0)
    {
        /* A, B & C direct */
        #ifdef VERBOSE
        if (verbose_flag > 2)
        {
            printf("ODE_newtonian.cpp -- compute_EOM_binary_triplets -- A %d, B %d, C %d direct\n",binary_A->index,binary_B->index,binary_C->index);
        }
        #endif

        double *r_A_vec = binary_A->r_vec;
        double *r_B_vec = binary_B->r_vec;
        double *r_C_vec = binary_C->r_vec;
        
        double r_A_P1 = norm3(r_A_vec);
        double r_A_P2 = r_A_P1 * r_A_P1;
        double r_C_P1 = norm3(r_C_vec);
        double r_C_Pminus_1 = 1.0/r_C_P1;
        double r_C_Pminus_2 = r_C_Pminus_1 * r_C_Pminus_1;
        double r_C_Pminus_4 = r_C_Pminus_2 * r_C_Pminus_2;
        double r_C_Pminus_5 = r_C_Pminus_1 * r_C_Pminus_4;
        double r_C_Pminus_7 = r_C_Pminus_2 * r_C_Pminus_5;
        double r_C_Pminus_9 = r_C_Pminus_2 * r_C_Pminus_7;
        
        double r_A_vec_dot_r_B_vec = dot3(r_A_vec,r_B_vec);
        double r_A_vec_dot_r_C_vec = dot3(r_A_vec,r_C_vec);
        double r_A_vec_dot_r_C_vec_P2 = r_A_vec_dot_r_C_vec * r_A_vec_dot_r_C_vec;
        double r_B_vec_dot_r_C_vec = dot3(r_B_vec,r_C_vec);
        
        H = c_3div2 * C_tr * ( r_A_vec_dot_r_C_vec_P2 * r_B_vec_dot_r_C_vec * r_C_Pminus_7 - 2.0 * r_A_vec_dot_r_C_vec * r_A_vec_dot_r_B_vec * r_C_Pminus_5 - r_B_vec_dot_r_C_vec * r_A_P2 * r_C_Pminus_5 );
        *hamiltonian += H;
        *KS_V += H;
        
        if (compute_hamiltonian_only == true)
        {
            return;
        }
        
        double C_a_A = -c_3div2 * C_tr / mu_A;
        double C_a_B = -c_3div2 * C_tr / mu_B;
        double C_a_C = -c_3div2 * C_tr / mu_C;
        
        for (i=0; i<3; i++)
        {
            binary_A->a_vec[i] += C_a_A * ( 2.0 * r_A_vec_dot_r_C_vec * r_B_vec_dot_r_C_vec * r_C_Pminus_7 * r_C_vec[i] - 2.0 * ( r_A_vec_dot_r_B_vec * r_C_vec[i] + r_A_vec_dot_r_C_vec * r_B_vec[i] ) * r_C_Pminus_5 - 2.0 * r_B_vec_dot_r_C_vec * r_C_Pminus_5 * r_A_vec[i] );
            binary_B->a_vec[i] += C_a_B * ( r_A_vec_dot_r_C_vec_P2 * r_C_Pminus_7 * r_C_vec[i] - 2.0 * r_A_vec_dot_r_C_vec * r_C_Pminus_5 * r_A_vec[i] - r_A_P2 * r_C_Pminus_5 * r_C_vec[i] );
            binary_C->a_vec[i] += C_a_C * ( ( 2.0 * r_A_vec_dot_r_C_vec * r_B_vec_dot_r_C_vec * r_A_vec[i] + r_A_vec_dot_r_C_vec_P2 * r_B_vec[i] ) * r_C_Pminus_7 - 7.0 * r_A_vec_dot_r_C_vec_P2 * r_B_vec_dot_r_C_vec * r_C_Pminus_9 * r_C_vec[i] - 2.0 * r_A_vec_dot_r_B_vec * r_C_Pminus_5 * r_A_vec[i] + 10.0 * r_A_vec_dot_r_C_vec * r_A_vec_dot_r_B_vec * r_C_Pminus_7 * r_C_vec[i] - r_A_P2 * r_C_Pminus_5 * r_B_vec[i] + 5.0 * r_B_vec_dot_r_C_vec * r_A_P2 * r_C_Pminus_7 * r_C_vec[i] );
        }
    }
    else
    {
        
        printf("FATAL ERROR in newtonian.cpp: invalid combination of integration methods %d %d %d for binary triplet combination with indices %d %d %d (note: an orbit that is averaged cannot have children orbits that are directly integrated. \n",binary_A->integration_method,binary_B->integration_method,binary_C->integration_method,binary_A->index,binary_B->index,binary_C->index);
        exit(-1);
    }

}

double compute_order_of_magnitude_secular_timescale_for_pair(ParticlesMap *particlesMap, int inner_binary_index, int outer_binary_index, int connecting_child_in_outer_binary)
{

    Particle *inner_binary = (*particlesMap)[inner_binary_index];
    Particle *outer_binary = (*particlesMap)[outer_binary_index];
    
    Particle *P_sibling;
    if (connecting_child_in_outer_binary==1)
    {
        P_sibling = (*particlesMap)[outer_binary->child2];
    }
    else if (connecting_child_in_outer_binary==2)
    {
        P_sibling = (*particlesMap)[outer_binary->child1];
    }

    double M_in = inner_binary->mass;
    double m3 = P_sibling->mass;
    double M_out = M_in + m3;

    double a_in = inner_binary->a;
    double a_out = outer_binary->a;
    
    double e_out = outer_binary->e;

    double P_in = compute_orbital_period_from_semimajor_axis(M_in,a_in);
    double P_out = compute_orbital_period_from_semimajor_axis(M_out,a_out);
    double t_sec = (P_out*P_out/P_in)*(M_out/m3)*pow(1.0 - e_out*e_out,3.0/2.0);
        
    return t_sec;
}

}
