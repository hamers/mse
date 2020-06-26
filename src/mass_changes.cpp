/* MSE */

#include "evolve.h"
#include "mass_changes.h"

extern "C"
{

int ODE_handle_stellar_winds(Particle *p)
{
    /* p is assumed to be a binary */
    
//    ParticlesMapIterator it_p;
    //std::vector<int>::iterator it_parent_p,it_parent_q;

//    int seed = orbital_phases_random_seed;
//    int index=0;
    
    //int is_binary;
    double factor_h_vec,factor_spin_vec;


    //for (it_p = particlesMap->begin(); it_p != particlesMap->end(); it_p++)
    //{
        //Particle *p = (*it_p).second;
        
        /* assuming constant semimajor axis with mass loss */
        //double factor_h_vec = (child1_mass_dot/(2.0*child1_mass))*(1.0 + child2_mass/child1_mass_plus_child2_mass) + (child2_mass_dot/(2.0*child2_mass))*(1.0 + child1_mass/child1_mass_plus_child2_mass); 
        /* assuming constant SPECIFIC orbital angular momentum with mass loss, i.e. a(m1+m2)=const. */
        //is_binary = p->is_binary;
        //printf("??? %d %d\n",is_binary,p->is_binary);
        
        //if (p->is_binary == 1)
        //{
    double m1 = p->child1_mass;
    double m2 = p->child2_mass;

    double m1dot = p->child1_mass_dot_wind + p->child1_mass_dot_wind_accretion + p->child1_mass_dot_adiabatic_ejection;
    double m2dot = p->child2_mass_dot_wind + p->child2_mass_dot_wind_accretion + p->child2_mass_dot_adiabatic_ejection;;

    factor_h_vec = m1dot/m1 + m2dot/m2 - (m1dot + m2dot)/(m1+m2);
    //factor_h_vec = p->child1_mass_dot_wind/p->child1_mass + p->child2_mass_dot_wind/p->child2_mass - (p->child1_mass_dot_wind + p->child2_mass_dot_wind)/p->child1_mass_plus_child2_mass; /* assuming constant SPECIFIC orbital angular momentum with mass loss, i.e. a(m1+m2)=const. */
            //printf("factor_h_vec %g\n",factor_h_vec);
//        }
        //else
        //{
        
            /* assuming constant spin angular momentum of the body */
            //factor_spin_vec = - (p->mass_dot/p->mass + 2.0*p->radius_dot/p->radius);
//        }
//
        /* set external time derivatives if appropriate */
        //double h_vec_dot[3] = {h_vec_x_dot,h_vec_y_dot,h_vec_z_dot};
        //double e_vec_dot[3] = {e_vec_x_dot,e_vec_y_dot,e_vec_z_dot};
        //double spin_vec_dot[3] = {spin_vec_x_dot,spin_vec_y_dot,spin_vec_z_dot};

    for (int i=0; i<3; i++)
    {
    //        if (p->is_binary == 1)
      //      {
        p->dh_vec_dt[i] += p->h_vec[i]*factor_h_vec; //+ h_vec_dot[i];
    }
            //else
          //  {
                //dspin_vec_dt[i] += spin_vec[i]*factor_spin_vec; // + spin_vec_dot[i];
        //    }
      //  }
    //}
    
    return 0;

}

int ODE_handle_RLOF(ParticlesMap *particlesMap, Particle *p)
{
    Particle *child1 = (*particlesMap)[p->child1];
    Particle *child2 = (*particlesMap)[p->child2];
    
    if (child1->is_binary == false and child2->is_binary == false)
    /* For now, only allow mass transfer between two single stars in a binary (e.g., no transfer from tertiary to inner binary */
    {
     //   printf("ok %d %d %d\n",p->index,child1->index,child2->index);
        ODE_handle_RLOF_emt(p, child1, child2);
    }
    else if (child1->is_binary == false and child2->is_binary == true)
    {
        ODE_handle_RLOF_triple_mass_transfer(particlesMap, p, child1, child2);
    }
    else if (child1->is_binary == true and child2->is_binary == false)
    {
        ODE_handle_RLOF_triple_mass_transfer(particlesMap, p, child2, child1);
    }
    
    return 0;
}

int ODE_handle_RLOF_triple_mass_transfer(ParticlesMap *particlesMap, Particle *outer_binary, Particle *donor, Particle *inner_binary)
{
    /* Mass transfer from tertiary star to inner binary. */
    
    if (donor->is_binary == true or inner_binary->is_binary == false)
    {
        printf("mass_changes.cpp -- ODE_handle_RLOF_triple_mass_transfer -- ERROR: donor %d should be star; inner binary %d should be a binary!\n",donor->index,inner_binary->index);
        exit(-1);
    }

    if (donor->include_mass_transfer_terms==false)
    {
        return 0;
    }
    
    Particle *child1 = (*particlesMap)[inner_binary->child1];
    Particle *child2 = (*particlesMap)[inner_binary->child2];
    
    double m_donor = donor->mass;
    double m_inner_binary = inner_binary->mass;

    double R_donor = donor->radius;
       
    double a_in = inner_binary->a;
    double a_out = outer_binary->a;

    double q = m_donor/m_inner_binary;
    double R_Lc = roche_radius_pericenter_eggleton(a_out,q); /* with argument "a", actually computes circular Roche lobe radius */
    if (R_donor < R_Lc)
    {
        /* There is no RLOF of tertiary star; do nothing. */
        return 0;
    }

    printf("mass_changes.cpp -- ODE_handle_RLOF_triple_mass_transfer -- outer_binary %d donor %d inner_binary %d\n",outer_binary->index,donor->index,inner_binary->index);
    
    double m_C1 = child1->mass;
    double m_C2 = child2->mass;

    double e_in = inner_binary->e;
    double e_out = outer_binary->e;

    double m_donor_dot = donor->mass_dot_RLOF_triple;
    double m_C1_dot = child1->mass_dot_RLOF_triple;
    double m_C2_dot = child2->mass_dot_RLOF_triple;

        
    /* Effect on inner orbit -- CE-like behavior. */
    double a_in_dot = inner_binary->triple_mass_transfer_a_in_dot; /* This was computed beforehand in binary_evolution.cpp -- triple_stable_mass_transfer_evolution() */
    double e_in_dot = 0.0;
    
    /* Effect on outer orbit -- non-conservative mass transfer. */
    double beta = -(m_C1_dot+m_C2_dot)/m_donor_dot;
    double gamma = m_donor/m_inner_binary; /* Isotropic re-emission (Pols, lecture notes on binary evolution, Chapter 7) */
    double a_out_dot = compute_a_dot_circular_non_conservative_mass_transfer(a_out,m_donor,m_donor_dot,m_inner_binary,beta,gamma);
    double e_out_dot = 0.0;

    /* Effect on spins. */
    double Omega_donor = donor->spin_vec_norm;
    double J_spin_donor_dot = m_donor*R_donor*R_donor*Omega_donor;
    double Omega_donor_dot = compute_Omega_dot_from_J_dot_mass_transfer(J_spin_donor_dot, Omega_donor, m_donor, m_donor_dot, donor->core_mass, R_donor, donor->radius_dot, donor->core_radius, donor->sse_k2, donor->sse_k3);
    
    /* Assume the spins in the inner binary are unaffected. */

    if (Omega_donor_dot != Omega_donor_dot)
    {
        printf("mass_changes.cpp -- ODE_handle_RLOF_triple_mass_transfer -- Omega_donor_dot %g\n",Omega_donor_dot);
        exit(-1);
    }


    /* Update dmass_dt. */
    donor->dmass_dt += m_donor_dot;
    child1->dmass_dt += m_C1_dot;
    child2->dmass_dt += m_C2_dot;
    
    /* Update h,e, and spin vectors. */
    double h_in_dot_div_h_in = compute_h_dot_div_h(m_C1, m_C1_dot, m_C2, m_C2_dot, a_in, a_in_dot, e_in, e_in_dot);
    double h_out_dot_div_h_out = compute_h_dot_div_h(m_donor, m_donor_dot, m_inner_binary, m_C1_dot+m_C2_dot, a_out, a_out_dot, e_out, e_out_dot);
 
    for (int i=0; i<3; i++)
    {
        inner_binary->dh_vec_dt[i] += inner_binary->h_vec[i] * h_in_dot_div_h_in;
        outer_binary->dh_vec_dt[i] += outer_binary->h_vec[i] * h_out_dot_div_h_out;

        donor->dspin_vec_dt[i] += donor->spin_vec_unit[i] * Omega_donor_dot;
        
        //p->de_vec_dt[i] += p->e_vec_unit[i] * de_dt;
    }
    
    
    
    return 0;
}

double compute_a_dot_circular_non_conservative_mass_transfer(double a, double m_donor, double m_donor_dot, double m_accretor, double beta, double gamma)
{
    /* Pols, lecture notes on binary evolution, Chapter 7, eq. (7.14). */
    return -2.0*a*(m_donor_dot/m_donor) * ( 1.0 - beta*(m_donor/m_accretor) - (1.0 - beta)*(gamma + c_1div2)*m_donor/(m_donor+m_accretor) );
}

double compute_Omega_dot_from_J_dot_mass_transfer(double J_spin_dot, double Omega, double mass, double mass_dot, double core_mass, double radius, double radius_dot, double core_radius, double sse_k2, double sse_k3)
{
    double I = compute_moment_of_inertia(mass, core_mass, radius, core_radius, sse_k2, sse_k3);
    double I_dot = sse_k2*mass_dot*radius*radius + 2.0*sse_k2*(mass - core_mass)*radius*radius_dot; /* Approximate; neglects changes in core masses and radii and k2 and k3  */
    double Omega_dot = (J_spin_dot - I_dot*Omega)/I;
    return Omega_dot;
}

int ODE_handle_RLOF_emt(Particle *p, Particle *child1, Particle *child2)
{
    double a = p->a;
    double e = p->e;
    double q,R_Lc,x,E_0;
    bool in_RLOF;
    int flag;
    
    /* child1 -> child2 */
    if (child1->include_mass_transfer_terms==true)
    {
        q = child1->mass/child2->mass;
        R_Lc = roche_radius_pericenter_eggleton(a,q); /* with argument "a", actually computes circular Roche lobe radius */
        //printf("R_lc %g a %g q %g M %g R %g\n",R_Lc,a,q,child1->mass,child1->radius);
        x = R_Lc/child1->radius;
        flag = determine_E_0(e, x, &E_0, &in_RLOF);
        //printf("1 x %g q %g E_0 %g\n",x,q,E_0);
        if (in_RLOF == true and flag == 0)
        {
            compute_RLOF_emt_model(p,child1,child2,x,E_0);
            //printf("Delta R/R %g\n",(child1->radius-R_Lc)/child1->radius);
            #ifdef DEBUG
            printf("mass_changes.cpp -- IN RLOF *1 index %d x %g E_0 %g\n",child1->index,x,E_0);
            #endif
        }
    }

    /* child2 -> child1 */
    if (child2->include_mass_transfer_terms==true)
    {
        q = child2->mass/child1->mass;
        R_Lc = roche_radius_pericenter_eggleton(a,q);
        x = R_Lc/child2->radius;
        flag = determine_E_0(e, x, &E_0, &in_RLOF);
        //printf("2 x %g q %g E_0 %g\n",x,q,E_0);
        if (in_RLOF == true and flag == 0)
        {
            compute_RLOF_emt_model(p,child2,child1,x,E_0);
            #ifdef DEBUG
            printf("mass_changes.cpp -- IN RLOF *2 index %d x %g E_0 %g\n",child2->index,x,E_0);
            #endif
        }
    }
    
    return 0;
}

int compute_RLOF_emt_model(Particle *p, Particle *donor, Particle *accretor, double x, double E_0)
{
    /* Implementation of the analytic model of https://ui.adsabs.harvard.edu/abs/2019ApJ...872..119H/abstract,
     * with one important modification: instead of conservative transfer, adapt it to the non-conservative case. 
     * This is achieved by setting M_a_dot = - beta * M_d_dot, where beta=1 in the conservative case.
     * This modifies all the equations in a very simple way: make the replacement q -> q * beta. 
     * For reference: the "Sepinsky" model (https://ui.adsabs.harvard.edu/abs/2007ApJ...667.1170S/abstract) has
     * fm = 1.0, fa = \sqrt(1-e^2), fe = fa*(1-e), and fomega = 0.0. */
    int i;
    
    double M_d = donor->mass;
    double M_a = accretor->mass;
    double M = M_d + M_a;
    double q = M_d/M_a;
    double tau = donor->emt_tau;
    double R_d = donor->radius;
    double R_a = accretor->radius;
    double R_d_p2 = R_d*R_d;
    double R_a_p2 = R_a*R_a;
    
    double Omega_d = donor->spin_vec_norm;
    double Omega_a = accretor->spin_vec_norm;
    
    double a = p->a;
    double e = p->e;
    
    if (e < epsilon)
    {
        e = epsilon;
    }

    if (M_d <= epsilon or M_a <= epsilon) /* No more mass left; skip. */
    {
        return 0;
    }
    
    double fm,fa,fe,fomega,ga,ge,ha,he;
    fm=fa=fe=fomega=ga=ge=ha=he=0.0;

    double finite_size_term_a = 0.0;
    double finite_size_term_e = 0.0;

    double M_d_dot_av = donor->mass_dot_RLOF;
    double M_a_dot_av = accretor->mass_dot_RLOF;

    if (fabs(M_d_dot_av) <= epsilon or fabs(M_a_dot_av) <= epsilon) /* Effectively zero mass transfer rate. */
    {
        return 0;
    }

    double beta = -M_a_dot_av/M_d_dot_av; /* Quantifies non-conservativeness. */

    //printf("B %g\n",beta);
    
    /* TO DO: allow for user-specified MA_tau */
    double n = sqrt(CONST_G*M/(a*a*a));
    double MA_tau = n*tau; 

    double cos_eccentric_anomaly,sin_eccentric_anomaly;
    compute_eccentric_anomaly_from_mean_anomaly(MA_tau, e, &cos_eccentric_anomaly, &sin_eccentric_anomaly);

    double E_tau = atan2(sin_eccentric_anomaly,cos_eccentric_anomaly);
    if (E_tau > M_PI) /* E_tau cannot not be larger than Pi */
    {
        E_tau = M_PI; 
    }
    
    fm = fm_function(e,x,E_0,E_tau);
    fa = fa_function(e,x,E_0,E_tau);
    fe = fe_function(e,x,E_0,E_tau);
    fomega = fomega_function(e,x,E_0,E_tau);

    donor->emt_fm = fm;

    if (fabs(fm) <= epsilon)
    {
        fm = epsilon;
    }

    if (donor->emt_ejection_radius_mode > 0 or accretor->emt_accretion_radius > 0.0)
    {
       
        ha = ha_function(e,x,E_0);
        he = he_function(e,x,E_0);

        finite_size_term_a = - q * beta * (accretor->emt_accretion_radius/a) * ha;
        finite_size_term_e = - q * beta * (accretor->emt_accretion_radius/a) * he;

        if (donor->emt_ejection_radius_mode == 1) /*  limit of small donor spin */
        {
            ga = ga_function(e,x,E_0);
            ge = ge_function(e,x,E_0);
            double XL0 = XL0_q_function(q);
            
            finite_size_term_a += XL0*ga;
            finite_size_term_e += XL0*ge;
        }
        else if (donor->emt_ejection_radius_mode == 2) /* limit of large mass ratio q */
        {
            double temp = 1.0 - e;
            double n_peri = n*sqrt( (1.0 + e)/(temp*temp*temp) );
            double ospin_hat = norm3(donor->spin_vec) / n_peri;
            double XL0 = pow(ospin_hat,-2.0/3.0) * temp / pow( 1.0+e, 1.0/3.0);

            finite_size_term_a += XL0*ha;
            finite_size_term_e += XL0*he;
                
        }
    }

    //double P_orb = compute_orbital_period(p);
    //double M_d_dot_av = -(M_d/P_orb)*fm;
    //double M_a_dot_av = -M_d_dot_av;
    
    if (fabs(M_d_dot_av) > 1.0)
    {
        printf("mass_changes.cpp -- changing M_d_MT %g to -1.0\n",M_d_dot_av);
        M_d_dot_av = -1.0;
        M_a_dot_av = -M_d_dot_av;
    }
    
    double common_factor = -2.0*(M_d_dot_av/M_d)*(1.0/fm);
    //printf("CF %g\n",common_factor);
    if (x<=0.0)
    {
        printf("mass_changes.cpp -- ERROR: x<=0 (x=%g)\n",x);
        M_d_dot_av=0.0;
        M_a_dot_av=0.0;
    }

    //printf("mass_changes.cpp -- e %g x %g E_0 %g m_d %g m_a %g fm %g md %g\n",e,x,E_0,M_d,M_a,fm,M_d_dot_av);
    

    
    double da_dt = common_factor*a*( fa*(1.0 - q*beta) + finite_size_term_a );
    double de_dt = common_factor*( fe*(1.0 - q*beta) + finite_size_term_e );
    double domega_dt = common_factor*fomega*( 1.0 - q*beta );

    /* Update dmass_dt. */
    donor->dmass_dt += M_d_dot_av;
    accretor->dmass_dt += M_a_dot_av;

    //double factor_h_vec = M_d_dot_av/M_d + M_a_dot_av/M_a - c_1div2*(M_d_dot_av + M_a_dot_av)/M + c_1div2*(da_dt/a) - e*de_dt/(1.0 - e*e);
    double factor_h_vec = compute_h_dot_div_h(M_d, M_d_dot_av, M_a, M_a_dot_av, a, da_dt, e, de_dt);
    
    if (da_dt != da_dt || de_dt != de_dt || domega_dt != domega_dt)
    {
        printf("mass_changes.cpp -- ERROR: nans in dots of a/e/omega, %g %g %g beta %g M_a_dot_av %g M_a_dot_av %g\n",da_dt,de_dt,domega_dt,beta,M_a_dot_av,M_a_dot_av);
    }

    /* Compute spin changes due to RLOF */
    double J_spin_donor_dot = M_d_dot_av*R_d_p2*Omega_d;
    //double I_donor = compute_moment_of_inertia(M_d, donor->core_mass, R_d, donor->core_radius, donor->sse_k2, donor->sse_k3);
    //double I_donor_dot = donor->sse_k2*donor->dmass_dt*R_d_p2 + 2.0*donor->sse_k2*(M_d - donor->core_mass)*R_d*donor->radius_dot; /* Approximate; neglects changes in core masses and radii and k2 and k3  */
    //double Omega_d_dot = (J_spin_donor_dot - I_donor_dot*Omega_d)/I_donor;
    double Omega_d_dot = compute_Omega_dot_from_J_dot_mass_transfer(J_spin_donor_dot, Omega_d, M_d, M_d_dot_av, donor->core_mass, R_d, donor->radius_dot, donor->core_radius, donor->sse_k2, donor->sse_k3);
    
    double J_spin_accretor_dot;
    if (accretor->accretion_disk_is_present == true)
    {
        J_spin_accretor_dot = M_a_dot_av*sqrt(CONST_G * M_a * R_a);
    }
    else
    {
        double r_disk = 1.7*accretor->accretion_disk_r_min;
        J_spin_accretor_dot = M_a_dot_av*sqrt(CONST_G * M_a * r_disk);
    }
    //double I_accretor = compute_moment_of_inertia(M_a, accretor->core_mass, R_a, accretor->core_radius, accretor->sse_k2, accretor->sse_k3);
    //double I_accretor_dot = accretor->sse_k2*accretor->dmass_dt*R_a_p2 + 2.0*accretor->sse_k2*(M_a - accretor->core_mass)*R_a*accretor->radius_dot; /* Approximate; neglects changes in core masses and radii and k2 and k3  */
    //double Omega_a_dot = (J_spin_accretor_dot - I_accretor_dot*Omega_a)/I_accretor;
    double Omega_a_dot = compute_Omega_dot_from_J_dot_mass_transfer(J_spin_accretor_dot, Omega_a, M_a, M_a_dot_av, accretor->core_mass, R_a, accretor->radius_dot, accretor->core_radius, accretor->sse_k2, accretor->sse_k3);
    
    if (Omega_d_dot != Omega_d_dot or Omega_a_dot != Omega_a_dot or factor_h_vec!=factor_h_vec or de_dt!=de_dt or domega_dt!=domega_dt)
    {
        printf("mass_changes.cpp -- compute_RLOF_emt_model -- Omega_d_dot %g Omega_a_dot %g factor_h_vec %g de_dt %g domega_dt %g\n",Omega_d_dot,Omega_a_dot,factor_h_vec,de_dt,domega_dt);
        exit(-1);
    }
    

    for (i=0; i<3; i++)
    {
        p->dh_vec_dt[i] += p->h_vec[i] * factor_h_vec;
        p->de_vec_dt[i] += p->e_vec_unit[i] * de_dt  +  e * p->q_vec_unit[i] * domega_dt;

        /* Assume mass transfer does not affect the directions of the spins. */
        donor->dspin_vec_dt[i] += donor->spin_vec_unit[i] * Omega_d_dot;
        accretor->dspin_vec_dt[i] += accretor->spin_vec_unit[i] * Omega_a_dot;
        
        //#ifdef VERBOSE
        //printf("mass_changes.cpp -- p->a %g p->dh_vec_dt[i] %g p->de_vec_dt[i] %g da_dt %g de_dt %g domega_dt %g \n",p->a,p->dh_vec_dt[i],p->de_vec_dt[i],da_dt,de_dt,domega_dt);
        //#endif
    }


    return 0;
}

int determine_E_0(double e, double x, double *E_0, bool *in_RLOF)
{
    if (x >= 1.0/(1.0-e))
    {
        *E_0 = 0.0;
        *in_RLOF = false;
    }
    else if (x <= 1.0/(1.0+e))
    {
        *E_0 = M_PI;
        *in_RLOF = true;
    }
    else
    {
        *E_0 = acos( (1.0/e)*(1.0 - 1.0/x) );
        *in_RLOF = true;
    }
    if (*E_0 != *E_0)
    {
        printf("mass_changes.cpp -- ERROR in determine_E_0; E_0 = %g x  = %g e = %g \n",*E_0,x,e);
        return -1;
    }
    //printf("x %g e %g E_0... %g\n",x,e,*E_0);
    return 0;
}





double fm_function(double e, double x, double E0, double Etau)
{
    return -(-192*E0 + 576*E0*x - 576*E0*pow(x,2) - 288*pow(e,2)*E0*pow(x,2) + 192*E0*pow(x,3) + 288*pow(e,2)*E0*pow(x,3) + 72*pow(e,2)*E0*x*(4 - 8*x + (4 + pow(e,2))*pow(x,2))*cos(Etau) - 
      96*e*(-1 + x)*(2 - 4*x + (2 + 3*pow(e,2))*pow(x,2))*sin(E0) + 6*pow(e,4)*pow(x,3)*sin(2*E0 - 3*Etau) - 8*pow(e,3)*pow(x,3)*sin(3*E0 - 3*Etau) + 
      3*pow(e,4)*pow(x,3)*sin(4*E0 - 3*Etau) + 72*pow(e,3)*pow(x,2)*sin(E0 - 2*Etau) - 72*pow(e,3)*pow(x,3)*sin(E0 - 2*Etau) - 72*pow(e,2)*pow(x,2)*sin(2*E0 - 2*Etau) + 
      72*pow(e,2)*pow(x,3)*sin(2*E0 - 2*Etau) + 24*pow(e,3)*pow(x,2)*sin(3*E0 - 2*Etau) - 24*pow(e,3)*pow(x,3)*sin(3*E0 - 2*Etau) - 288*e*x*sin(E0 - Etau) + 
      576*e*pow(x,2)*sin(E0 - Etau) - 288*e*pow(x,3)*sin(E0 - Etau) - 72*pow(e,3)*pow(x,3)*sin(E0 - Etau) + 72*pow(e,2)*x*sin(2*E0 - Etau) - 144*pow(e,2)*pow(x,2)*sin(2*E0 - Etau) + 
      72*pow(e,2)*pow(x,3)*sin(2*E0 - Etau) + 18*pow(e,4)*pow(x,3)*sin(2*E0 - Etau) - 288*e*x*sin(E0 + Etau) + 576*e*pow(x,2)*sin(E0 + Etau) - 288*e*pow(x,3)*sin(E0 + Etau) - 
      72*pow(e,3)*pow(x,3)*sin(E0 + Etau) - 72*pow(e,2)*pow(x,2)*sin(2*(E0 + Etau)) + 72*pow(e,2)*pow(x,3)*sin(2*(E0 + Etau)) - 8*pow(e,3)*pow(x,3)*sin(3*(E0 + Etau)) + 
      72*pow(e,2)*x*sin(2*E0 + Etau) - 144*pow(e,2)*pow(x,2)*sin(2*E0 + Etau) + 72*pow(e,2)*pow(x,3)*sin(2*E0 + Etau) + 18*pow(e,4)*pow(x,3)*sin(2*E0 + Etau) + 
      72*pow(e,3)*pow(x,2)*sin(E0 + 2*Etau) - 72*pow(e,3)*pow(x,3)*sin(E0 + 2*Etau) + 24*pow(e,3)*pow(x,2)*sin(3*E0 + 2*Etau) - 24*pow(e,3)*pow(x,3)*sin(3*E0 + 2*Etau) + 
      6*pow(e,4)*pow(x,3)*sin(2*E0 + 3*Etau) + 3*pow(e,4)*pow(x,3)*sin(4*E0 + 3*Etau))/(192.*M_PI);
}
double fa_function(double e, double x, double E0, double Etau)
{
    return (192*E0 - 576*E0*x + 576*E0*pow(x,2) + 288*pow(e,2)*E0*pow(x,2) - 192*E0*pow(x,3) - 288*pow(e,2)*E0*pow(x,3) + 72*pow(e,2)*E0*x*(4 - 8*x + (4 + pow(e,2))*pow(x,2))*cos(Etau) - 
     96*e*(-1 + x)*(2 - 4*x + (2 + 3*pow(e,2))*pow(x,2))*sin(E0) + 6*pow(e,4)*pow(x,3)*sin(2*E0 - 3*Etau) + 8*pow(e,3)*pow(x,3)*sin(3*E0 - 3*Etau) + 
     3*pow(e,4)*pow(x,3)*sin(4*E0 - 3*Etau) + 72*pow(e,3)*pow(x,2)*sin(E0 - 2*Etau) - 72*pow(e,3)*pow(x,3)*sin(E0 - 2*Etau) + 72*pow(e,2)*pow(x,2)*sin(2*E0 - 2*Etau) - 
     72*pow(e,2)*pow(x,3)*sin(2*E0 - 2*Etau) + 24*pow(e,3)*pow(x,2)*sin(3*E0 - 2*Etau) - 24*pow(e,3)*pow(x,3)*sin(3*E0 - 2*Etau) + 288*e*x*sin(E0 - Etau) - 576*e*pow(x,2)*sin(E0 - Etau) + 
     288*e*pow(x,3)*sin(E0 - Etau) + 72*pow(e,3)*pow(x,3)*sin(E0 - Etau) + 72*pow(e,2)*x*sin(2*E0 - Etau) - 144*pow(e,2)*pow(x,2)*sin(2*E0 - Etau) + 
     72*pow(e,2)*pow(x,3)*sin(2*E0 - Etau) + 18*pow(e,4)*pow(x,3)*sin(2*E0 - Etau) + 288*e*x*sin(E0 + Etau) - 576*e*pow(x,2)*sin(E0 + Etau) + 288*e*pow(x,3)*sin(E0 + Etau) + 
     72*pow(e,3)*pow(x,3)*sin(E0 + Etau) + 72*pow(e,2)*pow(x,2)*sin(2*(E0 + Etau)) - 72*pow(e,2)*pow(x,3)*sin(2*(E0 + Etau)) + 8*pow(e,3)*pow(x,3)*sin(3*(E0 + Etau)) + 
     72*pow(e,2)*x*sin(2*E0 + Etau) - 144*pow(e,2)*pow(x,2)*sin(2*E0 + Etau) + 72*pow(e,2)*pow(x,3)*sin(2*E0 + Etau) + 18*pow(e,4)*pow(x,3)*sin(2*E0 + Etau) + 
     72*pow(e,3)*pow(x,2)*sin(E0 + 2*Etau) - 72*pow(e,3)*pow(x,3)*sin(E0 + 2*Etau) + 24*pow(e,3)*pow(x,2)*sin(3*E0 + 2*Etau) - 24*pow(e,3)*pow(x,3)*sin(3*E0 + 2*Etau) + 
     6*pow(e,4)*pow(x,3)*sin(2*E0 + 3*Etau) + 3*pow(e,4)*pow(x,3)*sin(4*E0 + 3*Etau))/(192.*M_PI);
}
double fe_function(double e, double x, double E0, double Etau)
{
    return -((-1 + pow(e,2))*(12*e*E0*x*(4 - 8*x + (4 + pow(e,2))*pow(x,2))*cos(Etau) + 
        (32 - 96*x + 96*pow(x,2) + 48*pow(e,2)*pow(x,2) - 32*pow(x,3) - 48*pow(e,2)*pow(x,3) + 3*pow(e,3)*pow(x,3)*cos(E0 - 3*Etau) + pow(e,3)*pow(x,3)*cos(3*E0 - 3*Etau) + 
           8*pow(e,2)*pow(x,2)*cos(2*E0 - 2*Etau) - 8*pow(e,2)*pow(x,3)*cos(2*E0 - 2*Etau) + 24*e*x*cos(E0 - Etau) - 48*e*pow(x,2)*cos(E0 - Etau) + 24*e*pow(x,3)*cos(E0 - Etau) + 
           6*pow(e,3)*pow(x,3)*cos(E0 - Etau) + 32*pow(e,2)*pow(x,2)*cos(2*Etau) - 32*pow(e,2)*pow(x,3)*cos(2*Etau) + 24*e*x*cos(E0 + Etau) - 48*e*pow(x,2)*cos(E0 + Etau) + 
           24*e*pow(x,3)*cos(E0 + Etau) + 6*pow(e,3)*pow(x,3)*cos(E0 + Etau) + 8*pow(e,2)*pow(x,2)*cos(2*(E0 + Etau)) - 8*pow(e,2)*pow(x,3)*cos(2*(E0 + Etau)) + 
           pow(e,3)*pow(x,3)*cos(3*(E0 + Etau)) + 3*pow(e,3)*pow(x,3)*cos(E0 + 3*Etau))*sin(E0)))/(32.*M_PI);
}
double fomega_function(double e, double x, double E0, double Etau)
{
    return -(sqrt(1 - pow(e,2))*x*sin(Etau)*(-48*E0 + 96*E0*x - 48*E0*pow(x,2) - 12*pow(e,2)*E0*pow(x,2) + 4*(6 - 12*x + (6 + pow(e,2))*pow(x,2))*sin(2*E0) + pow(e,2)*pow(x,2)*sin(4*E0) - 
        2*pow(e,2)*pow(x,2)*sin(2*E0 - 2*Etau) + pow(e,2)*pow(x,2)*sin(4*E0 - 2*Etau) - 24*e*x*sin(E0 - Etau) + 24*e*pow(x,2)*sin(E0 - Etau) + 8*e*x*sin(3*E0 - Etau) - 
        8*e*pow(x,2)*sin(3*E0 - Etau) - 24*e*x*sin(E0 + Etau) + 24*e*pow(x,2)*sin(E0 + Etau) - 2*pow(e,2)*pow(x,2)*sin(2*(E0 + Etau)) + 8*e*x*sin(3*E0 + Etau) - 
        8*e*pow(x,2)*sin(3*E0 + Etau) + pow(e,2)*pow(x,2)*sin(4*E0 + 2*Etau)))/(32.*M_PI);
}

/* Finite-size-related */
double ga_function(double e, double x, double E0)
{
    return (4*E0*x*(-8*(3 + (-3 + x)*x) + pow(e,2)*(12 + (-8 + pow(e,2))*pow(x,2))) + 64*sqrt(1 - pow(e,2))*atan(((1 + e)*tan(E0/2.))/sqrt(1 - pow(e,2))) + 
     e*x*(-16*(6 + pow(e,2)*(-3 + x) - 4*x)*x*sin(E0) + e*(8*(3 + x*(-6 + (2 + pow(e,2))*x))*sin(2*E0) + e*x*(-16*(-1 + x)*sin(3*E0) + 3*e*x*sin(4*E0)))))/(32.*M_PI);
}
double ge_function(double e, double x, double E0)
{
    return -((-1 + pow(e,2))*(12*E0*(-2 + pow(e,2)*x*(6 + x*(-9 + (4 + pow(e,2))*x))) + 48*sqrt(1 - pow(e,2))*atan(((1 + e)*tan(E0/2.))/sqrt(1 - pow(e,2))) - 
        6*e*(-4 + x*(24 + 8*(-3 + x)*x + pow(e,2)*x*(-15 + 14*x)))*sin(E0) + pow(e,2)*x*(6*(6 + x*(-15 + 2*(4 + pow(e,2))*x))*sin(2*E0) + e*x*(2*(9 - 10*x)*sin(3*E0) + 3*e*x*sin(4*E0))))
      )/(48.*e*M_PI);
}
double ha_function(double e, double x, double E0)
{
    return ((8*atan(sqrt(-((1 + e)/(-1 + e)))*tan(E0/2.)))/sqrt(1 - pow(e,2)) + e*(-12*x - (-4 + pow(e,2))*pow(x,3) - 4/(-1 + e*cos(E0)))*sin(E0) + 
     x*(-2*E0*(6 + x*(-6 + (2 + pow(e,2))*x)) - pow(e,2)*x*(-3*(-2 + x)*sin(2*E0) + e*x*sin(3*E0))))/(4.*M_PI);
}
double he_function(double e, double x, double E0)
{
    return (144*pow(1 - pow(e,2),1.5)*x*atan(sqrt(-((1 + e)/(-1 + e)))*tan(E0/2.)) + 48*sqrt(1 - pow(e,2))*(1 + 3*(-1 + pow(e,2))*x)*atan(((1 + e)*tan(E0/2.))/sqrt(1 - pow(e,2))) + 
     ((-1 + pow(e,2))*(-24*E0 - 36*pow(e,2)*E0*pow(x,2) + 24*pow(e,2)*E0*pow(x,3) - 12*e*E0*(-2 + pow(e,2)*pow(x,2)*(-3 + 2*x))*cos(E0) - 
          3*e*(-8 + 48*x - 3*(16 + 3*pow(e,2))*pow(x,2) + 2*(8 + 7*pow(e,2))*pow(x,3))*sin(E0) + 72*pow(e,2)*x*sin(2*E0) - 126*pow(e,2)*pow(x,2)*sin(2*E0) + 
          60*pow(e,2)*pow(x,3)*sin(2*E0) + 16*pow(e,4)*pow(x,3)*sin(2*E0) + 27*pow(e,3)*pow(x,2)*sin(3*E0) - 26*pow(e,3)*pow(x,3)*sin(3*E0) + 4*pow(e,4)*pow(x,3)*sin(4*E0))
        )/(-1 + e*cos(E0)))/(48.*e*M_PI);
}

double XL0_q_function(double q)
{
    double Aplus = pow(q*(54 + pow(q,2) + 6*sqrt(3)*sqrt(27 + pow(q,2))),0.3333333333333333);
    double Aminus = pow(q*(54 + pow(q,2) - 6*sqrt(3)*sqrt(27 + pow(q,2))),0.3333333333333333);
    return (3 + sqrt(3)*sqrt(3 + Aminus + Aplus - 2*q) - sqrt(3)*sqrt(6 - Aplus - 4*q - pow(q,2)/Aplus + (6*sqrt(3)*(1 + q))/sqrt(3 + Aminus + Aplus - 2*q)))/6.;
}


}
