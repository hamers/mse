/* MSE */

#include "evolve.h"
#include "binary_evolution.h"

extern "C"
{

//struct ktype__ ktype_;

int handle_binary_evolution(ParticlesMap *particlesMap, double t_old, double t, double *dt_binary_evolution, int *integration_flag)
{
    if (*integration_flag != 0) /* Only do binary evolution if we can be sure the system is dynamically stable and we have well-defined orbits. */
    {
        return 0;
    }
    update_structure(particlesMap);
    handle_mass_transfer(particlesMap,t_old,t,dt_binary_evolution,integration_flag);
    
    //*dt_binary_evolution = 1.0e4;
    
    return 0;
}



int handle_mass_transfer(ParticlesMap *particlesMap, double t_old, double t, double *dt_binary_evolution, int *integration_flag)
{
    *dt_binary_evolution = 1.0e100;
    double min_dt = 1.0e0;
    double dt = t - t_old; /* Imposed timestep which will be used for secular ODE integration. */

    double epsilon_MT = 0.01;
    double dt_temp;

    //int kw1=8;
    //int kw2=11;
    //int knew = ktype_(&kw1,&kw2) - 100;
    //int knew = ktype_(&kw1,&kw2) - 100;
    //int knew = determine_merger_type(kw1,kw2);
    //printf("knew %d\n",);

    ParticlesMapIterator it_p;
    for (it_p = particlesMap->begin(); it_p != particlesMap->end(); it_p++)
    {

        Particle *donor = (*it_p).second;
        printf("i %d %d %d\n",donor->index,donor->is_binary,donor->is_bound);
        if (donor->is_binary == false and donor->is_bound == true)
        {
            
            if (donor->RLOF_flag == 1)
            {
                int flag = handle_mass_transfer_cases(particlesMap, donor->parent, donor->index, donor->sibling, integration_flag, t_old, t, dt_binary_evolution);

                #ifdef IGNORE
                if (flag == 1)
                {
                    dt_temp = compute_orbital_period(parent);
                    common_envelope_evolution(particlesMap, parent->index, donor->index, accretor->index, integration_flag);
                }
                else /* TO DO: include low-mass dynamical mass transfer case (evolv2.f: lines 1150 - 1233 */
                {
                   dt_temp = compute_orbital_period(parent); /* TO DO: change this to something more appropriate */
                    compute_stable_mass_transfer_properties(particlesMap, parent->index, donor->index, accretor->index, integration_flag);
                }
                
                //dt_temp = epsilon_MT*fabs(donor->mass/donor->mass_dot_RLOF);
                //printf("0 dt_temp %g *dt_binary_evolution %g\n",dt_temp,*dt_binary_evolution);
                *dt_binary_evolution = min(*dt_binary_evolution, dt_temp);
                //printf("1 dt_temp %g *dt_binary_evolution %g\n",dt_temp,*dt_binary_evolution);
                *dt_binary_evolution = max(*dt_binary_evolution, min_dt);
                //printf("2 dt_temp %g *dt_binary_evolution %g\n",dt_temp,*dt_binary_evolution);
                //printf("fm %g M_dot %g dt_MT %g \n",fm,donor->mass_dot_RLOF,*dt_binary_evolution);
                #ifdef VERBOSE
                printf("fm %g M_dot %g dt_MT %g \n",fm,donor->mass_dot_RLOF,*dt_binary_evolution);
                #endif
    //printf("e %g x %g E_0 %g m_d %g m_a %g fm %g md %g\n",e,x,E_0,M_d,M_a,fm,M_d_dot_av);
                
                #endif
                
            }
        }

    }

    update_structure(particlesMap);
    
    return 0;
    
}

double compute_q_crit_for_common_envelope_evolution(int kw, double mass, double core_mass)
{
    double q_crit;
    
    if (kw == 2)
    {
        q_crit = 4.0;
    }
    else if (kw == 3 or kw == 5 or kw == 6) // Hjellming & Webbink, 1987, ApJ, 318, 794. 
    {
        q_crit = 0.362 + 1.0/(3.0 * (1.0 - core_mass/mass));
    }
    else if (kw == 8 or kw == 9)
    {
        q_crit = 0.784;
    }
    else
    {
        q_crit = 3.0;
    }
    
    return q_crit;
}


int handle_mass_transfer_cases(ParticlesMap *particlesMap, int parent_index, int donor_index, int accretor_index, int *integration_flag, double t_old, double t, double *dt_binary_evolution)
{
    Particle *parent = (*particlesMap)[parent_index];
    Particle *donor = (*particlesMap)[donor_index];
    Particle *accretor = (*particlesMap)[accretor_index];
    
    int flag = -1;
    
    double P_orb = compute_orbital_period(parent);
    double M_donor = donor->mass;
    double M_accretor = accretor->mass;

    /* MT & dynamical timescales */
    double fm = donor->fm;
    double m_dot = compute_orbit_averaged_mass_transfer_rate_emt_model(M_donor,fm,P_orb);
    double t_MT = M_donor/fabs(m_dot);
    double R_donor = donor->radius;
    double t_dyn = compute_stellar_dynamical_timescale(M_donor, R_donor);

    /* Critical mass ratios. */
    double q = M_donor / M_accretor;
    double q_crit = compute_q_crit_for_common_envelope_evolution(donor->stellar_type, M_donor, donor->core_mass);
    double q_crit_low_mass_donor = 0.695;
    double q_crit_WD_donor = 0.628;
    
    int kw = donor->stellar_type;

    /* Different cases */
    if (kw == 0 and q > q_crit_low_mass_donor)
    {
        /* `CE'-like evolution from a low-mass MS star to any secondary. */
        flag = 1;
        dynamical_mass_transfer_low_mass_donor(particlesMap, parent->index, donor->index, accretor->index, t_old, t, integration_flag);
    }
    else if ( (kw >= 2 and kw <= 9 and kw != 7) and ((t_MT <= t_dyn or q > q_crit) or t_MT < P_orb) ) /* add donor convective envelope mass? */
    {
        /* `Standard CE evolution. */
        flag = 2;
        common_envelope_evolution(particlesMap, parent->index, donor->index, accretor->index, integration_flag);
    }
    else if (kw >= 10 and kw <= 12 and q > q_crit_WD_donor)
    {
        /* Dynamical transfer from WD */
        flag = 3;
        dynamical_mass_transfer_WD_donor(particlesMap, parent->index, donor->index, accretor->index, t_old, t, integration_flag);
    }   
    else if (kw == 13)
    {
        flag = 4;
        mass_transfer_NS_donor(particlesMap, parent->index, donor->index, accretor->index, t_old, t, integration_flag);
    }
    else if (kw == 14)
    {
        flag = 5;
        mass_transfer_BH_donor(particlesMap, parent->index, donor->index, accretor->index, t_old, t, integration_flag);
    }
    else
    {
        /* Other cases: stable mass transfer. */
        flag = 6;
        stable_mass_transfer_evolution(particlesMap, parent->index, donor->index, accretor->index, t_old, t, integration_flag);
    }
    /* TO DO: contact cases! */
    
    if (flag == -1)
    {
        printf("binary_evolution.cpp -- fatal error in handle_mass_transfer_cases, parent %d donor %d accretor %d\n",parent->index,donor->index,accretor->index);
        exit(-1);
    }
    
    return flag;
}

double compute_orbit_averaged_mass_transfer_rate_emt_model(double M_donor, double fm, double P_orb)
{
    return -(M_donor/P_orb) * fm;
}

double compute_stellar_dynamical_timescale(double M, double R)
{
    return sqrt( ( R*R*R )/(CONST_G * M) );
}


int dynamical_mass_transfer_low_mass_donor(ParticlesMap *particlesMap, int parent_index, int donor_index, int accretor_index, double t_old, double t, int *integration_flag)
{
    /* Low-mass dynamical mass transfer case (BSE evolv2.f: lines 1150 - 1233) */
   return 0;
}

int dynamical_mass_transfer_WD_donor(ParticlesMap *particlesMap, int parent_index, int donor_index, int accretor_index, double t_old, double t, int *integration_flag)
{
    /* Dynamical mass transfer case WD donor (BSE evolv2.f: lines 1290 - 1342) */
    
    Particle *parent = (*particlesMap)[parent_index];
    Particle *donor = (*particlesMap)[donor_index];
    Particle *accretor = (*particlesMap)[accretor_index];
    
    int flag = -1;
    
    double P_orb = compute_orbital_period(parent);
    double M_donor = donor->mass;
    double M_accretor = accretor->mass;

    double fm = donor->fm;
    double m_dot = compute_orbit_averaged_mass_transfer_rate_emt_model(M_donor,fm,P_orb);
    
    double dm1; /* Mass lost from companion due to RLOF */
    double dm2; /* Mass accreted by companion. 
    * The `emt' model does not support non-conservative transfer.
    * Instead, any mass not accreted will be assumed to be lost in an isotropic wind. 
    * The latter is computed by adjusting the wind mass loss of the donor (effectively, it loses more mass in the wind). */
    
    double dt = t - t_old;
    dm1 = m_dot * dt;
    
    //double dme = 2.08d-03*eddfac*(1.d0/(1.d0 + zpars(11)))*rad(j2)*tb // TO DO: include case for eddfac<10
    
    dm2 = dm1;

    bool destroyed = false;
    int kw1 = donor->stellar_type;
    int kw2 = accretor->stellar_type;
    if (kw1 == 10 and kw2 == 10)
    {
        /* Assume the energy released by ignition of the triple-alpha reaction is enough to destroy both stars. */
        destroyed = true;
    }
    else if (kw1 == 10 or kw2 == 10)
    {
        /* Should be helium overflowing onto a CO or ONe core in which case the 
        *  helium swells up to form a giant envelope so a HeGB star is formed. 
        *  Allowance for the rare case of CO or ONe flowing onto He is made. */

        int kst = 9;
        double accretor_core_mass = accretor->core_mass; 
        if(kw2 == 10)
        {
            accretor_core_mass = dm2;
        }
        double accretor_age = accretor->age*yr_to_Myr;
        gntage_(&accretor_core_mass,&M_accretor,&kst,accretor->zpars,&accretor->sse_initial_mass,&accretor_age);
        accretor->stellar_type = kst;
        accretor->epoch = t - accretor_age * Myr_to_yr;
    }
    else if (kw2 <= 12)
    {
        accretor->sse_initial_mass = M_accretor;
        if (kw1 == 12 and kw2 == 11) /* Mixture of ONe and CO will result in an ONe product. */
        {
            accretor->stellar_type = 12;
        }
    }

    /* TO DO: handle the destruction... */
#ifdef IGNORE
         kstar(j1) = 15
         mass(j1) = 0.d0
*
* Might be a supernova that destroys the system.
*
         if(kstar(j2).le.11.and.mass(j2).gt.mch)then
            kstar(j2) = 15
            mass(j2) = 0.d0
         endif
         coel = .true.
         goto 135
#endif
    
    
   return 0;
}

int mass_transfer_NS_donor(ParticlesMap *particlesMap, int parent_index, int donor_index, int accretor_index, double t_old, double t, int *integration_flag)
{
    /* Mass transfer case NS donor (BSE evolv2.f: lines 1343 - 1354) */
   return 0;
}

int mass_transfer_BH_donor(ParticlesMap *particlesMap, int parent_index, int donor_index, int accretor_index, double t_old, double t, int *integration_flag)
{
    /* Mass transfer case BH donor (BSE evolv2.f: lines 1355 - 1365) */
   return 0;
}


int common_envelope_evolution(ParticlesMap *particlesMap, int binary_index, int index1, int index2, int *integration_flag)
{
    /*
     * **
      SUBROUTINE COMENV(M01,M1,MC1,AJ1,JSPIN1,KW1,
     &                  M02,M2,MC2,AJ2,JSPIN2,KW2,
     &                  ZPARS,ECC,SEP,JORB,COEL)
*
* Common Envelope Evolution.
* Mostly copied from BSE.
*
*     Author : C. A. Tout
*     Date :   18th September 1996
*
*     Redone : J. R. Hurley
*     Date :   7th July 1998
*
* Note units in SSE/BSE: length in RSUN, time in Myr, luminosity in LSun
*/

    printf("CE binary_index %d index1 %d index2 %d\n",binary_index,index1,index2);
    int i;

    Particle *star1 = (*particlesMap)[index1];
    Particle *star2 = (*particlesMap)[index2];
    double M1 = star1->mass;
    double M2 = star2->mass;

    /* CE evolution is not handled by the ODE integrator; correspondingly, mass_dot_RLOF should be set to zero. */
    star1->mass_dot_RLOF = 0.0;
    star2->mass_dot_RLOF = 0.0;

    /* Orbital properties */
    double SEP;
    double SEPF; // final separation
    double SEPL;
    double ECC;

    Particle *binary;
    double h_vec_unit[3],e_vec_unit[3];
    double initial_momentum[3];    
    if (*integration_flag == 0) /* secular mode */
    {
        set_up_derived_ODE_quantities(particlesMap); /* for setting a, e, etc. */
        binary = (*particlesMap)[binary_index];
        SEP = binary->a/CONST_R_SUN; // initial separation
        ECC = binary->e;
        get_unit_vector(binary->h_vec,h_vec_unit);
        get_unit_vector(binary->e_vec,e_vec_unit);
    }
    else /* Direct N-body mode */
    {
        double h_vec[3],e_vec[3],r[3],v[3];
        double true_anomaly;
        
        for (i=0; i<3; i++)
        {
            r[i] = star1->R_vec[i] - star2->R_vec[i];
            v[i] = star1->V_vec[i] - star2->V_vec[i];
            initial_momentum[i] = M1 * star1->V_vec[i] + M2 * star2->V_vec[i];
        }
        from_cartesian_to_orbital_vectors(M1,M2,r,v,e_vec,h_vec,&true_anomaly);
        ECC = norm3(e_vec);
        SEP = compute_a_from_h(M1,M2,norm3(h_vec),ECC);
        get_unit_vector(h_vec,h_vec_unit);
        get_unit_vector(e_vec,e_vec_unit);
        printf("CE Integration_flag>0 a %g e %g\n",SEP,ECC);
    }

    int CEFLAG = binary_evolution_CE_energy_flag;
    double ALPHA1 = binary->common_envelope_alpha;
    
    double fac = 1.0;
    
    bool COEL = false;
    
    /* Get current state of star 1 (primary) */
    double TM1,TN1;
    double *GB1,*TSCLS1,*LUMS1;
    GB1 = new double[10];
    TSCLS1 = new double[20];
    LUMS1 = new double[10];    
    
    int KW1 = star1->stellar_type;
    double M01 = star1->sse_initial_mass;
    double MC1 = star1->core_mass;
    double R1 = star1->radius/CONST_R_SUN;
    double RC1 = star1->core_radius/CONST_R_SUN;
    double L1 = star1->luminosity/CONST_L_SUN;
    double MENV1,RENV1,K21;
    //double K21;
    
    double AJ1 = star1->age * yr_to_Myr;
    double *ZPARS1 = star1->zpars;    
    //printf("t1 %g %d %g %g %g %g\n",ZPARS1[0],KW1,M01,M1,TM1,TN1);
    star_(&KW1, &M01, &M1, &TM1, &TN1, TSCLS1, LUMS1, GB1, ZPARS1);
    //printf("t2 %g %d %g %g %g %g\n",ZPARS1[0],KW1,M01,M1,TM1,TN1);
    //printf("arg hrd %g %g %g %g %g \n",M01,AJ1,M1,TM1,TN1);
    hrdiag_(&M01,&AJ1,&M1,&TM1,&TN1,TSCLS1,LUMS1,GB1,ZPARS1, \
        &R1,&L1,&KW1,&MC1,&RC1,&MENV1,&RENV1,&K21);

//      OSPIN1 = JSPIN1/(K21*R1*R1*(M1-MC1)+K3*RC1*RC1*MC1)
    double MENVD1 = MENV1 / (M1 - MC1);
    double RZAMS1 = rzamsf_(&M01);
    double LAMB1 = celamf_(&KW1,&M01,&L1,&R1,&RZAMS1,&MENVD1,&fac);

    /* Get current state of star 2 (secondary) */
    double TM2,TN2;
    double *GB2,*TSCLS2,*LUMS2;
    GB2 = new double[10];
    TSCLS2 = new double[20];
    LUMS2 = new double[10];    
    
    int KW2 = star2->stellar_type;
    double M02 = star2->sse_initial_mass;
    double MC2 = star2->core_mass;
    double R2 = star2->radius/CONST_R_SUN;
    double RC2 = star2->core_radius/CONST_R_SUN;
    double L2 = star2->luminosity/CONST_L_SUN;
    double MENV2,RENV2,K22;
    //double K22;
    
    double AJ2 = star2->age * yr_to_Myr;
    double *ZPARS2 = star2->zpars;    
    //printf("t3 %g %d %g %g %g %g\n",ZPARS2[0],KW2,M02,M2,TM2,TN2);
    star_(&KW2, &M02, &M2, &TM2, &TN2, TSCLS2, LUMS2, GB2, ZPARS2);
    //printf("t4 %g %d %g %g %g %g\n",ZPARS2[0],KW2,M02,M2,TM2,TN2);    
    
    hrdiag_(&M02,&AJ2,&M2,&TM2,&TN2,TSCLS2,LUMS2,GB2,ZPARS2, \
        &R2,&L2,&KW2,&MC2,&RC2,&MENV2,&RENV2,&K22);
//      OSPIN1 = JSPIN1/(K21*R1*R1*(M1-MC1)+K3*RC1*RC1*MC1)
    double MENVD2 = MENV2 / (M2 - MC2);
    double RZAMS2 = rzamsf_(&M02);
    double LAMB2 = celamf_(&KW2,&M02,&L2,&R2,&RZAMS2,&MENVD2,&fac);


    /*
    * Calculate the binding energy of the primary star's envelope (multiplied by lambda).
    * Also, calculate initial orbital energy.
    * Energies are divided by -G.
    */

    double EBINDI = M1 * (M1 - MC1)/(LAMB1 * R1);
    double EORBI = M1 * M2/(2.0 * SEP); // use primary total mass and secondary total mass (HTP02 eq. 69)

    /*
    * If the secondary star is also giant-like, add its envelope's energy.
    */
    
    if (KW2 >= 2 and KW2 <= 9 and KW2 != 7)
    {
        EBINDI += M2 * (M2 - MC2) / (LAMB2 * R2);

        if(CEFLAG != 3) // use primary core mass instead of total mass
        {
            EORBI = MC1 * MC2/(2.0 * SEP); // secondary is giant; use its core mass
        }
    }
    else
    {
        if(CEFLAG !=  3) // use primary core mass instead of total mass
        {
            EORBI = MC1 * M2/(2.0 * SEP); // secondary is not a giant; use its total mass
        }
    }
    
    double ECIRC = EORBI/(1.0 - ECC*ECC); // Allow for an eccentric orbit.

    double EORBF = EORBI + EBINDI/ALPHA1; // Calculate the final orbital energy without coalescence (HTP02 eq. 71).
    double EBINDF;
    
    /* Generic variables */

    double TN;
    double *GB,*TSCLS,*LUMS;
    GB = new double[10];
    TSCLS = new double[20];
    LUMS = new double[10];    
    double MENV,RENV;

    
    /* Variables for potentially merged stars */
    int KW;
    double MC3;
    double MF;
    
    double TB,OORB;
    double XX,CONST;
    double DELY,DERI,DELMF; // used for the NR iteration procedure
    bool iterate = true;

    
    /* By default, do not apply any kicks; can be changed below depending on situation 
     * Note: apply_kick will be set back to the default true at the end below */
    star1->apply_kick = false;
    star1->instantaneous_perturbation_delta_mass = 0.0;
    
    star2->apply_kick = false;
    star2->instantaneous_perturbation_delta_mass = 0.0;

    
    if(KW2 <= 1 or KW2 == 7)
    {
        /*
        * If the secondary is on the main sequence, see if it fills its Roche lobe.
        */
        
        SEPF = MC1 * M2/(2.0 * EORBF); // secondary does not have a core; assume it did not lose mass
        double Q1 = MC1 / M2;
        double Q2 = 1.0 / Q1;
         
        double RL1 = roche_radius_pericenter_eggleton(SEP,Q1) / SEP; // dimensionless RL radius (SEP cancels out)
        double RL2 = roche_radius_pericenter_eggleton(SEP,Q2) / SEP; // dimensionless RL radius (SEP cancels out)

        if (RC1/RL1 >= R2/RL2)
        {

            /* The helium core of a very massive star of type 4 may actually fill
            * its Roche lobe in a wider orbit with a very low-mass secondary. */
            
            if (RC1 > RL1 * SEPF)
            {
               COEL = true;
               SEPL = RC1 / RL1; // separation at which the stars coelesced
            }
        }
        else
        {
            if (R2 > RL2 * SEPF)
            {
               COEL = true;
               SEPL = R2 / RL2;
            }
        }

        if (COEL == true) /* Coalescence - calculate final binding energy. */
        {
            //KW = KTYPE(KW1,KW2) - 100
            KW = determine_merger_type(KW1,KW2);
            MC3 = MC1;
            if (KW2 == 7 and KW == 4)
            {
                MC3 = MC1 + M2; // note ASH: replaced `MC3' in BSE code to `MC1' -- does not change anything functionally, but is more clear
            }
            
            EORBF = max(MC1 * M2/(2.0 * SEPL) , EORBI);
            EBINDF = EBINDI - ALPHA1 * (EORBF - EORBI);
        }
        else // Primary becomes a black hole, neutron star, white dwarf or helium star.
        {
            MF = M1;
            M1 = MC1;
            star_(&KW1,&M01,&M1,&TM1,&TN,TSCLS1,LUMS,GB,ZPARS1);
            hrdiag_(&M01,&AJ1,&M1,&TM1,&TN,TSCLS1,LUMS,GB,ZPARS1, \
                &R1,&L1,&KW1,&MC1,&RC1,&MENV1,&RENV1,&K21);
            
            if (KW1 >= 13)
            {
                star1->apply_kick = true;
                star1->kick_distribution = 1;
                //star1->instantaneous_perturbation_delta_mass = M1 - MF;
            }
            
            // TO DO: include kicks 
            //IF(KW1.GE.13)THEN
               //CALL kick(KW1,MF,M1,M2,ECC,SEPF,JORB,VS)
               //IF(ECC.GT.1.D0) GOTO 30
            //ENDIF
         }
    }
    else
    {
        /*
        * Degenerate or giant secondary. Check if the least massive core fills its Roche lobe.
        */

        SEPF = MC1 * MC2/(2.0 * EORBF);
        double Q1 = MC1 / MC2;
        double Q2 = 1.0 / Q1;
        
        double RL1 = roche_radius_pericenter_eggleton(SEP,Q1) / SEP; // dimensionless RL radius (SEP cancels out)
        double RL2 = roche_radius_pericenter_eggleton(SEP,Q2) / SEP; // dimensionless RL radius (SEP cancels out)
        
        if (RC1/RL1 >= RC2/RL2)
        {
            if (RC1 > RL1 * SEPF)
            {
                COEL = true;
                SEPL = RC1 / RL1;
            }
        }
        else
        {
            if (RC2 > RL2 * SEPF)
            {
                COEL = true;
                SEPL = RC2 / RL2;
            }
        }

        if (COEL == true)
        {
            
            SEPF = 0.0; // ASH: redundant?
            if(KW2 >= 13)
            {
            /* If the secondary was a neutron star or black hole the outcome
            * is an unstable Thorne-Zytkow object that leaves only the core. */

               MC1 = MC2;
               M1 = MC1;
               MC2 = 0.0;
               M2 = 0.0;
               KW1 = KW2;
               KW2 = 15;
               AJ1 = 0.0;
        //* The envelope mass is not required in this case.
        //*
               goto label30; // TO DO: proper handling
            }
            //ENDIF
            //else // secondary was NOT an NS/BH
                
            //KW = KTYPE(KW1,KW2) - 100
            KW = determine_merger_type(KW1,KW2);
            MC3 = MC1 + MC2;

            /*
            * Calculate the final envelope binding energy.
            */
            
            EORBF = max(MC1 * MC2/(2.0 * SEPL), EORBI);
            EBINDF = EBINDI - ALPHA1*(EORBF - EORBI);

            /*
            * Check if we have the merging of two degenerate cores and if so
            * then see if the resulting core will survive or change form.
            */
            
            if (KW1 == 6 and (KW2 == 6 or KW2 >= 11))
            {
               dgcore_(&KW1,&KW2,&KW,&MC1,&MC2,&MC3,&EBINDF);
            }

            if(KW1 <= 3 and M01 <= ZPARS1[1])
            {
                if ( (KW2 >= 2 and KW2 <= 3 and M02 <= ZPARS2[1]) or KW2 == 10)
                {
                    dgcore_(&KW1,&KW2,&KW,&MC1,&MC2,&MC3,&EBINDF);
                    if(KW >= 10) // new star is compact object
                    {
                        KW1 = KW;
                        M1 = MC3;
                        MC1 = MC3;
                        if (KW < 15)
                        {
                            M01 = MC3;
                        }
                        AJ1 = 0.0;
                        MC2 = 0.0;
                        M2 = 0.0;
                        KW2 = 15;
                        
                        goto label30;
                        //GOTO 30 // TO DO: proper handling
                    }
                }
            }
        }
        else
        {
            
            /*
            * The cores do not coalesce - assign the correct masses and ages.
            */
            
            MF = M1;
            M1 = MC1;
            star_(&KW1,&M01,&M1,&TM1,&TN,TSCLS1,LUMS,GB,ZPARS1);
            hrdiag_(&M01,&AJ1,&M1,&TM1,&TN,TSCLS1,LUMS,GB,ZPARS1, \
                &R1,&L1,&KW1,&MC1,&RC1,&MENV1,&RENV1,&K21);
            if (KW1 >= 13)
            {
                star1->apply_kick = true;
                star1->kick_distribution = 1;
                // TO DO: handle kick
               //CALL kick(KW1,MF,M1,M2,ECC,SEPF,JORB,VS)
               //IF(ECC.GT.1.D0) GOTO 30
            }
            
            MF = M2;
            KW = KW2;
            M2 = MC2;
            star_(&KW2,&M02,&M2,&TM2,&TN,TSCLS2,LUMS,GB,ZPARS2);
            hrdiag_(&M02,&AJ2,&M2,&TM2,&TN,TSCLS2,LUMS,GB,ZPARS2, \
                &R2,&L2,&KW2,&MC2,&RC2,&MENV2,&RENV2,&K22);
            if(KW2 >= 13 and KW < 13) /* secondary became an NS */
            {
                star2->apply_kick = true;
                star2->kick_distribution = 1;
                // TO DO: handle kick
               //CALL kick(KW2,MF,M2,M1,ECC,SEPF,JORB,VS)
               //IF(ECC.GT.1.D0) GOTO 30
            }
            
        }
    } // back at main level
    
    double MC22; // effective core mass of the secondary
    double FAGE1,FAGE2;


    /* If making a single helium burning star, calculate the fractional age 
     * depending on the amount of helium that has burnt. */
    if (COEL == true)
    {
        MC22 = MC2;
        if(KW == 4 or KW == 7)
        {
            if (KW1 <= 3)
            {
               FAGE1 = 0.0;
            }
            else if (KW1 >= 6)
            {
               FAGE1 = 1.0;
            }
            else
            {
               FAGE1 = (AJ1 - TSCLS1[1])/(TSCLS1[12] - TSCLS1[1]);
            }

            if (KW2 <= 3 or KW2 == 10)
            {
               FAGE2 = 0.0;
            }
            else if (KW2 == 7)
            {
               FAGE2 = AJ2/TM2;
               MC22 = M2;
            }
            else if (KW2 >= 6)
            {
               FAGE2 = 1.0;
            }
            else
            {
               FAGE2 = (AJ2 - TSCLS2[1])/(TSCLS2[12] - TSCLS2[1]);
            }
        }
        //printf("COEL True FAGE1 %g FAGE 2 %g\n",FAGE1,FAGE2);
    }

    /* Now calculate the final mass following coelescence.  This requires a Newton-Raphson iteration. 
    * Note ASH: as a temporary measure, taking ZPARS2 for the coelesced object */

   
    if (COEL == true)
    {

        /* Calculate the orbital spin just before coalescence. */
         //TB = (SEPL/AURSUN)*SQRT(SEPL/(AURSUN*(MC1+MC2)))
        TB = TWOPI * (SEPL * CONST_R_SUN) * sqrt( SEPL * CONST_R_SUN / (CONST_G * (MC1 + MC2) ) );
        OORB = TWOPI/TB;

        /* Set up Newton-Raphson iteration */
        XX = 1.0 + ZPARS2[6];

        if(EBINDF <= 0.0)
        {
           MF = MC3;
           iterate = false;
        }
        else
        {
           CONST = pow(M1 + M2, XX) * (M1 - MC1 + M2 - MC22) * (EBINDF/EBINDI);
        }
    
        if (iterate == true)
        {
            /* Initial Guess. */
            MF = max( MC1 + MC22, (M1 + M2) * pow(EBINDF/EBINDI, 1.0/XX) );

            /* Start iteration */
            while(true)
            {
                DELY = pow(MF, XX) * (MF - MC1 - MC22) - CONST;
                
    //            if (fabs( pow(DELY/MF, 1.D0 + XX) <= 1.0e-2)) /* Note ASH: this case is commented out in BSE */
    //            {
    //                break;
    //            }
                if (fabs(DELY/MF) <= 1.0e-3)
                {
                    break;
                }

                DERI = pow(MF, ZPARS2[6]) * ( (1.0 + XX) * MF - XX * (MC1 + MC22) );
                DELMF = DELY/DERI;
                MF = MF - DELMF;
            }
        }
        
        /* Set the masses and separation. */
        if (MC22 == 0.0)
        {
            MF = max(MF, MC1 + M2);
        }
        M2 = 0.0;
        M1 = MF;
        KW2 = 15;

        /* Combine the core masses. */

        if (KW == 2)
        {
            star_(&KW,&M1,&M1,&TM2,&TN,TSCLS2,LUMS,GB,ZPARS2);
            if(GB[8] >= MC1)
            {
               M01 = M1;
               AJ1 = TM2 + (TSCLS2[0] - TM2) * (AJ1 - TM1) / (TSCLS1[0] - TM1);
               star_(&KW,&M01,&M1,&TM1,&TN,TSCLS1,LUMS,GB,ZPARS2);
            }
        }
        else if (KW == 7)
        {
            M01 = M1;
            star_(&KW,&M01,&M1,&TM1,&TN,TSCLS1,LUMS,GB,ZPARS2);
            AJ1 = TM1 * (FAGE1 * MC1 + FAGE2 * MC22) / (MC1 + MC22);
        }
        else if (KW == 4 or MC2 > 0.0 or KW != KW1)
        {
            if (KW == 4)
            {
                AJ1 = (FAGE1 * MC1 + FAGE2 * MC22) / (MC1 + MC22); /* Note ASH: this is actually an age factor (dimensionless); AJ1 will be changed below in gntage_ */
            }
            MC1 = MC1 + MC2;
            MC2 = 0.0;

            /* Obtain a new age for the giant. */

            gntage_(&MC1,&M1,&KW,ZPARS2,&M01,&AJ1);
            star_(&KW,&M01,&M1,&TM1,&TN,TSCLS1,LUMS,GB,ZPARS2);
         }
         
         hrdiag_(&M01,&AJ1,&M1,&TM1,&TN,TSCLS1,LUMS,GB,ZPARS2, \
            &R1,&L1,&KW,&MC1,&RC1,&MENV1,&RENV1,&K21);
         //JSPIN1 = OORB*(K21*R1*R1*(M1-MC1)+K3*RC1*RC1*MC1)
         KW1 = KW;
         ECC = 0.0;
    }
    else if (COEL == false) /* No coalescence */
    {
        
        /* Check if any eccentricity remains in the orbit by first using 
        * energy to circularise the orbit before removing angular momentum. 
        * (note this should not be done in case of CE SN ... fix).  */

         if(EORBF < ECIRC)
         {
            ECC = sqrt(1.0 - EORBF/ECIRC);
         }
         else
         {
            ECC = 0.0;
         }

         /* Set both cores in co-rotation with the orbit on exit of CE */

         //TB = (SEPF/AURSUN)*SQRT(SEPF/(AURSUN*(M1+M2)))
         TB = TWOPI * (SEPF * CONST_R_SUN) * sqrt( SEPF * CONST_R_SUN / (CONST_G * (M1 + M2) ) );
         OORB = TWOPI/TB;
         //double JORB = M1 * M2 / (M1 + M2) * sqrt(1.0 - ECC * ECC) * SEPF * SEPF * OORB;
//*        JSPIN1 = OORB*(K21*R1*R1*(M1-MC1)+K3*RC1*RC1*MC1)
//*        JSPIN2 = OORB*(K22*R2*R2*(M2-MC2)+K3*RC2*RC2*MC2)

        /* or, leave the spins of the cores as they were on entry.
        * Tides will deal with any synchronization later. */

//         JSPIN1 = OSPIN1 * (K21*R1*R1*(M1-MC1)+K3*RC1*RC1*MC1);
//         JSPIN2 = OSPIN2 * (K22*R2*R2*(M2-MC2)+K3*RC2*RC2*MC2);

    }



    /* Handle orbital changes / merge binaries in case of coalescence */
    
    label30:
    {
    
    bool unbound_orbits = false;
    double spin_vec_unit[3];
    printf("COEL %d\n",COEL);
    if (COEL == false)
    {
        
        /* Handle effect of mass loss on the rest of the system */
        printf("no coel M1 %g M2 %g kw1 %d kw2 %d\n",M1,M2,KW1,KW2);

        /* First, handle effects of possible SNe kicks. Mass loss is not yet taken into account */
        star1->instantaneous_perturbation_delta_mass = 0.0;
        star2->instantaneous_perturbation_delta_mass = 0.0;
        
        handle_SNe_in_system(particlesMap, &unbound_orbits, integration_flag);

        set_old_parameters_for_adiabatic_mass_loss(particlesMap); /* copy the old masses and h & e vectors for use of adiabatic mass loss below */
        
        /* Next, assume instantaneous mass loss for ALL orbits. 
         * Note: `binary' will be overwritten at the very end below. 
         * Other orbits, for which mass loss is expected to be adiabatic,
         * are adjusted below in compute_new_orbits_assuming_adiabatic_mass_loss. */
        star1->apply_kick = false; /* The kicks were already taken into account above */
        star2->apply_kick = false; /* The kicks were already taken into account above */

        double Delta_M1 = M1 - star1->mass;
        double Delta_M2 = M2 - star2->mass;
        star1->instantaneous_perturbation_delta_mass = Delta_M1; // final minus initial mass
        star2->instantaneous_perturbation_delta_mass = Delta_M2; // final minus initial mass        
        
        handle_SNe_in_system(particlesMap, &unbound_orbits, integration_flag); /* this will update the masses of stars 1 and 2 */
        
        /* Lastly, `correct' the orbits (h & e vectors) for which the effect of mass loss was actually expected to be adiabatic,
         * i.e., time of CE >> t_orb */
        star1->delta_m_adiabatic_mass_loss = Delta_M1;
        star2->delta_m_adiabatic_mass_loss = Delta_M2;
        compute_new_orbits_assuming_adiabatic_mass_loss(particlesMap, binary->common_envelope_timescale);
        
        //apply_instantaneous_mass_changes_and_kicks(particlesMap, integration_flag);
        printf("Ce -- post mass changes\n");
        print_system(particlesMap);
        /* Update the stars and orbit */
        star1->stellar_type = KW1;
        
        //star1->mass = M1; // handled by SNe function
        star1->core_mass = MC1;
        star1->sse_initial_mass = M01;
        star1->convective_envelope_mass = MENV1;
        
        star1->epoch = 0.0; /* is this correct? */
        star1->age = AJ1 * Myr_to_yr;

        star1->radius = R1 * CONST_R_SUN;
        star1->core_radius = RC1 * CONST_R_SUN;
        star1->convective_envelope_radius = RENV1 * CONST_R_SUN;
        
        star1->luminosity = L1 * CONST_L_SUN;
        
        /* Spins: either assume the spins do not change (binary_evolution_CE_spin_flag=0) */
        /* Alternatively: spins corotate with final orbit (binary_evolution_CE_spin_flag=1) */
        if (binary_evolution_CE_spin_flag == 1)
        {
            get_unit_vector(star1->spin_vec,spin_vec_unit);
            for (int i=0; i<3; i++)
            {
                star1->spin_vec[i] = OORB * spin_vec_unit[i];
            }
        }

        star2->stellar_type = KW2;
        
        //star2->mass = M2; // handled by SNe function
        star2->core_mass = MC2;
        star2->sse_initial_mass = M02;
        star2->convective_envelope_mass = MENV2;
        
        star2->epoch = 0.0; /* is this correct? */
        star2->age = AJ2 * Myr_to_yr;

        star2->radius = R2 * CONST_R_SUN;
        star2->core_radius = RC2 * CONST_R_SUN;
        star2->convective_envelope_radius = RENV2 * CONST_R_SUN;
        
        star2->luminosity = L2 * CONST_L_SUN;

        if (binary_evolution_CE_spin_flag == 1)
        {
            get_unit_vector(star2->spin_vec,spin_vec_unit);
            for (int i=0; i<3; i++)
            {
                star2->spin_vec[i] = OORB * spin_vec_unit[i];
            }
        }

        /* Set new e & h vectors for the binary.
         * Assume they only change in magnitude.
         * Note: this overrides the changes made to the binary's orbit after calling handle_SNe_in_system above */

        double a = SEPF * CONST_R_SUN;
        printf("a %g\n",a);
        double e = ECC;
        double h = compute_h_from_a(M1,M2,a,e);

        double e_vec_unit[3],h_vec_unit[3];
        if (*integration_flag == 0) /* secular */
        {
            for (int i=0; i<3; i++)
            {
                binary->e_vec[i] = e * e_vec_unit[i];
                binary->h_vec[i] = h * h_vec_unit[i];
            }

            /* Below: not strictly neccesary */
            binary->a = a;
            binary->e = ECC;
        
            /* Assume the binary true anomaly is not affected */
            set_binary_masses_from_body_masses(particlesMap);
            set_positions_and_velocities(particlesMap); /* update positions and velocities of all bodies based on updated binary e and h vectors */
        }
    
    }
    else if (COEL == true)
    {
        if (*integration_flag == 0) /* secular */
        {
            binary->is_binary = false; /* The binary becomes a body */
            particlesMap->erase(star1->index);
            particlesMap->erase(star2->index);


            /* Handle effect of mass loss on the rest of the system */

            set_old_parameters_for_adiabatic_mass_loss(particlesMap); /* copy the old masses and h & e vectors for use of adiabatic mass loss below */
            
            /* Next, assume instantaneous mass loss for ALL orbits. 
             * Orbits, for which mass loss is expected to be adiabatic,
             * are adjusted below in compute_new_orbits_assuming_adiabatic_mass_loss. */
            binary->apply_kick = false;

            binary->mass = star1->mass + star2->mass; // the old mass
            double Delta_M = M1 - binary->mass;
            binary->instantaneous_perturbation_delta_mass = M1 - binary->mass;

            handle_SNe_in_system(particlesMap, &unbound_orbits, integration_flag); /* mass of `binary' will be updated */

            //apply_instantaneous_mass_changes_and_kicks(particlesMap, integration_flag);
            
            //apply_user_specified_instantaneous_perturbation(particlesMap);
            //reset_instantaneous_perturbation_quantities(particlesMap);
            
            //printf("mold %g mnew %g\n",binary->mass,M1);
            
            /* Lastly, `correct' the orbits (h & e vectors) for which the effect of mass loss was actually expected to be adiabatic,
             * i.e., time of CE >> t_orb */
            binary->delta_m_adiabatic_mass_loss = Delta_M;
            compute_new_orbits_assuming_adiabatic_mass_loss(particlesMap, binary->common_envelope_timescale);


            //printf("post merge\n");
            //print_system(particlesMap);
            /* Should coelesced objects get a kick? 
             * For now, assume no kick
             * Note: mergers between BHs/NSs should never occur here (CE evolution), rather as collisions */

            binary->apply_kick = false;
            //printf("done SNe\n");
            
            /* Update the merged star */
            binary->stellar_type = KW1;
            
            //star1->mass = M1; // handled by SNe function
            binary->core_mass = MC1;
            binary->sse_initial_mass = M01;
            binary->convective_envelope_mass = MENV1;
            
            binary->epoch = 0.0; /* is this correct? */
            binary->age = AJ1 * Myr_to_yr;

            binary->radius = R1 * CONST_R_SUN;
            binary->core_radius = RC1 * CONST_R_SUN;
            binary->convective_envelope_radius = RENV1 * CONST_R_SUN;
            
            binary->luminosity = L1 * CONST_L_SUN;
            
            /* Set the spin equal to the orbital frequency just before coalescence (previously calculated as OORB).
             * Assume the direction is equal to the previous orbital orientation. */

            for (int i=0; i<3; i++)
            {
                binary->spin_vec[i] = OORB * spin_vec_unit[i];
            }
            
            binary->apply_kick = true; /* Default value for (stellar evolution) bodies */
        }
        else
        {
            /* Assign merger remnant to star1; remove star2 */
            particlesMap->erase(star2->index);
            star1->stellar_type = KW1;

            for (i=0; i<3; i++)
            {
                star1->V_vec[i] = initial_momentum[i]/M1; /* set new velocity according to linear momentum conservation */
            }

            star1->mass = M1;
            star1->core_mass = MC1;
            star1->sse_initial_mass = M01;
            star1->convective_envelope_mass = MENV1;
            
            star1->epoch = 0.0; /* is this correct? */
            star1->age = AJ1 * Myr_to_yr;

            star1->radius = R1 * CONST_R_SUN;
            star1->core_radius = RC1 * CONST_R_SUN;
            star1->convective_envelope_radius = RENV1 * CONST_R_SUN;
            
            star1->luminosity = L1 * CONST_L_SUN;
            
            /* Set the spin equal to the orbital frequency just before coalescence (previously calculated as OORB).
             * Assume the direction is equal to the previous orbital orientation. */

            for (int i=0; i<3; i++)
            {
                star1->spin_vec[i] = OORB * spin_vec_unit[i];
            }
            
            star1->apply_kick = true; /* Default value for (stellar evolution) bodies */
        }
    }
    
    /* Reset kick parameters if the binary remains */
    if (COEL == false)
    {
        star1->apply_kick = true;
        star2->apply_kick = true;
    }

    } // end of label30

    return 1;
}

void set_old_parameters_for_adiabatic_mass_loss(ParticlesMap *particlesMap)
{
    int i;
    ParticlesMapIterator it;
    for (it = particlesMap->begin(); it != particlesMap->end(); it++)
    {
        Particle *b = (*it).second;
        b->delta_m_adiabatic_mass_loss = 0.0;
        if (b->is_binary == true)
        {
            b->P_orb_adiabatic_mass_loss = compute_orbital_period(b);
            b->child1_mass_adiabatic_mass_loss = b->child1_mass;
            b->child2_mass_adiabatic_mass_loss = b->child2_mass;
            for (i=0; i<3; i++)
            {
                b->h_vec_old_adiabatic_mass_loss[i] = b->h_vec[i];
                b->e_vec_old_adiabatic_mass_loss[i] = b->e_vec[i];
            }
        }
    }
}

void compute_new_orbits_assuming_adiabatic_mass_loss(ParticlesMap *particlesMap, double mass_loss_timescale)
{
    set_binary_masses_from_body_masses(particlesMap); // to set the correct delta_m_adiabatic for each orbit

    int i;    
    ParticlesMapIterator it;
    for (it = particlesMap->begin(); it != particlesMap->end(); it++)
    {
        Particle *b = (*it).second;
        if (b->is_binary == true)
        {
            if (b->P_orb_adiabatic_mass_loss < mass_loss_timescale) // criterion for adiabatic mass loss
            {
                double m1_old = b->child1_mass_adiabatic_mass_loss;
                double m2_old = b->child2_mass_adiabatic_mass_loss;
                double delta_m1 = b->delta_child1_mass_adiabatic_mass_loss;
                double delta_m2 = b->delta_child2_mass_adiabatic_mass_loss;
                //double a_old = compute_a_from_h(m1_old,m2_old, norm3(b->h_vec_old_adiabatic_mass_loss), norm3(b->e_vec_old_adiabatic_mass_loss));
                //double delta_a = - (b->delta_m_adiabatic_mass_loss/(m1_old + m2_old)) * a_old;
                //double a = a_old + delta_a;
                //double h_old = compute_h_from_a(m1_old,m2_old,a,e);
                double h_old = norm3(b->h_vec_old_adiabatic_mass_loss);
                double h_new = h_old * ( (delta_m1/m1_old) + (delta_m2/m2_old) - (delta_m1 + delta_m2)/(m1_old + m2_old) );
                for (i=0; i<3; i++)
                {
                    /* Directions of h and e and magnitude of e do not change for adiabatic mass loss. */
                    b->h_vec[i] = (b->h_vec_old_adiabatic_mass_loss[i] / h_old) * h_new;
                    b->e_vec[i] = b->e_vec_old_adiabatic_mass_loss[i];
                }
            }
        }
    }
}

int common_envelope_evolution2(ParticlesMap *particlesMap, int binary_index, int donor_index, int accretor_index)
{
    /* Old function; can be ignored */
    Particle *binary = (*particlesMap)[binary_index];
    Particle *donor = (*particlesMap)[donor_index];
    Particle *accretor = (*particlesMap)[accretor_index];
    
    double M_donor = donor->mass;
    double M_core_donor = donor->core_mass;
    double M_envelope_donor = M_donor - M_core_donor; // do not take convective envelope mass
    double M_accretor = accretor->mass;
    
    double alpha = binary->common_envelope_alpha;
    double lambda = donor->common_envelope_lambda;
    
    double a_i = binary->a;
    double e = binary->e;
    double rp_i = a_i * (1.0 - e);
    
    double q = M_donor / M_accretor;
    double r_L_d = roche_radius_pericenter_eggleton(rp_i,q) / a_i; // take R_L at periapsis
    
    double a_f = a_i * (M_core_donor / M_donor) * M_accretor / ( M_accretor + 2.0 * M_envelope_donor / (alpha * lambda * r_L_d) );

    binary->a = a_f;
    binary->e = epsilon; // assume circularisation after CE
    donor->mass = M_core_donor; // what about stellar type?
    
    /* calc new e & h vectors; update R&V */
    return 0; 
}

int stable_mass_transfer_evolution(ParticlesMap *particlesMap, int parent_index, int donor_index, int accretor_index, double t_old, double t, int *integration_flag)
{
    /* Stable mass transfer case (BSE evolv2.f: lines 1370 - 1585) */
    

   //                 donor->mass_dot_RLOF = m_dot;
//                    accretor->mass_dot_RLOF = -m_dot;

    return 0;
}

}
