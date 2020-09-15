/* MSE */

#include "evolve.h"
#include "common_envelope_evolution.h"

extern "C"
{

int common_envelope_evolution(ParticlesMap *particlesMap, int binary_index, int index1, int index2, double t, int *integration_flag)//, ParticlesMapIterator &it_p)
{
    //printf("binary_evolution.cpp -- common_envelope_evolution -- start; binary_index %d index1 %d index2 %d integration_flag %d\n",binary_index,index1,index2,*integration_flag);
    //int i;

    Particle *star1 = (*particlesMap)[index1];
    Particle *star2 = (*particlesMap)[index2];
    
    if (star2->is_binary == false)
    {
        /* Both objects are physical stars. Apply the "standard" CE scheme from HTP02. */
        binary_common_envelope_evolution(particlesMap, binary_index, index1, index2, t, integration_flag);
    }
    else
    {
        /* A tertiary star is going into CE with a companion binary. */
        triple_common_envelope_evolution(particlesMap, binary_index, index1, index2, t, integration_flag);
    }
    
    return 0;
}

void binary_common_envelope_evolution(ParticlesMap *particlesMap, int binary_index, int index1, int index2, double t, int *integration_flag)//, ParticlesMapIterator &it_p)
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



    #ifdef LOGGING
    Log_info_type log_info;
    log_info.binary_index = binary_index;
    log_info.index1 = index1;
    log_info.index2 = index2;
    update_log_data(particlesMap, t, *integration_flag, 5, log_info);
    #endif

    printf("common_envelope_evolution.cpp -- binary_common_envelope_evolution() -- start; binary_index %d index1 %d index2 %d integration_flag %d\n",binary_index,index1,index2,*integration_flag);
    int i;

    Particle *star1 = (*particlesMap)[index1];
    Particle *star2 = (*particlesMap)[index2];

    if (star1->evolve_as_star == true and star2->evolve_as_star == false)
    {
        collision_product_star_planet(particlesMap, binary_index, index1, index2,t,integration_flag);
        return;
    }
    if (star1->evolve_as_star == false and star2->evolve_as_star == true)
    {
        collision_product_star_planet(particlesMap, binary_index, index2, index1,t,integration_flag);
        return;
    }
    if (star1->evolve_as_star == false and star2->evolve_as_star == false)
    {
        collision_product_planet_planet(particlesMap, binary_index, index1, index2,t,integration_flag);
        return;
    }
    
    double M1_old,M1;
    M1_old = M1 = star1->mass;
    double M2_old,M2;
    M2_old = M2 = star2->mass;

    /* Orbital properties */
    double SEP;
    double SEPF; // final separation
    double SEPL;
    double ECC;

    Particle *binary;
    double h_vec_unit[3],e_vec_unit[3];
    double initial_momentum[3],initial_R_CM[3];
    if (*integration_flag == 0) /* secular mode */
    {
        set_up_derived_quantities(particlesMap); /* for setting a, e, etc. */
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
            initial_R_CM[i] = (M1 * star1->R_vec[i] + M2 * star2->R_vec[i]) / (M1 + M2);
        }
        from_cartesian_to_orbital_vectors(M1,M2,r,v,e_vec,h_vec,&true_anomaly);
        ECC = norm3(e_vec);
        SEP = compute_a_from_h(M1,M2,norm3(h_vec),ECC)/CONST_R_SUN;
        get_unit_vector(h_vec,h_vec_unit);
        get_unit_vector(e_vec,e_vec_unit);
        printf("CE Integration_flag>0 a %g e %g\n",SEP,ECC);
    }

    int CEFLAG = binary_evolution_CE_energy_flag;
    double ALPHA1 = star1->common_envelope_alpha;
    
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

    
    /* By default, do not apply any kicks; can be changed below depending on new stellar types. */
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
               SEPL = RC1 / RL1; // separation at which the stars coalesced
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
            
            EORBF = CV_max(MC1 * M2/(2.0 * SEPL) , EORBI);
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
                //star1->kick_distribution = 1;
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
               goto label30;
            }
            //ENDIF
            //else // secondary was NOT an NS/BH
                
            //KW = KTYPE(KW1,KW2) - 100
            KW = determine_merger_type(KW1,KW2);
            MC3 = MC1 + MC2;

            /*
            * Calculate the final envelope binding energy.
            */
            
            EORBF = CV_max(MC1 * MC2/(2.0 * SEPL), EORBI);
            EBINDF = EBINDI - ALPHA1*(EORBF - EORBI);

            /*
            * Check if we have the merging of two degenerate cores and if so
            * then see if the resulting core will survive or change form.
            */
            
            if (KW1 == 6 and (KW2 == 6 or KW2 >= 11))
            {
               dgcore_(&KW1,&KW2,&KW,&MC1,&MC2,&MC3,&EBINDF); /* ASH: dgcore_ only potentially changes KW (third argument) */
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
                //star1->kick_distribution = 1;
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
                //star2->kick_distribution = 1;
            }
            
        }
    } // back at main level
    
    double MC22; // Effective core mass of the secondary
    double FAGE1,FAGE2; // Age factors


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
    * Note ASH: as a temporary measure, taking metallicity of star2 for the coelesced object */

    double z_new;
    double *zpars_new;
    zpars_new = new double[20];

    if (COEL == true)
    {
        z_new = star2->metallicity;
        zcnsts_(&z_new,zpars_new);

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
            MF = CV_max( MC1 + MC22, (M1 + M2) * pow(EBINDF/EBINDI, 1.0/XX) );

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
            MF = CV_max(MF, MC1 + M2);
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

            /* Obtain a new age and initial mass for the giant. */

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
            ECC = epsilon;
            //ECC = 1.0e-8;
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
    //double spin_vec_unit[3];
    printf("COEL %d\n",COEL);
    if (COEL == false)
    {
        
        /* Handle effect of mass loss on the rest of the system */
        printf("no coel M1 %g M2 %g kw1 %d kw2 %d SEPF %g ECC %g\n",M1,M2,KW1,KW2,SEPF,ECC);

        /* First, handle effects of possible SNe kicks (apply_kick = true for star1 and/or star2). Mass loss is not yet taken into account */
        star1->instantaneous_perturbation_delta_mass = 0.0;
        star2->instantaneous_perturbation_delta_mass = 0.0;
        
        int integration_flag_dummy;
        handle_SNe_in_system(particlesMap, &unbound_orbits, &integration_flag_dummy);

        printf("CE -- post handle_SNe_in_system\n");
        print_system(particlesMap, *integration_flag);

        /* Next, handle effects of mass loss */
        double Delta_M1 = M1 - star1->mass;
        double Delta_M2 = M2 - star2->mass;

        if (*integration_flag == 0) /* Take into account either instantaneous or adiabatic mass loss, depending on the time-scale */
        {

            //binary->mass = star1->mass + star2->mass; // the old mass -- is this necessary to set?
            //binary->apply_kick = false;
                
            handle_instantaneous_and_adiabatic_mass_changes_in_orbit(particlesMap, star1, star2, Delta_M1, Delta_M2, binary->common_envelope_timescale, &integration_flag_dummy);
        }
        else /* Always assume instantaneous mass loss in the N-body case */
        {
            star1->instantaneous_perturbation_delta_mass = Delta_M1;
            star2->instantaneous_perturbation_delta_mass = Delta_M2;;
            
            handle_SNe_in_system(particlesMap, &unbound_orbits, &integration_flag_dummy);
        }
        
        //apply_instantaneous_mass_changes_and_kicks(particlesMap, integration_flag);
        printf("CE -- post mass changes\n");
        print_system(particlesMap,*integration_flag);
        /* Update the stars and orbit */
        star1->stellar_type = KW1;
        
        //star1->mass = M1; // handled by SNe function
        star1->core_mass = MC1;
        star1->sse_initial_mass = M01;
        star1->convective_envelope_mass = MENV1;
        
        star1->age = AJ1 * Myr_to_yr;
        //star1->epoch = 0.0; /* is this correct? */
        star1->epoch = t - star1->age;
        
        star1->radius = R1 * CONST_R_SUN;
        star1->core_radius = RC1 * CONST_R_SUN;
        star1->convective_envelope_radius = RENV1 * CONST_R_SUN;
        
        star1->luminosity = L1 * CONST_L_SUN;
        
        /* Spins: either assume the spins do not change (binary_evolution_CE_spin_flag=0) */
        /* Alternatively: spins corotate and align with final orbit (binary_evolution_CE_spin_flag=1) */
        if (binary_evolution_CE_spin_flag == 1)
        {
            //get_unit_vector(star1->spin_vec,spin_vec_unit);
            for (int i=0; i<3; i++)
            {
                //star1->spin_vec[i] = OORB * spin_vec_unit[i];
                star1->spin_vec[i] = OORB * h_vec_unit[i];
            }
        }

        star2->stellar_type = KW2;
        
        //star2->mass = M2; // handled by SNe function
        star2->core_mass = MC2;
        star2->sse_initial_mass = M02;
        star2->convective_envelope_mass = MENV2;
        
        //star2->epoch = 0.0; /* is this correct? */
        star2->age = AJ2 * Myr_to_yr;
        star2->epoch = t - star2->age;

        star2->radius = R2 * CONST_R_SUN;
        star2->core_radius = RC2 * CONST_R_SUN;
        star2->convective_envelope_radius = RENV2 * CONST_R_SUN;
        
        star2->luminosity = L2 * CONST_L_SUN;

        if (binary_evolution_CE_spin_flag == 1)
        {
            //get_unit_vector(star2->spin_vec,spin_vec_unit);
            for (int i=0; i<3; i++)
            {
                //star2->spin_vec[i] = OORB * spin_vec_unit[i];
                star2->spin_vec[i] = OORB * h_vec_unit[i];
            }
        }

        /* Set new e & h vectors for the binary.
         * Assume they only change in magnitude.
         * Note: this overrides the changes made to the binary's orbit after calling handle_SNe_in_system above */

        double a = SEPF * CONST_R_SUN;
        double e = ECC;
        double h = compute_h_from_a(M1,M2,a,e);

        printf("a %g e %g h %g\n",a,e,h);
        print_system(particlesMap,*integration_flag);
        
        if (*integration_flag == 0) /* secular */
        {
            for (int i=0; i<3; i++)
            {
                binary->e_vec[i] = e * binary->e_vec_unit[i];
                binary->h_vec[i] = h * binary->h_vec_unit[i];
            }

            /* Below: not strictly neccesary */
            binary->a = a;
            binary->e = ECC;
        
            /* Assume the binary true anomaly is not affected */
            
            set_binary_masses_from_body_masses(particlesMap);
            set_positions_and_velocities(particlesMap); /* update positions and velocities of all bodies based on updated binary e and h vectors */
        }        
        else /* N-body; need to adjust positions and velocities of the two stars corresponding to the new a and e */
        {
            double R1_vec_new[3],R2_vec_new[3],V1_vec_new[3],V2_vec_new[3];
            compute_new_positions_and_velocities_given_new_semimajor_axis_and_eccentricity(M1_old,star1->R_vec,star1->V_vec,M2_old,star2->R_vec,star2->V_vec,M1,R1_vec_new,V1_vec_new,M2,R2_vec_new,V2_vec_new,a,e);
            for (int i=0; i<3; i++)
            {
                star1->R_vec[i] = R1_vec_new[i];
                star1->V_vec[i] = V1_vec_new[i];
                star2->R_vec[i] = R2_vec_new[i];
                star2->V_vec[i] = V2_vec_new[i];
            }
        }
        
        /* Override the integration flag (cf. handle_instantaneous_and_adiabatic_mass_changes_in_orbit above) */
        /* Is this necessary? Done at the end of binary_evolution() */
        bool unbound_orbits = check_for_unbound_orbits(particlesMap);
        if (unbound_orbits == true)
        {
            *integration_flag = 5;
            printf("binary_evolution.cpp -- CE -- Unbound orbits in system; switching to integration_flag %d\n",*integration_flag);
        }

        star1->RLOF_flag = 0;
        star2->RLOF_flag = 0;
        
        printf("end\n");
        print_system(particlesMap,*integration_flag);
   
    }
    else if (COEL == true)
    {
        if (*integration_flag == 0) /* secular */
        {
            #ifdef IGNORE
            /* Handle effect of mass loss on the rest of the system */

            set_old_parameters_for_adiabatic_mass_loss(particlesMap); /* copy the old masses and h & e vectors for use of adiabatic mass loss below */
            
            /* Next, assume instantaneous mass loss for ALL orbits. 
             * Orbits, for which mass loss is expected to be adiabatic,
             * are adjusted below in compute_new_orbits_assuming_adiabatic_mass_loss. */

            binary->mass = star1->mass + star2->mass; // the old mass
            double Delta_M = M1 - binary->mass;
            binary->instantaneous_perturbation_delta_mass = M1 - binary->mass;

            int integration_flag_dummy;
            handle_SNe_in_system(particlesMap, &unbound_orbits, &integration_flag_dummy); /* mass of `binary' will be updated */

            //apply_instantaneous_mass_changes_and_kicks(particlesMap, integration_flag);
            
            //apply_user_specified_instantaneous_perturbation(particlesMap);
            //reset_instantaneous_perturbation_quantities(particlesMap);
            
            //printf("mold %g mnew %g\n",binary->mass,M1);
            
            /* Lastly, `correct' the orbits (h & e vectors) for which the effect of mass loss was actually expected to be adiabatic,
             * i.e., time of CE >> t_orb */
            binary->delta_m_adiabatic_mass_loss = Delta_M;
            compute_new_orbits_assuming_adiabatic_mass_loss(particlesMap, binary->common_envelope_timescale);

            #endif

            /* Should coelesced objects get a kick? 
             * For now, assume no kick
             * Note: mergers between BHs/NSs (implying possible recoil) should never occur here (CE evolution), rather as collisions */
            star1->apply_kick = false;
            star2->apply_kick = false;
            
            //binary->mass = star1->mass + star2->mass; // the old mass -- is this necessary to set again here?
            //binary->apply_kick = false;
            
            double Delta_M1 = M1 - binary->mass;
            double Delta_M2 = 0.0;
            int integration_flag_dummy;
            handle_instantaneous_and_adiabatic_mass_changes_in_orbit(particlesMap, star1, star2, Delta_M1, Delta_M2, binary->common_envelope_timescale, &integration_flag_dummy);

            binary->is_binary = false; /* The binary becomes a body */
            //it_p = particlesMap->erase(particlesMap->find(star1->index));
            //it_p = particlesMap->erase(particlesMap->find(star2->index));

            particlesMap->erase(star1->index);
            particlesMap->erase(star2->index);

            //printf("post merge\n");
            //print_system(particlesMap);

            binary->apply_kick = false;
            binary->merged = false;
            //printf("done SNe\n");
            
            /* Update the merged star */
            binary->stellar_type = KW1;
            
            //star1->mass = M1; // handled by SNe function
            binary->core_mass = MC1;
            binary->sse_initial_mass = M01;
            binary->convective_envelope_mass = MENV1;
            
            //binary->epoch = 0.0; /* is this correct? */
            binary->age = AJ1 * Myr_to_yr;
            binary->epoch = t - binary->age;

            binary->radius = R1 * CONST_R_SUN;
            binary->core_radius = RC1 * CONST_R_SUN;
            binary->convective_envelope_radius = RENV1 * CONST_R_SUN;
            
            binary->luminosity = L1 * CONST_L_SUN;
            
            binary->metallicity = z_new;
            binary->zpars = zpars_new;

            /* Set the spin equal to the orbital frequency just before coalescence (previously calculated as OORB).
             * Assume the direction is equal to the orbital direction before coalescence. */

            for (int i=0; i<3; i++)
            {
                binary->spin_vec[i] = OORB * h_vec_unit[i];
            }
            printf("OORB %g S %g %g %g\n",OORB,binary->spin_vec[0],binary->spin_vec[1],binary->spin_vec[2]);
            binary->apply_kick = false;
            binary->RLOF_flag = 0;
            
            /* Override the integration flag (cf. handle_instantaneous_and_adiabatic_mass_changes_in_orbit above) */
            bool unbound_orbits = check_for_unbound_orbits(particlesMap);
            if (unbound_orbits == true)
            {
                *integration_flag = 3;
                printf("binary_evolution.cpp -- common_envelope_evolution -- Unbound orbits in system; switching to integration_flag %d\n",*integration_flag);
            }
        }
        else
        {
            /* Assign merger remnant to star1; remove star2 */
            //it_p = particlesMap->erase(particlesMap->find(star2->index));
            particlesMap->erase(star2->index);
            star1->stellar_type = KW1;

            for (i=0; i<3; i++)
            {
                star1->R_vec[i] = initial_R_CM[i];
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
                star1->spin_vec[i] = OORB * h_vec_unit[i];
            }

            star1->apply_kick = false; 
            star1->RLOF_flag = 0;
            //double t = 0.0;
            //double t_old = 0.0;
            //double t_out,dt_nbody;
            //integrate_nbody_system(particlesMap,integration_flag, 0.0, 0.0, &t_out, &dt_nbody);
            
        }
    }
    
    /* Reset kick parameters if the binary remains */
//    if (COEL == false)
//    {
//        star1->apply_kick = true;
//        star2->apply_kick = true;
//    }

    } // end of label30

    printf("binary_evolution.cpp -- common_envelope_evolution -- end\n");
    print_system(particlesMap,*integration_flag);

    #ifdef LOGGING
    //Log_info_type log_info;
    log_info.binary_index = binary_index;
    log_info.index1 = index1;
    log_info.index2 = index2;
    update_log_data(particlesMap, t, *integration_flag, 6, log_info);
    #endif


    return;
}


void triple_common_envelope_evolution(ParticlesMap *particlesMap, int binary_index, int index1, int index2, double t, int *integration_flag)//, ParticlesMapIterator &it_p)
{
    /* An extremely simplified treatment of triple CE evolution (an inner binary entering a tertiary star's envelope). 
     * Compute the new outer orbit assuming the binding energy of the tertiary star's envelope is converted into 
     * orbital energy of the outer orbit (for a given alpha_CE parameter of the tertiary star). 
     * Assume that the inner binary orbit does not change due to the CE (a big assumption!). 
     * Also, the tertiary star always simply loses its envelope, and the inner binary stars are not affected in their evolution. 
     * If the new outer orbital separation is too small for dynamical stability, envoke N-body integration in the future.
     * This may lead to (two-body) collisions later due to multi-body dynamical interactions. 
     * Otherwise, adjust the outer orbit, and take into account mass changes on the entire system.
     * SNe from the tertiary star (if it becomes a neutron star) might still disrupt the system.
    /* Units used here: length: RSun; time: Myr */
    
    printf("common_envelope_evolution.cpp -- triple_common_envelope_evolution() -- start; binary_index %d index1 %d index2 %d integration_flag %d\n",binary_index,index1,index2,*integration_flag);

    Particle *star3 = (*particlesMap)[index1];
    Particle *inner_binary = (*particlesMap)[index2];
    Particle *outer_binary = (*particlesMap)[binary_index];

    /* Check whether all prerequisites are met. */
    if (star3->is_binary == true or inner_binary->is_binary == false)
    {
        printf("common_envelope_evolution.cpp -- triple_common_envelope_evolution() -- ERROR: tertiary should be a star, and companion a binary\n");
        exit(-1);
    }

    Particle *c1 = (*particlesMap)[inner_binary->child1];
    Particle *c2 = (*particlesMap)[inner_binary->child2];

    if (c1->is_binary == true or c2->is_binary == true)
    {
        printf("common_envelope_evolution.cpp -- triple_common_envelope_evolution() -- at least of the components in the inner binary is itself a binary -- skipping. \n");
        return;
    }
    
    if (*integration_flag != 0)
    {
        printf("common_envelope_evolution.cpp -- triple_common_envelope_evolution() -- integration flag is nonzero; skipping. \n");
        return;
    }

    /* Define primary (star1) as the most massive star in the inner binary */
    Particle *star1, *star2;
    if (c1->mass > c2->mass)
    {
        star1 = (*particlesMap)[inner_binary->child1];
        star2 = (*particlesMap)[inner_binary->child2];
    }
    else
    {
        star1 = (*particlesMap)[inner_binary->child2];
        star2 = (*particlesMap)[inner_binary->child1];
    }

    set_up_derived_quantities(particlesMap); /* for setting a, e, etc. for all orbits */
    
    
    /* Get initial state of tertiary star */
    double tm3,tn3;
    double *GB3,*tscls3,*lums3;
    GB3 = new double[10];
    tscls3 = new double[20];
    lums3 = new double[10];    
    
    int kw3 = star3->stellar_type;
    double M3 = star3->mass;
    double M3_sse_init = star3->sse_initial_mass;
    double MC3 = star3->core_mass;
    double R3 = star3->radius/CONST_R_SUN;
    double RC3 = star3->core_radius/CONST_R_SUN;
    double L3 = star3->luminosity/CONST_L_SUN;
    double M_env3,R_env3,k2_3;
    
    double age3 = star3->age * yr_to_Myr;
    double *zpars3 = star3->zpars;    
    //printf("t1 %g %d %g %g %g %g\n",ZPARS1[0],KW1,M01,M1,TM1,TN1);
    star_(&kw3, &M3_sse_init, &M3, &tm3, &tn3, tscls3, lums3, GB3, zpars3);
    //printf("t2 %g %d %g %g %g %g\n",ZPARS1[0],KW1,M01,M1,TM1,TN1);
    //printf("arg hrd %g %g %g %g %g \n",M01,AJ1,M1,TM1,TN1);
    hrdiag_(&M3_sse_init,&age3,&M3,&tm3,&tm3,tscls3,lums3,GB3,zpars3, \
        &R3,&L3,&kw3,&MC3,&RC3,&M_env3,&R_env3,&k2_3);

//      OSPIN1 = JSPIN1/(K21*R1*R1*(M1-MC1)+K3*RC1*RC1*MC1)

    double fac = 1.0;
    double M_envd3 = M_env3 / (M3 - MC3);
    double R_ZAMS3 = rzamsf_(&M3_sse_init);
    double lambda3 = celamf_(&kw3,&M3_sse_init,&L3,&R3,&R_ZAMS3,&M_envd3,&fac);
    //double alpha3 = star3->common_envelope_alpha;
    double alpha3 = star3->triple_common_envelope_alpha;
    
    double M_inner_binary = inner_binary->mass;
    double a_in_i = inner_binary->a/CONST_R_SUN;
    double a_out_i = outer_binary->a/CONST_R_SUN;
        
    //double e_in = inner_binary->e;
    double e_out_i = outer_binary->e;


    /* Calculate the binding energy of the tertiary star's envelope (multiplied by lambda).
    * Also, calculate initial orbital energy.
    * Energies are divided by -G. */

    double E_bind3_i = M3 * (M3 - MC3)/(lambda3 * R3);
    
    double M3_CE_eff;
    if (binary_evolution_CE_energy_flag == 3)
    {
        M3_CE_eff = M3;
    }
    else
    {
        M3_CE_eff = MC3;
    }
    double E_orb_out_i = M3_CE_eff * M_inner_binary/(2.0 * a_out_i);
    
    double E_orb_out_circ_i = E_orb_out_i/(1.0 - e_out_i*e_out_i); // Allow for an eccentric orbit.

    /* Calculate the final orbital energy/orbit ignoring possible coalescence (HTP02 eq. 71). */
    double E_orb_out_f = E_orb_out_i + E_bind3_i/alpha3; 
    
    double a_out_f = MC3 * M_inner_binary/(2.0 * E_orb_out_f);

    double e_out_f;
    if (E_orb_out_f < E_orb_out_circ_i)
    {
        e_out_f = sqrt(1.0 - E_orb_out_f/E_orb_out_circ_i);
    }
    else
    {
        e_out_f = epsilon;
    }
    double rp_out_f = a_out_f*(1.0 - e_out_f);

    /* Check if the new outer orbit would be wide enough to accommodate the inner binary. */
    double rel_INCL;
    get_inclination_relative_to_parent(particlesMap,inner_binary->index,&rel_INCL); /* Assume the mutual inclination does not change during the CE */
    double rp_out_f_crit = compute_rp_out_crit_MA01(a_in_i, MC3/M_inner_binary, e_out_f, rel_INCL);
    
    bool stable = false;
    if (rp_out_f > rp_out_f_crit)
    {
        stable = true;
    }


    /* Regardless of the outcome, assume the tertiary star is always stripped of its envelope.
     * Get properties of stripped tertiary star.
     * Note: the tertiary star might get a kick. */
    double M3_f = MC3;
    int kw3_f = kw3; // could be changed below
    double M3_sse_init_f = M3_sse_init; // could be changed below
    double R3_f,L3_f,MC3_f,RC3_f,M_env3_f,R_env3_f,k2_3_f,age3_f;
    star_(&kw3_f,&M3_sse_init_f,&M3_f,&tm3,&tn3,tscls3,lums3,GB3,zpars3);
    hrdiag_(&M3_sse_init_f,&age3_f,&M3_f,&tm3,&tn3,tscls3,lums3,GB3,zpars3, \
        &R3_f,&L3_f,&kw3_f,&MC3_f,&RC3_f,&M_env3_f,&R_env3_f,&k2_3_f);

    if (kw3_f >= 13)
    {
        star3->apply_kick = true;
        //star3->kick_distribution = 1;
    }

    /* Update the tertiary star. 
     * Its mass will be updated later. */
    star3->stellar_type = kw3_f;
    
    star3->core_mass = MC3_f;
    star3->sse_initial_mass = M3_sse_init_f;
    star3->convective_envelope_mass = M_env3_f;
    
    star3->age = age3_f * Myr_to_yr;
    star3->epoch = t - star3->age;
    
    star3->radius = R3_f * CONST_R_SUN;
    star3->core_radius = RC3_f * CONST_R_SUN;
    star3->convective_envelope_radius = R_env3_f * CONST_R_SUN;
    
    star3->luminosity = L3_f * CONST_L_SUN;
    
    
    /* Adjust e & h vectors of the outer binary. 
     * Note the different units that were used in the above. */

    double a,e; /* The new outer orbit a & e */
    
    if (stable == true)
    {
        /* Final orbit is determined by available orbital energy */
        a = a_out_f * CONST_R_SUN;
        e = e_out_f;
    }
    else
    {
        /* Unstable; for the the N-body integration, park the inner binary
         * at the unstable boundary in a circular orbit. */
        a = rp_out_f_crit * CONST_R_SUN;
        e = epsilon;
    }
    
    double h = compute_h_from_a(M3_f,M_inner_binary,a,e);
    for (int i=0; i<3; i++)
    {
        outer_binary->e_vec[i] = e * outer_binary->e_vec_unit[i];
        outer_binary->h_vec[i] = h * outer_binary->h_vec_unit[i];
    }


    if (stable == true)
    {

        /* Assume the binary true anomaly is not affected */
        
        set_binary_masses_from_body_masses(particlesMap);
        set_positions_and_velocities(particlesMap); /* update positions and velocities of all bodies based on updated binary e and h vectors */

        /* Set spin of tertiary star. 
         * Assume it is either unaffected (binary_evolution_CE_spin_flag=0, or, alternatively,
         * it corotates and aligns with outer orbit (binary_evolution_CE_spin_flag=1).
         * Assume the spins of the inner binary stars are unaffected. */

        double P_orb_out_f = compute_orbital_period_from_semimajor_axis(M3_f+M_inner_binary,a);
         //double P_orb_out_f = TWOPI * (a_out_f * CONST_R_SUN) * sqrt( a_out_f * CONST_R_SUN / (CONST_G * (MC3 + M_inner_binary) ) );
        double Omega_orb_out_f = TWOPI/P_orb_out_f;

        if (binary_evolution_CE_spin_flag == 1)
        {
            //get_unit_vector(star1->spin_vec,spin_vec_unit);
            for (int i=0; i<3; i++)
            {
                star3->spin_vec[i] = Omega_orb_out_f * outer_binary->h_vec_unit[i];
            }
        }

        /* Ab initio, the new triple subsystem is dynamically stable. */
        /* However, mass loss from the CE might affect orbits exterior to the triple subsystem. 
         * Take this into account depending on the mass loss times-scale. 
         * This could change the integration flag to non-zero. */
         
        double Delta_M_in = 0.0;
        double Delta_M3 = M3_f - M3;

        outer_binary->apply_kick = false;
        handle_instantaneous_and_adiabatic_mass_changes_in_orbit(particlesMap, inner_binary, star3, Delta_M_in, Delta_M3, outer_binary->common_envelope_timescale, integration_flag); /* Will update mass of tertiary star; inner binary should be unaffected. */
    }
    else
    {
        *integration_flag = 1;
        star3->mass = M3_f;
    }

    /* Finally, apply a possible kick of the tertiary star (which could unbind the system if it were stable before). */
    star1->apply_kick = false;
    star2->apply_kick = false;
    star3->instantaneous_perturbation_delta_mass = 0.0;
//    int integration_flag_dummy;
    bool unbound_orbits;
    handle_SNe_in_system(particlesMap, &unbound_orbits, integration_flag);

    return;
}

int common_envelope_evolution_old(ParticlesMap *particlesMap, int binary_index, int donor_index, int accretor_index)
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

}
