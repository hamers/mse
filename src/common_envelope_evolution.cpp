/* MSE */

#include "evolve.h"
#include "common_envelope_evolution.h"

extern "C"
{

int common_envelope_evolution(ParticlesMap *particlesMap, int binary_index, int index1, int index2, double t, int *integration_flag)//, ParticlesMapIterator &it_p)
{

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
    /*************
    * Common Envelope Evolution.
    * Mostly adopted from BSE.
    * Note units in SSE/BSE: length in RSUN, time in Myr, luminosity in LSun */

    #ifdef LOGGING
    Log_info_type log_info;
    log_info.binary_index = binary_index;
    log_info.index1 = index1;
    log_info.index2 = index2;
    update_log_data(particlesMap, t, *integration_flag, LOG_CE_START, log_info);
    #endif

    #ifdef VERBOSE
    if (verbose_flag > 0)
    {
        printf("common_envelope_evolution.cpp -- binary_common_envelope_evolution() -- start; binary_index %d index1 %d index2 %d integration_flag %d\n",binary_index,index1,index2,*integration_flag);
    }
    #endif

    int i;

    Particle *star1 = (*particlesMap)[index1];
    Particle *star2 = (*particlesMap)[index2];

    if (star1->object_type == 1 and star2->object_type == 2)
    {
        collision_product_star_planet(particlesMap, binary_index, index1, index2,t,integration_flag);
        return;
    }
    if (star1->object_type == 2 and star2->object_type == 1)
    {
        collision_product_star_planet(particlesMap, binary_index, index2, index1,t,integration_flag);
        return;
    }
    if (star1->object_type == 2 and star2->object_type == 2)
    {
        collision_product_planet_planet(particlesMap, binary_index, index1, index2,t,integration_flag);
        return;
    }
    
    double M1_old,M1;
    M1_old = M1 = star1->mass;
    double M2_old,M2;
    M2_old = M2 = star2->mass;

    /* Orbital properties */
    double SEP; // unit: RSun
    double SEPF; // final separation; unit: RSun
    double SEPL; // unit: RSun
    double ECC;

    double initial_momentum[3],initial_R_CM[3],initial_V_CM[3];
    double h_vec[3],h_vec_unit[3],e_vec[3],e_vec_unit[3],r_vec[3],v_vec[3];

    get_initial_binary_orbital_properties_from_position_and_velocity(star1->R_vec, star1->V_vec, star2->R_vec, star2->V_vec, M1, M2, \
        r_vec, v_vec, initial_momentum, initial_R_CM, initial_V_CM, h_vec, e_vec);
    
    ECC = norm3(e_vec);
    if (ECC <= epsilon)
    {
        ECC = epsilon;
    }
    
    SEP = compute_a_from_h(M1,M2,norm3(h_vec),ECC)/CONST_R_SUN;

    if (SEP <= epsilon) /* Take instantenous separation if semimajor axis not well defined */
    {
        double r[3];
        for (int i=0; i<3; i++)
        {
            r[i] = star1->R_vec[i] - star2->R_vec[i];
        }
        SEP = norm3(r);
    }

    get_unit_vector(h_vec,h_vec_unit);
    get_unit_vector(e_vec,e_vec_unit);

    #ifdef VERBOSE
    if (verbose_flag > 1)
    {
        printf("common_envelope_evolution.cpp -- binary_common_envelope_evolution() -- init a/au %g e %g\n",SEP*CONST_R_SUN,ECC);
    }
    #endif
    
    int CEFLAG = binary_evolution_CE_energy_flag;
    double ALPHA1 = star1->common_envelope_alpha;
    
    double fac = binary_evolution_CE_recombination_fraction;
    bool COEL = false;
    
    /* Get current state of star 1 (primary) */
    double TM1,TN1;
    double *GB1,*TSCLS1,*LUMS1;
    GB1 = new double[10];
    TSCLS1 = new double[20];
    LUMS1 = new double[10];    
    
    int KW1 = star1->stellar_type;
    int KW1_old = KW1;
    double M01 = star1->sse_initial_mass;
    double MC1 = star1->core_mass;
    double MC1_old = MC1;
    double R1 = star1->radius/CONST_R_SUN;
    double RC1 = star1->core_radius/CONST_R_SUN;
    double L1 = star1->luminosity/CONST_L_SUN;
    double MENV1,RENV1,K21;
    
    double AJ1 = star1->age * yr_to_Myr;
    double *ZPARS1 = star1->zpars;   
    //printf("pre k %d AJ1 %g m01 %g m1 %g tm1 %g tn1 %g r1 %g l1 %g mc1 %g rc1 %g menv1 %g renv %g\n",KW1,AJ1,M01,M1,TM1,TN1,R1,L1,MC1,RC1,MENV1,RENV1);
    star_(&KW1, &M01, &M1, &TM1, &TN1, TSCLS1, LUMS1, GB1, ZPARS1);
    //printf("TM1 %g TN1 %g TSCLS1 %g\n",TM1,TN1,TSCLS1[0]);
    hrdiag_(&M01,&AJ1,&M1,&TM1,&TN1,TSCLS1,LUMS1,GB1,ZPARS1, \
        &R1,&L1,&KW1,&MC1,&RC1,&MENV1,&RENV1,&K21);
    //printf("post k %d AJ1 %g m01 %g m1 %g tm1 %g tn1 %g r1 %g l1 %g mc1 %g rc1 %g menv1 %g renv %g\n",KW1,AJ1,M01,M1,TM1,TN1,R1,L1,MC1,RC1,MENV1,RENV1);
    if (MC1 <= epsilon)
    {
        printf("common_envelope_evolution.cpp -- binary_common_envelope_evolution -- no core present (MC1 = %g)! M1 = %g; skipping\n",MC1,M1);
        //error_code = 39;
        //longjmp(jump_buf,1);
        
        delete[] GB1;
        delete[] TSCLS1;
        delete[] LUMS1;

        return;
            
        //MC1 = M1;
    }

    double MENVD1 = MENV1 / (M1 - MC1);
    double RZAMS1 = rzamsf_(&M01);
    double LAMB1 = celamf_(&KW1,&M01,&L1,&R1,&RZAMS1,&MENVD1,&fac);
    double spin_vec_1_norm = norm3(star1->spin_vec);
    
    /* Get current state of star 2 (secondary) */
    double TM2,TN2;
    double *GB2,*TSCLS2,*LUMS2;
    GB2 = new double[10];
    TSCLS2 = new double[20];
    LUMS2 = new double[10];    
    
    int KW2 = star2->stellar_type;
    int KW2_old = KW2;
    double M02 = star2->sse_initial_mass;
    double MC2 = star2->core_mass;
    double MC2_old = MC2;
    double R2 = star2->radius/CONST_R_SUN;
    double RC2 = star2->core_radius/CONST_R_SUN;
    double L2 = star2->luminosity/CONST_L_SUN;
    double MENV2,RENV2,K22;
    
    double AJ2 = star2->age * yr_to_Myr;
    double *ZPARS2 = star2->zpars;    

    star_(&KW2, &M02, &M2, &TM2, &TN2, TSCLS2, LUMS2, GB2, ZPARS2);
    hrdiag_(&M02,&AJ2,&M2,&TM2,&TN2,TSCLS2,LUMS2,GB2,ZPARS2, \
        &R2,&L2,&KW2,&MC2,&RC2,&MENV2,&RENV2,&K22);

    double MENVD2 = MENV2 / (M2 - MC2);
    double RZAMS2 = rzamsf_(&M02);
    double LAMB2 = celamf_(&KW2,&M02,&L2,&R2,&RZAMS2,&MENVD2,&fac);
    double spin_vec_2_norm = norm3(star2->spin_vec);

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

    bool skip_final_age_and_mass_determination = false;

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

            if (KW1 != KW1_old)
            {
                if (KW1 >= 13 or (KW1 >= 10 and KW1 <= 12 and star1->include_WD_kicks == true))
                {
                    star1->apply_kick = true;
                }
                
                if (NS_model == 1)
                {
                    if (KW1 == 13)
                    {
                        double ospin;
                        compute_NS_formation_properties_Ye19_model(false, &ospin, &star1->magnetic_field_strength_gauss);
                        rescale_vector(star1->spin_vec, ospin/norm3(star1->spin_vec));
                        
                        star1->time_of_NS_formation = t;
                        star1->initial_NS_period_s = compute_spin_period_from_spin_angular_frequency(ospin) * yr_to_s;
                        star1->initial_magnetic_field_strength_gauss = star1->magnetic_field_strength_gauss;
                    }
                }
            }
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

            KW = determine_merger_type(KW1,KW2);
            MC3 = MC1 + MC2;

            /*
            * Calculate the final envelope binding energy.
            */
            
            EORBF = CV_max(MC1 * MC2/(2.0 * SEPL), EORBI);
            EBINDF = EBINDI - ALPHA1*(EORBF - EORBI);

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

                if (KW1 != 15)
                {
                    R1 = determine_sse_compact_object_radius_RSun(KW1,M1);
                }
                if (KW1 == 15) /* SNe occurred */
                {
                    log_info.index1 = star1->index;
                    log_info.index2 = star2->index;
                    update_log_data(particlesMap, t, *integration_flag, LOG_SNE_START, log_info);
                }
                skip_final_age_and_mass_determination = true;
               //goto label30;
            }
                
            /*
            * Check if we have the merging of two degenerate cores and if so
            * then see if the resulting core will survive or change form.
            */
            
            if (skip_final_age_and_mass_determination == false)
            {
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
                            
                            if (KW1 != 15)
                            {
                                R1 = determine_sse_compact_object_radius_RSun(KW1,M1);
                            }
                            if (KW1 == 15) /* SNe occurred */
                            {
                                log_info.index1 = star1->index;
                                log_info.index2 = star2->index;
                                update_log_data(particlesMap, t, *integration_flag, LOG_SNE_START, log_info);
                            }

                            skip_final_age_and_mass_determination = true;
                            //goto label30;
                        }
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
            if (KW1 != KW1_old)
            {
                if (KW1 >= 13 or (KW1 >= 10 and KW1 <= 12 and star1->include_WD_kicks == true))
                {
                    star1->apply_kick = true;
                }
                
                if (NS_model == 1)
                {
                    if (KW1 == 13)
                    {
                        double ospin;
                        compute_NS_formation_properties_Ye19_model(false, &ospin, &star1->magnetic_field_strength_gauss);
                        rescale_vector(star1->spin_vec, ospin/norm3(star1->spin_vec));
                        
                        star1->time_of_NS_formation = t;
                        star1->initial_NS_period_s = compute_spin_period_from_spin_angular_frequency(ospin) * yr_to_s;
                        star1->initial_magnetic_field_strength_gauss = star1->magnetic_field_strength_gauss;
                    }
                }
            }
            
            MF = M2;
            KW = KW2;
            M2 = MC2;
            star_(&KW2,&M02,&M2,&TM2,&TN,TSCLS2,LUMS,GB,ZPARS2);
            hrdiag_(&M02,&AJ2,&M2,&TM2,&TN,TSCLS2,LUMS,GB,ZPARS2, \
                &R2,&L2,&KW2,&MC2,&RC2,&MENV2,&RENV2,&K22);
            if (KW2 != KW2_old)
            {
                if((KW2 >= 13 and KW < 13) or (KW2 >= 10 and KW <= 12 and star2->include_WD_kicks == true)) /* secondary became an NS/WD */
                {
                    star2->apply_kick = true;
                }

                if (NS_model == 1)
                {
                    if (KW2 == 13)
                    {
                        double ospin;
                        compute_NS_formation_properties_Ye19_model(false, &ospin, &star2->magnetic_field_strength_gauss);
                        rescale_vector(star2->spin_vec, ospin/norm3(star2->spin_vec));
                        
                        star2->time_of_NS_formation = t;
                        star2->initial_NS_period_s = compute_spin_period_from_spin_angular_frequency(ospin) * yr_to_s;
                        star2->initial_magnetic_field_strength_gauss = star2->magnetic_field_strength_gauss;
                    }
                }

            }
        }
    } // back at main level
    
    double MC22; // Effective core mass of the secondary
    double FAGE1,FAGE2; // Age factors


    /* If making a single helium burning star, calculate the fractional age 
     * depending on the amount of helium that has burnt. */
    if (COEL == true and skip_final_age_and_mass_determination == false)
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
    }

    /* Now calculate the final mass following coelescence.  This requires a Newton-Raphson iteration. 
    * Note ASH: for now, take metallicity of star2 for the coelesced object */

    double z_new;
    double *zpars_new;
    zpars_new = new double[20];

    if (COEL == true)
    {
        z_new = star2->metallicity;
        zcnsts_(&z_new,zpars_new);

        /* Calculate the orbital spin just before coalescence. */

        TB = TWOPI * (SEPL * CONST_R_SUN) * sqrt( SEPL * CONST_R_SUN / (CONST_G * (MC1 + MC2) ) );
        OORB = TWOPI/TB;

        if (skip_final_age_and_mass_determination == false)
        {
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
            else
            {
                star_(&KW,&M01,&M1,&TM1,&TN,TSCLS1,LUMS,GB,ZPARS2);
            }
            
            //printf("HRDIAG0 M01 %g M1 %g AJ1 %g KW %d TN %g TM1 %g\n",M01,M1,AJ1,KW1,TN,TM1);
            hrdiag_(&M01,&AJ1,&M1,&TM1,&TN,TSCLS1,LUMS,GB,ZPARS2, \
                &R1,&L1,&KW,&MC1,&RC1,&MENV1,&RENV1,&K21);
            //printf("HRDIAG1 M01 %g M1 %g AJ1 %g KW %d R1 %g L1 %g\n",M01,M1,AJ1,KW1,R1,L1);
            KW1 = KW;
            ECC = 0.0;
        }
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
         }

         /* Set both cores in co-rotation with the orbit on exit of CE */

         TB = TWOPI * (SEPF * CONST_R_SUN) * sqrt( SEPF * CONST_R_SUN / (CONST_G * (M1 + M2) ) );
         OORB = TWOPI/TB;
         
         check_number(OORB,"common_envelope_evolution.cpp -- binary_common_envelope_evolution()","OORB", true);
    }

    /* Handle orbital changes / merge binaries in case of coalescence */
    
    //label30:
    {

    *integration_flag = 1;

    #ifdef VERBOSE
    if (verbose_flag > 1)
    {
        printf("common_envelope_evolution.cpp -- binary_common_envelope_evolution() -- pre handle mass loss -- COEL %d OORB %g\n",COEL,OORB);
        print_system(particlesMap,*integration_flag);
    }
    #endif

    
    bool unbound_orbits = false;

    double final_R_CM[3],final_V_CM[3],final_momentum[3];
    handle_gradual_mass_loss_event_in_system(particlesMap, star1, star2, M1, M1_old, M2, M2_old, star1->common_envelope_timescale, \
        r_vec, v_vec, initial_R_CM, initial_V_CM, final_R_CM, final_V_CM, final_momentum);
    double Omega_crit,Omega;

    #ifdef VERBOSE
    if (verbose_flag > 1)
    {
        printf("common_envelope_evolution.cpp -- binary_common_envelope_evolution() -- post handle mass loss -- COEL %d OORB %g\n",COEL,OORB);
        print_system(particlesMap,*integration_flag);
    }
    #endif


    if (COEL == false)
    {
        
        /* Take into account kicks */
        double V_kick_vec[3];
        if (star1->apply_kick == true)
        {
            #ifdef LOGGING
            log_info.index1 = star1->index;
            if (KW1 >= 10 and KW1 <= 12 and star1->include_WD_kicks == true)
            {
                update_log_data(particlesMap, t, *integration_flag, LOG_WD_KICK_START, log_info);
            }
            else
            {
                update_log_data(particlesMap, t, *integration_flag, LOG_SNE_START, log_info);
            }
            #endif

            /* Temporarily set some properties for the purposes of kick velocity sampling */
            star1->stellar_type = KW1;
            star1->mass = M1_old;
            star1->instantaneous_perturbation_delta_mass = M1 - M1_old;
            star1->core_mass_old = MC1_old;
            sample_kick_velocity(star1, &V_kick_vec[0], &V_kick_vec[1], &V_kick_vec[2]);
            for (int i=0; i<3; i++)
            {
                star1->V_vec[i] += V_kick_vec[i];
            }
        }
        if (star2->apply_kick == true)
        {
            #ifdef LOGGING
            log_info.index1 = star2->index;
            if (KW2 >= 10 and KW2 <= 12 and star2->include_WD_kicks == true)
            {
                update_log_data(particlesMap, t, *integration_flag, LOG_WD_KICK_START, log_info);
            }
            else
            {
                update_log_data(particlesMap, t, *integration_flag, LOG_SNE_START, log_info);
            }
            #endif

            /* Temporarily set some properties for the purposes of kick velocity sampling */
            star2->stellar_type = KW2;
            star2->mass = M2_old;
            star2->instantaneous_perturbation_delta_mass = M2 - M2_old;
            star2->core_mass_old = MC2_old;
            sample_kick_velocity(star2, &V_kick_vec[0], &V_kick_vec[1], &V_kick_vec[2]);
            for (int i=0; i<3; i++)
            {
                star2->V_vec[i] += V_kick_vec[i];
            }
        }

        /* Update all properties of the two stars */
        star1->stellar_type = KW1;
        star1->mass = M1;        
        star1->core_mass = MC1;
        star1->sse_initial_mass = M01;
        star1->convective_envelope_mass = MENV1;
        
        star1->age = AJ1 * Myr_to_yr;
        star1->epoch = t - star1->age;
        
        star1->radius = R1 * CONST_R_SUN;
        star1->core_radius = RC1 * CONST_R_SUN;
        star1->convective_envelope_radius = RENV1 * CONST_R_SUN;
        star1->luminosity = L1 * CONST_L_SUN;
        
        star2->stellar_type = KW2;
        star2->mass = M2;
        star2->core_mass = MC2;
        star2->sse_initial_mass = M02;
        star2->convective_envelope_mass = MENV2;
        
        star2->age = AJ2 * Myr_to_yr;
        star2->epoch = t - star2->age;

        star2->radius = R2 * CONST_R_SUN;
        star2->core_radius = RC2 * CONST_R_SUN;
        star2->convective_envelope_radius = RENV2 * CONST_R_SUN;
        star2->luminosity = L2 * CONST_L_SUN;

        /* Spins: either assume the spins do not change (binary_evolution_CE_spin_flag=0) */
        /* Alternatively: spins corotate and align with final orbit (binary_evolution_CE_spin_flag=1) */
        /* Limit spins to critical rotation */
        if (binary_evolution_CE_spin_flag == 1)
        {
            Omega_crit = compute_breakup_angular_frequency(star1->mass,star1->radius);
            Omega = OORB;
            if (OORB >= Omega_crit)
            {
                #ifdef VERBOSE
                if (verbose_flag > 0)
                {
                    printf("common_envelope_evolution.cpp -- binary_common_envelope_evolution() -- Limiting spin frequency of star %d from %g to breakup rate %g\n",star1->index,OORB,Omega_crit);
                    print_system(particlesMap,*integration_flag);
                }
                #endif
                
                Omega = Omega_crit;
            }
            
            for (int i=0; i<3; i++)
            {
                star1->spin_vec[i] = Omega * h_vec_unit[i];
            }
        }

        if (binary_evolution_CE_spin_flag == 1)
        {
            Omega_crit = compute_breakup_angular_frequency(star2->mass,star2->radius);
            Omega = OORB;
            if (OORB >= Omega_crit)
            {
                #ifdef VERBOSE
                if (verbose_flag > 0)
                {
                    printf("common_envelope_evolution.cpp -- binary_common_envelope_evolution() -- Limiting spin frequency of star %d from %g to breakup rate %g\n",star2->index,OORB,Omega_crit);
                    print_system(particlesMap,*integration_flag);
                }
                #endif

                Omega = Omega_crit;
            }

            for (int i=0; i<3; i++)
            {
                star2->spin_vec[i] = Omega * h_vec_unit[i];
            }
        }

        /* Adjust positions and velocities of the two stars to be consistent with the post-CE orbit */
        /* Set new e & h vectors for the binary.
         * Assume they only change in magnitude.
         * Note: this overrides the changes made to the binary's orbit after calling handle_SNe_in_system above */

        double a = SEPF * CONST_R_SUN;
        double e = ECC;
        double h = compute_h_from_a(M1,M2,a,e);

        #ifdef VERBOSE
        if (verbose_flag > 1)
        {
            printf("common_envelope_evolution.cpp -- binary_common_envelope_evolution() -- no COEL post a %g e %g h %g\n",a,e,h);
            print_system(particlesMap,*integration_flag);
        }
        #endif

        double R1_vec_new[3],R2_vec_new[3],V1_vec_new[3],V2_vec_new[3];

        compute_new_positions_and_velocities_given_new_semimajor_axis_and_eccentricity(M1_old,star1->R_vec,star1->V_vec,M2_old,star2->R_vec,star2->V_vec,M1,R1_vec_new,V1_vec_new,M2,R2_vec_new,V2_vec_new,a,e);
        for (int i=0; i<3; i++)
        {
            star1->R_vec[i] = R1_vec_new[i];
            star1->V_vec[i] = V1_vec_new[i];
            star2->R_vec[i] = R2_vec_new[i];
            star2->V_vec[i] = V2_vec_new[i];
        }

        /* Reset some parameters */
        reset_ODE_mass_dot_quantities(star1);
        reset_ODE_mass_dot_quantities(star2);
    }
    else if (COEL == true)
    {
        /* Assign merger remnant to star1; remove star2 */
        particlesMap->erase(index2);

        /* Take into account kicks */
        double V_kick_vec[3];
        if (star1->apply_kick == true)
        {
            #ifdef LOGGING
            log_info.index1 = star1->index;
            if (KW1 >= 10 and KW1 <= 12 and star1->include_WD_kicks == true)
            {
                update_log_data(particlesMap, t, *integration_flag, LOG_WD_KICK_START, log_info);
            }
            else
            {
                update_log_data(particlesMap, t, *integration_flag, LOG_SNE_START, log_info);
            }
            #endif
            
            /* Temporarily set some properties for the purposes of kick velocity sampling */
            star1->stellar_type = KW1;
            star1->mass = M1_old;
            star1->instantaneous_perturbation_delta_mass = M1 - M1_old;
            star1->core_mass_old = MC1_old;
            
            sample_kick_velocity(star1, &V_kick_vec[0], &V_kick_vec[1], &V_kick_vec[2]);
            for (int i=0; i<3; i++)
            {
                star1->V_vec[i] += V_kick_vec[i];
            }
        }

        /* Merged star properties */
        star1->stellar_type = KW1;
        star1->mass = M1;
        star1->core_mass = MC1;
        star1->sse_initial_mass = M01;
        star1->convective_envelope_mass = MENV1;
        
        star1->epoch = t - star1->age;
        star1->age = AJ1 * Myr_to_yr;

        star1->radius = R1 * CONST_R_SUN;
        star1->core_radius = RC1 * CONST_R_SUN;
        star1->convective_envelope_radius = RENV1 * CONST_R_SUN;
        star1->luminosity = L1 * CONST_L_SUN;
        
        /* Position/velocity 
         * New position is CM position of original two stars
         * Set new velocity according to linear momentum conservation */
        for (i=0; i<3; i++)
        {
            star1->R_vec[i] = final_R_CM[i];
            star1->V_vec[i] = final_momentum[i]/M1; 
        }

            
        /* Set the spin equal to the orbital frequency just before coalescence (previously calculated as OORB).
         * Assume the direction is equal to the previous orbital orientation. 
         * Limit to breakup rotation */

        if (binary_evolution_CE_spin_flag == 1)
        {
            Omega_crit = compute_breakup_angular_frequency(star1->mass,star1->radius);
            Omega = OORB;
            if (OORB >= Omega_crit)
            {
                #ifdef VERBOSE
                if (verbose_flag > 0)
                {
                    printf("common_envelope_evolution.cpp -- binary_common_envelope_evolution() -- Limiting spin frequency of star %d from %g to breakup rate %g\n",star1->index,OORB,Omega_crit);
                    print_system(particlesMap,*integration_flag);
                }
                #endif
                
                Omega = Omega_crit;
            }

            for (int i=0; i<3; i++)
            {
                star1->spin_vec[i] = Omega * h_vec_unit[i];
            }
        }
        //printf("COEL %g Omega %g %d OORB %g Omega_crit %g\n",norm3(star1->spin_vec),Omega,binary_evolution_CE_spin_flag,OORB,Omega_crit);
        
        /* Reset some parameters */
        reset_ODE_mass_dot_quantities(star1);
    }

    /* Special treatment for NSs */
    if (NS_model == 1)
    {
        if (COEL == true)
        {
            if (KW1 == 13)
            {
                if (KW2_old == 13 and KW1_old < 13) /* star2 becomes star1; copy B field */
                {
                    star1->magnetic_field_strength_gauss = star2->magnetic_field_strength_gauss;
                }

                if ((KW1_old == 13 and KW2_old < 13) or (KW2_old == 13 and KW1_old < 13))
                {
                    if ((KW1_old == 13 and (determine_if_NS_is_MSP(spin_vec_1_norm,star1->magnetic_field_strength_gauss) == true)) or (KW2_old == 13 and (determine_if_NS_is_MSP(spin_vec_2_norm,star2->magnetic_field_strength_gauss) == true)))
                    {
                        /* This block is entered if either of the original objects was a MSP, and the new object is also a NS */
                        double ospin;
                        compute_NS_formation_properties_Ye19_model(true, &ospin, &star1->magnetic_field_strength_gauss);
                        rescale_vector(star1->spin_vec,ospin/norm3(star1->spin_vec));

                        star1->time_of_NS_formation = t;
                        star1->initial_NS_period_s = compute_spin_period_from_spin_angular_frequency(ospin) * yr_to_s;
                        star1->initial_magnetic_field_strength_gauss = star1->magnetic_field_strength_gauss;
                    }
                }
            }
        }
    }

    } // end of label30

    //*integration_flag = determine_orbits_in_system_using_nbody(particlesMap);

    remove_massless_remnants_from_system(particlesMap, integration_flag);
    reset_RLOF_flags(particlesMap);

    #ifdef LOGGING
    //Log_info_type log_info;
    log_info.binary_index = binary_index;
    log_info.index1 = index1;
    log_info.index2 = index2;
    update_log_data(particlesMap, t, *integration_flag, LOG_CE_END, log_info);
    #endif

    #ifdef VERBOSE
    if (verbose_flag > 0)
    {
        printf("common_envelope_evolution.cpp -- binary_common_envelope_evolution() -- end\n");
        print_system(particlesMap,*integration_flag);
    }
    #endif

    delete[] GB1;
    delete[] TSCLS1;
    delete[] LUMS1;

    delete[] GB2;
    delete[] TSCLS2;
    delete[] LUMS2;
    
    delete[] GB;
    delete[] TSCLS;
    delete[] LUMS;

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

    #ifdef VERBOSE
    if (verbose_flag > 0)
    {
        printf("common_envelope_evolution.cpp -- triple_common_envelope_evolution() -- start; binary_index %d index1 %d index2 %d integration_flag %d\n",binary_index,index1,index2,*integration_flag);
        print_system(particlesMap,*integration_flag);
    }
    #endif
    
    #ifdef LOGGING
    Log_info_type log_info;
    log_info.binary_index = binary_index;
    log_info.index1 = index1;
    log_info.index2 = index2;
    update_log_data(particlesMap, t, *integration_flag, LOG_TRIPLE_CE_START, log_info);
    #endif

    Particle *star3 = (*particlesMap)[index1];
    Particle *inner_binary = (*particlesMap)[index2];
    Particle *outer_binary = (*particlesMap)[binary_index];

    /* Check whether all prerequisites are met. */
    if (star3->is_binary == true or inner_binary->is_binary == false)
    {
        printf("common_envelope_evolution.cpp -- triple_common_envelope_evolution() -- ERROR: tertiary should be a star, and companion a binary\n");
        error_code = 8;
        longjmp(jump_buf,1);
        //exit(-1);
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

    star_(&kw3, &M3_sse_init, &M3, &tm3, &tn3, tscls3, lums3, GB3, zpars3);
    hrdiag_(&M3_sse_init,&age3,&M3,&tm3,&tm3,tscls3,lums3,GB3,zpars3, \
        &R3,&L3,&kw3,&MC3,&RC3,&M_env3,&R_env3,&k2_3);

    double fac = binary_evolution_CE_recombination_fraction;
    double M_envd3 = M_env3 / (M3 - MC3);
    double R_ZAMS3 = rzamsf_(&M3_sse_init);
    double lambda3 = celamf_(&kw3,&M3_sse_init,&L3,&R3,&R_ZAMS3,&M_envd3,&fac);
    double alpha3 = star3->triple_common_envelope_alpha;
    
    double M_inner_binary = inner_binary->mass;
    double a_in_i = inner_binary->a/CONST_R_SUN;
    double a_out_i = outer_binary->a/CONST_R_SUN;
        
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
    int kw3_f = kw3; // could be changed by hrdiag() below
    double M3_sse_init_f = M3_sse_init; // could be changed by hrdiag() below
    double age3_f = age3; // could be changed by hrdiag() below
    double R3_f,L3_f,MC3_f,RC3_f,M_env3_f,R_env3_f,k2_3_f;
    
    star_(&kw3_f,&M3_sse_init_f,&M3_f,&tm3,&tn3,tscls3,lums3,GB3,zpars3);
    hrdiag_(&M3_sse_init_f,&age3_f,&M3_f,&tm3,&tn3,tscls3,lums3,GB3,zpars3, \
        &R3_f,&L3_f,&kw3_f,&MC3_f,&RC3_f,&M_env3_f,&R_env3_f,&k2_3_f);

    if (kw3_f >= 13 or (kw3_f >= 10 and kw3_f <= 12 and star3->include_WD_kicks == true))
    {
        star3->apply_kick = true;
        
        #ifdef LOGGING
        log_info.index1 = star1->index;
        update_log_data(particlesMap, t, *integration_flag, LOG_WD_KICK_START, log_info);
        #endif
        
        if (NS_model == 1)
        {
            if (kw3_f != kw3 and kw3_f == 13)
            {
                double ospin;
                compute_NS_formation_properties_Ye19_model(false, &ospin, &star3->magnetic_field_strength_gauss);
                rescale_vector(star3->spin_vec, ospin/norm3(star3->spin_vec));
                
                star3->time_of_NS_formation = t;
                star3->initial_NS_period_s = compute_spin_period_from_spin_angular_frequency(ospin) * yr_to_s;
                star3->initial_magnetic_field_strength_gauss = star3->magnetic_field_strength_gauss;
            }
        }

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


    set_binary_masses_from_body_masses(particlesMap);
    set_positions_and_velocities(particlesMap); /* update positions and velocities of all bodies based on updated binary e and h vectors */


    if (stable == true)
    {

        /* Assume the binary true anomaly is not affected */
        
//        set_binary_masses_from_body_masses(particlesMap);
//        set_positions_and_velocities(particlesMap); /* update positions and velocities of all bodies based on updated binary e and h vectors */

        /* Set spin of tertiary star. 
         * Assume it is either unaffected (binary_evolution_CE_spin_flag=0, or, alternatively,
         * it corotates and aligns with outer orbit (binary_evolution_CE_spin_flag=1).
         * Assume the spins of the inner binary stars are unaffected. */

        double P_orb_out_f = compute_orbital_period_from_semimajor_axis(M3_f+M_inner_binary,a);
         //double P_orb_out_f = TWOPI * (a_out_f * CONST_R_SUN) * sqrt( a_out_f * CONST_R_SUN / (CONST_G * (MC3 + M_inner_binary) ) );
        double Omega_orb_out_f = TWOPI/P_orb_out_f;

        if (binary_evolution_CE_spin_flag == 1)
        {
            double Omega_crit = compute_breakup_angular_frequency(star1->mass,star1->radius);
            double Omega = Omega_orb_out_f;
            if (Omega_orb_out_f >= Omega_crit)
            {
                #ifdef VERBOSE
                if (verbose_flag > 0)
                {
                    printf("common_envelope_evolution.cpp -- triple_common_envelope_evolution() -- Limiting spin frequency of star %d from %g to breakup rate %g\n",star3->index,Omega_orb_out_f,Omega_crit);
                    print_system(particlesMap,*integration_flag);
                }
                #endif
                
                Omega = Omega_crit;
            }
            
            for (int i=0; i<3; i++)
            {
                star3->spin_vec[i] = Omega * outer_binary->h_vec_unit[i];
            }
        }

        double final_R_CM[3], final_V_CM[3], final_momentum[3];
        handle_gradual_mass_loss_event_in_system_triple_CE(particlesMap, star3, c1, c2, M3_f, M3, star3->common_envelope_timescale, \
            star3->R_vec,star3->V_vec,final_R_CM, final_V_CM, final_momentum);

        /* Ab initio, the new triple subsystem is dynamically stable. */
        /* However, mass loss from the CE might affect orbits exterior to the triple subsystem. 
         * Take this into account depending on the mass loss times-scale. 
         * This could change the integration flag to non-zero. */

        #ifdef IGNORE
        double Delta_M_in = 0.0;
        double Delta_M3 = M3_f - M3;

        outer_binary->apply_kick = false;
        handle_instantaneous_and_adiabatic_mass_changes_in_orbit(particlesMap, inner_binary, star3, Delta_M_in, Delta_M3, outer_binary->common_envelope_timescale, integration_flag); /* Will update mass of tertiary star; inner binary should be unaffected. */
        #endif
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

    bool unbound_orbits;
    handle_SNe_in_system(particlesMap, &unbound_orbits, integration_flag);

    /* Reset some simulation parameters for the tertiary star */
    reset_ODE_mass_dot_quantities(star3);

    delete[] GB3;
    delete[] tscls3;
    delete[] lums3;
    
    #ifdef LOGGING
    log_info.binary_index = binary_index;
    log_info.index1 = index1;
    log_info.index2 = index2;
    update_log_data(particlesMap, t, *integration_flag, LOG_TRIPLE_CE_END, log_info);
    #endif
    
    

    #ifdef VERBOSE
    if (verbose_flag > 0)
    {
        printf("common_envelope_evolution.cpp -- triple_common_envelope_evolution() -- end; binary_index %d index1 %d index2 %d integration_flag %d\n",binary_index,index1,index2,*integration_flag);
        print_system(particlesMap,*integration_flag);
    }
    #endif
    
    return;
}

}
