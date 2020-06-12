/* MSE */

#include "evolve.h"
#include "stellar_evolution.h"

extern "C"
{

struct value1__ value1_;
struct value2__ value2_;
struct value3__ value3_;
struct value4__ value4_;
struct value5__ value5_;
struct flags__ flags_;
struct points__ points_;

int initialize_stars(ParticlesMap *particlesMap)
{

    value1_.neta = 0.5;
    value1_.bwind = 0.0;
    value1_.hewind = 0.5;
    value1_.mxns = 3.0;
    value2_.alpha1 = 1.0;
    value2_.lambda = 1.0;
    value3_.idum = 0;
    value4_.sigma = 190.0;
    value4_.bhflag = 1;
    //        value5_.beta = 0.0;
    //        value5_.xi = 0.0;
    //        value5_.acc2 = 0.0;
    //        value5_.epsnov = 0.0;
    //        value5_.eddfac = 0.0;
    //        value5_.gamma = 0.0;
    flags_.ceflag = binary_evolution_CE_energy_flag;
    flags_.tflag = 0;
    flags_.ifflag = 0;
    flags_.nsflag = 1;
    flags_.wdflag = 1;
    points_.pts1 = 0.05;
    points_.pts2 = 0.01;
    points_.pts3 = 0.02;
    
    double sse_initial_mass,mass,mt;
    double z;

    set_up_derived_quantities(particlesMap); /* compute a & e for all orbits; used at the end below to adjust h vectors such that possible initial stellar evolution would not change the orbital elements  */

    int i;
    int kw,kw_desired;
    double tm,tn;
    double age;
    double *GB,*tscls,*lums;
    double r,lum,mc,rc,menv,renv,k2;
            
           
    ParticlesMapIterator it_p;
    for (it_p = particlesMap->begin(); it_p != particlesMap->end(); it_p++)
    {
        Particle *p = (*it_p).second;
        if (p->is_binary == false and p->evolve_as_star == true)
        {
            z = p->metallicity;

            double *zpars;
            zpars = new double[20];
            zcnsts_(&z,zpars);
            p->zpars = zpars;
            
            kw_desired = p->stellar_type;
            sse_initial_mass = p->sse_initial_mass;
            mt = p->mass;
            //age = p->age*yr_to_Myr;
            age = 0.0;

            GB = new double[10];
            tscls = new double[20];
            lums = new double[10];    
            
//            star_(&kw, &mass, &mt, &tm, &tn, tscls, lums, GB, zpars);
//            hrdiag_(&mass,&age,&mt,&tm,&tn,tscls,lums,GB,zpars,
//                &r,&lum,&kw,&mc,&rc,&menv,&renv,&k2);
            
            double dtp=0.0;
            double ospin=0.0;
            double epoch = 0.0;
            double tms;
            double tphys = 0.0;
            double tphysf = 0.0;;
            double dt;

            if (kw_desired == 0 or kw_desired == 1) /* ZAMS star; no aging */
            {
                kw = kw_desired;
                evolv1_(&kw,&sse_initial_mass,&mt,&r,&lum,&mc,&rc,&menv,&renv,&ospin,&epoch,&tms,&tphys,&tphysf,&dtp,&z,zpars,&k2);
            }
            else /* age star until reaching desired stellar type */
            {
                //printf("STEV PRE kw_desired %d ev m %g\n",kw_desired,mass);
                double dtm,dtr;
                kw = 1;
                int j=0;
                while (kw < kw_desired)
                {
                    //r=lum=mc=rc=menv=renv=tms=0.0;
                    //printf("input %d kw minit %g mt %g tphys %g tphysf %g\n",kw,p->sse_initial_mass,mt,tphys,tphysf);
                    evolv1_(&kw,&sse_initial_mass,&mt,&r,&lum,&mc,&rc,&menv,&renv,&ospin,&epoch,&tms,&tphys,&tphysf,&dtp,&z,zpars,&k2);
                    age = tphysf - epoch;

                    //star_(&kw, &mass, &mt, &tm, &tn, tscls, lums, GB, zpars);
                    //deltat_(&kw,&age,&tm,&tn,tscls,&dtm,&dtr);
                    dt = get_new_dt(kw,sse_initial_mass,mt,age,dt,zpars);
                    //dt = max(dtm,dtr);
//                    dt = min(dtr, dtm);
//                    dt = max(dt,1.0d-07*age);
//                    dt = max(dt,1.0e-15);

                    //dtm = min(dtr, dt);
                    //dtm = max(dtm,1.0d-07*age);
    
                    //printf("&kw %d, &mass %g, &mt %g, &tm %g, &tn %g dtm %g dtr %g epoch %g\n",kw, mass, mt, tm, tn, dtm, dtr,epoch);
                    tphys = tphysf;
                    tphysf += dt;
                    //printf("kw %d tphys %g tphysf %g dt %g\n",kw,tphys,tphysf,dt);
                    j++;
                    
                    if (j>10000)
                    {
                        printf("stellar_evolution.cpp -- unable to evolve star with index %d to desired stellar type %d; current stellar type %d; stopping initial evolution lum %g\n",p->index,kw_desired,kw,lum);
                        break;
                    }
                }
            }
            //printf("pre S %d\n",kw);
  //          star_(&kw, &sse_initial_mass, &mt, &TM1, &TN1, TSCLS1, LUMS1, GB1, ZPARS1);
    //printf("t2 %g %d %g %g %g %g\n",ZPARS1[0],KW1,M01,M1,TM1,TN1);

    //        hrdiag_(&M01,&AJ1,&M1,&TM1,&TN1,TSCLS1,LUMS1,GB1,ZPARS1, \
                &R1,&L1,&KW1,&MC1,&RC1,&MENV1,&RENV1,&K21);
            #ifdef IGNORE
            star_(&kw, &sse_initial_mass, &mt, &tm, &tn, tscls, lums, GB, zpars);
            printf("pre SS %d age %g sse_initial_mass %g mt %g\n",kw,age,sse_initial_mass,mt);
            hrdiag_(&sse_initial_mass,&age,&mt,&tm,&tn,tscls,lums,GB,zpars,
                &r,&lum,&kw,&mc,&rc,&menv,&renv,&k2);
            hrdiag_(&sse_initial_mass,&age,&mt,&tm,&tn,tscls,lums,GB,zpars,
                &r,&lum,&kw,&mc,&rc,&menv,&renv,&k2);
            hrdiag_(&sse_initial_mass,&age,&mt,&tm,&tn,tscls,lums,GB,zpars,
                &r,&lum,&kw,&mc,&rc,&menv,&renv,&k2);
            
            printf("arg hrd %g %g %g %g %g \n",sse_initial_mass,age,mt,tm,tn);
            printf("post S %d\n",kw);
            #endif
            
            //printf("stellar_evolution.cpp -- initialize_star -- done evolve; kw %d sse_initial_mass %g mt %g\n",kw,sse_initial_mass,mt);
            p->stellar_type = kw;

            p->sse_initial_mass = sse_initial_mass;
            p->mass = mt;
            p->radius = r*CONST_R_SUN;

            p->age = age*Myr_to_yr;
            p->sse_main_sequence_timescale = tms*Myr_to_yr;
            p->epoch = epoch*Myr_to_yr;

            p->luminosity = lum*CONST_L_SUN;
            p->core_mass = mc;
            p->core_radius = rc*CONST_R_SUN;
            p->convective_envelope_mass = menv;
            p->convective_envelope_radius = renv*CONST_R_SUN;

            p->sse_k2 = k2;
            p->sse_k3 = 0.21;

            /* Set up spins parallel with parent orbit */
            Particle *parent = (*particlesMap)[p->parent];
            //double h_vec[3] = {parent->h_vec_x,parent->h_vec_y,parent->h_vec_z};
            double *h_vec = parent->h_vec;
            double h = norm3(h_vec);
            
            for (i=0; i<3; i++)
            {
                p->spin_vec[i] = ospin*h_vec[i]/h;
            }
            //p->spin_vec_x = ospin*h_vec[0]/h;
            //p->spin_vec_y = ospin*h_vec[1]/h;
            //p->spin_vec_z = ospin*h_vec[2]/h;
            //printf("initialize %g %d %g %g %g\n",ospin,p->index,p->spin_vec_x,p->spin_vec_y,p->spin_vec_z);
            
            p->check_for_RLOF_at_pericentre = true;
        }
    }

    update_structure(particlesMap); /* set new binary masses */

    /* Adjust h vectors such that possible initial stellar evolution would not change the original orbital elements  */
    for (it_p = particlesMap->begin(); it_p != particlesMap->end(); it_p++)
    {
        Particle *p = (*it_p).second;
        if (p->is_binary == true)
        {
            double h = (p->child1_mass*p->child2_mass*sqrt(CONST_G*p->a/(p->child1_mass + p->child2_mass)))*sqrt(1.0 - p->e * p->e); /* new h based on new masses, but old a & e */
            for (i=0; i<3; i++)
            {
                p->h_vec[i] = p->h_vec_unit[i] * h;
            }
        }
    }

    return 0;
}

#ifdef IGNORE
int initialize_spins(ParticlesMap *particlesMap)
{
    ParticlesMapIterator it_p;
    for (it_p = particlesMap->begin(); it_p != particlesMap->end(); it_p++)
    {
        Particle *p = (*it_p).second;
        if (p->is_binary == 0 and p->evolve_as_star == 1)
        {
            Particle *parent = ParticlesMap[p->parent];
#endif            

int evolve_stars(ParticlesMap *particlesMap, double start_time, double end_time, double *stellar_evolution_timestep, bool get_timestep_only, bool *apply_SNe_effects)
{
    /* TO DO: apsidal_motion_constant? */
    //double mass = 5.0;

    value1_.neta = 0.5;
    value1_.bwind = 0.0;
    value1_.hewind = 0.5;
    value1_.mxns = 3.0;
    value2_.alpha1 = 1.0;
    value2_.lambda = 1.0;
    value3_.idum = 0;
    value4_.sigma = 190.0;
    value4_.bhflag = 1;
    //        value5_.beta = 0.0;
    //        value5_.xi = 0.0;
    //        value5_.acc2 = 0.0;
    //        value5_.epsnov = 0.0;
    //        value5_.eddfac = 0.0;
    //        value5_.gamma = 0.0;
    flags_.ceflag = binary_evolution_CE_energy_flag;
    flags_.tflag = 0;
    flags_.ifflag = 0;
    flags_.nsflag = 1;
    flags_.wdflag = 1;
    points_.pts1 = 0.05;
    points_.pts2 = 0.01;
    points_.pts3 = 0.02;

    ParticlesMapIterator it_p;
    std::vector<int>::iterator it_parent_p,it_parent_q;


    /* Go through all stars in the system and evolve them */
    double sse_initial_mass,sse_time_step;
    double mt,tphys,tphysf,dtp,epoch,ospin,ospin_old,age;
    double dt,dt_min;
    double r,lum,mc,rc,menv,renv,tms,k2;
    double r_old,mt_old;
    double rzams,fac,menv_fraction;
    int kw,kw_old;
    double spin_vec[3];
    dt_min = 1.0e100;
    
    double z;
    double *zpars;
    
    bool reset_h_vectors = false;
    *apply_SNe_effects = false;

    for (it_p = particlesMap->begin(); it_p != particlesMap->end(); it_p++)
    {
        Particle *p = (*it_p).second;
        if (p->is_binary == 1)
        {
            //double h_vec[3];
            //h_vec[0] = p->h_vec_x;
            //h_vec[1] = p->h_vec_y;
            //h_vec[2] = p->h_vec_z;
            //printf("old h %g \n",norm3(h_vec));
            p->child1_mass_old = p->child1_mass;
            p->child2_mass_old = p->child2_mass;
        }
    }
            
    for (it_p = particlesMap->begin(); it_p != particlesMap->end(); it_p++)
    {
        Particle *p = (*it_p).second;
        if (p->is_binary == false and p->evolve_as_star == true)
        {
            /* Extract stellar evolution parameters */
            mt = mt_old = p->mass;
            sse_initial_mass = p->sse_initial_mass;
            kw = kw_old = p->stellar_type;
            epoch = p->epoch*yr_to_Myr;
            r = r_old = p->radius/CONST_R_SUN;
            
            z = p->metallicity;
            zpars = p->zpars;
            //printf("SE %g\n",zpars[0]);
            //spin_vec[0] = p->spin_vec_x;
            //spin_vec[1] = p->spin_vec_y;
            //spin_vec[2] = p->spin_vec_z;
            //ospin_old = ospin = norm3(spin_vec);
            ospin_old = ospin = norm3(p->spin_vec);
            //printf("extract ? %g %g %g\n",spin_vec[0],spin_vec[1],spin_vec[2]);

            sse_time_step = p->sse_time_step*yr_to_Myr;
            age = p->age*yr_to_Myr;

            tphys = start_time*yr_to_Myr;
            tphysf = end_time*yr_to_Myr;
            double desired_tphysf = tphysf;
            
            dtp = 0.0;
            #ifdef DEBUG
            printf("stellar_evolution.cpp -- sse1 kw %d mt %g r %g sse_time_step %g epoch %g age %g tphys %g tphysf %g ospin %g\n",kw,mt,r,sse_time_step,epoch,age,tphys,tphysf,ospin);
            #endif
            //printf("Sx %g Sy %g Sz %g\n",p->spin_vec_x,p->spin_vec_y,p->spin_vec_z);
            if (get_timestep_only == true)
            {
                //printf("old dt %g kw %d\n",sse_time_step,kw);
                age = tphys - epoch;
                sse_time_step = get_new_dt(kw,sse_initial_mass,mt,age,sse_time_step,zpars);
                //printf("new dt %g kw %d\n",sse_time_step,kw);
            }
            else
            {
                #ifdef DEBUG
                printf("stellar_evolution.cpp -- INPUT evolv1_ kw %d m0 %.15f mt %.15f epoch %.15f tphys %.15f tphysf %.15f dtp %.15f z %.15f zpars[0] %.15f \n",kw,sse_initial_mass,mt,epoch,tphys,tphysf,dtp,z,zpars[0]);
                #endif
                evolv1_(&kw,&sse_initial_mass,&mt,&r,&lum,&mc,&rc,&menv,&renv,&ospin,&epoch,&tms,&tphys,&tphysf,&dtp,&z,zpars,&k2);

                if (tphysf != desired_tphysf)
                {
                    printf("ERROR tphysf != desired_tphysf ");
                }

                age = tphysf - epoch;
                sse_time_step = get_new_dt(kw,sse_initial_mass,mt,age,sse_time_step,zpars);
                #ifdef DEBUG
                printf("stellar_evolution.cpp -- sse2 kw %d mt %g r %g lum %g sse_time_step %g epoch %g age %g  tphys %g tphysf %g ospin %g\n",kw,mt,r,lum,sse_time_step,epoch,age,tphys,tphysf,ospin);
                #endif

                /* Write new stellar evolution parameters */
                if (kw_old==kw or kw < 13)
                {
                    dt = end_time - start_time;
                    p->mass_dot_wind = (mt - mt_old)/dt;
                    //printf("SE %.15f \n",p->mass_dot_wind);
                    //p->radius_dot = (r*CONST_R_SUN - p->radius)/dt;
                    p->radius_dot = CONST_R_SUN*(r - r_old)/dt;
                    //printf("p->radius_dot %g\n",p->radius_dot);
                    
                    p->ospin_dot = (ospin - ospin_old)/dt;
                    //printf("p->ospin_dot %g p->radius_dot %g\n",p->ospin_dot,p->radius_dot);
                    //p->ospin_dot = 0.0;
                    //printf("Md %g Rd %g\n",p->mass_dot,p->radius_dot);
                    //p->radius_dot = 0.1/dt;
                    //printf("old R %g m %g\n",p->radius,p->mass);
                    //printf("predicted R %g m %g dR/dt %g\n",p->radius + p->radius_dot*dt,p->mass + p->mass_dot*dt,p->radius_dot);
                    p->instantaneous_perturbation_delta_mass = 0.0;
                }
                else /* kw change with new kw=13 or kw=14 */
                { /* CAREFULL: the semimajor axis might still change in this case in secular code */
                    p->mass_dot_wind = 0.0;
                    p->radius_dot = 0.0;
                    p->ospin_dot = 0.0;
                    
                    /* New mass will be handled by handle_SNe_in_system(), but radius has to be adjusted here */
                    p->radius = r*CONST_R_SUN;
                    
                    //printf("switch kw_old %d kw_new %d t %g\n",kw_old,kw,start_time);
                    //reset_h_vectors = true;
                    *apply_SNe_effects = true;
                    
                    p->instantaneous_perturbation_delta_mass = mt - mt_old; /* Will be used by handle_SNe_in_system() */
                    //printf("dm %g\n",p->instantaneous_perturbation_delta_mass);
                    //p->instantaneous_perturbation_delta_mass = 0.1;

                }
                //printf("dots1 % g %g %g %g\n",p->mass_dot,p->radius_dot,r*CONST_R_SUN,p->radius);
                //printf("dots2 % g %g %g %g\n",p->mass_dot,p->radius_dot,mt,p->mass);

                //p->mass_dot = 0.0;
                //p->radius_dot = 0.0;

                //p->mass = mt;
                p->sse_initial_mass = sse_initial_mass;
                p->stellar_type = kw;
                p->epoch = epoch*Myr_to_yr;
                //printf("ospin %g ospin_old %g\n",ospin,ospin_old);
                //if (ospin_old != 0.0)
                if (1==0)
                {
                    //p->spin_vec_x *= ospin/ospin_old;
                    //p->spin_vec_y *= ospin/ospin_old;
                    //p->spin_vec_z *= ospin/ospin_old;
                }
                
                p->age = age*Myr_to_yr;
                p->sse_main_sequence_timescale = tms*Myr_to_yr;
                //p->radius = r*CONST_R_SUN;
                p->luminosity = lum*CONST_L_SUN;
                p->core_mass = mc;
                p->core_radius = rc*CONST_R_SUN;
                p->convective_envelope_mass = menv;
                p->convective_envelope_radius = renv*CONST_R_SUN;

                p->sse_k2 = k2;
                p->sse_k3 = 0.21;
                //printf("ks %g %g\n",p->sse_k2,p->sse_k3);
                //rzams = rzamsf_(&sse_initial_mass);
                //fac = value2_.lambda;
                //menv_fraction = menv/(mt - mc);
                //p->common_envelope_lambda = celamf_(&kw,&sse_initial_mass,&lum,&r,&rzams,&menv_fraction,&fac);
                //printf("kw %d lambda %g menv_fraction %g\n",kw,p->common_envelope_lambda,menv_fraction);
                
            }
            p->sse_time_step = sse_time_step*Myr_to_yr;

            /* Common time step */
            if (p->sse_time_step < dt_min)
            {
                dt_min = p->sse_time_step;
            }

        }
    }
    
    //printf("sse dt %g\n",dt_min);
    *stellar_evolution_timestep = dt_min;
    
    #ifdef IGNORE
    if (reset_h_vectors == true)
    {
        int N_bodies, N_binaries, N_root_finding;
        determine_binary_parents_and_levels(particlesMap, &N_bodies, &N_binaries, &N_root_finding);
        set_binary_masses_from_body_masses(particlesMap);

        for (it_p = particlesMap->begin(); it_p != particlesMap->end(); it_p++)
        {
            Particle *p = (*it_p).second;
            if (p->is_binary == 1)
            {

                double factor_h = ( (p->child1_mass * p->child2_mass) / (p->child1_mass_old * p->child2_mass_old) ) * sqrt( (p->child1_mass_old + p->child2_mass_old) / (p->child1_mass + p->child2_mass) );
                printf("resetting h vectors f %g\n",factor_h);
                p->h_vec_x *= factor_h;
                p->h_vec_y *= factor_h;
                p->h_vec_z *= factor_h;
                
                double h_vec[3];
                h_vec[0] = p->h_vec_x;
                h_vec[1] = p->h_vec_y;
                h_vec[2] = p->h_vec_z;
                printf("new h %g \n",norm3(h_vec));
            }

        }
    }
    #endif
    
    return 0;
}

double get_new_dt(int kw, double mass, double mt, double age, double dt, double *zpars)
{
    //printf("get_new_dt %d %g %g %g %g\n",kw,mass,mt,age,dt);
    double dtm,dtr;
    double tm,tn;

    double *GB,*tscls,*lums;
    GB = new double[10];
    tscls = new double[20];
    lums = new double[10];    
    //printf("ASDASDSAD %g\n",dt);
    star_(&kw, &mass, &mt, &tm, &tn, tscls, lums, GB, zpars);
    //printf("tm %g age %g\n",tm,age);
    deltat_(&kw,&age,&tm,&tn,tscls,&dt,&dtr);
    //printf("qweqwe %g %g\n",dt,dtr);
    //dtm = std::min(dtr, dt);
    //dtm = std::max(dtm,1.0d-07*age);

    dtm = CV_min(dtr, dt);
    dtm = CV_max(dtm,1.0d-07*age);
    dtm = CV_max(dtm,1.0e-15);
    //printf("gdt %g %g %g\n",dtr,dtm,dt);
    
    return dtm;
}


void update_stellar_evolution_properties(Particle *p)
{
    /* hrdiag: 
     * Computes the new mass, luminosity, radius & stellar type.
    *  Input (MASS, AJ, TM, TN, LUMS & TSCLS) supplied by routine STAR. */
    
    double tm,tn;

    double *GB,*tscls,*lums;
    GB = new double[10];
    tscls = new double[20];
    lums = new double[10];  
    
    double age = p->age*yr_to_Myr;
    double r,rc,lum,menv,renv,k2;
    star_(&p->stellar_type, &p->sse_initial_mass, &p->mass, &tm, &tn, tscls, lums, GB, p->zpars);
    hrdiag_(&p->sse_initial_mass,&age,&p->mass,&tm,&tn,tscls,lums,GB,p->zpars, \
            &r,&lum,&p->stellar_type,&p->core_mass,&rc,&menv,&renv,&k2);

    p->radius = r*CONST_R_SUN;

    p->age = age*Myr_to_yr;

    p->luminosity = lum*CONST_L_SUN;
    p->core_radius = rc*CONST_R_SUN;
    p->convective_envelope_mass = menv;
    p->convective_envelope_radius = renv*CONST_R_SUN;
}

double compute_moment_of_inertia(double mass, double core_mass, double radius, double core_radius, double k2, double k3)
{
    return k2*(mass - core_mass)*radius*radius + k3*core_mass*core_radius*core_radius;
}

}
