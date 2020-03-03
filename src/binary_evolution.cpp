/* SecularMultiple */
/* Adrian Hamers November 2019 */

#include "evolve.h"
#include "binary_evolution.h"

extern "C"
{


int handle_binary_evolution(ParticlesMap *particlesMap, double *dt_binary_evolution)
{
    
    compute_mass_transfer_rates(particlesMap,dt_binary_evolution);
    
    //*dt_binary_evolution = 1.0e4;
    
    return 0;
}



int compute_mass_transfer_rates(ParticlesMap *particlesMap, double *dt_binary_evolution)
{
    *dt_binary_evolution = 1.0e100;
    double min_dt = 1.0e0;

    double epsilon = 0.01;
    double dt_temp;

    ParticlesMapIterator it_p;
    for (it_p = particlesMap->begin(); it_p != particlesMap->end(); it_p++)
    {
        Particle *donor = (*it_p).second;
        if (donor->is_binary == 0)
        {
            if (donor->RLOF_flag == 1)
            {
                
                Particle *parent = (*particlesMap)[donor->parent];
                Particle *accretor = (*particlesMap)[donor->sibling];
                double P_orb = compute_orbital_period(parent);
            
                double fm = donor->fm;
    /* child1 -> child2 */
    //q = child1->mass/child2->mass;
    //R_Lc = roche_radius_pericenter_eggleton(a,q); /* with argument "a", actually computes circular Roche lobe radius */
    //printf("R_lc %g a %g q %g M %g R %g\n",R_Lc,a,q,child1->mass,child1->radius);
    //x = R_Lc/child1->radius;
    //determine_E_0(e, x, &E_0, &in_RLOF);
    
                donor->mass_dot_RLOF = -(donor->mass/P_orb)*fm;
                accretor->mass_dot_RLOF = -donor->mass_dot_RLOF;
                
                dt_temp = epsilon*fabs(donor->mass/donor->mass_dot_RLOF);
                //printf("0 dt_temp %g *dt_binary_evolution %g\n",dt_temp,*dt_binary_evolution);
                *dt_binary_evolution = min(*dt_binary_evolution, dt_temp);
                //printf("1 dt_temp %g *dt_binary_evolution %g\n",dt_temp,*dt_binary_evolution);
                *dt_binary_evolution = max(*dt_binary_evolution, min_dt);
                //printf("2 dt_temp %g *dt_binary_evolution %g\n",dt_temp,*dt_binary_evolution);
                
                #ifdef VERBOSE
                printf("fm %g M_dot %g dt_MT %g \n",fm,donor->mass_dot_RLOF,*dt_binary_evolution);
                #endif
    //printf("e %g x %g E_0 %g m_d %g m_a %g fm %g md %g\n",e,x,E_0,M_d,M_a,fm,M_d_dot_av);
                
                
            }
        }

    }
    
    return 0;
    
}

int determine_mass_transfer_stability(ParticlesMap *particlesMap)
{
    
    
    return 0;
}

}



#ifdef IGNORE
 8    dme = 2.08d-03*eddfac*(1.d0/(1.d0 + zpars(11)))*rad(j2)*tb
      supedd = .false.
      novae = .false.
      disk = .false.
*
* Determine whether the transferred material forms an accretion 
* disk around the secondary or hits the secondary in a direct 
* stream, by using eq.(1) of Ulrich & Burger (1976, ApJ, 206, 509) 
* fitted to the calculations of Lubow & Shu (1974, ApJ, 198, 383). 
*
*     if(kstar(j2).ge.10) disk = .true.
      rmin = 0.0425d0*sep*(q(j2)*(1.d0+q(j2)))**(1.d0/4.d0)
      if(rmin.gt.rad(j2)) disk = .true.
*
* Kelvin-Helmholtz time from the modified classical expression.
*
      do 13 , k = 1,2
         tkh(k) = 1.0d+07*mass(k)/(rad(k)*lumin(k))
         if(kstar(k).le.1.or.kstar(k).eq.7.or.kstar(k).ge.10)then
            tkh(k) = tkh(k)*mass(k)
         else
            tkh(k) = tkh(k)*(mass(k) - massc(k))
         endif
 13   continue
*
* Dynamical timescale for the primary.
*
      tdyn = 5.05d-05*SQRT(rad(j1)**3/mass(j1))
*
* Identify special cases.
*
      if(kstar(j1).eq.2)then
         qc = 4.d0
      elseif(kstar(j1).eq.3.or.kstar(j1).eq.5.or.kstar(j1).eq.6)then
*        qc = (1.67d0-zpars(7)+2.d0*(massc(j1)/mass(j1))**5)/2.13d0
* Alternatively use condition of Hjellming & Webbink, 1987, ApJ, 318, 794. 
         qc = 0.362 + 1.0/(3.0*(1.0 - massc(j1)/mass(j1)))
* Or allow all cases to avoid common-envelope.  
*        qc = 100.d0
      elseif(kstar(j1).eq.8.or.kstar(j1).eq.9)then
         qc = 0.784d0
      else
         qc = 3.d0
      endif
*
      if(kstar(j1).eq.0.and.q(j1).gt.0.695d0)then
*
* This will be dynamical mass transfer of a similar nature to
* common-envelope evolution.  The result is always a single
* star placed in *2.
*
         taum = SQRT(tkh(j1)*tdyn)
         dm1 = mass(j1)
         if(kstar(j2).le.1)then
*
* Restrict accretion to thermal timescale of secondary.
*
            dm2 = taum/tkh(j2)*dm1
            mass(j2) = mass(j2) + dm2
*
* Rejuvenate if the star is still on the main sequence.
*
            mass0(j2) = mass(j2)
            CALL star(kstar(j2),mass0(j2),mass(j2),tmsnew,tn,
     &                tscls,lums,GB,zpars)
* If the star has no convective core then the effective age decreases,
* otherwise it will become younger still.
            if(mass(j2).lt.0.35d0.or.mass(j2).gt.1.25d0)then
               aj(j2) = tmsnew/tms(j2)*aj(j2)*(mass(j2) - dm2)/mass(j2)
            else
               aj(j2) = tmsnew/tms(j2)*aj(j2)
            endif
            epoch(j2) = tphys - aj(j2)
         elseif(kstar(j2).le.6)then
*
* Add all the material to the giant's envelope.
*
            dm2 = dm1
            mass(j2) = mass(j2) + dm2
            if(kstar(j2).eq.2)then
               mass0(j2) = mass(j2)
               CALL star(kstar(j2),mass0(j2),mass(j2),tmsnew,tn,tscls,
     &                   lums,GB,zpars)
               aj(j2) = tmsnew + tscls(1)*(aj(j2)-tms(j2))/tbgb(j2)
               epoch(j2) = tphys - aj(j2)
            endif
         elseif(kstar(j2).le.12)then
*
* Form a new giant envelope.
*
            dm2 = dm1
            kst = ktype(kstar(j1),kstar(j2))
            if(kst.gt.100) kst = kst - 100
            if(kst.eq.4)then
               aj(j2) = aj(j2)/tms(j2)
               massc(j2) = mass(j2)
            endif
*
* Check for planets or low-mass WDs.
*
            if((kstar(j2).eq.10.and.mass(j2).lt.0.05d0).or.
     &         (kstar(j2).ge.11.and.mass(j2).lt.0.5d0))then
               kst = kstar(j1)
               mass(j1) = mass(j2) + dm2
               mass(j2) = 0.d0
            else
               mass(j2) = mass(j2) + dm2
               CALL gntage(massc(j2),mass(j2),kst,zpars,
     &                     mass0(j2),aj(j2))
               epoch(j2) = tphys - aj(j2)
            endif
            kstar(j2) = kst
         else
*
* The neutron star or black hole simply accretes at the Eddington rate.
*
            dm2 = MIN(dme*taum/tb,dm1)
            if(dm2.lt.dm1) supedd = .true. 
            mass(j2) = mass(j2) + dm2
         endif
         coel = .true.
         if(mass(j2).gt.0.d0)then
            mass(j1) = 0.d0
            kstar(j1) = 15
         else
            kstar(j1) = kstar(j2)
            kstar(j2) = 15
         endif
         goto 135
      elseif(((ABS(ABS(2*kstar(j1)-11)-3).eq.2.or.kstar(j1).eq.9).
     &        and.(q(j1).gt.qc.or.radx(j1).le.radc(j1))).or.
     &        (kstar(j1).eq.2.and.q(j1).gt.qc).or.
     &        (kstar(j1).eq.4.and.q(j1).gt.qc))then
*
* Common-envelope evolution.
*
         m1ce = mass(j1)
         m2ce = mass(j2)
         CALL comenv(mass0(j1),mass(j1),massc(j1),aj(j1),jspin(j1),
     &               kstar(j1),mass0(j2),mass(j2),massc(j2),aj(j2),
     &               jspin(j2),kstar(j2),zpars,ecc,sep,jorb,coel)
*
         jp = MIN(80,jp + 1)
         bpp(jp,1) = tphys
         bpp(jp,2) = mass(1)
         if(kstar(1).eq.15) bpp(jp,2) = mass0(1)
         bpp(jp,3) = mass(2)
         if(kstar(2).eq.15) bpp(jp,3) = mass0(2)
         bpp(jp,4) = float(kstar(1))
         bpp(jp,5) = float(kstar(2))
         bpp(jp,6) = sep
         bpp(jp,7) = ecc
         bpp(jp,8) = rad(1)/rol(1)
         bpp(jp,9) = rad(2)/rol(2)
         bpp(jp,10) = 7.0
*
         epoch(j1) = tphys - aj(j1)
         if(coel)then
            com = .true.
            goto 135
         endif
         epoch(j2) = tphys - aj(j2)
         if(ecc.gt.1.d0)then
            if(kstar(1).ge.13)then
               rc = corerd(kstar(1),mass(1),mass(1),zpars(2))
               ospin(1) = jspin(1)/(k3*rc*rc*mass(1))
            endif
            if(kstar(2).ge.13)then
               rc = corerd(kstar(2),mass(2),mass(2),zpars(2))
               ospin(2) = jspin(2)/(k3*rc*rc*mass(2))
            endif
            goto 135
         endif
*
* Next step should be made without changing the time.
*
         dm1 = m1ce - mass(j1)
         dm2 = mass(j2) - m2ce
         dm22 = dm2
         dtm = 0.d0
*
* Reset orbital parameters as separation may have changed.
*
         tb = (sep/aursun)*SQRT(sep/(aursun*(mass(1)+mass(2))))
         oorb = twopi/tb
      elseif(kstar(j1).ge.10.and.kstar(j1).le.12.and.
     &       q(j1).gt.0.628d0)then
*
* Dynamic transfer from a white dwarf.  Secondary will have KW > 9.
*
         taum = SQRT(tkh(j1)*tdyn)
         dm1 = mass(j1)
         if(eddfac.lt.10.d0)then
            dm2 = MIN(dme*taum/tb,dm1)
            if(dm2.lt.dm1) supedd = .true. 
         else
            dm2 = dm1
         endif
         mass(j2) = mass(j2) + dm2
*
         if(kstar(j1).eq.10.and.kstar(j2).eq.10)then
*
* Assume the energy released by ignition of the triple-alpha reaction 
* is enough to destroy the star. 
*
            kstar(j2) = 15
            mass(j2) = 0.d0
         elseif(kstar(j1).eq.10.or.kstar(j2).eq.10)then
*
* Should be helium overflowing onto a CO or ONe core in which case the 
* helium swells up to form a giant envelope so a HeGB star is formed. 
* Allowance for the rare case of CO or ONe flowing onto He is made. 
*
            kst = 9
            if(kstar(j2).eq.10) massc(j2) = dm2
            CALL gntage(massc(j2),mass(j2),kst,zpars,mass0(j2),aj(j2))
            kstar(j2) = kst
            epoch(j2) = tphys - aj(j2)
         elseif(kstar(j2).le.12)then
            mass0(j2) = mass(j2)
            if(kstar(j1).eq.12.and.kstar(j2).eq.11)then
*
* Mixture of ONe and CO will result in an ONe product. 
*
               kstar(j2) = 12
            endif
         endif
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
      elseif(kstar(j1).eq.13)then
*
* Gamma ray burster?
*
         dm1 = mass(j1)
         mass(j1) = 0.d0
         kstar(j1) = 15
         dm2 = dm1
         mass(j2) = mass(j2) + dm2
         kstar(j2) = 14
         coel = .true.
         goto 135
      elseif(kstar(j1).eq.14)then
*
* Both stars are black holes.  Let them merge quietly.
*
         dm1 = mass(j1)
         mass(j1) = 0.d0
         kstar(j1) = 15
         dm2 = dm1
         mass(j2) = mass(j2) + dm2
         coel = .true.
         goto 135
      else
*
* Mass transfer in one Kepler orbit.
*
         dm1 = 3.0d-06*tb*(LOG(rad(j1)/rol(j1))**3)*
     &         MIN(mass(j1),5.d0)**2
         if(kstar(j1).eq.2)then
            mew = (mass(j1) - massc(j1))/mass(j1)
            dm1 = MAX(mew,0.01d0)*dm1
         elseif(kstar(j1).ge.10)then
*           dm1 = dm1*1.0d+03/MAX(rad(j1),1.0d-04)
            dm1 = dm1*1.0d+03*mass(j1)/MAX(rad(j1),1.0d-04)
         endif
         kst = kstar(j2)
*
#endif
