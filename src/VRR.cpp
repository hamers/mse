/*
*/

#include "types.h"
#include "VRR.h"
#include "../interface.h" /* for parameters */
#include "structure.h" /* for determine_binary_parents_and_levels */
#include <stdio.h>

extern "C"
{

double compute_VRR_perturbations(ParticlesMap *particlesMap, int index, double time)
{
    
    int c;
    
    Particle *p = (*particlesMap)[index];

    double a = p->a;
    double e = p->e;
    double h = p->h;
    double *e_vec_unit = p->e_vec_unit;
    double *h_vec_unit = p->h_vec_unit;
    double q_vec_unit[3];
    cross3(h_vec_unit,e_vec_unit,q_vec_unit);
    
    double hamiltonian = 0.0;


    /* Mass precesssion (for all models; in case of model 2: make sure not to include twice) */
    if (p->VRR_include_mass_precession == 1)
    {
        double VRR_mass_precession_rate = p->VRR_mass_precession_rate;
        //printf("VRR_mass_precession_rate %g\n",VRR_mass_precession_rate);

        for (c=0; c<3; c++)
        {
            p->de_vec_dt[c] += e*VRR_mass_precession_rate*q_vec_unit[c];
        }
    }
    
    if (p->VRR_model == 1)
    {
        double Omega_vec[3] = {p->VRR_Omega_vec_x,p->VRR_Omega_vec_y,p->VRR_Omega_vec_z};
        double Omega_vec_cross_h_vec_unit[3];
        
        cross3(Omega_vec,h_vec_unit,Omega_vec_cross_h_vec_unit);

        //printf("Omega %g %g %g\n",Omega_vec[0],Omega_vec[1],Omega_vec[2]);
        //printf("O %g %g %g\n",Omega_top_vec_cross_h_vec_unit[0],Omega_top_vec_cross_h_vec_unit[1],Omega_top_vec_cross_h_vec_unit[2]);

        for (c=0; c<3; c++)
        {
            p->dh_vec_dt[c] += h*Omega_vec_cross_h_vec_unit[c];
            //printf("add top %g \n",h*Omega_top_vec_cross_h_vec_unit[c]);
        }
    }
    
    if (p->VRR_model == 2)
    {
        printf("ERROR: VRR Model 2 not supported!\n");
        exit(-1);
        
#ifdef IGNORE
        /* Distortion model */
//        printf("Distortion model\n");
       
        /* extract current orbital elements from orbital vectors */
        double x_vec[3] = {1.0,0.0,0.0};
        double y_vec[3] = {0.0,1.0,0.0};
        double z_vec[3] = {0.0,0.0,1.0};

        /* `i' = `inclination'; `c' = `cos'; `s' = `sin' */
        double ci = dot3(z_vec,h_vec_unit);
        double si = sqrt(1.0 - ci*ci);

        double Omega_vec[3],Omega_vec_unit[3];
        cross3( z_vec,h_vec_unit, Omega_vec);
        double Omega_vec_norm = norm3(Omega_vec);
        for (c=0; c<3; c++)
        {
            Omega_vec_unit[c] = Omega_vec[c]/Omega_vec_norm;
        }
        /* `obig' = `Omega' = `longitude of ascending node' */
        double cobig = dot3(x_vec,Omega_vec_unit);
        double sobig = dot3(y_vec,Omega_vec_unit);
        
        /* `omega' = `argument of periapsis' */
        double comega = dot3(e_vec_unit,Omega_vec_unit);
        double somega = -dot3(q_vec_unit,Omega_vec_unit);

        //double tempe[3] = {cobig*comega - sobig*somega*ci, sobig*comega + cobig*somega*ci,somega*si};
        //double temph[3] = {sobig*si, -cobig*si, ci};
        
        //printf("test e 1 %g %g 2 %g %g 3 %g %g\n",e_vec_unit[0],tempe[0],e_vec_unit[1],tempe[1],e_vec_unit[2],tempe[2]);
        //printf("test h 1 %g %g 2 %g %g 3 %g %g\n",h_vec_unit[0],temph[0],h_vec_unit[1],temph[1],h_vec_unit[2],temph[2]);
        
        /* model-related parameters */
        /* distortion vector */
        double S_init[3] = {p->VRR_initial_distortion_axis_vec_x,p->VRR_initial_distortion_axis_vec_y,p->VRR_initial_distortion_axis_vec_z};
        double S_final[3] = {p->VRR_final_distortion_axis_vec_x,p->VRR_final_distortion_axis_vec_y,p->VRR_final_distortion_axis_vec_z};
        double VRR_initial_time = p->VRR_initial_time;
        double VRR_final_time = p->VRR_final_time;
        double factor_t = (time - VRR_initial_time)/(VRR_final_time - VRR_initial_time);
        
        if (factor_t < 0.0)
        {
            factor_t = 0.0;
        }
        if (factor_t >= 1.0)
        {
            factor_t = 1.0;
        }
        double S[3];
        
        for (c=0; c<3; c++)
        {
            S[c] = S_init[c] + (S_final[c] - S_init[c])*factor_t;
        }
        double S_norm = norm3(S);
        for (c=0; c<3; c++)
        {
            S[c] /= S_norm;
        }
        double Sx = S[0];
        double Sy = S[1];
        double Sz = S[2];
        
        //printf("DM Sx %g Sy %g Sz %g\n",Sx,Sy,Sz);
        double fracdip = p->VRR_dipole_fraction;
        double fracquad = 1.0 - fracdip;
        
        if (fracdip < 0.0 || fracdip > 1.0)
        {
            printf("Invalid dipole fraction -- exiting\n");
            exit(-1);
        }

        double Ad = p->VRR_AD_parameter;
        double Am = p->VRR_AM_parameter;
        double Aj = p->VRR_AJ_parameter;
        
        
        //printf("Distortion model parameters Ad %g Am %g Aj %g S[0] %g S[1] %g S[2] %g fracdip %g fracquad %g\n",Ad,Am,Aj,S[0],S[1],S[2],fracdip,fracquad);

        /* time unit in the ODE's */
        double P_orb = compute_orbital_period(p);
        double nu_radial = 2.0*M_PI/P_orb;
        
        double rg = CONST_G*p->child1_mass_plus_child2_mass/(CONST_C_LIGHT_P2); // gravitational radius
        double nu_0 = nu_radial*3.0*(rg/a); // GR precession frequency

        double ell = sqrt(1.0 - e*e);
        double ell_p2 = ell*ell;
        double ell_p3 = ell*ell_p2;
        
        double e_p2 = e*e;
        double e_p4 = e_p2*e_p2;
        
        //printf("Distortion model test P %g nu_0 %g nu_rad %g ell %g rg %g a %g\n",P_orb,nu_0,nu_radial,ell,rg,a);
        
        double Acoef = comega*cobig - ci*somega*sobig;
        double Bcoef = -somega*cobig - ci*comega*sobig;
        double Ccoef = comega*sobig + ci*somega*cobig;
        double Dcoef = -somega*sobig + ci*comega*cobig;
        double Fcoef = si*somega;
        double Gcoef = si*comega;

        /* dipole */
        double hd = Ad*e*(Sx*Acoef + Sy*Ccoef + Sz*Fcoef);
        hamiltonian += hd;

        double dip1 = -Ad*Sz*e*Gcoef - Ad*Sx*e*Bcoef - Ad*Sy*e*Dcoef; // associated with dell/dtau
        double dip2 = Ad*e/ell*somega*(-Sx*sobig + Sy*cobig - Sz*ci/si); // associated with dOmega/dtau
        double dip3 = -(ell/e)*(hd/e) + Ad*e/ell*ci*somega*(Sx*sobig - Sy*cobig + Sz*ci/si); // associated with domega/dtau
        double dip4 = Ad*Sx*e*Ccoef - Ad*Sy*e*Acoef; // associated with dellz/dtau
        
//      dip2 = Ad*e/ell*somega*(-Sx*sobig + Sy*cobig)
//      dip3 = -(ell/e)*(hd/e) + Ad*e/ell*ci*somega*(Sx*sobig - Sy*cobig)
//      else
//      dCFdlz =  Fcoef*cobig*somega/ell
//      dDGdlz =  Gcoef*cobig*comega/ell
//      dAFdlz = -Fcoef*sobig*somega/ell
//      dBGdlz = -Gcoef*sobig*comega/ell

        /* quadrupole */
        double gfunc = ((1.0-2.0*e_p2)*(1.0-ell)+e_p4)/e_p2/ell;
        double hfunc = 1.0 - gfunc;
        double Sx2 = Sx*Sx;
        double Sy2 = Sy*Sy;
        double Sz2 = Sz*Sz;
        double Sxy = Sx*Sy;
        double Syz = Sy*Sz;
        double Sxz = Sx*Sz;
      
        double quad1 = 0.0;
        double quad2 = 0.0;
        double quad3 = 0.0;
        double quad4 = 0.0;
      
        double dA2do = 2.0*Acoef*Bcoef;
        double dB2do = -dA2do;
        double dC2do = 2.0*Ccoef*Dcoef;
        double dD2do = -dC2do;
        double dF2do = 2.0*Fcoef*Gcoef;
        double dG2do = -dF2do;
        double dACdo = Acoef*Dcoef + Bcoef*Ccoef;
        double dBDdo = -dACdo;
        double dCFdo = Ccoef*Gcoef + Dcoef*Fcoef;
        double dDGdo = -dCFdo;
        double dAFdo = Acoef*Gcoef + Bcoef*Fcoef;
        double dBGdo = -dAFdo;
        
        quad1 = Sx2*(dA2do*gfunc + dB2do*hfunc) \
            + Sy2*(dC2do*gfunc + dD2do*hfunc) \
            + Sz2*(dF2do*gfunc + dG2do*hfunc) \
            + 2.0*Sxy*(dACdo*gfunc + dBDdo*hfunc) \
            + 2.0*Syz*(dCFdo*gfunc + dDGdo*hfunc) \
            + 2.0*Sxz*(dAFdo*gfunc + dBGdo*hfunc);
        quad1 = -Ad*quad1;
        
        double dA2dlz = -2.0*Acoef/ell*sobig*somega;
        double dB2dlz = -2.0*Bcoef/ell*sobig*comega;
        double dC2dlz =  2.0*Ccoef/ell*cobig*somega;
        double dD2dlz =  2.0*Dcoef/ell*cobig*comega;
        double dF2dlz = -2.0/ell*somega*somega*ci;
        double dG2dlz = -2.0/ell*comega*comega*ci;
        double dACdlz = somega/ell*(Acoef*cobig - Ccoef*sobig);
        double dBDdlz = comega/ell*(Bcoef*cobig - Dcoef*sobig);

        double dCFdlz = somega/ell*(-Ccoef*ci/si + Fcoef*cobig);
        double dDGdlz = comega/ell*(-Dcoef*ci/si + Gcoef*cobig);
        double dAFdlz = -somega/ell*(Acoef*ci/si + Fcoef*sobig);
        double dBGdlz = -comega/ell*(Bcoef*ci/si + Gcoef*sobig);

        quad2 = Sx2*(dA2dlz*gfunc + dB2dlz*hfunc) \
            + Sy2*(dC2dlz*gfunc + dD2dlz*hfunc) \
            + Sz2*(dF2dlz*gfunc + dG2dlz*hfunc) \
            + 2.0*Sxy*(dACdlz*gfunc + dBDdlz*hfunc) \
            + 2.0*Syz*(dCFdlz*gfunc + dDGdlz*hfunc) \
            + 2.0*Sxz*(dAFdlz*gfunc + dBGdlz*hfunc);
        quad2 = Ad*quad2;

        quad3 = Ad*quad3;
        
        double dA2dobig = -2.0*Ccoef*Acoef;
        double dB2dobig = -2.0*Bcoef*Dcoef;
        double dC2dobig = -dA2dobig;
        double dD2dobig = -dB2dobig;
        double dF2dobig = 0.0;
        double dG2dobig = 0.0;
        double dACdobig = -Ccoef*Ccoef + Acoef*Acoef;
        double dBDdobig =  Bcoef*Bcoef - Dcoef*Dcoef;
        double dCFdobig = Fcoef*Acoef;
        double dDGdobig = Gcoef*Bcoef;
        double dAFdobig = -Fcoef*Ccoef;
        double dBGdobig = -Gcoef*Dcoef;
        
        quad4 = Sx2*(dA2dobig*gfunc + dB2dobig*hfunc) \
            + Sy2*(dC2dobig*gfunc + dD2dobig*hfunc) \
            + 2.0*Sxy*(dACdobig*gfunc + dBDdobig*hfunc) \
            + 2.0*Syz*(dCFdobig*gfunc + dDGdobig*hfunc) \
            + 2.0*Sxz*(dAFdobig*gfunc + dBGdobig*hfunc);
        quad4 = -Ad*quad4;
        
        double dell_dtau = fracdip*dip1 + fracquad*quad1;
        double dobig_dtau = fracdip*dip2 + fracquad*quad2;
        double domega_dtau = fracdip*dip3 + fracquad*quad3;
        double dellz_dtau = fracdip*dip4 + fracquad*quad4;
        
//        dell_dtau = 0.0;
//        dobig_dtau = 0.0;
//        domega_dtau = 0.0;
//        dellz_dtau = 0.0;
        
        if (p->VRR_include_GR_precession == 1)
        {
            domega_dtau += 1.0/(ell_p2);
        }
        if (p->VRR_include_mass_precession == 1)
        {
            domega_dtau -= Am*ell/(1.0+ell);
        }
        if (p->VRR_include_frame_dragging == 1)
        {
            dobig_dtau += Aj/(3.0*ell_p3);
            domega_dtau -= Aj*ci/ell_p3;
        }

        double dell_dt = nu_0*dell_dtau;
        double dobig_dt = nu_0*dobig_dtau;
        double domega_dt = nu_0*domega_dtau;
        double dellz_dt = nu_0*dellz_dtau;
        
        double di_dt = (1.0/(ell*si))*(dell_dt*ci - dellz_dt);
        
        double dh_dt = (h/ell)*dell_dt;
        double dh_unit_dt[3];
        dh_unit_dt[0] = cobig*si*dobig_dt + sobig*ci*di_dt;
        dh_unit_dt[1] = sobig*si*dobig_dt - cobig*ci*di_dt;
        dh_unit_dt[2] = -si*di_dt;
        
        double de_dt = (-ell/e)*dell_dt;
        double de_unit_dt[3];
        de_unit_dt[0] = -sobig*comega*dobig_dt - cobig*somega*domega_dt - cobig*somega*ci*dobig_dt - sobig*comega*ci*domega_dt + sobig*somega*si*di_dt;
        de_unit_dt[1] = cobig*comega*dobig_dt - sobig*somega*domega_dt - sobig*somega*ci*dobig_dt + cobig*comega*ci*domega_dt - cobig*somega*si*di_dt;
        de_unit_dt[2] = comega*si*domega_dt + somega*ci*di_dt;
        
        //printf("test e %g\n",dot3(e_vec_unit,de_unit_dt));
        //printf("test h %g\n",dot3(h_vec_unit,dh_unit_dt));
        
        for (c=0; c<3; c++)
        {
            p->de_vec_dt[c] += e*de_unit_dt[c] + e_vec_unit[c]*de_dt;
            p->dh_vec_dt[c] += h*dh_unit_dt[c] + h_vec_unit[c]*dh_dt;
            
            //printf("DM e sc %g\n",de_dt);
            //printf("DM h sc %g\n",dh_dt);

            //printf("DM e unit %g %g %g\n",de_unit_dt[c]);
            //printf("DM h unit %g %g %g\n",dh_unit_dt[c]);
        }
#endif
        
    }
        
        
    if (p->VRR_model == 3) 
    {
        /* Bar-Or model */
        //printf("Bar-Or model\n");

        double eta_20_init = p->VRR_eta_20_init;
        double eta_a_22_init = p->VRR_eta_a_22_init;
        double eta_b_22_init = p->VRR_eta_b_22_init;
        double eta_a_21_init = p->VRR_eta_a_21_init;
        double eta_b_21_init = p->VRR_eta_b_21_init;

        double eta_20_final = p->VRR_eta_20_final;
        double eta_a_22_final = p->VRR_eta_a_22_final;
        double eta_b_22_final = p->VRR_eta_b_22_final;
        double eta_a_21_final = p->VRR_eta_a_21_final;
        double eta_b_21_final = p->VRR_eta_b_21_final;

        //printf("test etas init %g %g %g %g %g\n",eta_20_init,eta_a_22_init,eta_b_22_init,eta_a_21_init,eta_b_21_init);
        //printf("test etas final %g %g %g %g %g\n",eta_20_final,eta_a_22_final,eta_b_22_final,eta_a_21_final,eta_b_21_final);

        double Omega_vec[3];
        //double common_factor = (3.0/2.0)*sqrt(5.0/(3.0*M_PI));
        
        double common_factor = sqrt(c_3div4); /* changed Feb 5 2018 */
        
        int N_x = 3;
        int N_y = 3;
        double M_init[3][3] = {{eta_a_22_init,eta_b_22_init,-eta_a_21_init},{eta_b_22_init,-eta_a_22_init,-eta_b_21_init},{-eta_a_21_init,-eta_b_21_init,sqrt(3.0)*eta_20_init}};
        double M_final[3][3] = {{eta_a_22_final,eta_b_22_final,-eta_a_21_final},{eta_b_22_final,-eta_a_22_final,-eta_b_21_final},{-eta_a_21_final,-eta_b_21_final,sqrt(3.0)*eta_20_final}};
        
        //printf("test 1 M_init %g %g %g\n",M_init[0][0],M_init[1][1],M_init[2][2]);
        scalar_times_matrix(N_x,N_y,M_init,common_factor);
        scalar_times_matrix(N_x,N_y,M_final,common_factor);
        //printf("test 2 M_init %g %g %g\n",M_init[0][0],M_init[1][1],M_init[2][2]);
        
        double VRR_initial_time = p->VRR_initial_time;
        double VRR_final_time = p->VRR_final_time;
        double factor_t = (time - VRR_initial_time)/(VRR_final_time - VRR_initial_time);

        if (factor_t < 0.0)
        {
            factor_t = 0.0;
        }
        if (factor_t >= 1.0)
        {
            factor_t = 1.0;
        }
        double M[3][3];
        for (int i_x=0; i_x<N_x; i_x++)
        {
            for (int i_y=0; i_y<N_y; i_y++)
            {
                M[i_x][i_y] = M_init[i_x][i_y] + (M_final[i_x][i_y] - M_init[i_x][i_y])*factor_t;
            }
        }
        //printf("ts %g %g %g\n",time,VRR_initial_time,VRR_final_time);
        //printf("test 2 M %g %g %g t %g \n",M[0][0],M[1][1],M[2][2],factor_t);
        matrix_on_vector(N_x,N_y,M,h_vec_unit,Omega_vec);

        //printf("test x M %g %g %g\n",M[0][0],M[0][1],M[0][2]);
        //printf("test y M %g %g %g\n",M[1][0],M[1][1],M[1][2]);
        //printf("test z M %g %g %g\n",M[2][0],M[2][1],M[2][2]);
        //printf("h_vec_unit %g %g %g\n",h_vec_unit[0],h_vec_unit[1],h_vec_unit[2]);
        //printf("Omega_vec %g %g %g\n",Omega_vec[0],Omega_vec[1],Omega_vec[2]);

        //Omega_vec[0] = common_factor*( eta_a_22*h_vec_unit
        double Omega_cross_h_vec_unit[3];
        cross3(Omega_vec,h_vec_unit,Omega_cross_h_vec_unit);

        //printf("test O %g %g %g Oxh %g %g %g\n",Omega_vec[0],Omega_vec[1],Omega_vec[2],Omega_cross_h_vec_unit[0],Omega_cross_h_vec_unit[1],Omega_cross_h_vec_unit[2]);
        for (c=0; c<3; c++)
        {
            p->dh_vec_dt[c] += h*Omega_cross_h_vec_unit[c];
        }
    }


    return hamiltonian;
}



}
