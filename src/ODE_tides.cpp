/* MSE */

#include "evolve.h"
#include "ODE_tides.h"
#include "tools.h" /* for orbital_period */
#include <stdio.h>

extern "C"
{
double const MINIMUM_MASS_FOR_RADIATIVE_DAMPING_MSUN = 1.2; // in future: make user-adjustable
int const MAIN_SEQUENCE = 1;
int const CHeB = 4;
int const HeMS = 7;
int const HeWD = 10;


void handle_secular_tidal_evolution(ParticlesMap *particlesMap, Particle *p)
{

    Particle *child1 = (*particlesMap)[p->child1];
    Particle *child2 = (*particlesMap)[p->child2];
    
    if (child1->include_tidal_friction_terms == true || child1->include_tidal_bulges_precession_terms == true || child1->include_rotation_precession_terms == true)
    {
        if (child1->tides_method == 0 || child1->tides_method == 1)
        {
            compute_EOM_equilibrium_tide_BO_full(particlesMap,p->index,child1->index,child2->index,child1->include_tidal_friction_terms,child1->include_tidal_bulges_precession_terms,child1->include_rotation_precession_terms,child1->minimum_eccentricity_for_tidal_precession,child1->tides_method);
        }
        else if (child1->tides_method == 2)
        {
            compute_EOM_equilibrium_tide(particlesMap,p->index,child1->index,child2->index,child1->include_tidal_friction_terms,child1->include_tidal_bulges_precession_terms,child1->include_rotation_precession_terms,child1->minimum_eccentricity_for_tidal_precession);
        }
    }
    if (child2->include_tidal_friction_terms == true || child2->include_tidal_bulges_precession_terms == true || child2->include_rotation_precession_terms == true)
    {
        if (child2->tides_method == 0 || child2->tides_method == 1)
        {
            compute_EOM_equilibrium_tide_BO_full(particlesMap,p->index,child2->index,child1->index,child2->include_tidal_friction_terms,child2->include_tidal_bulges_precession_terms,child2->include_rotation_precession_terms,child2->minimum_eccentricity_for_tidal_precession,child2->tides_method);
        }
        else if (child2->tides_method == 2)
        {
            compute_EOM_equilibrium_tide(particlesMap,p->index,child2->index,child1->index,child2->include_tidal_friction_terms,child2->include_tidal_bulges_precession_terms,child2->include_rotation_precession_terms,child2->minimum_eccentricity_for_tidal_precession);
        }

    }
}


int const PreMS = 17;

bool check_for_radiative_damping(int stellar_type, double mass, double convective_envelope_mass, double convective_envelope_radius)
{
    if ((stellar_type == MAIN_SEQUENCE) && (mass/CONST_MSUN >= MINIMUM_MASS_FOR_RADIATIVE_DAMPING_MSUN))
    {
        return true;
    }
    else if ((stellar_type == CHeB) || (stellar_type == HeMS))
    {
        return true;
    }
    else
    {
        return false;
    }
}   
bool check_for_convective_damping(int stellar_type)
{
    if (stellar_type == PreMS)
    {
        return true;
    }
    
    if (stellar_type < HeWD)
    {
        return true;
    }
    else
    {
        return false;
    }
}

double from_k_AM_div_T_to_t_V(double k_AM_div_T, double apsidal_motion_constant)
{
    double f = 2.0*apsidal_motion_constant + 1.0;
    return 3.0 * f * f / k_AM_div_T;
}

double compute_t_V(Particle *star, Particle *companion, double semimajor_axis)
{
    int tides_viscous_time_scale_prescription = star->tides_viscous_time_scale_prescription;
    double t_V = 1.0e10; /* large value by default, i.e. weak tides if tides_viscous_time_scale_prescription is not given the correct value */
    
    if (tides_viscous_time_scale_prescription == 0 or star->object_type != 1)
    {
        t_V = star->tides_viscous_time_scale;
    }
    else if (tides_viscous_time_scale_prescription == 1)
    {
        t_V = compute_t_V_hurley
        (
            star->stellar_type,
            star->mass,
            star->convective_envelope_mass,
            companion->mass,
            semimajor_axis,
            star->radius,
            star->convective_envelope_radius,
            star->luminosity,
            star->spin_vec_norm,
            star->gyration_radius,
            star->apsidal_motion_constant
        ); 
    }

    #ifdef VERBOSE
    if (verbose_flag > 2)
    {
        printf("ODE_tides.cpp -- compute_t_V star %d companion %d t_V %g\n",star->index,companion->index,t_V);
    }
    #endif

    return t_V;
}

double compute_t_V_hurley
(
    int stellar_type,
    double mass,
    double convective_envelope_mass,
    double companion_mass,
    double semimajor_axis,
    double radius,
    double convective_envelope_radius,
    double luminosity,
    double spin_angular_frequency,
    double gyration_radius,
    double apsidal_motion_constant
)
{
    #ifdef VERBOSE
    if (verbose_flag > 2)
    {
        printf("ODE_tides.cpp -- compute_t_V_hurley kw %d m %g menv %g mcomp %g a %g r %g renv %g lum %g omega %g rg %g kam %g \n",stellar_type,mass,convective_envelope_mass,companion_mass,semimajor_axis,radius,convective_envelope_radius,luminosity,spin_angular_frequency,gyration_radius,apsidal_motion_constant);
    }
    #endif
    
    bool USE_RADIATIVE_DAMPING = check_for_radiative_damping(stellar_type,mass,convective_envelope_mass,convective_envelope_radius);
    bool USE_CONVECTIVE_DAMPING = check_for_convective_damping(stellar_type);
    double k_AM_div_T,t_V;
    
    if (USE_CONVECTIVE_DAMPING == true && ((convective_envelope_mass <= 0.0) || (convective_envelope_radius <= 0.0)))
    {
        USE_RADIATIVE_DAMPING = true;
    }
    if (radius <= 0.0)
    {
        return 1.0e100;
    }
    
    if (USE_RADIATIVE_DAMPING == true) // radiative damping
    {
        double E2 = 1.592e-09*pow(mass/CONST_MSUN,2.84); // Hurley prescription; Zahn, 1977, A&A, 57, 383 and 1975, A&A, 41, 329
        
#ifdef IGNORE
        // Izzard's prescription not yet implemented
		else if(stardata->preferences->E2_prescription==E2_IZZARD)
		{
			if(stardata->star[star_number].stellar_type<HERTZSPRUNG_GAP)
			{
				double fburn=1.0;
				double x=stardata->star[star_number].aj/stardata->star[star_number].tms;
				if(x<fburn)
				{
					/* log mass and Z */
					double logm=log10(stardata->star[star_number].mass);
					double logz=log10(stardata->common.metallicity);
 
					/* fits for am and E20 */
					double am = 0.15*sin(3.2*logm)+0.31*logm;
					double E20 = -1.23250e+01+1.04550e+01*logm-4.15420e-01*logz-7.18650e+00*logm*logm+1.97940e+00*logm*logm*logm;
			  		E20=pow(10.0,E20);

					/* calc log(E2/E20) */
					E2 = -pow(x+am,4.0)*pow(MAX(1,x/0.95),30.0);
			  
					/* hence E2 */
					E2 = E20 * pow(10.0,E2);
			  
					/* verbosity */
					/*
					if(star->starnum==1)
					{
						printf("E2 kw=%d I=%g (fburn=%g x=%g E20=%g) H=%g\n",
						star->stellar_type,
						E2,
						fburn,x,
						E20,E2_Hurley);
					}
					*/
				}
				else
				{
					/* no conv core */
					E2=0.0;
				}
			}
			else
			{
				E2=0.0;
			}
		}
#endif
        
        k_AM_div_T = E2*pow(1.0 + companion_mass/mass,5.0/6.0)*radius*sqrt(CONST_G*mass/(pow(semimajor_axis,5.0)));
        t_V = from_k_AM_div_T_to_t_V(k_AM_div_T,apsidal_motion_constant);
        return t_V;
        
    }
    else if (USE_CONVECTIVE_DAMPING == true) // convective damping
    {
        double P_orb = compute_orbital_period_from_semimajor_axis(mass + companion_mass, semimajor_axis);
        double P_spin,P_tid;
        
        if (spin_angular_frequency < epsilon)
        {
            P_tid = P_orb;
        }
        else
        {
            P_spin = TWOPI/spin_angular_frequency;
            P_tid = 1.0/( 1e-10 + fabs( 1.0/P_orb - 1.0/P_spin) );
        }
        double radius_arg = (radius - (1.0/2.0)*convective_envelope_radius);
        if (radius_arg <= 0.0)
        {
            radius_arg = convective_envelope_radius;
        }
            
        double tau_convective = pow( (convective_envelope_mass*convective_envelope_radius*radius_arg)/(3.0*luminosity), 1.0/3.0);

        double f_convective = pow(P_tid/(2.0*tau_convective),2.0);
        f_convective = CV_min(1.0,f_convective);

        k_AM_div_T = (2.0/21.0)*(f_convective/tau_convective)*(convective_envelope_mass/mass);
        t_V = from_k_AM_div_T_to_t_V(k_AM_div_T,apsidal_motion_constant);
        
        return t_V;

    }
    else // degenerate damping -- 1984MNRAS.207..433C
    {
        k_AM_div_T = 2.564e-8 * gyration_radius*gyration_radius * pow((luminosity/CONST_L_SUN)/(mass/CONST_MSUN),5.0/7.0);

        t_V = from_k_AM_div_T_to_t_V(k_AM_div_T,apsidal_motion_constant);
        
        return t_V;
    }
}

double compute_EOM_equilibrium_tide_BO_full(ParticlesMap *particlesMap, int binary_index, int star_index, int companion_index, bool include_tidal_friction_terms, bool include_tidal_bulges_precession_terms, bool include_rotation_precession_terms, double minimum_eccentricity_for_tidal_precession, int tides_method)
/* Barker & Ogilvie (2009; http://adsabs.harvard.edu/abs/2009MNRAS.395.2268B) */

/* NOTE: in SecularMultiple, the h-vector is defined as the orbital angular momentum vector,
 * NOT the SPECIFIC orbital angular momentum vector. Compared to the notation used by Eggleton,
 * h_vec_SecularMultiple = mu*h_vec_Eggleton where mu = m*M/(m+M) is the reduced mass.
 * In particular, note the line `star->dspin_vec_dt[i] += -dh_vec_dt_star[i]/I;' */
{
   
    Particle *binary = (*particlesMap)[binary_index];
    Particle *star = (*particlesMap)[star_index];
    Particle *companion = (*particlesMap)[companion_index];
    
    if (star->is_binary == true)
    {
        #ifdef VERBOSE
        if (verbose_flag > 2)
        {
            printf("ODE_tides.cpp -- compute_EOM_equilibrium_tide_BO_full -- star is actually a binary -- not doing tides; binary_index %d star_index %d companion_index %d\n",binary_index,star_index,companion_index);
        }
        #endif

        return 0;
    }
    if (star->is_bound == false or companion->is_bound == false)
    {
        #ifdef VERBOSE
        if (verbose_flag > 2)
        {
            printf("ODE_tides.cpp -- compute_EOM_equilibrium_tide_BO_full -- tar/companion is not bound -- not doing tides; binary_index %d star_index %d companion_index %d\n",binary_index,star_index,companion_index);
        }
        #endif

        return 0;
    }

    if (star->stellar_type == 14) /* No tides for BHs */
    {
        return 0;
    }
    
    /* orbit quantities */
    double e = binary->e;
    double e_p2 = binary->e_p2;
    double a = binary->a;

    if (e >= 1.0 or a <= 0.0)
    {
        //printf("ODE_tides.cpp -- compute_EOM_equilibrium_tide_BO_full -- invalid e %g and/or a %g; skipping\n",e,a);
        //return 0;
    }
    
    double h = binary->h;
    double *e_vec = binary->e_vec;
    double *h_vec = binary->h_vec;
    double *e_vec_unit = binary->e_vec_unit;
    double *h_vec_unit = binary->h_vec_unit;
    double *q_vec_unit = binary->q_vec_unit;
    double j = binary->j;
    double j_p2 = binary->j_p2;
    double j_p3 = binary->j_p3;
    double j_p4 = binary->j_p4;
    double j_p8 = j_p4*j_p4;
    double j_p10 = j_p2*j_p8;  
    double j_p13 = j_p3*j_p10; 
    double j_p4_inv = 1.0/j_p4; 
    double j_p10_inv = 1.0/j_p10; 
    double j_p13_inv = 1.0/j_p13;
    double P_orb = compute_orbital_period_from_semimajor_axis(binary->mass,a);
    double n = 2.0*M_PI/P_orb; /* mean motion */

    /* stellar properties */
    double *spin_vec = star->spin_vec;
    double M = star->mass;
    double m = companion->mass;
    double R = star->radius;
    
    double k_AM;
    double I; // moment of intertia
    if (star->object_type == 1) // star 
    {
        I = compute_moment_of_inertia(star->stellar_type, M, star->core_mass, R, star->core_radius, star->sse_k2, star->sse_k3);
        k_AM = compute_apsidal_motion_constant(star);
    }
    else
    {
        I = star->gyration_radius*M*R*R;
        k_AM = star->apsidal_motion_constant;
    }

    double t_V = compute_t_V(star,companion,a);
    star->tides_viscous_time_scale = t_V;

    #ifdef VERBOSE
    if (verbose_flag > 2)
    {
        printf("ODE_tides.cpp -- compute_EOM_equilibrium_tide_BO_full -- binary_index %d star_index %d companion_index %d t_V %.15g rg %.15g kAM %.15g R %.15g I %.15g\n",binary_index,star_index,companion_index,t_V,star->gyration_radius,k_AM,star->radius,I);
    }
    #endif

    
    if (t_V!=t_V)
    {
        printf("ODE_tides.cpp -- FATAL ERROR: t_V is NaN\n");
        printf("binary_index %d star_index %d companion_index %d\n",binary_index,star_index,companion_index);
        printf("t_V %g \n",t_V);
        printf("st %d\n",star->stellar_type);
        printf("pr %d\n",star->tides_viscous_time_scale_prescription);
        printf("M %g\n",M);
        printf("star->convective_envelope_mass %g\n",star->convective_envelope_mass);
        printf("m %g\n",m);
        printf("a %.15g\n",a);
        printf("spin %g %g %g norm %g\n",spin_vec[0],spin_vec[1],spin_vec[2],star->spin_vec_norm);
        printf("e %.15g\n",e);
        printf("rp %g\n",a*(1.0-e));
        printf("R %g\n",R);
        printf("rg %g\n",star->gyration_radius);
        printf("star->convective_envelope_radius %g\n",star->convective_envelope_radius);
        printf("star->luminosity %g\n",star->luminosity);
        printf("k_AM %g\n",k_AM);
        print_system(particlesMap,0);
        //exit(-1);
        error_code = 22;
    }

    double R_div_a = R/a;
	double R_div_a_p5 = pow(R_div_a,5.0);
    double R_div_a_p8 = pow(R_div_a,8.0);

    double f = 2.0 * k_AM + 1.0;
    double t_f_inv = (9.0 / t_V) * R_div_a_p8 * (m * (m + M) / (M * M)) * f * f;

    double f_tides1 = f_tides1_function_BO(e_p2,j_p10_inv,j_p13_inv);
    double f_tides2 = f_tides2_function_BO(e_p2,j_p10_inv,j_p13_inv);
    double f_tides3 = f_tides3_function_BO(e_p2,j_p10_inv,j_p13_inv);
    double f_tides4 = f_tides4_function_BO(e_p2,j_p10_inv,j_p13_inv);
    double f_tides5 = f_tides5_function_BO(e_p2,j_p10_inv,j_p13_inv);

    double spin_vec_dot_e_vec = dot3(spin_vec,e_vec);
    double spin_vec_dot_h_vec = dot3(spin_vec,h_vec);

    double spin_vec_dot_e_vec_unit = dot3(spin_vec,e_vec_unit);
    double spin_vec_dot_h_vec_unit = dot3(spin_vec,h_vec_unit);
    double spin_vec_dot_q_vec_unit = dot3(spin_vec,q_vec_unit);
    double mu = m*M/(m+M);
    double C = m*k_AM*R_div_a_p5/(mu*n);
    
    check_number(spin_vec[0],"tides.cpp -- compute_EOM_equilibrium_tide_BO_full","spin_vec[0]", true);
    check_number(spin_vec[1],"tides.cpp -- compute_EOM_equilibrium_tide_BO_full","spin_vec[1]", true);
    check_number(spin_vec[2],"tides.cpp -- compute_EOM_equilibrium_tide_BO_full","spin_vec[2]", true);
    check_number(h_vec_unit[0],"tides.cpp -- compute_EOM_equilibrium_tide_BO_full","h_vec_unit[0]", true);
    check_number(h_vec_unit[1],"tides.cpp -- compute_EOM_equilibrium_tide_BO_full","h_vec_unit[1]", true);
    check_number(h_vec_unit[2],"tides.cpp -- compute_EOM_equilibrium_tide_BO_full","h_vec_unit[2]", true);
    check_number(e_vec_unit[0],"tides.cpp -- compute_EOM_equilibrium_tide_BO_full","e_vec_unit[0]", true);
    check_number(e_vec_unit[1],"tides.cpp -- compute_EOM_equilibrium_tide_BO_full","e_vec_unit[1]", true);
    check_number(e_vec_unit[2],"tides.cpp -- compute_EOM_equilibrium_tide_BO_full","e_vec_unit[2]", true);
    check_number(q_vec_unit[0],"tides.cpp -- compute_EOM_equilibrium_tide_BO_full","q_vec_unit[0]", true);
    check_number(q_vec_unit[1],"tides.cpp -- compute_EOM_equilibrium_tide_BO_full","q_vec_unit[1]", true);
    check_number(q_vec_unit[2],"tides.cpp -- compute_EOM_equilibrium_tide_BO_full","q_vec_unit[2]", true);
    check_number(spin_vec_dot_e_vec_unit,"tides.cpp -- compute_EOM_equilibrium_tide_BO_full","spin_vec_dot_e_vec_unit", true);
    check_number(spin_vec_dot_h_vec_unit,"tides.cpp -- compute_EOM_equilibrium_tide_BO_full","spin_vec_dot_h_vec_unit", true);
    check_number(spin_vec_dot_q_vec_unit,"tides.cpp -- compute_EOM_equilibrium_tide_BO_full","spin_vec_dot_q_vec_unit", true);
    check_number(C,"tides.cpp -- compute_EOM_equilibrium_tide_BO_full","C", true);
    check_number(j_p4_inv,"tides.cpp -- compute_EOM_equilibrium_tide_BO_full","j_p4_inv", true);

    double C_rot;
    double X_rot = 0.0;
    double Y_rot = 0.0;
    double Z_rot = 0.0;
    double Z_TB = 0.0;
    if (include_rotation_precession_terms == true)
    {
        if (tides_method == 1)
        {
            C_rot = C*j_p4_inv*spin_vec_dot_h_vec_unit;
            X_rot = -C_rot*spin_vec_dot_e_vec_unit;
            Y_rot = -C_rot*spin_vec_dot_q_vec_unit;
        }

        Z_rot = C*c_1div2*j_p4_inv*(2.0*spin_vec_dot_h_vec_unit*spin_vec_dot_h_vec_unit - spin_vec_dot_q_vec_unit*spin_vec_dot_q_vec_unit - spin_vec_dot_e_vec_unit*spin_vec_dot_e_vec_unit);
    }    
    if (include_tidal_bulges_precession_terms == true)
    {
        Z_TB = C*15.0*n*n*(mu/M)*f_tides2;
    }

    double X = X_rot;
    double Y = Y_rot;
    double Z = Z_rot + Z_TB;

    check_number(X_rot,"tides.cpp -- compute_EOM_equilibrium_tide_BO_full","X_rot", true);
    check_number(Y_rot,"tides.cpp -- compute_EOM_equilibrium_tide_BO_full","Y_rot", true);
    check_number(Z_rot,"tides.cpp -- compute_EOM_equilibrium_tide_BO_full","Z_rot", true);
    check_number(Z_TB,"tides.cpp -- compute_EOM_equilibrium_tide_BO_full","Z_TB", true);
    check_number(h,"tides.cpp -- compute_EOM_equilibrium_tide_BO_full","h", true);
    check_number(n,"tides.cpp -- compute_EOM_equilibrium_tide_BO_full","n", true);
    double dh_vec_dt_star_i;

    if (binary->exclude_for_secular_integration == true)// and binary->exclude_rotation_and_bulges_precession_in_case_of_isolated_binary == false)
    {
        X = Y = Z = 0.0;
    }
    //printf("X %g Y %g Z %g\n",X,Y,Z);
    for (int i=0; i<3; i++)
    {

        if (include_tidal_friction_terms == true)
        {
            dh_vec_dt_star_i = -t_f_inv*( (h/(2.0*n))*(spin_vec_dot_e_vec*f_tides5*e_vec[i] - spin_vec[i]*f_tides3) \
                + h_vec[i]*(f_tides4 - spin_vec_dot_h_vec*f_tides2/(2.0*n*h)) );
            binary->dh_vec_dt[i] += dh_vec_dt_star_i;

            binary->de_vec_dt[i] += -(t_f_inv/h)*( spin_vec_dot_e_vec*f_tides2*h_vec[i]/(2.0*n) \
                + 9.0*e_vec[i]*(f_tides1*h - c_11div18*spin_vec_dot_h_vec*f_tides2/n) );
            star->dspin_vec_dt[i] += -dh_vec_dt_star_i/I;

            #ifdef VERBOSE
            if (verbose_flag > 2)
            {
                printf("ODE_tides.cpp -- compute_EOM_equilibrium_tide_BO_full -- star %d bin %d t_f_inv %g binary->dh_vec_dt[i] %g binary->de_vec_dt[i] %g star->dspin_vec_dt[i] %g\n",star->index,binary->index,t_f_inv,binary->dh_vec_dt[i],binary->de_vec_dt[i],star->dspin_vec_dt[i]);
            }
            #endif
            
        }
        if (include_rotation_precession_terms == true || include_tidal_bulges_precession_terms == true)
        {
            if (e >= minimum_eccentricity_for_tidal_precession and (binary->parent != -1 or binary->exclude_rotation_and_bulges_precession_in_case_of_isolated_binary == false)) /* Do not include orbital precession if e is very small, or if the binary is isolated */
            {
                binary->de_vec_dt[i] += e*(Z*q_vec_unit[i] - Y*h_vec_unit[i]);

                dh_vec_dt_star_i = h*(-X*q_vec_unit[i] + Y*e_vec_unit[i]);
                binary->dh_vec_dt[i] += dh_vec_dt_star_i;
                star->dspin_vec_dt[i] += -dh_vec_dt_star_i/I;
            }
        }
        
        check_number(binary->dh_vec_dt[i],"tides.cpp -- compute_EOM_equilibrium_tide_BO_full","binary->dh_vec_dt[i]", true);
        check_number(binary->de_vec_dt[i],"tides.cpp -- compute_EOM_equilibrium_tide_BO_full","binary->de_vec_dt[i]", true);
        check_number(star->dspin_vec_dt[i],"tides.cpp -- compute_EOM_equilibrium_tide_BO_full","star->dspin_vec_dt[i]", true);
    }

    #ifdef VERBOSE
    if (verbose_flag > 1)
    {
        for (int i=0; i<3; i++)
        {
            printf("ODE_tides.cpp -- compute_EOM_equilibrium_tide_BO_full star %d comp %d i %d spin %.15f de/dt %.15g dh/dt %.15g dspin/dt %.15g a %.15g e %.15g\n",star->index,companion->index,i,spin_vec[i],binary->de_vec_dt[i],binary->dh_vec_dt[i],star->dspin_vec_dt[i],a,e);
        }
    }
    #endif

    return 0;
}


double f_tides1_function_BO(double e_p2, double j_p10_inv, double j_p13_inv)
{
    return j_p13_inv*(1.0 + e_p2*(c_15div4 + e_p2*(c_15div8 + e_p2*c_5div64)));
}
double f_tides2_function_BO(double e_p2, double j_p10_inv, double j_p13_inv)
{
    return j_p10_inv*(1.0 + e_p2*(c_3div2 + e_p2*c_1div8));
}
double f_tides3_function_BO(double e_p2, double j_p10_inv, double j_p13_inv)
{
    return j_p10_inv*(1.0 + e_p2*(c_9div2 + e_p2*c_5div8));
}
double f_tides4_function_BO(double e_p2, double j_p10_inv, double j_p13_inv)
{
    return j_p13_inv*(1.0 + e_p2*(c_15div2 + e_p2*(c_45div8 + e_p2*c_5div16)));
}
double f_tides5_function_BO(double e_p2, double j_p10_inv, double j_p13_inv)
{
    return j_p10_inv*(3.0 + c_1div2*e_p2);
}  

void compute_estimated_tidal_apsidal_motion_timescales(double M, double m, double R, double k_AM, double *spin_vec, double a, double e, double *e_vec, double *h_vec, double *t_apsidal_motion)
{
    /* This function estimates the timescale for tides to advance the apsidal line.
     * Is supposed to be called only occasionally, so not much attempt made at optimisation. 
     * M: subject; m: companion */
     
    double mu = m*M/(m+M);
    double P_orb = compute_orbital_period_from_semimajor_axis(m+M,a);
    double n = 2.0*M_PI/P_orb; /* mean motion */
    double C = m*k_AM*pow(R/a,5.0)/(mu*n);
    
    double e_p2 = e*e;
    double j = sqrt(1.0-e_p2);
    double j_p4_inv = pow(j,-4.0);
    double j_p10_inv = pow(j,-10.0);
    double j_p13_inv = pow(j,-13.0);
    
    double f_tides2 = f_tides2_function_BO(e_p2,j_p10_inv,j_p13_inv);

    double spin_vec_dot_e_vec = dot3(spin_vec,e_vec);
    double spin_vec_dot_h_vec = dot3(spin_vec,h_vec);

    double h = norm3(h_vec);
    double q_vec_unit[3];
    cross3(h_vec,e_vec,q_vec_unit); /* not yet normalised */
    double q_norm = norm3(q_vec_unit);
    for (int i=0; i<3; i++)
    {
        q_vec_unit[i] = q_vec_unit[i]/q_norm;
    }

    double spin_vec_dot_e_vec_unit = dot3(spin_vec,e_vec)/e;
    double spin_vec_dot_h_vec_unit = dot3(spin_vec,h_vec)/h;
    double spin_vec_dot_q_vec_unit = dot3(spin_vec,q_vec_unit);
    
    double Z_rot = C*c_1div2*j_p4_inv*(2.0*spin_vec_dot_h_vec_unit*spin_vec_dot_h_vec_unit - spin_vec_dot_q_vec_unit*spin_vec_dot_q_vec_unit - spin_vec_dot_e_vec_unit*spin_vec_dot_e_vec_unit);
    double Z_TB = C*15.0*n*n*(mu/M)*f_tides2;

    double Z = fabs(Z_rot + Z_TB);
    if (Z <= epsilon)
    {
        Z = epsilon;
    }
    *t_apsidal_motion = 1.0/Z;
}


double compute_EOM_equilibrium_tide(ParticlesMap *particlesMap, int binary_index, int star_index, int companion_index, bool include_tidal_friction_terms, bool include_tidal_bulges_precession_terms, bool include_rotation_precession_terms, double minimum_eccentricity_for_tidal_precession)

/* Equilibrium tide in vector form adopted from Eggleton & Kisseleva 1998 */
 
/* NOTE: in SecularMultiple, the h-vector is defined as the orbital angular momentum vector,
 * NOT the SPECIFIC orbital angular momentum vector. Compared to the notation used by Eggleton,
 * h_vec_SecularMultiple = mu*h_vec_Eggleton where mu = m*M/(m+M) is the reduced mass.
 * In particular, note the line `star->dspin_vec_dt[i] += -dh_vec_dt_star[i]/I;' */

{

    #ifdef VERBOSE
    if (verbose_flag > 2)
    {
        printf("ODE_tides.cpp -- compute_EOM_equilibrium_tide (Eggleton & Kisseleva 1998) -- bin %d star %d companion %d\n",binary_index,star_index,companion_index);
    }
    #endif

    Particle *binary = (*particlesMap)[binary_index];
    Particle *star = (*particlesMap)[star_index];
    Particle *companion = (*particlesMap)[companion_index];
    
    /* orbit quantities */
    double e = binary->e;
    double e_p2 = binary->e_p2;
    double a = binary->a;
    double h = binary->h;
    double *e_vec = binary->e_vec;
    double *h_vec = binary->h_vec;
    double *e_vec_unit = binary->e_vec_unit;
    double *h_vec_unit = binary->h_vec_unit;
    double *q_vec_unit = binary->q_vec_unit;    
    double j = binary->j;
    double j_p2 = binary->j_p2;
    double j_p3 = binary->j_p3;
    double j_p4 = binary->j_p4;
    double j_p8 = j_p4*j_p4;
    double j_p10 = j_p2*j_p8;  
    double j_p13 = j_p3*j_p10; 
    double j_p4_inv = 1.0/j_p4;
    double j_p10_inv = 1.0/j_p10; 
    double j_p13_inv = 1.0/j_p13;
    double P_orb = compute_orbital_period_from_semimajor_axis(binary->mass,binary->a);
    double n = 2.0*M_PI/P_orb; /* mean motion */
    
    /* stellar properties */
    double *spin_vec = star->spin_vec;
    double M = star->mass;
    double m = companion->mass;
    double mu = m*M/(m+M);
    double R = star->radius;
    double k_AM = star->apsidal_motion_constant;
    double t_V = compute_t_V(star,companion,a);
    star->tides_viscous_time_scale = t_V;
    double rg = star->gyration_radius;
    double I;
    if (star->object_type == 1) // star 
    {
        I = compute_moment_of_inertia(star->stellar_type, M, star->core_mass, R, star->core_radius, star->sse_k2, star->sse_k3);
        k_AM = compute_apsidal_motion_constant(star);
    }
    else
    {
        I = star->gyration_radius*M*R*R;
        k_AM = star->apsidal_motion_constant;
    }

    
    double R_div_a = R/a;
	double R_div_a_p5 = pow(R_div_a,5.0);
    double R_div_a_p8 = pow(R_div_a,8.0);
    double t_f_inv = (9.0/t_V)*R_div_a_p8*((M+m)*m/(M*M))*(1.0 + 2.0*k_AM)*(1.0 + 2.0*k_AM);

    double spin_vec_dot_e_vec_unit = dot3(spin_vec,e_vec_unit);
    double spin_vec_dot_h_vec_unit = dot3(spin_vec,h_vec_unit);
    double spin_vec_dot_q_vec_unit = dot3(spin_vec,q_vec_unit);
        
    double V,W,X,Y,Z;
    VWXYZ_tides_function(
        include_tidal_friction_terms, include_tidal_bulges_precession_terms, include_rotation_precession_terms, minimum_eccentricity_for_tidal_precession, \
        t_f_inv,k_AM, \
        m, M, mu, n, R_div_a_p5, \
        e, e_p2, j_p4_inv, j_p10_inv, j_p13_inv, \
        spin_vec_dot_h_vec_unit, spin_vec_dot_e_vec_unit, spin_vec_dot_q_vec_unit, \
        &V, &W, &X, &Y, &Z);

    #ifdef VERBOSE
    if (verbose_flag > 1)
    {
        printf("ODE_tides.cpp -- compute_EOM_equilibrium_tide (Eggleton & Kisseleva 1998) -- t_f_inv V %g W %g X %g Y %g Z %g\n",t_f_inv,V,W,X,Y,Z);
    }
    #endif

    double dh_vec_dt_star[3];
    
    for (int i=0; i<3; i++)
    {
        dh_vec_dt_star[i] = h*( Y*e_vec_unit[i] - W*h_vec_unit[i] - X*q_vec_unit[i] );
        binary->dh_vec_dt[i] += dh_vec_dt_star[i];

        binary->de_vec_dt[i] += e*( -V*e_vec_unit[i] - Y*h_vec_unit[i] + Z*q_vec_unit[i] );
        
        star->dspin_vec_dt[i] += -dh_vec_dt_star[i]/I; /* conservation of total angular momentum (orbit+spin) */
    }

    #ifdef VERBOSE
    if (verbose_flag > -1)
    {
        for (int i=0; i<3; i++)
        {
            printf("ODE_tides.cpp -- compute_EOM_equilibrium_tide (Eggleton & Kisseleva 1998) star %d comp %d i %d spin %.15f de/dt %.15g dh/dt %.15g dspin/dt %.15g a %.15g e %.15g\n",star->index,companion->index,i,spin_vec[i],binary->de_vec_dt[i],binary->dh_vec_dt[i],star->dspin_vec_dt[i],a,e);
        }
    }
    #endif

    return 0;
}

double f_tides1_function(double e_p2, double j_p10_inv, double j_p13_inv)
{
    return j_p13_inv*(1.0 + e_p2*(c_15div4 + e_p2*(c_15div8 + e_p2*c_5div64)));
}
double f_tides2_function(double e_p2, double j_p10_inv, double j_p13_inv)
{
    return j_p10_inv*(1.0 + e_p2*(c_3div2 + e_p2*c_1div8));
}
double f_tides3_function(double e_p2, double j_p10_inv, double j_p13_inv)
{
    return j_p13_inv*(1.0 + e_p2*(c_15div2 + e_p2*(c_45div8 + e_p2*c_5div16)));
}
double f_tides4_function(double e_p2, double j_p10_inv, double j_p13_inv)
{
    return j_p10_inv*(1.0 + e_p2*(3.0 + e_p2*c_3div8));
}
double f_tides5_function(double e_p2, double j_p10_inv, double j_p13_inv)
{
    return j_p10_inv*(1.0 + e_p2*(c_9div2 + e_p2*c_5div8));
}


int VWXYZ_tides_function
(
    bool include_tidal_friction_terms, bool include_tidal_bulges_precession_terms, bool include_rotation_precession_terms, double minimum_eccentricity_for_tidal_precession, \
    double t_f_inv, double k_AM, \
    double m, double M, double mu, double n, double R_div_a_p5, \
    double e, double e_p2, double j_p4_inv, double j_p10_inv,double j_p13_inv, \
    double spin_vec_dot_h_vec_unit, double spin_vec_dot_e_vec_unit,double spin_vec_dot_q_vec_unit, \
    double* V, double* W, double* X, double* Y, double* Z
)
{
   
    *V = 0.0;
    *W = 0.0;
    *X = 0.0;
    *Y = 0.0;
    *Z = 0.0;

    double f2 = f_tides2_function(e_p2,j_p10_inv,j_p13_inv); /* needed for both pure tidal dissipation and pure tidal bulges terms */

//    if (e < minimum_eccentricity_for_tidal_precession)
//    {
//        return 0;
//    }


    if (include_tidal_friction_terms == true)
    {
        double f1 = f_tides1_function(e_p2,j_p10_inv,j_p13_inv);
        double f3 = f_tides3_function(e_p2,j_p10_inv,j_p13_inv);
        double f4 = f_tides4_function(e_p2,j_p10_inv,j_p13_inv);
        double f5 = f_tides5_function(e_p2,j_p10_inv,j_p13_inv);
    
        *V += 9.0*t_f_inv*(f1    - c_11div18*(spin_vec_dot_h_vec_unit/n)*f2);
        *W += t_f_inv*(f3        - (spin_vec_dot_h_vec_unit/n)*f4);
        *X +=  -t_f_inv*spin_vec_dot_q_vec_unit*f5/(2.0*n);
        *Y += t_f_inv*spin_vec_dot_e_vec_unit*f2/(2.0*n);
//	printf("TF X %g Y %g V %g W  %g\n",*X,*Y,*Z,*V,*W);

    }

    if (e < minimum_eccentricity_for_tidal_precession)
    {
        return 0;
    }


    if ((include_tidal_bulges_precession_terms == true) || (include_rotation_precession_terms) == true)
    {
        double C = m*k_AM*R_div_a_p5/(mu*n);
    
        if (include_tidal_bulges_precession_terms == true)
        {
            *Z += C*15.0*n*n*(mu/M)*f2;
//            printf("include_tidal_bulges_precession_terms Z %g\n",*Z);
        }    
        if (include_rotation_precession_terms == true)
        {
            double C_XY = -C*spin_vec_dot_h_vec_unit*j_p4_inv;
            *X += C_XY*spin_vec_dot_e_vec_unit;
            *Y += C_XY*spin_vec_dot_q_vec_unit;            


//		printf("ROT X %g Y %g Z %g\n",C_XY*spin_vec_dot_e_vec_unit,C_XY*spin_vec_dot_q_vec_unit,C*c_1div2*j_p4_inv*(2.0*spin_vec_dot_h_vec_unit*spin_vec_dot_h_vec_unit - spin_vec_dot_q_vec_unit*spin_vec_dot_q_vec_unit - spin_vec_dot_e_vec_unit*spin_vec_dot_e_vec_unit));
//        
//            printf("Z1 %g\n",*Z);
            *Z += C*c_1div2*j_p4_inv*(2.0*spin_vec_dot_h_vec_unit*spin_vec_dot_h_vec_unit - spin_vec_dot_q_vec_unit*spin_vec_dot_q_vec_unit - spin_vec_dot_e_vec_unit*spin_vec_dot_e_vec_unit);
//            printf("Z2 %g\n",*Z);
//            printf("O %g %g %g\n",spin_vec_dot_h_vec_unit,spin_vec_dot_e_vec_unit,spin_vec_dot_q_vec_unit);
//            printf("include_rotation_precession_terms Z %g\n",*Z);
        }
    }
    
    return 0;
}

}
