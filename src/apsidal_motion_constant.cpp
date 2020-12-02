/* MSE */

#include "evolve.h"
#include "apsidal_motion_constant.h"

double compute_apsidal_motion_constant(Particle *star)
{
	static const double d_mass_list[] = {AMC_DATA_M_LIST};

	static const double d_HeMS[] = {AMC_DATA_HEMS};
	static const double d_WD[] = {AMC_DATA_WD};

	double val = 0.0;	/* Return value of the apsidal motion constant */
	
    int stellar_type = star->stellar_type;
    double mass = star->mass;

    if (log10(mass) > 2.10)
    {
        
        #ifdef VERBOSE
        if (verbose_flag > 0)
        {
            printf("apsidal_motion_constant.cpp -- WARNING: input mass %g exceeds maximum; setting effective mass to %g; index %d kw %d m %g age %g\n",mass,pow(10.0,2.10),star->index,stellar_type,star->mass,star->age);
        }
        #endif
        
        mass = pow(10.0,2.10);
    }

    if (stellar_type == 0)
    {
        val = 0.14327923; /* Approximate the star as being fully convective such that n = 3/2 */
    }
    
    else if ((stellar_type >= 1) and (stellar_type <= 6))
    {

        /*-----------------------------------------------------------------------------------------------------------------------------------*/
        /* Main sequence/Hertzsprung Gap/Giant/HeMS/AGB: determine from tabulated values from detailed models by Claret (2004), for Z = 0.02 */
        /*--------------------------------------------------------------------------------------------------------------------------*/

        /* Determine which two masses to use for which AMC is tabulated */
        double m_upper = 0.0;
        double m_lower = 0.0;
        int i = 0;

        while (m_upper < mass)
        {
            m_lower = m_upper;
            m_upper = pow(10.0,d_mass_list[i]);
            i++;
        }
        double k_AM_lower = AMC_data_function(d_mass_list[i-2],star);
        double k_AM_upper = AMC_data_function(d_mass_list[i-1],star);

        val = (k_AM_upper-k_AM_lower)/(m_upper-m_lower)*(mass-m_upper) + k_AM_upper;
        //printf("val %g k_AM_lower %g k_AM_upper %g %dml %g dmu %g \n",val,k_AM_lower,k_AM_upper,d_mass_list[i-2],d_mass_list[i-1]);
    }

    else if ((stellar_type >= 7) and (stellar_type <= 9))
    {

        /*--------------------------------------------------------------------------------------------------------------------------*/
        /* Stripped helium star: determine from fitted function from data of Vila (1977; https://ui.adsabs.harvard.edu/abs/1977ApJ...213..464V) (no age information used, just mass)
        /*--------------------------------------------------------------------------------------------------------------------------*/
    
        val = d_HeMS[0] + d_HeMS[1]*atan(d_HeMS[2]*pow(mass - d_HeMS[3],4.0));
    }

    else if ((stellar_type >= 10) and (stellar_type <= 12))
    {

        /*--------------------------------------------------------------------------------------------------------------------------*/
        /* White dwarf: determine from fitted function from data of Vila (1977; https://ui.adsabs.harvard.edu/abs/1977ApJ...213..464V) (no age information used, just mass)
        /*--------------------------------------------------------------------------------------------------------------------------*/
    
        val = d_WD[0] + d_WD[1]*(mass - d_WD[2]) + d_WD[3]*pow(mass - d_WD[4],2.0);
    }
    
    else if (stellar_type == 13) 
    {

        /* Take canonical constant value for n=1 polytrope (https://ui.adsabs.harvard.edu/abs/1955MNRAS.115..101B/abstract) */
        val = 0.25990728; // n=1
    }
    
    else if (stellar_type == 14)
    {

        /*--------------------------------------------------------------------------------------------------------------------------*/
        /* Black hole: see Binnington, 2009: k_AM = 0 for a non-rotating black hole
        /*--------------------------------------------------------------------------------------------------------------------------*/

        val = 0.0;
    }
    
    else if (stellar_type == 15) /* Massless remnant: k_AM = 0 */
    {
        val = 0.0;
    }
    
    else
    {
        /*--------------------------------------------------------------------------------------------------------------------------*/
        /* No clue what to do otherwise, so, use n=3/2 polytrope value as "best guess/typical value" (https://ui.adsabs.harvard.edu/abs/1955MNRAS.115..101B/abstract) */
        /*--------------------------------------------------------------------------------------------------------------------------*/

        val = 0.14327923;
    }
	
    check_number(val, "apsidal_motion_constant.cpp", "k_AM", true);
    
    if (val<0)
    {
        printf("apsidal_motion_constant.cpp -- ERROR: AMC is %g index %d kw %d m %g age %g\n",val,star->index,stellar_type,star->mass,star->age);
        exit(0);
    }

    return val;
}

void limit_tau(double *tau)
{
    if (*tau >= 1.0)
    {
        #ifdef VERBOSE
        if (verbose_flag > 1)
        {
            printf("apsidal_motion_constant.cpp -- limit_tau -- limiting tau from %g to 1.0\n",*tau);
        }
        #endif
        
        *tau = 1.0;
    }
}

double AMC_data_function(double log_m, Particle *star)
{
	int stellar_type = star->stellar_type;
	double age = yr_to_Myr*star->age; /* age in Myr since the timescales below are in Myr */
	double tau = 0.0;						/* Fractional age at specific evolutionary stage */

    double tms,tn;
    double *GB,*timescales,*lums;
    GB = new double[10];
    timescales = new double[20];
    lums = new double[10];  
    
    int stellar_type_old = stellar_type;
    star_(&stellar_type, &star->sse_initial_mass, &star->mass, &tms, &tn, timescales, lums, GB, star->zpars);

    if (stellar_type_old != stellar_type)
    {
        #ifdef VERBOSE
        if (verbose_flag > 0)
        {
            printf("apsidal_motion_constant.cpp -- WARNING: stellar type has changed from %d to %d!\n",stellar_type_old,stellar_type);
        }
        #endif
    }

    int T_BGB = 0;
    int T_HE_IGNITION = 1;
    int T_HE_BURNING = 2;

    #ifdef VERBOSE
    if (verbose_flag > 1)
    {
        printf("apsidal_motion_constant.cpp -- AMC_data_function -- input log_m %g age %g tms %g\n",log_m, age, tms);
    }
    #endif

	if (log_m == -0.0969)	/* 0.80 M_SUN */
	{
		if (stellar_type == 1)	/* MAIN SEQUENCE */
		{
			tau = age/tms;
            limit_tau(&tau);
			return d_m_m1[0] + d_m_m1[1]*tau;		
		}
		else if (stellar_type == 2)	/* HERTZSPRUNG GAP */
		{
			tau = (age - tms)/(timescales[T_BGB] - tms);
            limit_tau(&tau);
			if (tau < d_m_m1[4])
			{
				return d_m_m1[2] + d_m_m1[3]*tau;
			}
			if ((tau > d_m_m1[4]) and (tau < d_m_m1[7]))
			{
				return d_m_m1[5] + d_m_m1[6]*tau;
			}
			else
			{
				return d_m_m1[8] + d_m_m1[9]*tau;
			}
		}
		else if (stellar_type == 3) /* RED GIANT */
		{
			tau = (age - timescales[T_BGB])/(timescales[T_HE_IGNITION] - timescales[T_BGB]);
            limit_tau(&tau);
			return d_m_m1[10] + d_m_m1[11]*tau;
		}
		else if (stellar_type == 4) /* He MS */
		{
			return 0.04;
		}
		else if (stellar_type == 5 || stellar_type == 6) /* EAGB/TPAGB */
		{
			return 0.04;
		}
	}

	else if (log_m == 0) /* 1.00 M_SUN */
	{
		if (stellar_type == 1)	/* MAIN SEQUENCE */
		{
			tau = age/ tms;
            limit_tau(&tau);
//			printf("tau MS %g age %g tms %g\n",tau,age,*tms);
			return d_m_0000[0] + d_m_0000[1]*tau;		
		}
		else if (stellar_type == 2)	/* HERTZSPRUNG GAP */
		{
			tau = (age - tms)/(timescales[T_BGB] - tms);
            limit_tau(&tau);
//			printf("tau HG %g tms %g age %g T_BGB %g\n",tau,*tms,age,timescales[T_BGB]);
			if (tau < d_m_0000[4])
			{
				return d_m_0000[2] + d_m_0000[3]*tau;
			}
			if ((tau > d_m_0000[4]) and (tau < d_m_0000[7]))
			{
				return d_m_0000[5] + d_m_0000[6]*tau;
			}
			else
			{
				return d_m_0000[8] + d_m_0000[9]*tau;
			}
		}
		else if (stellar_type == 3) /* RED GIANT */
		{
			tau = (age - timescales[T_BGB])/(timescales[T_HE_IGNITION] - timescales[T_BGB]);
            limit_tau(&tau);
//			printf("tau RG %g age %g\n",tau,age);
			if (tau < d_m_0000[12])
			{
				return d_m_0000[10] + d_m_0000[11]*tau;
			}
			else
			{
				return d_m_0000[13] + d_m_0000[14]*tau;
			}
		}
		else if (stellar_type == 4) /* He MS */
		{
			return 0.04;
		}
		else if (stellar_type == 5 || stellar_type == 6) /* EAGB/TPAGB */
		{
			return 0.04;
		}
	}
	else if (log_m == 0.05)	/* 1.12 M_MSUN */
	{
		if (stellar_type == 1)	/* MAIN SEQUENCE */
		{
			tau = age/tms;
            limit_tau(&tau);
//			printf("tau MS %g age %g tms %g\n",tau,age,*tms);
			return d_m_0050[0] + d_m_0050[1]*tau;		
		}
		else if (stellar_type == 2)	/* HERTZSPRUNG GAP */
		{
			tau = (age - tms)/(timescales[T_BGB] - tms);
            limit_tau(&tau);
//			printf("tau HG %g tms %g age %g T_BGB %g\n",tau,*tms,age,timescales[T_BGB]);
			if (tau < d_m_0050[4])
			{
				return d_m_0050[2] + d_m_0050[3]*tau;
			}
			if ((tau > d_m_0050[4]) and (tau < d_m_0050[7]))
			{
				return d_m_0050[5] + d_m_0050[6]*tau;
			}
			else
			{
				return d_m_0050[8] + d_m_0050[9]*tau;
			}
		}
		else if (stellar_type == 3) /* RED GIANT */
		{
			tau = (age-timescales[T_BGB])/(timescales[T_HE_IGNITION]-timescales[T_BGB]);	/* See hrdiag & hrdiag_RG */
            limit_tau(&tau);
//			printf("tau RG %g age %g\n",tau,age);
			if (tau < d_m_0050[12])
			{
				return d_m_0050[10] + d_m_0050[11]*tau;
			}
			else
			{
				return d_m_0050[13] + d_m_0050[14]*tau;
			}
		}
		else if (stellar_type == 4) /* He MS */
		{
			return 0.04;
		}
		else if (stellar_type == 5 || stellar_type == 6) /* EAGB/TPAGB */
		{
			return 0.04;
		}
	}
	else if (log_m == 0.075)	/* 1.19 M_SUN */
	{
		if (stellar_type == 1)	/* MAIN SEQUENCE */
		{
			tau = age/tms;
            limit_tau(&tau);
//			printf("tau MS %g age %g tms %g\n",tau,age,*tms);
			return d_m_0075[0] + d_m_0075[1]*tau;		
		}
		else if (stellar_type == 2)	/* HERTZSPRUNG GAP */
		{
			tau = (age - tms)/(timescales[T_BGB] - tms);
            limit_tau(&tau);
//			printf("tau HG %g tms %g age %g T_BGB %g\n",tau,*tms,age,timescales[T_BGB]);
			if (tau < d_m_0075[4])
			{
				return d_m_0075[2] + d_m_0075[3]*tau;
			}
			if ((tau > d_m_0075[4]) and (tau < d_m_0075[7]))
			{
				return d_m_0075[5] + d_m_0075[6]*tau;
			}
			else
			{
				return d_m_0075[8] + d_m_0075[9]*tau;
			}
		}
		else if (stellar_type == 3) /* RED GIANT */
		{
			tau = (age-timescales[T_BGB])/(timescales[T_HE_IGNITION]-timescales[T_BGB]);	/* See hrdiag & hrdiag_RG */
            limit_tau(&tau);
//			printf("tau RG %g age %g\n",tau,age);
			if (tau < d_m_0075[12])
			{
				return d_m_0075[10] + d_m_0075[11]*tau;
			}
			else
			{
				return d_m_0075[13] + d_m_0075[14]*tau;
			}
		}
		else if (stellar_type == 4) /* He MS */
		{
			return 0.04;
		}
		else if (stellar_type == 5 || stellar_type == 6) /* EAGB/TPAGB */
		{
			return 0.04;
		}
	}
	else if (log_m == 0.10)	/* 1.26 M_SUN */
	{
		if (stellar_type == 1)	/* MAIN SEQUENCE */
		{
			tau = age/tms;
            limit_tau(&tau);
//			printf("tau MS %g age %g tms %g\n",tau,age,*tms);
			return d_m_0100[0] + d_m_0100[1]*tau;		
		}
		else if (stellar_type == 2)	/* HERTZSPRUNG GAP */
		{
			tau = (age - tms)/(timescales[T_BGB] - tms);
            limit_tau(&tau);
//			printf("tau HG %g tms %g age %g T_BGB %g\n",tau,*tms,age,timescales[T_BGB]);
			if (tau < d_m_0100[4])
			{
				return d_m_0100[2] + d_m_0100[3]*tau;
			}
			if ((tau > d_m_0100[4]) and (tau < d_m_0100[7]))
			{
				return d_m_0100[5] + d_m_0100[6]*tau;
			}
			if ((tau > d_m_0100[7]) and (tau < d_m_0100[10]))
			{
				return d_m_0100[8] + d_m_0100[9]*tau;
			}
			if ((tau > d_m_0100[10]) and (tau < d_m_0100[13]))
			{
				return d_m_0100[11] + d_m_0100[12]*tau;
			}
			else
			{
				return d_m_0100[14] + d_m_0100[15]*tau;
			}
		}
		else if (stellar_type == 3) /* RED GIANT */
		{
			tau = (age-timescales[T_BGB])/(timescales[T_HE_IGNITION]-timescales[T_BGB]);	/* See hrdiag & hrdiag_RG */
            limit_tau(&tau);
//			printf("tau RG %g age %g\n",tau,age);
			if (tau < d_m_0100[18])
			{
				return d_m_0100[16] + d_m_0100[17]*tau;
			}
			else
			{
				return d_m_0100[19] + d_m_0100[20]*tau;
			}
		}
		else if (stellar_type == 4) /* He MS */
		{
			return 0.04;
		}
		else if (stellar_type == 5 || stellar_type == 6) /* EAGB/TPAGB */
		{
			return 0.04;
		}
	}
	else if (log_m == 0.125)	/* 1.33 M_SUN */
	{
		if (stellar_type == 1)	/* MAIN SEQUENCE */
		{
			tau = age/tms;
            limit_tau(&tau);
//			printf("tau MS %g age %g tms %g\n",tau,age,*tms);
			return d_m_0125[0] + d_m_0125[1]*tau;
		}
		else if (stellar_type == 2)	/* HERTZSPRUNG GAP */
		{
			tau = (age - tms)/(timescales[T_BGB] - tms);
            limit_tau(&tau);
//			printf("tau HG %g tms %g age %g T_BGB %g\n",tau,*tms,age,timescales[T_BGB]);
			if (tau < d_m_0125[4])
			{
				return d_m_0125[2] + d_m_0125[3]*tau;
			}
			if ((tau > d_m_0125[4]) and (tau < d_m_0125[7]))
			{
				return d_m_0125[5] + d_m_0125[6]*tau;
			}
			if ((tau > d_m_0125[7]) and (tau < d_m_0125[10]))
			{
				return d_m_0125[8] + d_m_0125[9]*tau;
			}
			if ((tau > d_m_0125[10]) and (tau < d_m_0125[13]))
			{
				return d_m_0125[11] + d_m_0125[12]*tau;
			}
			else
			{
				return d_m_0125[14] + d_m_0125[15]*tau;
			}
		}
		else if (stellar_type == 3) /* RED GIANT */
		{
			tau = (age-timescales[T_BGB])/(timescales[T_HE_IGNITION]-timescales[T_BGB]);	/* See hrdiag & hrdiag_RG */
            limit_tau(&tau);
//			printf("tau RG %g age %g\n",tau,age);
			if (tau < d_m_0125[18])
			{
				return d_m_0125[16] + d_m_0125[17]*tau;
			}
			else
			{
				return d_m_0125[19] + d_m_0125[20]*tau;
			}
		}
		else if (stellar_type == 4) /* He MS */
		{
			return 0.04;
		}
		else if (stellar_type == 5 || stellar_type == 6) /* EAGB/TPAGB */
		{
			return 0.04;
		}
	}
	else if (log_m == 0.150)	/* 1.41 M_SUN */
	{
		if (stellar_type == 1)	/* MAIN SEQUENCE */
		{
			tau = age/tms;
            limit_tau(&tau);
//			printf("tau MS %g age %g tms %g\n",tau,age,*tms);
			if (tau < d_m_0150[2])
			{
				return d_m_0150[0] + d_m_0150[1]*tau;
			}
			else
			{
				return d_m_0150[3] + d_m_0150[4]*tau;
			}
		}
		else if (stellar_type == 2)	/* HERTZSPRUNG GAP */
		{
			tau = (age - tms)/(timescales[T_BGB] - tms);
            limit_tau(&tau);
//			printf("tau HG %g tms %g age %g T_BGB %g\n",tau,*tms,age,timescales[T_BGB]);
			if (tau < d_m_0150[7])
			{
				return d_m_0150[5] + d_m_0150[6]*tau;
			}
			if ((tau > d_m_0150[7]) and (tau < d_m_0150[10]))
			{
				return d_m_0150[8] + d_m_0150[9]*tau;
			}
			if ((tau > d_m_0150[10]) and (tau < d_m_0150[13]))
			{
				return d_m_0150[11] + d_m_0150[12]*tau;
			}
			if ((tau > d_m_0150[13]) and (tau < d_m_0150[16]))
			{
				return d_m_0150[14] + d_m_0150[15]*tau;
			}
			else
			{
				return d_m_0150[17] + d_m_0150[18]*tau;
			}
		}
		else if (stellar_type == 3) /* RED GIANT */
		{
			tau = (age-timescales[T_BGB])/(timescales[T_HE_IGNITION]-timescales[T_BGB]);	/* See hrdiag & hrdiag_RG */
            limit_tau(&tau);
//			printf("tau RG %g age %g\n",tau,age);
			if (tau < d_m_0150[21])
			{
				return d_m_0150[19] + d_m_0150[20]*tau;
			}
			else
			{
				return d_m_0150[22] + d_m_0150[23]*tau;
			}
		}
		else if (stellar_type == 4) /* He MS */
		{
			return 0.04;
		}
		else if (stellar_type == 5 || stellar_type == 6) /* EAGB/TPAGB */
		{
			return 0.04;
		}
	}
	else if (log_m == 0.173)	/* 1.49 M_SUN */
	{
		if (stellar_type == 1)	/* MAIN SEQUENCE */
		{
			tau = age/tms;
            limit_tau(&tau);
//			printf("tau MS %g age %g tms %g\n",tau,age,*tms);
			if (tau < d_m_0173[2])
			{
				return d_m_0173[0] + d_m_0173[1]*tau;
			}
			else
			{
				return d_m_0173[3] + d_m_0173[4]*tau;
			}
		}
		else if (stellar_type == 2)	/* HERTZSPRUNG GAP */
		{
			tau = (age - tms)/(timescales[T_BGB] - tms);
            limit_tau(&tau);
//			printf("tau HG %g tms %g age %g T_BGB %g\n",tau,*tms,age,timescales[T_BGB]);
			if (tau < d_m_0173[7])
			{
				return d_m_0173[5] + d_m_0173[6]*tau;
			}
			if ((tau > d_m_0173[7]) and (tau < d_m_0173[10]))
			{
				return d_m_0173[8] + d_m_0173[9]*tau;
			}
			if ((tau > d_m_0173[10]) and (tau < d_m_0173[13]))
			{
				return d_m_0173[11] + d_m_0173[12]*tau;
			}
			else
			{
				return d_m_0173[14] + d_m_0173[15]*tau;
			}
		}
		else if (stellar_type == 3) /* RED GIANT */
		{
			tau = (age-timescales[T_BGB])/(timescales[T_HE_IGNITION]-timescales[T_BGB]);	/* See hrdiag & hrdiag_RG */
            limit_tau(&tau);
//			printf("tau RG %g age %g\n",tau,age);
			if (tau < d_m_0173[18])
			{
				return d_m_0173[16] + d_m_0173[17]*tau;
			}
			else
			{
				return d_m_0173[19] + d_m_0173[20]*tau;
			}
		}
		else if (stellar_type == 4) /* He MS */
		{
			return 0.04;
		}
		else if (stellar_type == 5 || stellar_type == 6) /* EAGB/TPAGB */
		{
			return 0.04;
		}
	}
	else if (log_m == 0.200)	/* 1.58 M_SUN */
	{
		if (stellar_type == 1)	/* MAIN SEQUENCE */
		{
			tau = age/tms;
            limit_tau(&tau);
//			printf("tau MS %g age %g tms %g\n",tau,age,*tms);
			if (tau < d_m_0200[2])
			{
				return d_m_0200[0] + d_m_0200[1]*tau;
			}
			else
			{
				return d_m_0200[3] + d_m_0200[4]*tau;
			}
		}
		else if (stellar_type == 2)	/* HERTZSPRUNG GAP */
		{
			tau = (age - tms)/(timescales[T_BGB] - tms);
            limit_tau(&tau);
//			printf("tau HG %g tms %g age %g T_BGB %g\n",tau,*tms,age,timescales[T_BGB]);
			if (tau < d_m_0200[7])
			{
				return d_m_0200[5] + d_m_0200[6]*tau;
			}
			if ((tau > d_m_0200[7]) and (tau < d_m_0200[10]))
			{
				return d_m_0200[8] + d_m_0200[9]*tau;
			}
			if ((tau > d_m_0200[10]) and (tau < d_m_0200[13]))
			{
				return d_m_0200[11] + d_m_0200[12]*tau;
			}
			else
			{
				return d_m_0200[14] + d_m_0200[15]*tau;
			}
		}
		else if (stellar_type == 3) /* RED GIANT */
		{
			tau = (age-timescales[T_BGB])/(timescales[T_HE_IGNITION]-timescales[T_BGB]);	/* See hrdiag & hrdiag_RG */
            limit_tau(&tau);
//			printf("tau RG %g age %g\n",tau,age);
			if (tau < d_m_0200[18])
			{
				return d_m_0200[16] + d_m_0200[17]*tau;
			}
			else
			{
				return d_m_0200[19] + d_m_0200[20]*tau;
			}
		}
		else if (stellar_type == 4) /* He MS */
		{
			return 0.04;
		}
		else if (stellar_type == 5 || stellar_type == 6) /* EAGB/TPAGB */
		{
			return 0.04;
		}
	}
	else if (log_m == 0.250)	/* 1.78 M_SUN */
	{
		if (stellar_type == 1)	/* MAIN SEQUENCE */
		{
			tau = age/tms;
            limit_tau(&tau);
//			printf("tau MS %g age %g tms %g\n",tau,age,*tms);
			if (tau < d_m_0250[2])
			{
				return d_m_0250[0] + d_m_0250[1]*tau;
			}
			else
			{
				return d_m_0250[3] + d_m_0250[4]*tau;
			}
		}
		else if (stellar_type == 2)	/* HERTZSPRUNG GAP */
		{
			tau = (age - tms)/(timescales[T_BGB] - tms);
            limit_tau(&tau);
//			printf("tau HG %g tms %g age %g T_BGB %g\n",tau,*tms,age,timescales[T_BGB]);
			if (tau < d_m_0250[7])
			{
				return d_m_0250[5] + d_m_0250[6]*tau;
			}
			if ((tau > d_m_0250[7]) and (tau < d_m_0250[10]))
			{
				return d_m_0250[8] + d_m_0250[9]*tau;
			}
			if ((tau > d_m_0250[10]) and (tau < d_m_0250[13]))
			{
				return d_m_0250[11] + d_m_0250[12]*tau;
			}
			else
			{
				return d_m_0250[14] + d_m_0250[15]*tau;
			}
		}
		else if (stellar_type == 3) /* RED GIANT */
		{
			tau = (age-timescales[T_BGB])/(timescales[T_HE_IGNITION]-timescales[T_BGB]);	/* See hrdiag & hrdiag_RG */
            limit_tau(&tau);
//			printf("tau RG %g age %g\n",tau,age);
			if (tau < d_m_0250[18])
			{
				return d_m_0250[16] + d_m_0250[17]*tau;
			}
			else
			{
				return d_m_0250[19] + d_m_0250[20]*tau;
			}
		}
		else if (stellar_type == 4) /* He MS */
		{
			return 0.04;
		}
		else if (stellar_type == 5 || stellar_type == 6) /* EAGB/TPAGB */
		{
			return 0.04;
		}
	}
	else if (log_m == 0.300)	/* 2.00 M_SUN */
	{
		if (stellar_type == 1)	/* MAIN SEQUENCE */
		{
			tau = age/tms;
            limit_tau(&tau);
//			printf("tau MS %g age %g tms %g\n",tau,age,*tms);
			return d_m_0300[0] + d_m_0300[1]*tau;
		}
		else if (stellar_type == 2)	/* HERTZSPRUNG GAP */
		{
			tau = (age - tms)/(timescales[T_BGB] - tms);
//			printf("tau HG %g tms %g age %g T_BGB %g\n",tau,*tms,age,timescales[T_BGB]);
			if (tau < d_m_0300[4])
			{
				return d_m_0300[2] + d_m_0300[3]*tau;
			}
			if ((tau > d_m_0300[4]) and (tau < d_m_0300[7]))
			{
				return d_m_0300[5] + d_m_0300[6]*tau;
			}
			else
			{
				return d_m_0300[8] + d_m_0300[9]*tau;
			}
		}
		else if (stellar_type == 3) /* RED GIANT */
		{
			tau = (age-timescales[T_BGB])/(timescales[T_HE_IGNITION]-timescales[T_BGB]);	/* See hrdiag & hrdiag_RG */
            limit_tau(&tau);
//			printf("tau RG %g age %g\n",tau,age);
			if (tau < d_m_0300[12])
			{
				return d_m_0300[10] + d_m_0300[11]*tau;
			}
			if ((tau > d_m_0300[12]) and (tau < d_m_0300[15]))
			{
				return d_m_0300[13] + d_m_0300[14]*tau;
			}
			else
			{
				return d_m_0300[16] + d_m_0300[17]*tau;
			}
		}
		else if (stellar_type == 4) /* He MS */
		{
			return 0.04;
		}
		else if (stellar_type == 5 || stellar_type == 6) /* EAGB/TPAGB */
		{
			return 0.04;
		}
	}
	else if (log_m == 0.400)	/* 2.51 M_SUN */
	{
		if (stellar_type == 1)	/* MAIN SEQUENCE */
		{
			tau = age/tms;
            limit_tau(&tau);
//			printf("tau MS %g age %g tms %g\n",tau,age,*tms);
			return d_m_0400[0] + d_m_0400[1]*tau;
		}
		else if (stellar_type == 2)	/* HERTZSPRUNG GAP */
		{
			tau = (age - tms)/(timescales[T_BGB] - tms);
            limit_tau(&tau);
//			printf("tau HG %g tms %g age %g T_BGB %g\n",tau,*tms,age,timescales[T_BGB]);
			if (tau < d_m_0400[4])
			{
				return d_m_0400[2] + d_m_0400[3]*tau;
			}
			else
			{
				return d_m_0400[5] + d_m_0400[6]*tau;
			}
		}
		else if (stellar_type == 3) /* RED GIANT */
		{
			tau = (age-timescales[T_BGB])/(timescales[T_HE_IGNITION]-timescales[T_BGB]);	/* See hrdiag & hrdiag_RG */
            limit_tau(&tau);
//			printf("tau RG %g age %g\n",tau,age);
			if (tau < d_m_0400[9])
			{
				return d_m_0400[7] + d_m_0400[8]*tau;
			}
			if ((tau > d_m_0400[9]) and (tau < d_m_0400[12]))
			{
				return d_m_0400[10] + d_m_0400[11]*tau;
			}
			else
			{
				return d_m_0400[13] + d_m_0400[14]*tau;
			}
		}
		else if (stellar_type == 4) /* CHeB */
		{
			tau = (age - timescales[T_HE_IGNITION])/timescales[T_HE_BURNING];
            limit_tau(&tau);
//			printf("tau CHeB %g age %g\n",tau,age);
			if (tau < d_m_0400[17])
			{
				return d_m_0400[15] + d_m_0400[16]*tau;
			}
			else
			{
				return d_m_0400[18] + d_m_0400[19]*tau;
			}

		}
		else if (stellar_type == 5) /* EAGB */
		{
			return 0.05;
		}
		else if (stellar_type == 6) /* TPAGB */
		{
			return 0.044;
		}
	}
	else if (log_m == 0.500)	/* 3.16 M_SUN */
	{
		if (stellar_type == 1)	/* MAIN SEQUENCE */
		{
			tau = age/tms;
            limit_tau(&tau);
//			printf("tau MS %g age %g tms %g\n",tau,age,*tms);
			return d_m_0500[0] + d_m_0500[1]*tau;
		}
		else if (stellar_type == 2)	/* HERTZSPRUNG GAP */
		{
			tau = (age - tms)/(timescales[T_BGB] - tms);
            limit_tau(&tau);
//			printf("tau HG %g tms %g age %g T_BGB %g\n",tau,*tms,age,timescales[T_BGB]);
			if (tau < d_m_0500[4])
			{
				return d_m_0500[2] + d_m_0500[3]*tau;
			}
			else
			{
				return d_m_0500[5] + d_m_0500[6]*tau;
			}
		}
		else if (stellar_type == 3) /* RED GIANT */
		{
			tau = (age-timescales[T_BGB])/(timescales[T_HE_IGNITION]-timescales[T_BGB]);	/* See hrdiag & hrdiag_RG */
            limit_tau(&tau);
//			printf("tau RG %g age %g\n",tau,age);
			if (tau < d_m_0500[9])
			{
				return d_m_0500[7] + d_m_0500[8]*tau;
			}
			else
			{
				return d_m_0500[10] + d_m_0500[11]*tau;
			}
		}
		else if (stellar_type == 4) /* CHeB */
		{
			tau = (age - timescales[T_HE_IGNITION])/timescales[T_HE_BURNING];
            limit_tau(&tau);
			if (tau < d_m_0500[14])
			{
				return d_m_0500[12] + d_m_0500[13]*tau;
			}
			if ((tau > d_m_0500[14]) and (tau < d_m_0500[17]))
			{
				return d_m_0500[15] + d_m_0500[16]*tau;
			}
			else
			{
				return d_m_0500[18] + d_m_0500[19]*tau;
			}

		}
		else if (stellar_type == 5) /* EAGB */
		{
			return 0.056;

		}
		else if (stellar_type == 6) /* TPAGB */
		{
			return 0.07;
		}
	}
	else if (log_m == 0.600)	/* 3.98 M_SUN */
	{
		if (stellar_type == 1)	/* MAIN SEQUENCE */
		{
			tau = age/tms;
            limit_tau(&tau);
//			printf("tau MS %g age %g tms %g\n",tau,age,*tms);
			return d_m_0600[0] + d_m_0600[1]*tau;
		}
		else if (stellar_type == 2)	/* HERTZSPRUNG GAP */
		{
			tau = (age - tms)/(timescales[T_BGB] - tms);
            limit_tau(&tau);
//			printf("tau HG %g tms %g age %g T_BGB %g\n",tau,*tms,age,timescales[T_BGB]);
			if (tau < d_m_0600[4])
			{
				return d_m_0600[2] + d_m_0600[3]*tau;
			}
			else
			{
				return d_m_0600[5] + d_m_0600[6]*tau;
			}
		}
		else if (stellar_type == 3) /* RED GIANT */
		{
			tau = (age-timescales[T_BGB])/(timescales[T_HE_IGNITION]-timescales[T_BGB]);	/* See hrdiag & hrdiag_RG */
            limit_tau(&tau);
//			printf("tau RG %g age %g\n",tau,age);
			if (tau < d_m_0600[9])
			{
				return d_m_0600[7] + d_m_0600[8]*tau;
			}
			else
			{
				return d_m_0600[10] + d_m_0600[11]*tau;
			}
		}
		else if (stellar_type == 4) /* CHeB */
		{
			tau = (age - timescales[T_HE_IGNITION])/timescales[T_HE_BURNING];
            limit_tau(&tau);
			if (tau < d_m_0600[14])
			{
				return d_m_0600[12] + d_m_0600[13]*tau;
			}
			if ((tau > d_m_0600[14]) and (tau < d_m_0600[17]))
			{
				return d_m_0600[15] + d_m_0600[16]*tau;
			}
			else
			{
				return d_m_0600[18] + d_m_0600[19]*tau;
			}

		}
		else if (stellar_type == 5) /* EAGB */
		{
			return 0.050;

		}
		else if (stellar_type == 6) /* TPAGB */
		{
			return 0.041;
		}
	}
	else if (log_m == 0.700)	/* 5.01 M_SUN */
	{
		if (stellar_type == 1)	/* MAIN SEQUENCE */
		{
			tau = age/tms;
            limit_tau(&tau);
//			printf("tau MS %g age %g tms %g\n",tau,age,*tms);
			return d_m_0700[0] + d_m_0700[1]*tau;
		}
		else if (stellar_type == 2)	/* HERTZSPRUNG GAP */
		{
			tau = (age - tms)/(timescales[T_BGB] - tms);
            limit_tau(&tau);
//			printf("tau HG %g tms %g age %g T_BGB %g\n",tau,*tms,age,timescales[T_BGB]);
			if (tau < d_m_0700[4])
			{
				return d_m_0700[2] + d_m_0700[3]*tau;
			}
			else
			{
				return d_m_0700[5] + d_m_0700[6]*tau;
			}
		}
		else if (stellar_type == 3) /* RED GIANT */
		{
			tau = (age-timescales[T_BGB])/(timescales[T_HE_IGNITION]-timescales[T_BGB]);	/* See hrdiag & hrdiag_RG */
            limit_tau(&tau);
//			printf("tau RG %g age %g\n",tau,age);
			if (tau < d_m_0700[9])
			{
				return d_m_0700[7] + d_m_0700[8]*tau;
			}
			else
			{
				return d_m_0700[10] + d_m_0700[11]*tau;
			}
		}
		else if (stellar_type == 4) /* CHeB */
		{
			tau = (age - timescales[T_HE_IGNITION])/timescales[T_HE_BURNING];
            limit_tau(&tau);
			if (tau < d_m_0700[14])
			{
				return d_m_0700[12] + d_m_0700[13]*tau;
			}
			if ((tau > d_m_0700[14]) and (tau < d_m_0700[17]))
			{
				return d_m_0700[15] + d_m_0700[16]*tau;
			}
			if ((tau > d_m_0700[17]) and (tau < d_m_0700[20]))
			{
				return d_m_0700[18] + d_m_0700[19]*tau;
			}
			else
			{
				return d_m_0700[21] + d_m_0700[22]*tau;
			}

		}
		else if (stellar_type == 5) /* EAGB */
		{
			return 0.039;

		}
		else if (stellar_type == 6) /* TPAGB */
		{
			return 0.013;
		}
	}
	else if (log_m == 0.800)	/* 6.31 M_SUN */
	{
		if (stellar_type == 1)	/* MAIN SEQUENCE */
		{
			tau = age/tms;
            limit_tau(&tau);
//			printf("tau MS %g age %g tms %g\n",tau,age,*tms);
			return d_m_0800[0] + d_m_0800[1]*tau;
		}
		else if (stellar_type == 2)	/* HERTZSPRUNG GAP */
		{
			tau = (age - tms)/(timescales[T_BGB] - tms);
            limit_tau(&tau);
//			printf("tau HG %g tms %g age %g T_BGB %g\n",tau,*tms,age,timescales[T_BGB]);
			if (tau < d_m_0800[4])
			{
				return d_m_0800[2] + d_m_0800[3]*tau;
			}
			else
			{
				return d_m_0800[5] + d_m_0800[6]*tau;
			}
		}
		else if (stellar_type == 3) /* RED GIANT */
		{
			tau = (age-timescales[T_BGB])/(timescales[T_HE_IGNITION]-timescales[T_BGB]);	/* See hrdiag & hrdiag_RG */
            limit_tau(&tau);
//			printf("tau RG %g age %g\n",tau,age);
			if (tau < d_m_0800[9])
			{
				return d_m_0800[7] + d_m_0800[8]*tau;
			}
			else
			{
				return d_m_0800[10] + d_m_0800[11]*tau;
			}
		}
		else if (stellar_type == 4) /* CHeB */
		{
			tau = (age - timescales[T_HE_IGNITION])/timescales[T_HE_BURNING];
            limit_tau(&tau);
			if (tau < d_m_0800[14])
			{
				return d_m_0800[12] + d_m_0800[13]*tau;
			}
			if ((tau > d_m_0800[14]) and (tau < d_m_0800[17]))
			{
				return d_m_0800[15] + d_m_0800[16]*tau;
			}
			if ((tau > d_m_0800[17]) and (tau < d_m_0800[20]))
			{
				return d_m_0800[18] + d_m_0800[19]*tau;
			}
			else
			{
				return d_m_0800[21] + d_m_0800[22]*tau;
			}

		}
		else if (stellar_type == 5) /* EAGB */
		{
			return 0.038;

		}
		else if (stellar_type == 6) /* TPAGB */
		{
			return 0.047;
		}
	}
	else if (log_m == 0.900)	/* 7.94 M_SUN */
	{
		if (stellar_type == 1)	/* MAIN SEQUENCE */
		{
			tau = age/tms;
            limit_tau(&tau);
//			printf("tau MS %g age %g tms %g\n",tau,age,*tms);
			return d_m_0900[0] + d_m_0900[1]*tau;
		}
		else if (stellar_type == 2)	/* HERTZSPRUNG GAP */
		{
			tau = (age - tms)/(timescales[T_BGB] - tms);
            limit_tau(&tau);
//			printf("tau HG %g tms %g age %g T_BGB %g\n",tau,*tms,age,timescales[T_BGB]);
			if (tau < d_m_0900[4])
			{
				return d_m_0900[2] + d_m_0900[3]*tau;
			}
			else
			{
				return d_m_0900[5] + d_m_0900[6]*tau;
			}
		}
		else if (stellar_type == 3) /* RED GIANT */
		{
			tau = (age-timescales[T_BGB])/(timescales[T_HE_IGNITION]-timescales[T_BGB]);	/* See hrdiag & hrdiag_RG */
            limit_tau(&tau);
//			printf("tau RG %g age %g\n",tau,age);
			if (tau < d_m_0900[9])
			{
				return d_m_0900[7] + d_m_0900[8]*tau;
			}
			else
			{
				return d_m_0900[10] + d_m_0900[11]*tau;
			}
		}
		else if (stellar_type == 4) /* CHeB */
		{
			tau = (age - timescales[T_HE_IGNITION])/timescales[T_HE_BURNING];
            limit_tau(&tau);
			if (tau < d_m_0900[14])
			{
				return d_m_0900[12] + d_m_0900[13]*tau;
			}
			if ((tau > d_m_0900[14]) and (tau < d_m_0900[17]))
			{
				return d_m_0900[15] + d_m_0900[16]*tau;
			}
			if ((tau > d_m_0900[17]) and (tau < d_m_0900[20]))
			{
				return d_m_0900[18] + d_m_0900[19]*tau;
			}
			else
			{
				return d_m_0900[21] + d_m_0900[22]*tau;
			}

		}
		else if (stellar_type == 5) /* EAGB */
		{
			return 0.039;

		}
		else if (stellar_type == 6) /* TPAGB */
		{
			return 0.051;
		}
	}
	else if (log_m == 1.00)	/* 10.00 M_SUN */
	{
		if (stellar_type == 1)	/* MAIN SEQUENCE */
		{
			tau = age/tms;
            limit_tau(&tau);
//			printf("tau MS %g age %g tms %g\n",tau,age,*tms);
			return d_m_1000[0] + d_m_1000[1]*tau;
		}
		else if (stellar_type == 2)	/* HERTZSPRUNG GAP */
		{
			tau = (age - tms)/(timescales[T_BGB] - tms);
            limit_tau(&tau);
//			printf("tau HG %g tms %g age %g T_BGB %g\n",tau,*tms,age,timescales[T_BGB]);
			if (tau < d_m_1000[4])
			{
				return d_m_1000[2] + d_m_1000[3]*tau;
			}
			else
			{
				return d_m_1000[5] + d_m_1000[6]*tau;
			}
		}
		else if (stellar_type == 3) /* RED GIANT */
		{
			tau = (age-timescales[T_BGB])/(timescales[T_HE_IGNITION]-timescales[T_BGB]);	/* See hrdiag & hrdiag_RG */
            limit_tau(&tau);
			//printf("tau RG %g age %g timescales[T_BGB] %g\n",tau,age,timescales[T_BGB]);
			if (tau < d_m_1000[9])
			{
				return d_m_1000[7] + d_m_1000[8]*tau;
			}
			else
			{
				return d_m_1000[10] + d_m_1000[11]*tau;
			}
		}
		else if (stellar_type == 4) /* CHeB */
		{
			tau = (age - timescales[T_HE_IGNITION])/timescales[T_HE_BURNING];
            limit_tau(&tau);
			if (tau < d_m_1000[14])
			{
				return d_m_1000[12] + d_m_1000[13]*tau;
			}
			if ((tau > d_m_1000[14]) and (tau < d_m_1000[17]))
			{
				return d_m_1000[15] + d_m_1000[16]*tau;
			}
			if ((tau > d_m_1000[17]) and (tau < d_m_1000[20]))
			{
				return d_m_1000[18] + d_m_1000[19]*tau;
			}
			else
			{
				return d_m_1000[21] + d_m_1000[22]*tau;
			}

		}
		else if (stellar_type == 5) /* EAGB */
		{
			return 0.0054;

		}
		else if (stellar_type == 6) /* TPAGB */
		{
			return 0.051;
		}
	}
	else if (log_m == 1.10)	/* 12.59 M_SUN */
	{
		if (stellar_type == 1)	/* MAIN SEQUENCE */
		{
			tau = age/tms;
            limit_tau(&tau);
//			printf("tau MS %g age %g tms %g\n",tau,age,*tms);
			return d_m_1100[0] + d_m_1100[1]*tau;
		}
		else if (stellar_type == 2)	/* HERTZSPRUNG GAP */
		{
			tau = (age - tms)/(timescales[T_BGB] - tms);
            limit_tau(&tau);
//			printf("tau HG %g tms %g age %g T_BGB %g\n",tau,*tms,age,timescales[T_BGB]);
			if (tau < d_m_1100[4])
			{
				return d_m_1100[2] + d_m_1100[3]*tau;
			}
			else
			{
				return d_m_1100[5] + d_m_1100[6]*tau;
			}
		}
		else if (stellar_type == 3) /* RED GIANT */
		{
			tau = (age-timescales[T_BGB])/(timescales[T_HE_IGNITION]-timescales[T_BGB]);	/* See hrdiag & hrdiag_RG */
            limit_tau(&tau);
//			printf("tau RG %g age %g\n",tau,age);
			if (tau < d_m_1100[9])
			{
				return d_m_1100[7] + d_m_1100[8]*tau;
			}
			else
			{
				return d_m_1100[10] + d_m_1100[11]*tau;
			}
		}
		else if (stellar_type == 4) /* CHeB */
		{
			tau = (age - timescales[T_HE_IGNITION])/timescales[T_HE_BURNING];
            limit_tau(&tau);
			if (tau < d_m_1100[14])
			{
				return d_m_1100[12] + d_m_1100[13]*tau;
			}
			if ((tau > d_m_1100[14]) and (tau < d_m_1100[17]))
			{
				return d_m_1100[15] + d_m_1100[16]*tau;
			}
			if ((tau > d_m_1100[17]) and (tau < d_m_1100[20]))
			{
				return d_m_1100[18] + d_m_1100[19]*tau;
			}
			else
			{
				return d_m_1100[21] + d_m_1100[22]*tau;
			}

		}
		else if (stellar_type == 5) /* EAGB */
		{
			return 0.037;

		}
		else if (stellar_type == 6) /* TPAGB */
		{
			return 0.050;
		}
	}
	else if (log_m == 1.20)	/* 15.85 M_SUN */
	{
		if (stellar_type == 1)	/* MAIN SEQUENCE */
		{
			tau = age/tms;
            limit_tau(&tau);
//			printf("tau MS %g age %g tms %g\n",tau,age,*tms);
			return d_m_1200[0] + d_m_1200[1]*tau;
		}
		else if (stellar_type == 2)	/* HERTZSPRUNG GAP */
		{
			tau = (age - tms)/(timescales[T_BGB] - tms);
            limit_tau(&tau);
			//printf("tau HG %g tms %g age %g T_BGB %g\n",tau,*tms,age,timescales[T_BGB]);
			if (tau < d_m_1200[4])
			{
				return d_m_1200[2] + d_m_1200[3]*tau;
			}
			else
			{
				return d_m_1200[5] + d_m_1200[6]*tau;
			}
		}
		else if (stellar_type == 3) /* RED GIANT */
		{
			tau = (age-timescales[T_BGB])/(timescales[T_HE_IGNITION]-timescales[T_BGB]);	/* See hrdiag & hrdiag_RG */
            limit_tau(&tau);
//			printf("tau RG %g age %g\n",tau,age);
			if (tau < d_m_1200[9])
			{
				return d_m_1200[7] + d_m_1200[8]*tau;
			}
			else
			{
				return d_m_1200[10] + d_m_1200[11]*tau;
			}
		}
		else if (stellar_type == 4) /* CHeB */
		{
			tau = (age - timescales[T_HE_IGNITION])/timescales[T_HE_BURNING];
            limit_tau(&tau);
//			printf("tau CHeB %g \n",tau);
			if (tau < d_m_1200[14])
			{
				return d_m_1200[12] + d_m_1200[13]*tau;
			}
			if ((tau > d_m_1200[14]) and (tau < d_m_1200[17]))
			{
				return d_m_1200[15] + d_m_1200[16]*tau;
			}
			else
			{
				return d_m_1200[18] + d_m_1200[19]*tau;
			}

		}
		else
		{
			return 0.04;

		}
	}
	else if (log_m == 1.30)	/* 19.95 M_SUN */
	{
		if (stellar_type == 1)	/* MAIN SEQUENCE */
		{
			tau = age/tms;
            limit_tau(&tau);
//			printf("tau MS %g age %g tms %g\n",tau,age,*tms);
			return d_m_1300[0] + d_m_1300[1]*tau;
		}
		else if (stellar_type == 2)	/* HERTZSPRUNG GAP */
		{
			return 0.025;
		}
		else
		{
			return 0.04;
		}
	}
	else if (log_m == 1.40)	/* 25.12 M_SUN */
	{
		if (stellar_type == 1)	/* MAIN SEQUENCE */
		{
			tau = age/tms;
            limit_tau(&tau);
//			printf("tau MS %g age %g tms %g\n",tau,age,*tms);
			return d_m_1400[0] + d_m_1400[1]*tau;
		}
		else if (stellar_type == 2)	/* HERTZSPRUNG GAP */
		{
			return 0.02;
		}
		else
		{
			return 0.04;
		}
	}
	else if (log_m == 1.50)	/* 31.62 M_SUN */
	{
		if (stellar_type == 1)	/* MAIN SEQUENCE */
		{
			tau = age/tms;
            limit_tau(&tau);
			//printf("log_m 1.50 tau MS %g age %g tms %g %g %g\n",tau,age,tms,d_m_1500[0] , d_m_1500[1]);
			return d_m_1500[0] + d_m_1500[1]*tau;
		}
		else if (stellar_type == 2)	/* HERTZSPRUNG GAP */
		{
			return 0.03;
		}
		else
		{
			return 0.04;
		}
	}
	else if (log_m == 1.60)	/* 39.81 M_SUN */
	{
		if (stellar_type == 1)	/* MAIN SEQUENCE */
		{
			tau = age/tms;
            limit_tau(&tau);
			//printf("log_m 1.60 tau MS %g age %g tms %g %g %g\n",tau,age,tms,d_m_1600[0],d_m_1600[1]);
			return d_m_1600[0] + d_m_1600[1]*tau;
		}
		else if (stellar_type == 3)	/* Giant */
		{
			return 0.045;
		}
		else
		{
			return 0.001;
		}
	}
	else if (log_m == 1.70)	/* 50.12 M_SUN */
	{
		if (stellar_type == 1)	/* MAIN SEQUENCE */
		{
			tau = age/tms;
            limit_tau(&tau);
//			printf("tau MS %g age %g tms %g\n",tau,age,*tms);
			return d_m_1700[0] + d_m_1700[1]*tau;
		}
		else if (stellar_type == 3)	/* Giant */
		{
			return 0.038;
		}
		else
		{
			return 0.002;
		}
	}
	else if (log_m == 1.80)	/* 63.10 M_SUN */
	{
		if (stellar_type == 1)	/* MAIN SEQUENCE */
		{
			tau = age/tms;
            limit_tau(&tau);
//			printf("tau MS %g age %g tms %g\n",tau,age,*tms);
			if (tau < d_m_1800[2])
			{
				return d_m_1800[0] + d_m_1800[1]*tau;
			}
			else
			{
				return d_m_1800[3] + d_m_1800[4]*tau;
			}
		}
		else if (stellar_type == 3)	/* Giant */
		{
			return 0.01;
		}
		else
		{
			return 0.003;
		}
	}
	else if (log_m == 1.90)	/* 79.43 M_SUN */
	{
		if (stellar_type == 1)	/* MAIN SEQUENCE */
		{
			tau = age/tms;
            limit_tau(&tau);
//			printf("tau MS %g age %g tms %g\n",tau,age,*tms);
			if (tau < d_m_1900[2])
			{
				return d_m_1900[0] + d_m_1900[1]*tau;
			}
			else
			{
				return d_m_1900[3] + d_m_1900[4]*tau;
			}
		}
		else if (stellar_type == 3)	/* Giant */
		{
			return 0.005;
		}
		else
		{
			return 0.0002;
		}
	}
	else if (log_m == 2.00)	/* 100.00 M_SUN */
	{
		if (stellar_type == 1)	/* MAIN SEQUENCE */
		{
			tau = age/tms;
            limit_tau(&tau);
//			printf("tau MS %g age %g tms %g\n",tau,age,*tms);
			if (tau < d_m_2000[2])
			{
				return d_m_2000[0] + d_m_2000[1]*tau;
			}
			else
			{
				return d_m_2000[3] + d_m_2000[4]*tau;
			}
		}
		else if (stellar_type == 3)	/* Giant */
		{
			return 0.015;
		}
		else
		{
			return 0.00001;
		}
	}
	else if (log_m == 2.10)	/* 125.89 M_SUN */
	{
		if (stellar_type == 1)	/* MAIN SEQUENCE */
		{
			tau = age/tms;
            limit_tau(&tau);
//			printf("tau MS %g age %g tms %g\n",tau,age,*tms);
			if (tau < d_m_2100[2])
			{
				return d_m_2100[0] + d_m_2100[1]*tau;
			}
			else
			{
				return d_m_2100[3] + d_m_2100[4]*tau;
			}
		}
		else if (stellar_type == 3)	/* Giant */
		{
			return 0.015;
		}
		else
		{
			return 0.00002;
		}
	}

    return -1;
}					
