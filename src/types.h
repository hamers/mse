#include <math.h>
#include <cstdlib>
#include <map>
#include <vector>


extern "C"
{

/* constants */
/* units (cf. interface.py): 
 * unit_l = AU
 * unit_m = MSun
 * unit_t = yr
 */
 
//#define VERBOSE
 
#ifndef __FOUND_ROOT
#define ___FOUND_ROOT
#define FOUND_ROOT ((roots_found[i_root] == 1) || (roots_found[i_root] == -1))
#endif 

#ifndef __CONSTANTS
#define __CONSTANTS

// Default constants //
extern double CONST_G;
extern double CONST_G_P2;
extern double CONST_G_P3;
extern double CONST_C_LIGHT;
extern double CONST_C_LIGHT_P2;
extern double CONST_C_LIGHT_P4;
extern double CONST_C_LIGHT_P5;
extern double CONST_MSUN;
extern double CONST_R_SUN;
extern double CONST_L_SUN;
extern double CONST_KM_PER_S;
extern double CONST_PER_PC3;

// Default parameters //
extern double relative_tolerance;
extern double absolute_tolerance_eccentricity_vectors;
extern bool include_quadrupole_order_terms;
extern bool include_octupole_order_binary_pair_terms;
extern bool include_octupole_order_binary_triplet_terms;
extern bool include_hexadecupole_order_binary_pair_terms;
extern bool include_dotriacontupole_order_binary_pair_terms;
extern int random_seed;

extern bool include_flybys;
extern bool flybys_correct_for_gravitational_focussing;
extern int flybys_velocity_distribution;
extern int flybys_mass_distribution;
extern int flybys_N_enc;
extern int flybys_N_not_impulsive;
extern int flybys_reference_binary;
extern double flybys_t_next_encounter;
extern double flybys_mass_distribution_lower_value;
extern double flybys_mass_distribution_upper_value;
extern double flybys_encounter_sphere_radius;
extern double flybys_stellar_density;
extern double flybys_stellar_relative_velocity_dispersion;
extern double flybys_W_max;
extern double flybys_total_encounter_rate_at_R_enc;
extern double flybys_stellar_density_at_R_enc;
extern double flybys_internal_mass;
extern double flybys_internal_semimajor_axis;

#ifdef IGNORE
// Default constants //
double CONST_G = 4.0*M_PI*M_PI; 
double CONST_G_P2 = CONST_G*CONST_G;
double CONST_G_P3 = CONST_G_P2*CONST_G;
double CONST_C_LIGHT = 63239.72638679138;
double CONST_C_LIGHT_P2 = CONST_C_LIGHT*CONST_C_LIGHT;
double CONST_C_LIGHT_P4 = CONST_C_LIGHT_P2*CONST_C_LIGHT_P2;
double CONST_C_LIGHT_P5 = CONST_C_LIGHT_P4*CONST_C_LIGHT;
double CONST_MSUN = 1.0;
double CONST_R_SUN = 0.004649130343817401;
double CONST_L_SUN = 0.0002710404109745588;
double CONST_KM_PER_S = 0.210862;
double CONST_PER_PC3 = 1.14059e-16;

// Default parameters //
double relative_tolerance = 1.0e-12;
double absolute_tolerance_eccentricity_vectors = 1.0e-14;
bool include_quadrupole_order_terms = true;
bool include_octupole_order_binary_pair_terms = true;
bool include_octupole_order_binary_triplet_terms = true;
bool include_hexadecupole_order_binary_pair_terms = true;
bool include_dotriacontupole_order_binary_pair_terms = true;
int random_seed = 0;

bool include_flybys = true;
bool flybys_correct_for_gravitational_focussing = true;
int flybys_velocity_distribution = 0;
int flybys_mass_distribution = 0;
double flybys_mass_distribution_lower_value = 0.1;
double flybys_mass_distribution_upper_value = 100.0;
double flybys_encounter_sphere_radius = 1.0e4;
double flybys_stellar_density = 0.1*CONST_PER_PC3; /* density at infinity */
double flybys_stellar_relative_velocity_dispersion = 30.0*CONST_KM_PER_S;
double flybys_W_max = 0.0;
double flybys_total_encounter_rate_at_R_enc = 0.0;
double flybys_stellar_density_at_R_enc = 0.0; /* density at infinity */
double flybys_internal_mass = 0.0;
double flybys_internal_semimajor_axis = 0.0;
#endif

#ifdef IGNORE
#define CONST_C_LIGHT		(double)	63239726386.8
#define CONST_C_LIGHT_P2	(double)	CONST_C_LIGHT*CONST_C_LIGHT
#define CONST_C_LIGHT_P4	(double)	CONST_C_LIGHT_P2*CONST_C_LIGHT_P2
#define CONST_C_LIGHT_P5	(double)	CONST_C_LIGHT_P4*CONST_C_LIGHT
#define CONST_MSUN          (double)    1.0
#define CONST_R_SUN         (double)    0.00464913034382
#define CONST_L_SUN         (double)    2.71040410975e+14
#endif

#define yr_to_Myr           (double)    1.0e-6
#define Myr_to_yr           (double)    1.0e6
#define c_1div2             (double)    1.0/2.0
#define c_1div3             (double)    1.0/3.0
#define c_1div4             (double)    1.0/4.0
#define c_1div5             (double)    1.0/5.0
#define c_1div6             (double)    1.0/6.0
#define c_1div7             (double)    1.0/7.0
#define c_1div8             (double)    1.0/8.0
#define c_1div10            (double)    1.0/10.0
#define c_1div15            (double)    1.0/15.0
#define c_1div16            (double)    1.0/16.0
#define c_1div30            (double)    1.0/30.0
#define c_2div3             (double)    2.0/3.0
#define c_3div2             (double)    3.0/2.0
#define c_3div4             (double)    3.0/4.0
#define c_3div5             (double)    3.0/5.0
#define c_3div8             (double)    3.0/8.0
#define c_3div32            (double)    3.0/32.0
#define c_3div1024          (double)    3.0/1024.0
#define c_4div15            (double)    4.0/15.0
#define c_5div2             (double)    5.0/2.0
#define c_5div8             (double)    5.0/8.0
#define c_5div16            (double)    5.0/16.0
#define c_5div64            (double)    5.0/64.0
#define c_7div8             (double)    7.0/8.0
#define c_8div5             (double)    8.0/5.0
#define c_8div7             (double)    8.0/7.0
#define c_9div2             (double)    9.0/2.0
#define c_9div16            (double)    9.0/16.0
#define c_9div32            (double)    9.0/32.0
#define c_11div18           (double)    11.0/18.0
#define c_15div2            (double)    15.0/2.0
#define c_15div4            (double)    15.0/4.0
#define c_15div8            (double)    15.0/8.0
#define c_15div16           (double)    15.0/16.0
#define c_16div5            (double)    16.0/5.0
#define c_25div16           (double)    25.0/16.0
#define c_25div64           (double)    25.0/64.0
#define c_31div2            (double)    31.0/2.0
#define c_32div5            (double)    32.0/5.0
#define c_37div96           (double)    37.0/96.0
#define c_45div8            (double)    45.0/8.0
#define c_64div5            (double)    64.0/5.0
#define c_73div24           (double)    73.0/24.0
#define c_105div4096        (double)    105.0/4096.0
#define c_121div304         (double)    121.0/304.0
#define c_185div16          (double)    185.0/16.0
#define c_255div8           (double)    255.0/8.0
#define c_304div15          (double)    304.0/15.0
#endif


#ifndef __TABLES
#define __TABLES
#define MAX_ORDER (int) 5

#define TABLEWIDTH_A (int) 3
#define TABLELENGTH_A (int) 10
const double A_TABLE[TABLELENGTH_A][TABLEWIDTH_A] =
{{2, 0, -0.500000000000000}, {2, 2, 1.50000000000000}, {3, 
  1, -1.50000000000000}, {3, 3, 2.50000000000000}, {4, 0, 
  0.375000000000000}, {4, 2, -3.75000000000000}, {4, 4, 
  4.37500000000000}, {5, 1, 1.87500000000000}, {5, 
  3, -8.75000000000000}, {5, 5, 7.87500000000000}};
  
#define TABLEWIDTH_B (int) 8
#define TABLELENGTH_B (int) 47
#define HIGHEST_POWER_ECCP2_IN_B_TABLE (int) 3
const double B_TABLE[TABLELENGTH_B][TABLEWIDTH_B] =
{{2, 0, 0, 0, 1.00000000000000, 1.50000000000000, 0, 0}, {2, 1, 1, 
  0, -2.00000000000000, -0.500000000000000, 0, 0}, {2, 2, 0, 0, 
  0.500000000000000, -0.500000000000000, 0, 0}, {2, 2, 0, 
  2, -0.500000000000000, 0, 0, 0}, {2, 2, 2, 0, 2.50000000000000, 0, 
  0, 0}, {3, 0, 0, 0, 1.00000000000000, 3.00000000000000, 
  0.375000000000000, 0}, {3, 1, 1, 
  0, -2.50000000000000, -1.87500000000000, 0, 0}, {3, 2, 0, 0, 
  0.500000000000000, -0.375000000000000, -0.125000000000000, 0}, {3, 
  2, 0, 2, -0.500000000000000, -0.125000000000000, 0, 0}, {3, 2, 2, 0,
   3.75000000000000, 0.625000000000000, 0, 0}, {3, 3, 1, 
  0, -1.87500000000000, 1.87500000000000, 0, 0}, {3, 3, 1, 2, 
  1.87500000000000, 0, 0, 0}, {3, 3, 3, 0, -4.37500000000000, 0, 0, 
  0}, {4, 0, 0, 0, 1.00000000000000, 5.00000000000000, 
  1.87500000000000, 0}, {4, 1, 1, 
  0, -3.00000000000000, -4.50000000000000, -0.375000000000000, 0}, {4,
   2, 0, 0, 0.500000000000000, -0.125000000000000, -0.375000000000000,
   0}, {4, 2, 0, 2, -0.500000000000000, -0.375000000000000, 0, 0}, {4,
   2, 2, 0, 5.25000000000000, 2.62500000000000, 0, 0}, {4, 3, 1, 
  0, -2.25000000000000, 1.87500000000000, 0.375000000000000, 0}, {4, 
  3, 1, 2, 2.25000000000000, 0.375000000000000, 0, 0}, {4, 3, 3, 
  0, -7.00000000000000, -0.875000000000000, 0, 0}, {4, 4, 0, 0, 
  0.375000000000000, -0.750000000000000, 0.375000000000000, 0}, {4, 4,
   0, 2, -0.750000000000000, 0.750000000000000, 0, 0}, {4, 4, 0, 4, 
  0.375000000000000, 0, 0, 0}, {4, 4, 2, 0, 
  5.25000000000000, -5.25000000000000, 0, 0}, {4, 4, 2, 
  2, -5.25000000000000, 0, 0, 0}, {4, 4, 4, 0, 7.87500000000000, 0, 0,
   0}, {5, 0, 0, 0, 1.00000000000000, 7.50000000000000, 
  5.62500000000000, 0.312500000000000}, {5, 1, 1, 
  0, -3.50000000000000, -8.75000000000000, -2.18750000000000, 0}, {5, 
  2, 0, 0, 0.500000000000000, 
  0.250000000000000, -0.687500000000000, -0.0625000000000000}, {5, 2, 
  0, 2, -0.500000000000000, -0.750000000000000, -0.0625000000000000, 
  0}, {5, 2, 2, 0, 7.00000000000000, 7.00000000000000, 
  0.437500000000000, 0}, {5, 3, 1, 0, -2.62500000000000, 
  1.31250000000000, 1.31250000000000, 0}, {5, 3, 1, 2, 
  2.62500000000000, 1.31250000000000, 0, 0}, {5, 3, 3, 
  0, -10.5000000000000, -3.93750000000000, 0, 0}, {5, 4, 0, 0, 
  0.375000000000000, -0.687500000000000, 0.250000000000000, 
  0.0625000000000000}, {5, 4, 0, 2, -0.750000000000000, 
  0.625000000000000, 0.125000000000000, 0}, {5, 4, 0, 4, 
  0.375000000000000, 0.0625000000000000, 0, 0}, {5, 4, 2, 0, 
  7.00000000000000, -6.12500000000000, -0.875000000000000, 0}, {5, 4, 
  2, 2, -7.00000000000000, -0.875000000000000, 0, 0}, {5, 4, 4, 0, 
  13.1250000000000, 1.31250000000000, 0, 0}, {5, 5, 1, 
  0, -2.18750000000000, 4.37500000000000, -2.18750000000000, 0}, {5, 
  5, 1, 2, 4.37500000000000, -4.37500000000000, 0, 0}, {5, 5, 1, 
  4, -2.18750000000000, 0, 0, 0}, {5, 5, 3, 0, -13.1250000000000, 
  13.1250000000000, 0, 0}, {5, 5, 3, 2, 13.1250000000000, 0, 0, 
  0}, {5, 5, 5, 0, -14.4375000000000, 0, 0, 0}};

#define TABLEWIDTH_D (int) 8
#define TABLELENGTH_D (int) 15
const double D_TABLE[TABLELENGTH_D][TABLEWIDTH_D] = 
{
    {2, 0, 0, 0, 0, 0, 0, 1}, \
    {2, 0, 2, 0, 0, 0, 0, 2}, \
    {2, 0, 2, 0, 2, 0, 0, 3}, \
    {2, 0, 2, 0, 0, 0, 2, 4}, \
    {2, 2, 0, 0, 0, 0, 0, 5}, \
    {2, 2, 0, 2, 0, 0, 0, 6}, \
    {2, 2, 0, 0, 0, 2, 0, 7}, \
    {3, 1, 0, 1, 0, 0, 0, 8}, \
    {3, 1, 2, 0, 1, 1, 1, 9}, \
    {3, 1, 2, 1, 0, 0, 0, 10}, \
    {3, 1, 2, 1, 0, 0, 2, 11}, \
    {3, 1, 2, 1, 2, 0, 0, 12}, \
    {3, 3, 0, 1, 0, 0, 0, 13}, \
    {3, 3, 0, 1, 0, 2, 0, 14}, \
    {3, 3, 0, 3, 0, 0, 0, 15}
};

/* the numbers in each last entry of the D table refer to the functions defined below */

inline double D_TABLE_FUNC1(double ep, double ep_p2, double ef, double ef_p2, double ef_p4, double one_div_ef_p1, double one_div_ef_p2, double one_div_ef_p3, double ef_p2_minus_one, double sqrt_ef_p2_minus_one, double asec_minus_ef)
{
    return 2.0*(sqrt_ef_p2_minus_one + asec_minus_ef);
}
inline double D_TABLE_FUNC1_DER(double ep, double ep_p2, double ef, double ef_p2, double ef_p4, double one_div_ef_p1, double one_div_ef_p2, double one_div_ef_p3, double ef_p2_minus_one, double sqrt_ef_p2_minus_one, double asec_minus_ef)
{
    return 0.0;
}

inline double D_TABLE_FUNC2(double ep, double ep_p2, double ef, double ef_p2, double ef_p4, double one_div_ef_p1, double one_div_ef_p2, double one_div_ef_p3, double ef_p2_minus_one, double sqrt_ef_p2_minus_one, double asec_minus_ef)
{
    return (1.0 - ep_p2)*( c_1div3*one_div_ef_p2*sqrt_ef_p2_minus_one*(1.0 + 2.0*ef_p2) + asec_minus_ef);
}
inline double D_TABLE_FUNC2_DER(double ep, double ep_p2, double ef, double ef_p2, double ef_p4, double one_div_ef_p1, double one_div_ef_p2, double one_div_ef_p3, double ef_p2_minus_one, double sqrt_ef_p2_minus_one, double asec_minus_ef)
{
    return -2.0*ep*( c_1div3*one_div_ef_p2*sqrt_ef_p2_minus_one*(1.0 + 2.0*ef_p2) + asec_minus_ef);
}

inline double D_TABLE_FUNC3(double ep, double ep_p2, double ef, double ef_p2, double ef_p4, double one_div_ef_p1, double one_div_ef_p2, double one_div_ef_p3, double ef_p2_minus_one, double sqrt_ef_p2_minus_one, double asec_minus_ef)
{
    return -c_1div3*one_div_ef_p2*sqrt_ef_p2_minus_one*(1.0 + 2.0*ef_p2) - asec_minus_ef;
}
inline double D_TABLE_FUNC3_DER(double ep, double ep_p2, double ef, double ef_p2, double ef_p4, double one_div_ef_p1, double one_div_ef_p2, double one_div_ef_p3, double ef_p2_minus_one, double sqrt_ef_p2_minus_one, double asec_minus_ef)
{
    return 0.0;
}

inline double D_TABLE_FUNC4(double ep, double ep_p2, double ef, double ef_p2, double ef_p4, double one_div_ef_p1, double one_div_ef_p2, double one_div_ef_p3, double ef_p2_minus_one, double sqrt_ef_p2_minus_one, double asec_minus_ef)
{
    return c_2div3*one_div_ef_p2*ef_p2_minus_one*sqrt_ef_p2_minus_one;
}
inline double D_TABLE_FUNC4_DER(double ep, double ep_p2, double ef, double ef_p2, double ef_p4, double one_div_ef_p1, double one_div_ef_p2, double one_div_ef_p3, double ef_p2_minus_one, double sqrt_ef_p2_minus_one, double asec_minus_ef)
{
    return 0.0;
}

inline double D_TABLE_FUNC5(double ep, double ep_p2, double ef, double ef_p2, double ef_p4, double one_div_ef_p1, double one_div_ef_p2, double one_div_ef_p3, double ef_p2_minus_one, double sqrt_ef_p2_minus_one, double asec_minus_ef)
{
    return ep_p2*( c_1div3*one_div_ef_p2*sqrt_ef_p2_minus_one*(1.0 + 2.0*ef_p2) + asec_minus_ef);
}
inline double D_TABLE_FUNC5_DER(double ep, double ep_p2, double ef, double ef_p2, double ef_p4, double one_div_ef_p1, double one_div_ef_p2, double one_div_ef_p3, double ef_p2_minus_one, double sqrt_ef_p2_minus_one, double asec_minus_ef)
{
    return 2.0*ep*( c_1div3*one_div_ef_p2*sqrt_ef_p2_minus_one*(1.0 + 2.0*ef_p2) + asec_minus_ef);
}

inline double D_TABLE_FUNC6(double ep, double ep_p2, double ef, double ef_p2, double ef_p4, double one_div_ef_p1, double one_div_ef_p2, double one_div_ef_p3, double ef_p2_minus_one, double sqrt_ef_p2_minus_one, double asec_minus_ef)
{
    return c_2div3*one_div_ef_p2*ef_p2_minus_one*sqrt_ef_p2_minus_one;
}
inline double D_TABLE_FUNC6_DER(double ep, double ep_p2, double ef, double ef_p2, double ef_p4, double one_div_ef_p1, double one_div_ef_p2, double one_div_ef_p3, double ef_p2_minus_one, double sqrt_ef_p2_minus_one, double asec_minus_ef)
{
    return 0.0;
}

inline double D_TABLE_FUNC7(double ep, double ep_p2, double ef, double ef_p2, double ef_p4, double one_div_ef_p1, double one_div_ef_p2, double one_div_ef_p3, double ef_p2_minus_one, double sqrt_ef_p2_minus_one, double asec_minus_ef)
{
    return -c_1div3*one_div_ef_p2*sqrt_ef_p2_minus_one*(1.0 + 2.0*ef_p2) - asec_minus_ef;
}
inline double D_TABLE_FUNC7_DER(double ep, double ep_p2, double ef, double ef_p2, double ef_p4, double one_div_ef_p1, double one_div_ef_p2, double one_div_ef_p3, double ef_p2_minus_one, double sqrt_ef_p2_minus_one, double asec_minus_ef)
{
    return 0.0;
}
inline double D_TABLE_FUNC8(double ep, double ep_p2, double ef, double ef_p2, double ef_p4, double one_div_ef_p1, double one_div_ef_p2, double one_div_ef_p3, double ef_p2_minus_one, double sqrt_ef_p2_minus_one, double asec_minus_ef)
{
    return 2.0*(c_1div3*one_div_ef_p1*sqrt_ef_p2_minus_one*(1.0 + 2.0*ef_p2) + ef*asec_minus_ef);
}
inline double D_TABLE_FUNC8_DER(double ep, double ep_p2, double ef, double ef_p2, double ef_p4, double one_div_ef_p1, double one_div_ef_p2, double one_div_ef_p3, double ef_p2_minus_one, double sqrt_ef_p2_minus_one, double asec_minus_ef)
{
    return 0.0;
}
inline double D_TABLE_FUNC9(double ep, double ep_p2, double ef, double ef_p2, double ef_p4, double one_div_ef_p1, double one_div_ef_p2, double one_div_ef_p3, double ef_p2_minus_one, double sqrt_ef_p2_minus_one, double asec_minus_ef)
{
    return c_1div15*one_div_ef_p3*sqrt_ef_p2_minus_one*(2.0 - 9.0*ef_p2 - 8.0*ef_p4) - ef*asec_minus_ef;
}
inline double D_TABLE_FUNC9_DER(double ep, double ep_p2, double ef, double ef_p2, double ef_p4, double one_div_ef_p1, double one_div_ef_p2, double one_div_ef_p3, double ef_p2_minus_one, double sqrt_ef_p2_minus_one, double asec_minus_ef)
{
    return 0.0;
}
inline double D_TABLE_FUNC10(double ep, double ep_p2, double ef, double ef_p2, double ef_p4, double one_div_ef_p1, double one_div_ef_p2, double one_div_ef_p3, double ef_p2_minus_one, double sqrt_ef_p2_minus_one, double asec_minus_ef)
{
    return (1.0 - ep_p2)*( c_1div30*one_div_ef_p3*sqrt_ef_p2_minus_one*(-2.0 + 9.0*ef_p2 + 8.0*ef_p4) + c_1div2*ef*asec_minus_ef);
}
inline double D_TABLE_FUNC10_DER(double ep, double ep_p2, double ef, double ef_p2, double ef_p4, double one_div_ef_p1, double one_div_ef_p2, double one_div_ef_p3, double ef_p2_minus_one, double sqrt_ef_p2_minus_one, double asec_minus_ef)
{
    return -2.0*ep*( c_1div30*one_div_ef_p3*sqrt_ef_p2_minus_one*(-2.0 + 9.0*ef_p2 + 8.0*ef_p4) + c_1div2*ef*asec_minus_ef);
}
inline double D_TABLE_FUNC11(double ep, double ep_p2, double ef, double ef_p2, double ef_p4, double one_div_ef_p1, double one_div_ef_p2, double one_div_ef_p3, double ef_p2_minus_one, double sqrt_ef_p2_minus_one, double asec_minus_ef)
{
    return c_4div15*one_div_ef_p3*ef_p2_minus_one*ef_p2_minus_one*sqrt_ef_p2_minus_one;
}
inline double D_TABLE_FUNC11_DER(double ep, double ep_p2, double ef, double ef_p2, double ef_p4, double one_div_ef_p1, double one_div_ef_p2, double one_div_ef_p3, double ef_p2_minus_one, double sqrt_ef_p2_minus_one, double asec_minus_ef)
{
    return 0.0;
}
inline double D_TABLE_FUNC12(double ep, double ep_p2, double ef, double ef_p2, double ef_p4, double one_div_ef_p1, double one_div_ef_p2, double one_div_ef_p3, double ef_p2_minus_one, double sqrt_ef_p2_minus_one, double asec_minus_ef)
{
    return c_1div30*one_div_ef_p3*sqrt_ef_p2_minus_one*(2.0 - 9.0*ef_p2 - 8.0*ef_p4) - c_1div2*ef*asec_minus_ef;
}
inline double D_TABLE_FUNC12_DER(double ep, double ep_p2, double ef, double ef_p2, double ef_p4, double one_div_ef_p1, double one_div_ef_p2, double one_div_ef_p3, double ef_p2_minus_one, double sqrt_ef_p2_minus_one, double asec_minus_ef)
{
    return 0.0;
}
inline double D_TABLE_FUNC13(double ep, double ep_p2, double ef, double ef_p2, double ef_p4, double one_div_ef_p1, double one_div_ef_p2, double one_div_ef_p3, double ef_p2_minus_one, double sqrt_ef_p2_minus_one, double asec_minus_ef)
{
    return ep_p2*( c_1div10*one_div_ef_p3*sqrt_ef_p2_minus_one*(-2.0 + 9.0*ef_p2 + 8.0*ef_p4) + c_3div2*ef*asec_minus_ef );
}
inline double D_TABLE_FUNC13_DER(double ep, double ep_p2, double ef, double ef_p2, double ef_p4, double one_div_ef_p1, double one_div_ef_p2, double one_div_ef_p3, double ef_p2_minus_one, double sqrt_ef_p2_minus_one, double asec_minus_ef)
{
    return 2.0*ep*( c_1div10*one_div_ef_p3*sqrt_ef_p2_minus_one*(-2.0 + 9.0*ef_p2 + 8.0*ef_p4) + c_3div2*ef*asec_minus_ef );
}
inline double D_TABLE_FUNC14(double ep, double ep_p2, double ef, double ef_p2, double ef_p4, double one_div_ef_p1, double one_div_ef_p2, double one_div_ef_p3, double ef_p2_minus_one, double sqrt_ef_p2_minus_one, double asec_minus_ef)
{
    return c_1div10*one_div_ef_p3*sqrt_ef_p2_minus_one*(2.0 - 9.0*ef_p2 - 8.0*ef_p4) - c_3div2*ef*asec_minus_ef;
}
inline double D_TABLE_FUNC14_DER(double ep, double ep_p2, double ef, double ef_p2, double ef_p4, double one_div_ef_p1, double one_div_ef_p2, double one_div_ef_p3, double ef_p2_minus_one, double sqrt_ef_p2_minus_one, double asec_minus_ef)
{
    return 0.0;
}
inline double D_TABLE_FUNC15(double ep, double ep_p2, double ef, double ef_p2, double ef_p4, double one_div_ef_p1, double one_div_ef_p2, double one_div_ef_p3, double ef_p2_minus_one, double sqrt_ef_p2_minus_one, double asec_minus_ef)
{
    return c_4div15*one_div_ef_p3*ef_p2_minus_one*ef_p2_minus_one*sqrt_ef_p2_minus_one;
}
inline double D_TABLE_FUNC15_DER(double ep, double ep_p2, double ef, double ef_p2, double ef_p4, double one_div_ef_p1, double one_div_ef_p2, double one_div_ef_p3, double ef_p2_minus_one, double sqrt_ef_p2_minus_one, double asec_minus_ef)
{
    return 0.0;
}

#endif

/*	ODE solver macros	*/
#ifndef __ODE_MACROS
#define __ODE_MACROS
    #define Ith(v,i)    NV_Ith_S(v,i-1)       		/* Ith numbers components 1..NEQ */
    #define IJth(A,i,j) DENSE_ELEM(A,i-1,j-1) 		/* IJth numbers rows,cols 1..NEQ */
    #ifndef max
        #define max( a, b ) ( ((a) > (b)) ? (a) : (b) )
    #endif
    #ifndef min
        #define min(X,Y) ((X) < (Y) ? (X) : (Y))
    #endif
#endif

/* vector operators */
#ifndef __VECTOR_OPERATORS
#define __VECTOR_OPERATORS
inline void cross3(double a[3], double b[3], double result[3])
{
    result[0] = a[1]*b[2] - a[2]*b[1];
    result[1] = a[2]*b[0] - a[0]*b[2];
    result[2] = a[0]*b[1] - a[1]*b[0];
}
inline double norm3(double v[3])
{
    double result = sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
    return result;
}
inline double norm3_squared(double v[3])
{
    double result = v[0]*v[0] + v[1]*v[1] + v[2]*v[2];
    return result;
}
inline double dot3(double a[3], double b[3])
{
    double result = (a[0]*b[0] + a[1]*b[1] + a[2]*b[2]);
    return result;
}
inline void matrix_on_vector( int N_x, int N_y, double M[][3], double v[], double result[] )
{
    
    int i,j;

    for (j=0; j<N_y; j++)
    {
        result[j] = 0.0;
    }

    for (j=0; j<N_y; j++)
    {
        for (i=0; i<N_x; i++)
        {
            result[j] += M[i][j]*v[i];
        }
    }
}

inline void scalar_times_matrix(int N_x, int N_y, double M[][3], double common_factor)
{
    int i,j;
    for (i=0; i<N_x; i++)
    {
        for (j=0; j<N_y; j++)
        {
            M[i][j] *= common_factor;
        }
    }
}
#endif


/* sse interface */
#ifndef __STEV
#define __STEV

struct value1__ 
{
        double neta,bwind,hewind,mxns;
};

struct value2__ 
{
        double alpha1,lambda;
};

struct value3__ 
{
        int idum;
};

struct value4__ 
{
        double sigma;
        int bhflag;
};

struct value5__ 
{
        double beta,xi,acc2,epsnov,eddfac,gamma;
};

struct flags__ 
{
        int ceflag,tflag,ifflag,nsflag,wdflag;
};

struct points__ 
{
        double pts1,pts2,pts3;
};

extern "C" void evolv1_(int *kw, double *mass,double *mt, double *r, double *lum, double *mc, double *rc, double *menv, double *renv, double *ospin, double *epoch, double *tms, double *tphys, double *tphysf, double *dtp, double *z, double *zpars);
extern "C" void zcnsts_(double* z, double *zpars);
extern "C" void star_(int *kw, double *mass, double *mt, double *tm, double *tn, double *tscls, double *lums, double *GB, double *zpars);
extern "C" void deltat_(int *kw, double *age, double *tm, double *tn, double *tscls, double *dt, double *dtr);
extern "C" void hrdiag_(double *mass, double *aj, double *mt, double *tm, double *tn, double * tscls, double *lums, double *GB, double *zpars,
    double *r, double *lum, int *kw, double *mc, double *rc, double *menv, double *renv, double *k2);
extern "C" double vrotf_(double mt);    
#endif



/* classes */
#ifndef __Particle
#define __Particle
class Particle
{
    public:

    /* generic properties */
    int index,child1,child2;
    int parent;
    int sibling;
    std::vector<int> parents;
    std::vector<int> connecting_child_in_parents;
    int level,highest_level;
    int is_binary;
    double mass,child1_mass,child2_mass,total_system_mass;

    /*******************
    /* body properties *
     * ****************/

    /* general */
    double radius;
    double spin_vec_x,spin_vec_y,spin_vec_z;
    double spin_vec_x_dot,spin_vec_y_dot,spin_vec_z_dot;
    
    double x,y,z;
    double vx,vy,vz;

    /* stellar evolution */
    int stellar_type;
    int evolve_as_star;
    double sse_initial_mass,sse_time_step;
    double mass_dot_wind,radius_dot,radius_ddot,ospin_dot;
    double epoch,age;
    double metallicity;
    double *zpars;
    double convective_envelope_mass,convective_envelope_radius;
    double core_mass,core_radius;
    double luminosity;
    double child1_mass_old,child2_mass_old;
    double apsidal_motion_constant, gyration_radius;
    
    /* RLOF */
    int RLOF_flag; /* 0: not in RLOF; 1: in RLOF */
    double mass_dot_RLOF;
    int emt_ejection_radius_mode;
    double emt_accretion_radius;
    double emt_tau;
    int RLOF_at_pericentre_has_occurred_entering_RLOF;
    double fm;
    
    /* kicks */
    int kick_distribution;
    double kick_distribution_sigma;

    /* used in ODE solver only */
    double spin_vec[3],dspin_vec_dt[3]; 
    double spin_vec_norm;
    double dmass_dt,dradius_dt;    
    
    void set_ODE_quantities(double delta_time);
    void reset_ODE_quantities();
    
    /*********************
    /* binary properties *
     * ******************/

    /* general */
    double e_vec_x,e_vec_y,e_vec_z;
    double h_vec_x,h_vec_y,h_vec_z;

    double e_vec_x_dot,e_vec_y_dot,e_vec_z_dot;
    double h_vec_x_dot,h_vec_y_dot,h_vec_z_dot;
   
    double true_anomaly;
    
    /* PN terms */
    int include_pairwise_1PN_terms,include_pairwise_25PN_terms;
        
    /* tidal friction */
    int include_tidal_friction_terms,tides_method,include_tidal_bulges_precession_terms,include_rotation_precession_terms;
    double minimum_eccentricity_for_tidal_precession;
    
    double tides_viscous_time_scale;
    int tides_viscous_time_scale_prescription;

    /* root finding */
    int check_for_secular_breakdown,secular_breakdown_has_occurred;
    int check_for_dynamical_instability,dynamical_instability_has_occurred,dynamical_instability_criterion;
    int dynamical_instability_central_particle;
    double dynamical_instability_K_parameter;
    int check_for_physical_collision_or_orbit_crossing,physical_collision_or_orbit_crossing_has_occurred;
    int check_for_minimum_periapse_distance,minimum_periapse_distance_has_occurred;
    double check_for_minimum_periapse_distance_value;
    int check_for_RLOF_at_pericentre,check_for_RLOF_at_pericentre_use_sepinsky_fit,RLOF_at_pericentre_has_occurred;
    int check_for_GW_condition,GW_condition_has_occurred;    

    /* used in ODE solver only */
    double e_vec[3],h_vec[3];
    double e_vec_unit[3],h_vec_unit[3],q_vec_unit[3];
    double de_vec_dt[3],dh_vec_dt[3];
    double child1_mass_plus_child2_mass,child1_mass_minus_child2_mass,child1_mass_times_child2_mass;
    double child1_mass_dot_wind,child2_mass_dot_wind;
    double e,e_p2;
    double j,j_p2,j_p3,j_p4,j_p5; // j=sqrt(1-e^2)
    double h,a;
    
//    Particle(int index, int is_binary, int child1, int child2, double mass, double radius, double spin_vec_x, double spin_vec_y, double spin_vec_z, double e_vec_x, double e_vec_y, double e_vec_z, double h_vec_x, double h_vec_y, double h_vec_z, int include_pairwise_1PN_terms, int include_pairwise_25PN_terms, int include_tides_terms, double tides_Q_prime, double tides_gyration_radius, int check_for_secular_breakdown, int secular_breakdown_has_occurred, int check_for_dynamical_instability, int dynamical_instability_has_occurred, int dynamical_instability_criterion, int check_for_physical_collision, int physical_collision_has_occurred) : index(index), is_binary(is_binary), child1(child1), child2(child2), mass(mass), radius(radius), spin_vec_x(spin_vec_x), spin_vec_y(spin_vec_y), spin_vec_z(spin_vec_z), e_vec_x(e_vec_x), e_vec_y(e_vec_y), e_vec_z(e_vec_z), h_vec_x(h_vec_x), h_vec_y(h_vec_y), h_vec_z(h_vec_z), include_pairwise_1PN_terms(include_pairwise_1PN_terms), include_pairwise_25PN_terms(include_pairwise_25PN_terms), include_tides_terms(include_tides_terms), tides_Q_prime(tides_Q_prime), tides_gyration_radius(tides_gyration_radius), 

    /* user-specified instantaneous perturbations */
    int sample_orbital_phase_randomly;
    double instantaneous_perturbation_delta_mass;
    double instantaneous_perturbation_delta_x,instantaneous_perturbation_delta_y,instantaneous_perturbation_delta_z;
    double instantaneous_perturbation_delta_vx,instantaneous_perturbation_delta_vy,instantaneous_perturbation_delta_vz;

    int is_external;
    double external_t_ref,external_e,external_r_p;

    /* VRR */
    int VRR_model;
    //int VRR_include_GR_precession;
    int VRR_include_mass_precession;
    //int VRR_include_frame_dragging;
    double VRR_mass_precession_rate;
    double VRR_initial_time, VRR_final_time;
    
    /* Dipole model */
    double VRR_Omega_vec_x,VRR_Omega_vec_y,VRR_Omega_vec_z;
    
    /* Distortion model */
    //double VRR_dipole_fraction;
    //double VRR_initial_distortion_axis_vec_x,VRR_initial_distortion_axis_vec_y,VRR_initial_distortion_axis_vec_z;
    //double VRR_final_distortion_axis_vec_x,VRR_final_distortion_axis_vec_y,VRR_final_distortion_axis_vec_z;
    //double VRR_AD_parameter,VRR_AM_parameter,VRR_AJ_parameter;
    
    /* Bar-Or model */
    double VRR_eta_20_init,VRR_eta_a_22_init,VRR_eta_b_22_init,VRR_eta_a_21_init,VRR_eta_b_21_init;
    double VRR_eta_20_final,VRR_eta_a_22_final,VRR_eta_b_22_final,VRR_eta_a_21_final,VRR_eta_b_21_final;
    


    Particle(int index, int is_binary) : index(index), is_binary(is_binary)
    {
        /* default values */
        check_for_secular_breakdown = 0;
        check_for_dynamical_instability = 0;
        check_for_physical_collision_or_orbit_crossing = 1;
        check_for_minimum_periapse_distance = 0;
        check_for_RLOF_at_pericentre = 0;
        check_for_RLOF_at_pericentre_use_sepinsky_fit = 0; /* 0: use Eggleton formula for consistency with emt model */

        secular_breakdown_has_occurred = 0;
        dynamical_instability_has_occurred = 0;
        physical_collision_or_orbit_crossing_has_occurred = 0;
        minimum_periapse_distance_has_occurred = 0;
        RLOF_at_pericentre_has_occurred = 0;
       
        include_pairwise_1PN_terms = 1;
        include_pairwise_25PN_terms = 1;
        include_tidal_friction_terms = 1;
        include_tidal_bulges_precession_terms = 1;
        include_rotation_precession_terms = 1;

        radius = 1.0e-10; /* this must be set (to nonzero), otherwise ODE solver will have invalid ewt values */
        tides_method = 1; /* `full' tides equations of motion including spin-orbit terms for all orientations */
        tides_viscous_time_scale_prescription = 1; /* 0: constant, user-specified t_V; 1: using Hurley prescription */
        minimum_eccentricity_for_tidal_precession = 1.0e-3;
        spin_vec_x_dot = spin_vec_y_dot = spin_vec_z_dot = 0.0;

        /* stellar evolution */
        stellar_type = 1;
        evolve_as_star = 1;
        sse_initial_mass = 0.0;
        sse_time_step = 1.0;
        metallicity = 0.02;
        epoch = 0.0;
        age = 0.0;
        convective_envelope_mass = 0.0;
        convective_envelope_radius = 0.0;
        core_mass = 0.0;
        core_radius = 0.0;
        luminosity = 0.0;
        apsidal_motion_constant = 0.19;
        gyration_radius = 0.08;
        mass_dot_wind = radius_dot = radius_ddot = ospin_dot = 0.0;
                
        /* RLOF */
        RLOF_flag = 0;
        mass_dot_RLOF = 0.0;
        emt_ejection_radius_mode = 0;
        emt_accretion_radius = 0.0;
        emt_tau = 0.0;
        
        
        /* kicks */
        kick_distribution = 0;
        kick_distribution_sigma = 265.0; /* https://ui.adsabs.harvard.edu/abs/2005MNRAS.360..974H/abstract */
         
        sample_orbital_phase_randomly = 1;
        instantaneous_perturbation_delta_mass = 0.0;
        instantaneous_perturbation_delta_x = instantaneous_perturbation_delta_y = instantaneous_perturbation_delta_z = 0.0;
        instantaneous_perturbation_delta_vx = instantaneous_perturbation_delta_vy = instantaneous_perturbation_delta_vz = 0.0;
        
        spin_vec_x = 0.0;
        spin_vec_y = 0.0;
        spin_vec_z = 1.0e-10;
        
        e_vec_x_dot = e_vec_y_dot = e_vec_z_dot = 0.0;
        h_vec_x_dot = h_vec_y_dot = h_vec_z_dot = 0.0;
     
        is_external = 0;
        external_t_ref = 0.0;
        external_r_p = 1.0;
        external_e = 10.0;
        
        /* VRR */
        VRR_model = 0;
        //VRR_include_GR_precession = 0;
        VRR_include_mass_precession = 0;
        //VRR_include_frame_dragging = 0;
        VRR_mass_precession_rate = 0.0;
        VRR_initial_time = 0.0;
        VRR_final_time = 1.0;
        
        VRR_Omega_vec_x = VRR_Omega_vec_y = VRR_Omega_vec_z = 0.0;
        //VRR_dipole_fraction = 0.5;
        //VRR_initial_distortion_axis_vec_x = 0.0;
        //VRR_initial_distortion_axis_vec_y = 0.0;
        //VRR_initial_distortion_axis_vec_z = 0.0;
        //VRR_final_distortion_axis_vec_x = 0.0;
        //VRR_final_distortion_axis_vec_y = 0.0;
        //VRR_final_distortion_axis_vec_z = 0.0;
        //VRR_AD_parameter = 0.0;
        //VRR_AM_parameter = 0.0;
        //VRR_AJ_parameter = 0.0;
        
        VRR_eta_20_init = 0.0;
        VRR_eta_a_22_init = 0.0;
        VRR_eta_b_22_init = 0.0;
        VRR_eta_a_21_init = 0.0;
        VRR_eta_b_21_init = 0.0;

        VRR_eta_20_final = 0.0;
        VRR_eta_a_22_final = 0.0;
        VRR_eta_b_22_final = 0.0;
        VRR_eta_a_21_final = 0.0;
        VRR_eta_b_21_final = 0.0;

    }
};

typedef std::map<int, Particle *> ParticlesMap;
typedef std::map<int, Particle *>::iterator ParticlesMapIterator;

extern int highest_particle_index;
extern ParticlesMap particlesMap;
#endif

/* CVODE UserData */
#ifndef __UserData
#define __UserData
typedef struct {
	ParticlesMap *particlesMap;
    double hamiltonian;
    int N_root_finding;
    double start_time;
} *UserData;
#endif

}
