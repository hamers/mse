#include "constants.h"
#include <math.h>

double CONST_G = 4.0*M_PI*M_PI; 
double CONST_G_P2 = CONST_G*CONST_G;
double CONST_G_P3 = CONST_G_P2*CONST_G;
double CONST_C_LIGHT = 63239.72638679138;
double CONST_C_LIGHT_P2 = CONST_C_LIGHT*CONST_C_LIGHT;
double CONST_C_LIGHT_P3 = CONST_C_LIGHT_P2*CONST_C_LIGHT;
double CONST_C_LIGHT_P4 = CONST_C_LIGHT_P3*CONST_C_LIGHT;
double CONST_C_LIGHT_P5 = CONST_C_LIGHT_P4*CONST_C_LIGHT;
double CONST_C_LIGHT_P6 = CONST_C_LIGHT_P5*CONST_C_LIGHT;
double CONST_MSUN = 1.0;
double CONST_R_SUN = 0.004649130343817401;
double CONST_L_SUN = 0.0002710404109745588;
double CONST_KM_PER_S = 0.210862;
double CONST_PER_PC3 = 1.14059e-16;

/* Used in MSTAR only */
double SPEEDOFLIGHT = CONST_C_LIGHT;
double GCONST = CONST_G;

double MSTAR_maximum_separation_for_inclusion = 1.0e5;

bool MSTAR_include_PN_acc_10 = true;
bool MSTAR_include_PN_acc_20 = true;
bool MSTAR_include_PN_acc_25 = true;
bool MSTAR_include_PN_acc_30 = true;
bool MSTAR_include_PN_acc_35 = true;

bool MSTAR_include_PN_acc_SO = true;
bool MSTAR_include_PN_acc_SS = true;
bool MSTAR_include_PN_acc_Q = true;

bool MSTAR_include_PN_spin_SO = true;
bool MSTAR_include_PN_spin_SS = true;
bool MSTAR_include_PN_spin_Q = true;
