#ifndef __CONSTANTS
#define __CONSTANTS

extern double CONST_G;
extern double CONST_G_P2;
extern double CONST_G_P3;
extern double CONST_C_LIGHT;
extern double CONST_C_LIGHT_P2;
extern double CONST_C_LIGHT_P3;
extern double CONST_C_LIGHT_P4;
extern double CONST_C_LIGHT_P5;
extern double CONST_C_LIGHT_P6;
extern double CONST_MSUN;
extern double CONST_R_SUN;
extern double CONST_L_SUN;
extern double CONST_KM_PER_S;
extern double CONST_PER_PC3;

/* Used in MSTAR only */
extern double SPEEDOFLIGHT;
extern double GCONST;
extern double MSTAR_maximum_separation_for_inclusion;

#include <stdbool.h>
extern bool MSTAR_include_PN_acc_10;
extern bool MSTAR_include_PN_acc_20;
extern bool MSTAR_include_PN_acc_25;
extern bool MSTAR_include_PN_acc_30;
extern bool MSTAR_include_PN_acc_35;

extern bool MSTAR_include_PN_acc_SO;
extern bool MSTAR_include_PN_acc_SS;
extern bool MSTAR_include_PN_acc_Q;

extern bool MSTAR_include_PN_spin_SO;
extern bool MSTAR_include_PN_spin_SS;
extern bool MSTAR_include_PN_spin_Q;

#endif
