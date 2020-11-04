#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <float.h>
#include <limits.h>
#include "regularization.h"

#ifdef USE_PN

// vector products
double cross_product( double a[3], double b[3], double result[3] ){
        result[0] =    a[1]*b[2] - a[2]*b[1];
        result[1] = -( a[0]*b[2] - a[2]*b[0] );
        result[2] =    a[0]*b[1] - a[1]*b[0];
}

// relative velocity of particle 1 seen from particle 2 = v_1 - v_2. return v_rel^2.
inline static double relative_velocity(double *v1, double *v2,  double *v12) {
    double vrel2 = 0;
    for(int k=0;k<3;k++){
        v12[k] = v1[k] - v2[k];
        vrel2 += v12[k]*v12[k];
    }
    return vrel2;
}

#ifdef USE_PN_SPIN
void compute_PN_twobody_Acc(struct RegularizedRegion *R, double *Vel, double *AccPN, double *Spin, double *dSpin_PN, int i, int j ){
#else
void compute_PN_twobody_Acc(struct RegularizedRegion *R, double *Vel, double *AccPN, int i, int j ){
#endif

    // speed of light constants
    const double c  = SPEEDOFLIGHT;
    const double c2 = c*c;
    const double c3 = c*c2;
    const double c4 = c2*c2;
    const double c5 = c*c4;
    const double c6 = c3*c3;
    const double c7 = c*c6;
    const double pi2 = M_PI*M_PI;

    // mass constants
    const double m1   = R->Mass[i];
    const double m2   = R->Mass[j];
    const double m    = m1+m2 ;
    const double eta  = m1*m2/(m*m);
    const double eta2 = eta*eta;
    const double eta3 = eta2*eta;
    const double G    = GCONST;
    const double gm   = G*m;

    // relative separation & velocity & unit separation
    double velvec_i[3], velvec_j[3];
    for(int k=0;k<3;k++){
	velvec_i[k] = Vel[3*i+k];
	velvec_j[k] = Vel[3*j+k];
    }
    double rv[3] = {0};
    double vv[3] = {0};
    double nv[3] = {0};
    // rdot = (rv dot vv) / r
    double rdot = 0;
    // powers of r, v, rdot and G(m1+m2)/r
    double r2 = 0, r, r3, r4;
    double rdot2, rdot4, rdot6;
    double v2=0, v4=0, v6=0;
    double gmpr=0, gmpr2=0, gmpr3=0;

    // spin contributions to spin derivatives
    double Omega_SO[3] = {0};
    double Omega_SS[3] = {0};
    double Omega_Q[3]  = {0};
    double dSpin1[3]   = {0};
    double dSpin2[3]   = {0};

    // the PN acceleration coefficients
    double A1=0, A2=0, A2p5=0, A3=0, A3p5=0, Atot = 0;
    double B1=0, B2=0, B2p5=0, B3=0, B3p5=0, Btot = 0;

    // calculate relative quantities and derived quantities.
    // every relative quantity x = x_1 - x_2, where x_1 corresponds to
    // index i and x_2 to index j
    const int Nd = GLOBAL_ND;
    int path[Nd]; int sign[Nd]; int d;
    int proximity = check_relative_proximity(i,j,Nd,R,&d,path,sign);
    for(int k=0;k<3;k++){
    	if(proximity==1){
        	for(int l=0;l<d;l++){
                	int index = path[l];
                        rv[k] += sign[l]*R->State[3*index+k];
                }
        } else{
                rv[k] = R->Pos[3*i+k] - R->Pos[3*j+k];
        }
    }

    // relative velocity is also v_1 - v_2
    v2 = relative_velocity(velvec_i, velvec_j, vv);

    // expressions using relative separation, velocity and mass
    for (int k = 0; k < 3; k++) {
        rdot += rv[k]*vv[k];
        r2 += rv[k]*rv[k];
    }
    r = sqrt(r2);
    r3 = r2 * r;
    r4 = r2*r2;
    rdot /= r;
    rdot2 = rdot*rdot;
    rdot4 = rdot2*rdot2;
    rdot6 = rdot2*rdot4;
    v4 = v2*v2;
    v6 = v4*v2;
    gmpr = gm/r;
    gmpr2 = gmpr * gmpr;
    gmpr3 = gmpr * gmpr2;
    for(int k=0;k<3;k++){
        nv[k] = rv[k]/r;
    }

    double pn_acc[3] = {0}, pn_spin_acc[3] = {0};

    // spin-independent terms
    // conservative terms: 1.0 PN, 2.0 PN, 3.0 PN
    // radiative terms:    2.5 PN, 3.5 PN

    // 1.0 PN
    if( MSTAR_include_PN_acc_10 == true ){
    	A1 = 2.0*(2.0+eta)*gmpr - (1+3*eta)*v2 + 1.5e0*eta*rdot2;
    	B1 = 2.0*(2.0-eta)*rdot;
    }

    // 2.0 PN
    if( MSTAR_include_PN_acc_20 == true ){
    	A2 = -(3.0/4)*(12+29*eta)*gmpr2 - eta*(3-4*eta)*v4-(15./8)*eta*(1-3*eta)*rdot4 + (1./2)*eta*(13-4*eta)*gmpr*v2+(2+25*eta+2*eta2)*gmpr*rdot2 + (3./2)*eta*(3-4*eta)*v2*rdot2;
    	B2 = -(1.0/2)*rdot*((4+41*eta+8*eta2)*gmpr - eta*(15+4*eta)*v2 + 3*eta*(3+2*eta)*rdot2);
    }

    // 2.5 PN
    if( MSTAR_include_PN_acc_25 == true ){
    	A2p5=8.0/5.0*eta*gmpr*rdot*(17./3*gmpr+3*v2);
    	B2p5=-8.0/5.0*eta*gmpr*(3*gmpr+v2);
    }

    // 3.0 PN
    if( MSTAR_include_PN_acc_30 == true ){
    	A3=(16+(1399./12-41./16*pi2)*eta+71./2*eta2)*gmpr3+eta*(20827./840+123./64*pi2-eta2)*gmpr2*v2-(1+(22717./168+615./64*pi2)*eta+11./8*eta2-7*eta3)*gmpr2*rdot2-.25*eta*(11-49*eta+52*eta2)*v6+35./16*eta*(1-5*eta+5*eta2)*rdot6-.25*eta*(75+32*eta-40*eta2)*gmpr*v2*v2-.5*eta*(158-69*eta-60*eta2)*gmpr*rdot4+eta*(121-16*eta-20*eta2)*gmpr*v2*rdot2+3./8*eta*(20-79*eta+60*eta2)*v2*v2*rdot2-15./8*eta*(4-18*eta+17*eta2)*v2*rdot4;
    	B3=rdot*((4+(5849./840.+123./32.*pi2)*eta-25*eta2-8*eta3)*gmpr2+1./8.*eta*(65-152*eta-48*eta2)*v2*v2+15/8.*eta*(3-8*eta-2*eta2)*rdot4+eta*(15+27*eta+10*eta2)*gmpr*v2-1./6.*eta*(329+177*eta+108*eta2)*gmpr*rdot2-.75*eta*(16-37*eta-16*eta2)*v2*rdot2);
    }

    // 3.5 PN
    if( MSTAR_include_PN_acc_35 == true ){
    	A3p5=-8./5*eta*gmpr*rdot*(23./14*(43+14*eta)*gmpr2+3./28*(61+70*eta)*v2*v2+70*rdot4+1./42*(519-1267*eta)*gmpr*v2+.25*(147+188*eta)*gmpr*rdot2-15/4.*(19+2*eta)*v2*rdot2);
    	B3p5=8./5.*eta*gmpr*(1./42.*(1325+546*eta)*gmpr2+1./28.*(313+42*eta)*v2*v2+75*rdot4-1./42.*(205+777*eta)*gmpr*v2+1./12.*(205+424*eta)*gmpr*rdot2-.75*(113+2*eta)*v2*rdot2);
    }

    // sum all acceleration terms with correct powers of c
    Atot = A1/c2 + A2/c4 + A2p5/c5 + A3/c6 + A3p5/c7;
    Btot = B1/c2 + B2/c4 + B2p5/c5 + B3/c6 + B3p5/c7;
    for(int k=0;k<3;k++){
    	pn_acc[k] = gm/r2*(nv[k]*Atot  + vv[k]*Btot);
    }

#ifdef USE_PN_SPIN

    // spin-dependent terms
    // 1.0 PN: Omega_SO		// 1.5 PN: A_SO, Omega_SS, Omega_Q		// 2.0 PN: A_SS, A_Q
    // SO: spin-orbit, SS: spin-spin, Q: quadrupole

    // reduced mass, spin terms and products
    double mu           = (m1*m2)/(m1+m2);
    double sigma[3]     = {0};
    double S[3]         = {0};
    double n_cross_v[3] = {0};
    double S1_2         = 0;
    double S2_2         = 0;
    double S1_dot_S2    = 0;
    double n_dot_S1     = 0;
    double n_dot_S2     = 0;

    cross_product( nv, vv, n_cross_v );

    double *S_1 = &Spin[3*i];
    double *S_2 = &Spin[3*j];

    // spin supplementary quantities
    for(int k=0;k<3;k++){
	S[k]        = S_1[k] + S_2[k];
	sigma[k]    = m2/m1*S_1[k] + m1/m2*S_2[k];
	S1_dot_S2  += S_1[k]*S_2[k];
	n_dot_S1   += nv[k]*S_1[k];
	n_dot_S2   += nv[k]*S_2[k];
	S1_2       += S_1[k]*S_1[k];
        S2_2       += S_2[k]*S_2[k];

    }

    // spin contributions to acceleration: SS, S0 & Q
    double buffer1[3]    = {0};
    double buffer2[3]    = {0};
    double A_SO_term1    = {0};
    double A_SO_term2[3] = {0};
    double A_SO_term3[3] = {0};

    if( MSTAR_include_PN_acc_SO == true ){
    	for(int k=0;k<3;k++){
		A_SO_term1    += 6.0 * n_cross_v[k] * (S[k] + sigma[k]);
        	buffer1[k]     = 4.0 * S[k] + 3.0*sigma[k];
        	buffer2[k]     = 3.0 * rdot * (2.0*S[k] + sigma[k]);
    	}
    	cross_product( vv, buffer1, A_SO_term2 );
    	cross_product( nv, buffer2, A_SO_term3 );
    }

    if( MSTAR_include_PN_acc_SO == true ){
    	for(int k=0;k<3;k++){
    		pn_spin_acc[k] += (G/(c2 * r3)) * ( nv[k] * A_SO_term1 - A_SO_term2[k] + A_SO_term3[k] );
    	}
    }
    if( MSTAR_include_PN_acc_SS == true ){
    	for(int k=0; k<3;++k) {
		pn_spin_acc[k] += -3 * G / (mu * c2 * r4) * (S1_dot_S2 * nv[k] - 5 * n_dot_S1 * n_dot_S2 * nv[k] + n_dot_S1 * S_2[k] + n_dot_S2 * S_1[k]);
	}
    }
    if( MSTAR_include_PN_acc_Q == true ){
    	for(int k=0; k<3;++k) {
		pn_spin_acc[k] += -3./2*G*m/(r4 * c2) * ((1.0/(m1*m1)*((S1_2 - 5*n_dot_S1 * n_dot_S1) * nv[k] + 2*n_dot_S1 * S_1[k])) + (1.0/(m2*m2) * ((S2_2-5* n_dot_S2*n_dot_S2) * nv[k] + 2*n_dot_S2 * S_2[k])));
    	}
    }
    if( MSTAR_include_PN_acc_SO == true || MSTAR_include_PN_acc_SS == true || MSTAR_include_PN_acc_Q == true ){
    	for(int k=0;k<3;k++){
    		pn_acc[k] += pn_spin_acc[k];
    	}
    }

    // spin derivatives
    for(int k=0;k<3;k++){
	if( MSTAR_include_PN_spin_SO == true ){
		Omega_SO[k] = 2*G*mu/(c2*r2) * (1 + 3*m2/(4*m1)) * n_cross_v[k];
	}
	if( MSTAR_include_PN_spin_SS == true ){
		Omega_SS[k] = G/(c2*r3) * ( 3*n_dot_S2 * nv[k] - S_2[k] );
	}
	if( MSTAR_include_PN_spin_Q == true ){
		Omega_Q[k]  = G*m2/(c2*r3*m1) * ( 3*n_dot_S1 * nv[k] - S_1[k] );
	}
	buffer1[k] = Omega_SO[k] + Omega_SS[k] + Omega_Q[k];

	if( MSTAR_include_PN_spin_SO == true ){
		Omega_SO[k] = 2*G*mu/(c2*r2) * (1 + 3*m1/(4*m2)) * n_cross_v[k];
	}
	if( MSTAR_include_PN_spin_SS == true ){
		Omega_SS[k] = G/(c2*r3) * ( 3*n_dot_S1 * nv[k] - S_1[k] );
	}
	if( MSTAR_include_PN_spin_Q == true ){
		Omega_Q[k]  = G*m1/(c2*r3*m2) * ( 3*n_dot_S2 * nv[k] - S_2[k] );
	}
	buffer2[k] = Omega_SO[k] + Omega_SS[k] + Omega_Q[k];
    }

    cross_product( buffer1, S_1, dSpin1 );
    cross_product( buffer2, S_2, dSpin2 );

    // final spin derivatives
    if( MSTAR_include_PN_spin_SO == true || MSTAR_include_PN_spin_SS == true || MSTAR_include_PN_spin_Q == true ){
    	for(int k=0;k<3;k++){
		dSpin_PN[3*i+k] += dSpin1[k];
        	dSpin_PN[3*j+k] += dSpin2[k];
    	}
    }

#endif

    // final PN accelerations
    for(int k=0;k<3;k++){
        AccPN[3*i+k] += (m2/m)*pn_acc[k];
        AccPN[3*j+k] -= (m1/m)*pn_acc[k];
    }

}

#ifdef USE_PN_SPIN
void compute_Post_Newtonian_Acc(struct RegularizedRegion *R, double *Vel, double *Spin ){
#else
void compute_Post_Newtonian_Acc(struct RegularizedRegion *R, double *Vel ){
#endif

	for(int i=0;i<3*R->NumVertex;i++){
		R->AccPN[i]      = 0.0;
#ifdef USE_PN_SPIN
		R->Spin_dS_PN[i] = 0.0;
#endif
	}

	if(R->NumVertex<2){
		return;
	}

	for(int i=0;i<R->NumVertex;i++){
		for(int j=i+1;j<R->NumVertex;j++){
#ifdef USE_PN_SPIN
			compute_PN_twobody_Acc( R, Vel, R->AccPN, &Spin[0], &R->Spin_dS_PN[0], i, j );
#else
			compute_PN_twobody_Acc( R, Vel, R->AccPN, i, j );
#endif
		}
	}

	// compute the chained PN accelerations
	int hi,lo;
        for(int i=0;i<R->NumVertex-1;i++){

		int v1 = R->EdgeInMST[i].vertex1;
		int v2 = R->EdgeInMST[i].vertex2;

                if( R->Vertex[v1].level > R->Vertex[v2].level ){
                        hi = v1;
                        lo = v2;
                } else {
                        hi = v2;
                        lo = v1;
                }

                for(int k=0;k<3;k++){
                        R->MSTedgeAcc_PN[3*i+k] = R->AccPN[3*hi+k] - R->AccPN[3*lo+k];
                }
        }

}

#endif

