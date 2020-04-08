#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <float.h>
#include <limits.h>
#include "regularization.h"

#ifdef USE_PN

// relative velocity of particle 1 seen from particle 2 = v_1 - v_2. return v_rel^2.
inline static double relative_velocity(double *v1, double *v2,  double *v12) {
    double vrel2 = 0;
    for(int k=0;k<3;k++){
        v12[k] = v1[k] - v2[k];
        vrel2 += v12[k]*v12[k];
    }
    return vrel2;
}

// From Valtonen/Mikkola code, based on
// Mora & Will, Phys. Rev. D 69, 104021
// Original implementation by P. Pihajoki, modified by A. Rantala.

void compute_PN_twobody_Acc(struct RegularizedRegion *R, double *Vel, double *AccPN, int i, int j){

    const double c  = SPEEDOFLIGHT;
    const double c2 = c*c;
    const double c3 = c*c2;
    const double c4 = c2*c2;
    const double c5 = c*c4;
    const double c6 = c3*c3;
    const double c7 = c*c6;
    const double pi2 = M_PI*M_PI;

    double acc[3]      = {0};
    double acc_grad[3] = {0};

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

    // the PN acceleration coefficients
    double A1=0, A2=0, A2p5=0, A3=0, A3p5=0, Atot = 0, Agrad = 0;
    double B1=0, B2=0, B2p5=0, B3=0, B3p5=0, Btot = 0, Bgrad = 0;

    int k;

    // calculate relative quantities and derived quantities.
    // every relative quantity x = x_1 - x_2, where x_1 corresponds to
    // index i and x_2 to index j

    const int Nd = GLOBAL_ND;
    int path[Nd]; int sign[Nd]; int d;
    int proximity = check_relative_proximity(i,j,Nd,R,&d,path,sign);

    for(k=0;k<3;k++){
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

    for (k = 0; k < 3; k++) {
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

    for(k=0;k<3;k++){
        nv[k] = rv[k]/r;
    }

#if 0
    A1 = 2*(2+eta)*gmpr - (1+3*eta)*v2 + 1.5e0*eta*rdot2;
    B1 = 2*(2-eta)*rdot;
    A2 = -(3./4)*(12+29*eta)*gmpr2 - eta*(3-4*eta)*v4-(15./8)*eta*(1-3*eta)*rdot4 + (1./2)*eta*(13-4*eta)*gmpr*v2+(2+25*eta+2*eta2)*gmpr*rdot2 + (3./2)*eta*(3-4*eta)*v2*rdot2;
    B2 = -(1./2)*rdot*((4+41*eta+8*eta2)*gmpr - eta*(15+4*eta)*v2 + 3*eta*(3+2*eta)*rdot2);
#endif

    A2p5=8./5*eta*gmpr*rdot*(17./3*gmpr+3*v2);
    B2p5=-8./5.*eta*gmpr*(3*gmpr+v2);

#if 0
    A3=(16+(1399./12-41./16*pi2)*eta+71./2*eta2)*gmpr3+eta*(20827./840+123./64*pi2-eta2)*gmpr2*v2-(1+(22717./168+615./64*pi2)*eta+11./8*eta2-7*eta3)*gmpr2*rdot2-.25*eta*(11-49*eta+52*eta2)*v6+35./16*eta*(1-5*eta+5*eta2)*rdot6-.25*eta*(75+32*eta-40*eta2)*gmpr*v2*v2-.5*eta*(158-69*eta-60*eta2)*gmpr*rdot4+eta*(121-16*eta-20*eta2)*gmpr*v2*rdot2+3./8*eta*(20-79*eta+60*eta2)*v2*v2*rdot2-15./8*eta*(4-18*eta+17*eta2)*v2*rdot4;
    B3=rdot*((4+(5849./840.+123./32.*pi2)*eta-25*eta2-8*eta3)*gmpr2+1./8.*eta*(65-152*eta-48*eta2)*v2*v2+15/8.*eta*(3-8*eta-2*eta2)*rdot4+eta*(15+27*eta+10*eta2)*gmpr*v2-1./6.*eta*(329+177*eta+108*eta2)*gmpr*rdot2-.75*eta*(16-37*eta-16*eta2)*v2*rdot2);
    A3p5=-8./5*eta*gmpr*rdot*(23./14*(43+14*eta)*gmpr2+3./28*(61+70*eta)*v2*v2+70*rdot4+1./42*(519-1267*eta)*gmpr*v2+.25*(147+188*eta)*gmpr*rdot2-15/4.*(19+2*eta)*v2*rdot2);
    B3p5=8./5.*eta*gmpr*(1./42.*(1325+546*eta)*gmpr2+1./28.*(313+42*eta)*v2*v2+75*rdot4-1./42.*(205+777*eta)*gmpr*v2+1./12.*(205+424*eta)*gmpr*rdot2-.75*(113+2*eta)*v2*rdot2);
#endif

    // sum all terms with correct powers of c
    Atot = A1/c2 + A2/c4 + A2p5/c5 + A3/c6 + A3p5/c7;
    Btot = B1/c2 + B2/c4 + B2p5/c5 + B3/c6 + B3p5/c7;

    // collect gravitational radiation terms separately
    Agrad = A2p5/c5 + A3p5/c7;
    Bgrad = B2p5/c5 + B3p5/c7;

    for(k=0;k<3;k++){
        acc[k]      = gm/r2*(nv[k]*Atot  + vv[k]*Btot);
        acc_grad[k] = gm/r2*(nv[k]*Agrad + vv[k]*Bgrad);
    }

    for(k=0;k<3;k++){
        AccPN[3*i+k] += (m2/m)*acc[k];
        AccPN[3*j+k] -= (m1/m)*acc[k];
    }

}

void compute_Post_Newtonian_Acc(struct RegularizedRegion *R, double *Vel){

	// Init
	for(int i=0;i<3*R->NumVertex;i++){
		R->AccPN[i] = 0;
	}

	return;

	// Find the SMBHs and their indexes
	int Nbh=0;
	for(int i=0;i<R->NumVertex;i++){
		if(R->Vertex[i].type ==5){
			Nbh++;
		}
	}
	if(Nbh<2){
		return;
	}

	int c=0;
	int bhindex[Nbh];
	for(int i=0;i<R->NumVertex;i++){
                if(R->Vertex[i].type ==5){
			bhindex[c] = i;
                        c++;
                }
        }

	// Compute the PN accelerations for SMBHs

	for(int i=0;i<Nbh;i++){
		for(int j=i+1;j<Nbh;j++){
			compute_PN_twobody_Acc(R,Vel,R->AccPN,bhindex[i],bhindex[j]);
		}
	}

	// Finally the chained PN accelerations

	int hi,lo;
        for(int i=0;i<R->NumVertex-1;i++){
                int index = R->MSTedgeList[i];
                int v1 = R->Edge[index].vertex1;
                int v2 = R->Edge[index].vertex2;
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

