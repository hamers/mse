#include "call_sse.h"
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <math.h>

using namespace std;

extern "C" void evolv1_(int *kw, double *mass,double *mt, double *r, double *lum, double *mc, double *rc, double *menv, double *renv, double *ospin, double *epoch, double *tms, double *tphys, double *tphysf, double *dtp, double *z, double *zpars);
extern "C" void zcnsts_(double* z, double *zpars);
extern "C" void star_(int *kw, double *mass, double *mt, double *tm, double *tn, double *tscls, double *lums, double *GB, double *zpars);
extern "C" void deltat_(int *kw, double *age, double *tm, double *tn, double *tscls, double *dt, double *dtr);

struct value1__ value1_;
struct value2__ value2_;
struct value3__ value3_;
struct value4__ value4_;
struct value5__ value5_;
struct flags__ flags_;
struct points__ points_;

int main(int argc, char **argv)
{
        double z = 0.02;
        double mass = 5.0;

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
        flags_.ceflag = 0;
        flags_.tflag = 0;
        flags_.ifflag = 0;
        flags_.nsflag = 1;
        flags_.wdflag = 1;
        points_.pts1 = 0.05;
        points_.pts2 = 0.01;
        points_.pts3 = 0.02;
        
/*
      READ(22,*)mass,z,tphysf
      READ(22,*)neta,bwind,hewind,sigma
      READ(22,*)ifflag,wdflag,bhflag,nsflag,mxns
      READ(22,*)pts1,pts2,pts3
*/


        
//        double neta = 0.5;
//        double bwind = 0.0;
//        double hewind = 0.5;
//        double sigma = 190.0;

        double *zpars;
        zpars = new double[20];

        double mt = mass;
        int kw = 1;
        double tphys = 0.0;
        double epoch = 0.0;
        double ospin = 0.0;
        
        zcnsts_(&z,zpars);
        
        double r,lum,mc,rc,menv,renv,tms;

        
        //evolv1_(&kw,&mass,&mt,&r,&lum,&mc,&rc,&menv,&renv,&ospin,&epoch,&tms,&tphys,&tphysf,&dtp,&z,zpars);
        
        double dt = 0.0;
        dt = get_new_dt(kw,mass,mt,tphys,dt,epoch,zpars);
        //cout << "dt " << dt << endl;
        //exit(0);
        double tphysf = dt;

        //double dtp = 0.0;
        double dtp;
        
        if (1==0)
        {
            dtp = 0.0;
            tphysf = 12000.0;
            evolv1_(&kw,&mass,&mt,&r,&lum,&mc,&rc,&menv,&renv,&ospin,&epoch,&tms,&tphys,&tphysf,&dtp,&z,zpars);
        }
        else
        {
        while (tphysf <= 12000.0)
        //while (tphysf < 24.733452227138400)
        {
            dtp = tphysf;
            evolv1_(&kw,&mass,&mt,&r,&lum,&mc,&rc,&menv,&renv,&ospin,&epoch,&tms,&tphys,&tphysf,&dtp,&z,zpars);
            
            dt = get_new_dt(kw,mass,mt,tphysf,dt,epoch,zpars);
            //dt=10.0;

            cout << " t " << tphysf << " mt " << mt << " dt " << dt << " epoch " << epoch << endl;            
            //cout << " test epoch " << epoch << " " << tphys << " " << tm << " " << tn << " age " << age << " dtm " << dtm << " dtr " << dtr << endl;
            //cout << " test epoch " << epoch << " " << tphys << " " << tm << " " << tn << " age " << age << " dtm " << dtm << " dtr " << dtr << endl;

            tphys = tphysf;
            tphysf += dt;

        }
        }
        cout << "done kw " << kw << " mass " << mass << " mt " << mt << " r " << log10(r) << " lum " << log10(lum) << " tphys " << tphys << endl;

    return 1;
}

double get_new_dt(int kw, double mass, double mt, double tphys, double dt, double epoch, double *zpars)
{
    double dtm,dtr;
    double tm,tn,age;

    double *GB,*tscls,*lums;
    GB = new double[10];
    tscls = new double[20];
    lums = new double[10];    
    
    star_(&kw, &mass, &mt, &tm, &tn, tscls, lums, GB, zpars);

    age = tphys - epoch;
    
    deltat_(&kw,&age,&tm,&tn,tscls,&dt,&dtr);

    dtm = std::min(dtr, dt);
    dtm = std::max(dtm,1.0d-07*age);
    
    return dtm;
}
/*

      READ(22,*)mass,z,tphysf
      READ(22,*)neta,bwind,hewind,sigma
      READ(22,*)ifflag,wdflag,bhflag,nsflag,mxns
      READ(22,*)pts1,pts2,pts3


10.0 0.02 12000.
0.5 0.0 0.5 190.0
0 1 0 1 3.0 999
0.05 0.01 0.02

      integer kw,it,ip,jp,j,kwold,rflag
      integer nv
      parameter(nv=50000)

      real*8 mass,z,aj
      real*8 epoch,tphys,tphys2,tmold,tbgold
      real*8 mt,tm,tn,tphysf,dtp,tsave
      real*8 tscls(20),lums(10),GB(10),zpars(20)
      real*8 r,lum,mc,teff,rc,menv,renv,vs(3)
      real*8 ospin,jspin,djt,djmb,k2,k3
      parameter(k3=0.21d0)
      real*8 m0,r1,lum1,mc1,rc1,menv1,renv1,k21
      real*8 dt,dtm,dtr,dr,dtdr,dms,dml,mt2,rl
      real*8 tol,tiny
      parameter(tol=1.0d-10,tiny=1.0d-14)
      real*8 ajhold,rm0,eps,alpha2
      parameter(eps=1.0d-06,alpha2=0.09d0)
      real*8 mlwind,vrotf
      external mlwind,vrotf
      logical iplot,isave
      REAL*8 neta,bwind,hewind,mxns
      COMMON /VALUE1/ neta,bwind,hewind,mxns
      REAL*8 pts1,pts2,pts3
      COMMON /POINTS/ pts1,pts2,pts3
      REAL scm(50000,14),spp(20,3)
      COMMON /SINGLE/ scm,spp
*/
