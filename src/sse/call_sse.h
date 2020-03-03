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

/*
      READ(22,*)mass,z,tphysf
      READ(22,*)neta,bwind,hewind,sigma
      READ(22,*)ifflag,wdflag,bhflag,nsflag,mxns
      READ(22,*)pts1,pts2,pts3
*
*
* const_bse.h
*
      INTEGER idum
      COMMON /VALUE3/ idum
      INTEGER idum2,iy,ir(32)
      COMMON /RAND3/ idum2,iy,ir
      INTEGER ktype(0:14,0:14)
      COMMON /TYPES/ ktype
      INTEGER ceflag,tflag,ifflag,nsflag,wdflag
      COMMON /FLAGS/ ceflag,tflag,ifflag,nsflag,wdflag
      INTEGER bhflag
*
      REAL*8 neta,bwind,hewind,mxns,alpha1,lambda
      REAL*8 sigma,beta,xi,acc2,epsnov,eddfac,gamma
      COMMON /VALUE1/ neta,bwind,hewind,mxns
      COMMON /VALUE2/ alpha1,lambda
      COMMON /VALUE4/ sigma,bhflag
      COMMON /VALUE5/ beta,xi,acc2,epsnov,eddfac,gamma
      REAL*8 pts1,pts2,pts3
      COMMON /POINTS/ pts1,pts2,pts3
      REAL*8 dmmax,drmax
      COMMON /TSTEPC/ dmmax,drmax
      REAL scm(50000,14),spp(20,3)
      COMMON /SINGLE/ scm,spp
      REAL bcm(50000,34),bpp(80,10)
      COMMON /BINARY/ bcm,bpp
*
*/


extern struct value1__ value1_;
extern struct value2__ value2_;
extern struct value3__ value3_;
extern struct value4__ value4_;
extern struct value5__ value5_;
extern struct flags__ flags_;
extern struct points__ points_;

double get_new_dt(int kw, double mass, double mt, double tphys, double dt, double epoch, double *zpars);
