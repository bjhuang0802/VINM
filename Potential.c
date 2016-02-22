#include "Potential.h"
#include <math.h>
#define rcut  1.2

double PhiC(double qa,double qb,double r) {
   double Hee,krf,crf,er;
   Hee=KCoulomb*qa*qb/r;
   return(Hee);
}

double DPhiC(double qa,double qb,double r) {
   double DHee,krf,crf,er;
   DHee=KCoulomb*qa*qb*(-1.0/r/r);
   return(DHee);
}

double DDPhiC(double qa,double qb,double r) {
   double DDHee,krf,crf,er;
   DDHee=KCoulomb*qa*qb*(2.0/r/r/r);
   return(DDHee);
}

double PhiLJ(double r,double sig,double eps) {
   double r2,r4,r6,r12;
   double tmp;
   double s12,s6;
    
   s6=pow(sig,6.0);
   s12=s6*s6;
   r2=r*r;
   r4=r2*r2;
   r6=r4*r2;
   r12=r6*r6;
   tmp=4.0*eps*((s12/r12) - (s6/r6));
   return(tmp);
}

double DPhiLJ(double r,double sig,double eps) {
   double r2,r4,r6,r12,r7,r13;
   double tmp;
   double s12,s6;
   s6=pow(sig,6.0);
   s12=s6*s6;

   r2=r*r;
   r4=r2*r2;
   r6=r4*r2;
   r12=r6*r6;
   r7=r6*r;
   r13=r12*r;
   tmp=4.0*eps*((6.0*s6/r7) - (12.0*s12/r13));
   return(tmp);
}

double DDPhiLJ(double r,double sig,double eps) {
   double r2,r4,r6,r12,r8,r14;
   double tmp=0.0;
   double s12,s6;
   s6=pow(sig,6.0);
   s12=s6*s6;

   r2=r*r;
   r4=r2*r2;
   r6=r4*r2;
   r12=r6*r6;
   r8=r6*r2;
   r14=r12*r2;

   tmp=4.0*eps*((156.0*s12/r14) - (42.0*s6/r8));
   return(tmp);
}

double PhiIR(double r,double ***coe,long i) {
   double bond,Kb;
   bond = *(*(*(coe+0)+0)+i);
   Kb= *(*(*(coe+0)+1)+i);

   return(0.5*Kb*pow(r-bond,2.0));
}

double DPhiIR(double r,double ***coe,long i) {
   double bond,Kb;
   bond = *(*(*(coe+0)+0)+i);
   Kb= *(*(*(coe+0)+1)+i);
   
   return(Kb*(r-bond));
}

double DDPhiIR(double r,double ***coe,long i) {
   double bond,Kb;
   bond = *(*(*(coe+0)+0)+i);
   Kb= *(*(*(coe+0)+1)+i);
   
   return(Kb);
}
double PhiLJ14(double r,double ***coe,long i) {
   double C6,C12,r2,r4,r6,r12,tmp;
   C6 = *(*(*(coe+3)+0)+i);
   C12= *(*(*(coe+3)+1)+i);
   r2=r*r;
   r4=r2*r2;
   r6=r4*r2;
   r12=r6*r6;
   tmp=(C12/r12) - (C6/r6);
   return(tmp);
}

double DPhiLJ14(double r,double ***coe,long i) {
   double C6,C12,r2,r4,r6,r7,r12,r13,tmp;
   C6 = *(*(*(coe+3)+0)+i);
   C12= *(*(*(coe+3)+1)+i);
   r2=r*r;
   r4=r2*r2;
   r6=r4*r2;
   r7=r6*r;
   r13=r12*r;
   tmp=((6.0*C6/r7) - (12.0*C12/r13));
   
   return(tmp);
}

double DDPhiLJ14(double r,double ***coe,long i) {
   double C6,C12,r2,r4,r6,r8,r12,r14,tmp;
   C6 = *(*(*(coe+3)+0)+i);
   C12= *(*(*(coe+3)+1)+i);
   r2=r*r;
   r4=r2*r2;
   r6=r4*r2;
   r8=r6*r2;
   r12=r6*r6;
   r14=r12*r2;
   tmp=((156.0*C12/r14) - (42.0*C6/r8));
   
   return(tmp);
}
