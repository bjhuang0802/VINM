#include <IOUtility.h>
#include <math.h>
#include <DiagSM.h>
#include "Potential.h"

void GetParam(double *,double *,double *,double *,long **);
void GetConf(FILE *,long,double ***);
void GetCoeff(FILE *,long ***,double ***,long *);
void GetHessian(long,double ***,double **);
void SetZero(double **,long);
double Dt(long,long);
double DDt(long,long,long);
void SetHessianCNE(long,double ***,double **);
double GetR(double *,double *,double *);
void SetHessianCE(long,double ***,double **);
void SetHessianIR(long,double ***,double **);
void SetHessianLJ14(long,double ***,double **);
double CalY(double ***,long,long);
double CalDY(double ***,long,long,long,long);
double CalDDY(double ***,long,long,long,long,long,long);
double CalZ(long,long,long);
double CalDZ(long,long,long,long,long);
double CalDDZ(long,long,long,long,long,long,long);
double CalT(long ,long ) ;
double CalDT(long ,long ,long , long ) ;
double CalDDT(long ,long ,long , long ,long , long );
void SetHessianIA1(long,double ***,double **);
void SetHessianIA2(long,double ***,double **);
void SetHessianIA3(long,double ***,double **);
void SetHessianID(long,double ***,double **);
void DividMass(long,double **);
void SetSM(long,double **);
void SetHessianLJNE(long,double ***,double **);
void SetHessianLJE(long,double ***,double **);
double potLJSR(long NMolecules,double ***);

   double SMass[12],Charges[12],eps[12],sig[12];
   long NMolecules,NAtom,NConf,ff,***ato,np[4],**ecl;
   FILE *Conf,*Coeff,*NHBVector,*Output1,*Output2;
   double ***R,**H,*Omega,*Lambda,*e,***coe;
   double BL=2.92837,HBL,NHBL;

int main(void) {
   char filename[100];
   long i,j,k,l,x,y,z,index1,index2;
   long alpha,beta,gamma,delta,mu,nu;
   double density,BoxLength,HalfBoxLength,NegHalfBoxLength,all,sum[3];
   double DOSparallel,DOSperpendicular,backvector[3],parallel,perpendicular;
   
   Conf=getAfile("Configurations?\n",filename,"r");
   Coeff=getAfile("Coefficient?\n",filename,"r");
   Output1=getAfile("Output for eigenvalus?\n",filename,"w");
   Output2=getAfile("Output for eigenvectors\n",filename,"w");
   getAlong("Serial number of data set?\n",&ff);
   getAlong("How many configurations?\n",&NConf);
   getAlong("How many molecules?\n",&NMolecules);
   getAlong("How many atoms of a molecule?\n",&NAtom);
   ecl=(long **)AllocateMemory2(sizeof(long),NAtom,NAtom);
   GetParam(Charges,SMass,sig,eps,ecl);

   R=(double ***)AllocateMemory3(sizeof(double),NMolecules,NAtom,3);
   H=(double **)AllocateMemory2(sizeof(double),3*NAtom*NMolecules+1,3*NAtom*NMolecules+1);
   ato=(long ***)AllocateMemory3(sizeof(long),5,5,20);
   coe=(double ***)AllocateMemory3(sizeof(double),5,5,20);
   Omega=(double *)AllocateMemory1(sizeof(double),3*NAtom*NMolecules+1);
   Lambda=(double *)AllocateMemory1(sizeof(double),3*NAtom*NMolecules+1);
   e=(double *)AllocateMemory1(sizeof(double),3*NAtom*NMolecules+1);
   GetCoeff(Coeff,ato,coe,np);
   HBL=BL/2.0;
   NHBL=-HBL;
   printf("The system length= %12.7f\n",BL);
   for(i=1;i<NConf;i++) {
       GetConf(Conf,NMolecules,R);
       //printf("The LJ potential is : %12.8f\n", potLJSR( NMolecules,R));
       if(i<ff)continue;
       GetHessian(NMolecules,R,H);
       printf(" Hessian is done!\n");

       tred2(H,NMolecules*NAtom*3,Lambda,e);
       tqli(Lambda,e,NMolecules*NAtom*3,H);

       for(j=0;j<NMolecules*NAtom*3;j++) {
	   if(*(Lambda+j+1)>0)*(Omega+j+1)=5.30693*sqrt(*(Lambda+j+1) );
	   if(*(Lambda+j+1)<0)*(Omega+j+1)=-5.30693*sqrt(-*(Lambda+j+1) );
           fprintf(Output1,"%16.6e\n ",*(Omega+j+1));
//           fprintf(Output1,"%16.6e %5.2f\n ",*(Omega+j+1)+0.1,0.0);
           fprintf(Output2,"conf. %ld  eigenvalue= %10.5f\n",i,*(Omega+j+1));
           for(k=0;k<NMolecules;k++) {
	       for(alpha=0;alpha<NAtom;alpha++) {
	              for(l=0;l<3;l++) {
	                 fprintf(Output2,"%20.10e ", *(*(H+k*NAtom*3+alpha*3+l+1)+j+1));
	                 //printf("%4.1e ", *(*(H+k*3*NAtom+alpha*3+l+1)+j+1));
	              }
	                 fprintf(Output2,"\n");
               }
	       //printf("\n ");
	   }
       } 
   }

   fclose(Conf);
   fclose(Output1);
   free(e);
   free(Omega);
   Free2((void **)H,3*NAtom*NMolecules+1);
   Free3((void ***)R,NMolecules,3);
}

void GetHessian(long NMolecules,
   double ***R,double **H) {
   
   SetZero(H,NMolecules);
   SetHessianCNE(NMolecules,R,H);
   SetHessianCE(NMolecules,R,H);
   SetHessianLJNE(NMolecules,R,H);
   SetHessianLJE(NMolecules,R,H);
   SetHessianIR(NMolecules,R,H);
   SetHessianLJ14(NMolecules,R,H);
   SetHessianIA1(NMolecules,R,H);
   SetHessianIA2(NMolecules,R,H);
   SetHessianIA3(NMolecules,R,H);
   SetHessianID(NMolecules,R,H);
   SetSM(NMolecules,H);
   DividMass(NMolecules,H);
}

void SetSM(long NMolecules,double **H) {
   long i,j;
   
   for(i=0;i<3*NAtom*NMolecules;i++) {
      for(j=i+1;j<3*NAtom*NMolecules;j++) {
         *(*(H+j+1)+i+1)=*(*(H+i+1)+j+1);
      }
   }
}

void DividMass(long NMolecules,double **H) {
   long i,j,alpha,beta,mu,nu;
   
   for(i=0;i<NMolecules;i++) {
      for(j=0;j<NMolecules;j++) {
         for(alpha=0;alpha<NAtom;alpha++) {
            for(beta=0;beta<NAtom;beta++) {
               for(mu=0;mu<3;mu++) {
                  for(nu=0;nu<3;nu++) {
                     *(*(H+i*3*NAtom+alpha*3+mu+1)+j*3*NAtom+beta*3+nu+1) *=
                      SMass[alpha]*SMass[beta];
                  }
               }
            }
         }
      }
   }
}

void SetHessianID(long NMolecules,double ***R,double **H){
     long k,alpha,beta,gamma,delta,mu,nu,ib;
     double c1,c2,c3,c4;
     double diff1,diff2,the;
     double term1,term2,term3,T2;
     for(k=0;k<NMolecules;k++) {
        for(ib=0;ib<np[2];ib++){
	c1=*(*(*(coe+2)+0)+ib);
	c2=*(*(*(coe+2)+1)+ib);
	c3=*(*(*(coe+2)+2)+ib);
	c4=*(*(*(coe+2)+3)+ib);
	the=acos(CalT(k,ib));
	T2=(1-CalT(k,ib)*CalT(k,ib));
	for(alpha=0;alpha<4;alpha++){
	    for(beta=0;beta<4;beta++){
		gamma=*(*(*(ato+2)+alpha)+ib);
		delta=*(*(*(ato+2)+ beta)+ib);
		for(mu=0;mu<3;mu++){
		    for(nu=0;nu<3;nu++){
			diff1=0.5*(-c1*cos(the)+4*c2*cos(2*the)-9*c3*cos(3*the)+16*c4*cos(4*the));
			diff2=0.5*(-c1*sin(the)+2*c2*sin(2*the)-3*c3*sin(3*the)+4*c4*sin(4*the));
			term1=diff1*CalDT(k,ib,gamma,mu)*CalDT(k,ib,delta,nu)/T2;
			term2=diff2*CalT(k,ib)/pow(T2,1.5)*CalDT(k,ib,gamma,mu)*CalDT(k,ib,delta,nu);
			term3=diff2*CalDDT(k,ib,gamma,mu,delta,nu)/pow(T2,0.5);
		      *(*(H+k*3*NAtom+gamma*3+mu+1)+k*3*NAtom+delta*3+nu+1) += 
			  term1 + term2 +term3;
		    }
		}
	    }
	}
	}
     }
}
void SetHessianIA1(long NMolecules,
   double ***R,double **H) {

   long k,alpha,beta,gamma,delta,mu,nu,ib;
   double Yk,ACOS2,Kth;
   
   for(k=0;k<NMolecules;k++) {
      for(ib=0;ib<np[1];ib++){  
      Kth=*(*(*(coe+1)+0)+ib);
      Yk=CalY(R,k,ib);
      ACOS2=1.0/(1.0-Yk*Yk);
      for(alpha=0;alpha<3;alpha++) {
         for(beta=0;beta<3;beta++) {
	    gamma=*(*(*(ato+1)+alpha)+ib);
	    delta=*(*(*(ato+1)+beta )+ib);
            for(mu=0;mu<3;mu++) {
               for(nu=0;nu<3;nu++) {
                  *(*(H+k*3*NAtom+gamma*3+mu+1)+k*3*NAtom+delta*3+nu+1) +=
                     //2.0*Kth*ACOS2*CalDY(R,k,ib,gamma,mu)*
                     Kth*ACOS2*CalDY(R,k,ib,gamma,mu)*
                     CalDY(R,k,ib,delta,nu);
               }
            }
         }
      }
   }
   }
}

void SetHessianIA2(long NMolecules,double ***R,double **H) {
   long k,alpha,beta,gamma,delta,mu,nu,ib;
   double Yk,Kth,the;

   for(k=0;k<NMolecules;k++) {
      for(ib=0;ib<np[1];ib++){  
      Kth=*(*(*(coe+1)+0)+ib);
      the=*(*(*(coe+1)+1)+ib);
      Yk=CalY(R,k,ib);
      for(alpha=0;alpha<3;alpha++) {
         for(beta=0;beta<3;beta++) {
	    gamma=*(*(*(ato+1)+alpha)+ib);
	    delta=*(*(*(ato+1)+beta )+ib);
            for(mu=0;mu<3;mu++) {
               for(nu=0;nu<3;nu++) {
                  *(*(H+k*3*NAtom+gamma*3+mu+1)+k*3*NAtom+delta*3+nu+1) +=
                     //2.0*Kth*(acos(Yk)-the*3.141592654/180.0)*
                     Kth*(acos(Yk)-the*3.141592654/180.0)*
                     CalDDY(R,k,ib,gamma,mu,delta,nu)*
                     (0.0-1.0/sqrt(1.0-Yk*Yk));
               }
            }
         }
      }
   }
   }
}

void SetHessianIA3(long NMolecules,double ***R,double **H) {
   long k,alpha,beta,gamma,delta,mu,nu,ib;
   double Yk,Kth,the;

   for(k=0;k<NMolecules;k++) {
      for(ib=0;ib<np[1];ib++){  
      Kth=*(*(*(coe+1)+0)+ib);
      the=*(*(*(coe+1)+1)+ib);
      Yk=CalY(R,k,ib);
      for(alpha=0;alpha<3;alpha++) {
         for(beta=0;beta<3;beta++) {
	    gamma=*(*(*(ato+1)+alpha)+ib);
	    delta=*(*(*(ato+1)+beta )+ib);
            for(mu=0;mu<3;mu++) {
               for(nu=0;nu<3;nu++) {
                  *(*(H+k*3*NAtom+gamma*3+mu+1)+k*3*NAtom+delta*3+nu+1) +=
                     //2.0*Kth*(acos(Yk)-the*3.141592654/180.0)*
                     Kth*(acos(Yk)-the*3.141592654/180.0)*
                     CalDY(R,k,ib,gamma,mu)*
                     (0.0-1.0/pow(sqrt(1.0-Yk*Yk),3.0))*Yk*
                     CalDY(R,k,ib,delta,nu);
               }
            }
         }
      }
   }
   }
}

double CalDDY(double ***R,long k,long ib,long gamma,long mu,long delta,long nu){

   double rOH1,r2OH1,rOH2,r2OH2,ROH1[3],ROH2[3],Yk,DDY;
   double DgOgH1,DgOgH2,DdOdH1,DdOdH2;
   long i,ai,aj,ak;

   ai= *(*(*(ato+1)+0)+ib);
   aj= *(*(*(ato+1)+1)+ib);
   ak= *(*(*(ato+1)+2)+ib);
   
   rOH1=GetR(*(*(R+k)+aj),*(*(R+k)+ai),ROH1);
   r2OH1=rOH1*rOH1;
   rOH2=GetR(*(*(R+k)+aj),*(*(R+k)+ak),ROH2);
   r2OH2=rOH2*rOH2;
   Yk=CalY(R,k,ib);
   DgOgH1=Dt(gamma,aj)-Dt(gamma,ai);
   DgOgH2=Dt(gamma,aj)-Dt(gamma,ak);
   DdOdH1=Dt(delta,aj)-Dt(delta,ai);
   DdOdH2=Dt(delta,aj)-Dt(delta,ak);
   
   DDY=(ROH2[mu]*DgOgH1+ROH1[mu]*DgOgH2)*
      (0.0-ROH1[nu]*DdOdH1/r2OH1 -
      ROH2[nu]*DdOdH2/r2OH2)/(rOH1*rOH2);

   DDY=DDY+Dt(mu,nu)*(DgOgH1*DdOdH2+DgOgH2*DdOdH1)/(rOH1*rOH2);
   DDY=DDY-CalDY(R,k,ib,delta,nu)*
      (ROH1[mu]*DgOgH1/r2OH1 + ROH2[mu]*DgOgH2/r2OH2);
   DDY=DDY-Yk*Dt(mu,nu)*
      (DgOgH1*DdOdH1/r2OH1 + DgOgH2*DdOdH2/r2OH2);
   DDY=DDY-Yk*
      (-2.0*ROH1[mu]*DgOgH1*ROH1[nu]*DdOdH1/(r2OH1*r2OH1) -
      2.0*ROH2[mu]*DgOgH2*ROH2[nu]*DdOdH2/(r2OH2*r2OH2));
   return(DDY);
}

double CalDY(double ***R,long k,long ib,long gamma,long mu){
   
   double rOH1,rOH2,ROH1[3],ROH2[3],Yk,DY,r2OH1,r2OH2,DgOgH1,DgOgH2;
   long i,ai,aj,ak;

   ai= *(*(*(ato+1)+0)+ib);
   aj= *(*(*(ato+1)+1)+ib);
   ak= *(*(*(ato+1)+2)+ib);
   
   DgOgH1=Dt(gamma,aj)-Dt(gamma,ai);
   DgOgH2=Dt(gamma,aj)-Dt(gamma,ak);
   rOH1=GetR(*(*(R+k)+aj),*(*(R+k)+ai),ROH1);
   rOH2=GetR(*(*(R+k)+aj),*(*(R+k)+ak),ROH2);
   r2OH1=rOH1*rOH1;
   r2OH2=rOH2*rOH2;
   Yk=CalY(R,k,ib);
   DY=DgOgH1*ROH2[mu]/(rOH1*rOH2);
   DY=DY+DgOgH2*ROH1[mu]/(rOH1*rOH2);
   DY=DY-Yk*ROH1[mu]*DgOgH1/r2OH1;
   DY=DY-Yk*ROH2[mu]*DgOgH2/r2OH2;
   return(DY);
}

double CalY(double ***R,long k,long ib) {
   long i,ai,aj,ak;
   double ROH1[3],ROH2[3],rOH1,rOH2;
   double DOT;
   ai=*(*(*(ato+1)+0)+ib);
   aj=*(*(*(ato+1)+1)+ib);
   ak=*(*(*(ato+1)+2)+ib);
   rOH1=GetR(*(*(R+k)+aj),*(*(R+k)+ai),ROH1);
   rOH2=GetR(*(*(R+k)+aj),*(*(R+k)+ak),ROH2);
   DOT=0.0;
   for(i=0;i<3;i++)
      DOT=DOT+ROH1[i]*ROH2[i];
   return(DOT/(rOH1*rOH2));
}

double CalT(long k,long ib) {
    double tmp,cos1,cos2,cos3,sin1,sin2;
    cos1=CalZ(1,k,ib);
    cos2=CalZ(2,k,ib);
    cos3=CalZ(3,k,ib);
    sin1=sqrt(1-cos1*cos1);
    sin2=sqrt(1-cos2*cos2);
    tmp=(cos3+cos1*cos2)/(sin1*sin2);
    return(tmp);
}
double CalDT(long k,long ib,long gamma, long mu) {
    double tmp,cos1,cos2,cos3,sin1,sin2;
    double dcos3,dcos1,dcos2,bk;
    cos1=CalZ(1,k,ib);
    cos2=CalZ(2,k,ib);
    cos3=CalZ(3,k,ib);
    sin1=sqrt(1-cos1*cos1);
    sin2=sqrt(1-cos2*cos2);
    dcos1=CalDZ(1,k,ib,gamma,mu);
    dcos2=CalDZ(2,k,ib,gamma,mu);
    dcos3=CalDZ(3,k,ib,gamma,mu);
    bk=cos1/sin1/sin1*dcos1+cos2/sin2/sin2*dcos2;
    tmp=(1.0/sin1/sin2)*(dcos3+dcos1*cos2+cos1*dcos2+(cos3+cos1*cos2)*bk);
    return(tmp);
}
double CalDDT(long k,long ib,long gamma, long mu,long delta, long nu) {
    double tmp,cos1,cos2,cos3,sin1,sin2;
    double dcos3,dcos1,dcos2,ddcos3,ddcos2,ddcos1;
    double dsin1gm,dsin2gm,dsin1dn,dsin2dn,ddsin1,ddsin2;
    double bk1,bk2,bk3,bk4,bk4a,bk4b;

    cos1=CalZ(1,k,ib);
    cos2=CalZ(2,k,ib);
    cos3=CalZ(3,k,ib);
    sin1=sqrt(1-cos1*cos1);
    sin2=sqrt(1-cos2*cos2);
    dsin1gm=-cos1/sin1*CalDZ(1,k,ib,gamma,mu);
    dsin2gm=-cos2/sin2*CalDZ(2,k,ib,gamma,mu);
    dsin1dn=-cos1/sin1*CalDZ(1,k,ib,delta,nu);
    dsin2dn=-cos2/sin2*CalDZ(2,k,ib,delta,nu);
    ddcos1=CalDDZ(1,k,ib,gamma,mu,delta,nu);
    ddcos2=CalDDZ(2,k,ib,gamma,mu,delta,nu);
    ddcos3=CalDDZ(3,k,ib,gamma,mu,delta,nu);
    ddsin1=-CalDZ(1,k,ib,gamma,mu)*CalDZ(1,k,ib,delta,nu)/sin1/sin1/sin1-cos1/sin1*ddcos1;
    ddsin2=-CalDZ(2,k,ib,gamma,mu)*CalDZ(2,k,ib,delta,nu)/sin2/sin2/sin2-cos2/sin2*ddcos2;
    bk1=(1.0/sin1/sin2)*(ddcos3+CalDZ(2,k,ib,delta,nu)*CalDZ(1,k,ib,gamma,mu)
	    +cos2*ddcos1+cos1*ddcos2+CalDZ(1,k,ib,delta,nu)*CalDZ(2,k,ib,gamma,mu));

    bk2=(1.0/sin1/sin2)*(cos1/sin1/sin1*CalDZ(1,k,ib,delta,nu)+cos2/sin2/sin2*CalDZ(2,k,ib,delta,nu))
	*(cos2*CalDZ(1,k,ib,gamma,mu)+cos1*CalDZ(2,k,ib,gamma,mu)+CalDZ(3,k,ib,gamma,mu));

    bk3=-CalDT(k,ib,delta,nu)*(dsin1gm/sin1 + dsin2gm/sin2);
    bk4a=(ddsin1-dsin1gm*dsin1dn/sin1)/sin1;
    bk4b=(ddsin2-dsin2gm*dsin2dn/sin2)/sin2;
    bk4=-CalT(k,ib)*(bk4a+bk4b);
    tmp=bk1+bk2+bk3+bk4;
    return(tmp);
}


double CalZ(long ith,long k,long ib) {
   long i,ai,aj,ak,al;
   double ROH1[3],ROH2[3],rOH1,rOH2;
   double DOT;
   ai=*(*(*(ato+2)+0)+ib);
   aj=*(*(*(ato+2)+1)+ib);
   ak=*(*(*(ato+2)+2)+ib);
   al=*(*(*(ato+2)+3)+ib);
   if(ith==1){
   rOH1=GetR(*(*(R+k)+ai),*(*(R+k)+aj),ROH1);
   rOH2=GetR(*(*(R+k)+aj),*(*(R+k)+ak),ROH2);
   }
   if(ith==2){
   rOH1=GetR(*(*(R+k)+ak),*(*(R+k)+aj),ROH1);
   rOH2=GetR(*(*(R+k)+al),*(*(R+k)+ak),ROH2);
   }
   if(ith==3){
   rOH1=GetR(*(*(R+k)+ai),*(*(R+k)+aj),ROH1);
   rOH2=GetR(*(*(R+k)+al),*(*(R+k)+ak),ROH2);
   }
   DOT=0.0;
   for(i=0;i<3;i++)
      DOT=DOT+ROH1[i]*ROH2[i];
   return(DOT/(rOH1*rOH2));
}
double CalDZ(long ith,long k,long ib,long gamma,long mu){
   
   double rv1,rv2,v1[3],v2[3],Yk,DZ,r2v1,r2v2,DDt1,DDt2;
   long i,ai,aj,ak,al;

   ai= *(*(*(ato+2)+0)+ib);
   aj= *(*(*(ato+2)+1)+ib);
   ak= *(*(*(ato+2)+2)+ib);
   al= *(*(*(ato+2)+3)+ib);
   
   if(ith==1){
   rv1=GetR(*(*(R+k)+ai),*(*(R+k)+aj),v1);
   rv2=GetR(*(*(R+k)+aj),*(*(R+k)+ak),v2);
   DDt1=Dt(gamma,ai)-Dt(gamma,aj);
   DDt2=Dt(gamma,aj)-Dt(gamma,ak);
   }
   if(ith==2){
   rv1=GetR(*(*(R+k)+ak),*(*(R+k)+aj),v1);
   rv2=GetR(*(*(R+k)+al),*(*(R+k)+ak),v2);
   DDt1=Dt(gamma,ak)-Dt(gamma,aj);
   DDt2=Dt(gamma,al)-Dt(gamma,ak);
   }
   if(ith==3){
   rv1=GetR(*(*(R+k)+ai),*(*(R+k)+aj),v1);
   rv2=GetR(*(*(R+k)+al),*(*(R+k)+ak),v2);
   DDt1=Dt(gamma,ai)-Dt(gamma,aj);
   DDt2=Dt(gamma,al)-Dt(gamma,ak);
   }
   r2v1=rv1*rv1;
   r2v2=rv2*rv2;
   Yk= CalZ(ith,k,ib);
   DZ= (DDt1*v2[mu]+ DDt2*v1[mu]) /(rv1*rv2) 
       -Yk*(v1[mu]*DDt1/r2v1 + v2[mu]*DDt2/r2v2);
   return(DZ);
}
double CalDDZ(long ith,long k,long ib,long gamma,long mu,long delta,long nu){

   double rOH1,r2OH1,rOH2,r2OH2,v1[3],v2[3],Yk,DDZ;
   double rij,rjk,rkl,r2ij,r2jk,r2kl;
   double DgOgH1,DgOgH2,DdOdH1,DdOdH2;
   double A,DA,B1,B2,B2b,dB1,dB2,db2ba,db2bb;
   long i,ai,aj,ak,al;

   ai= *(*(*(ato+2)+0)+ib);
   aj= *(*(*(ato+2)+1)+ib);
   ak= *(*(*(ato+2)+2)+ib);
   al= *(*(*(ato+2)+3)+ib);
   
   Yk=CalZ(ith,k,ib);
   if(ith==1){
     rij=GetR(*(*(R+k)+ai),*(*(R+k)+aj),v1);
     rjk=GetR(*(*(R+k)+aj),*(*(R+k)+ak),v2);
     r2ij=rij*rij; r2jk=rjk*rjk;
     A=1.0/(rij*rjk);
     DA=(DDt(delta,aj,ai)*v1[nu]/r2ij + DDt(delta,ak,aj)*v2[nu]/r2jk)*A;
     B1=DDt(gamma,aj,ak)*v1[mu]+ DDt(gamma,ai,aj)*v2[mu]; 
     B2b=(DDt(gamma,ai,aj)*v1[mu]*rjk/rij + DDt(gamma,aj,ak)*v2[mu]*rij/rjk);
     B2=-Yk*B2b;
     dB1=Dt(mu,nu)*(DDt(gamma,ai,aj)*DDt(delta,aj,ak)+DDt(gamma,aj,ak)*DDt(delta,ai,aj));
     db2ba=rjk/rij*DDt(gamma,ai,aj)
         *(Dt(mu,nu)*DDt(delta,ai,aj) + v1[mu]*(DDt(delta,aj,ai)*v1[nu]/r2ij+DDt(delta,aj,ak)*v2[nu]/r2jk)); 
     db2bb=rij/rjk*DDt(gamma,aj,ak)
         *(Dt(mu,nu)*DDt(delta,aj,ak) + v2[mu]*(DDt(delta,ai,aj)*v1[nu]/r2ij+DDt(delta,ak,aj)*v2[nu]/r2jk)); 
     dB2=-CalDZ(ith,k,ib,delta,nu)*B2b - Yk*(db2ba+db2bb);
   }
   if(ith==2){
     rjk=GetR(*(*(R+k)+aj),*(*(R+k)+ak),v1);
     rkl=GetR(*(*(R+k)+ak),*(*(R+k)+al),v2);
     r2jk=rjk*rjk; r2kl=rkl*rkl;
     A=1.0/(rjk*rkl);
     DA=(DDt(delta,ak,aj)*v1[nu]/r2jk + DDt(delta,al,ak)*v2[nu]/r2kl)*A;
     B1=DDt(gamma,ak,al)*v1[mu]+ DDt(gamma,aj,ak)*v2[mu]; 
     B2b=(DDt(gamma,aj,ak)*v1[mu]*rkl/rjk + DDt(gamma,ak,al)*v2[mu]*rjk/rkl);
     B2=-Yk*B2b;
     dB1=Dt(mu,nu)*(DDt(gamma,aj,ak)*DDt(delta,ak,al)+DDt(gamma,ak,al)*DDt(delta,aj,ak));
     db2ba=rkl/rjk*DDt(gamma,aj,ak)
         *(Dt(mu,nu)*DDt(delta,aj,ak) + v1[mu]*(DDt(delta,ak,aj)*v1[nu]/r2jk+DDt(delta,ak,al)*v2[nu]/r2kl)); 
     db2bb=rjk/rkl*DDt(gamma,ak,al)
         *(Dt(mu,nu)*DDt(delta,ak,al) + v2[mu]*(DDt(delta,aj,ak)*v1[nu]/r2jk+DDt(delta,al,ak)*v2[nu]/r2kl)); 
     dB2=-CalDZ(ith,k,ib,delta,nu)*B2b - Yk*(db2ba+db2bb);
   }
   if(ith==3){
     rij=GetR(*(*(R+k)+aj),*(*(R+k)+ai),v1);
     rkl=GetR(*(*(R+k)+ak),*(*(R+k)+al),v2);
     r2ij=rij*rij; r2kl=rkl*rkl;
     A=1.0/(rij*rkl);
     DA=(DDt(delta,ai,aj)*v1[nu]/r2ij + DDt(delta,al,ak)*v2[nu]/r2kl)*A;
     B1=DDt(gamma,ak,al)*v1[mu]+ DDt(gamma,aj,ai)*v2[mu]; 
     B2b=(DDt(gamma,aj,ai)*v1[mu]*rkl/rij + DDt(gamma,ak,al)*v2[mu]*rij/rkl);
     B2=-Yk*B2b;
     dB1=Dt(mu,nu)*(DDt(gamma,aj,ai)*DDt(delta,ak,al)+DDt(gamma,ak,al)*DDt(delta,aj,ai));
     db2ba=rkl/rij*DDt(gamma,aj,ai)
         *(Dt(mu,nu)*DDt(delta,aj,ai) + v1[mu]*(DDt(delta,ai,aj)*v1[nu]/r2ij+DDt(delta,ak,al)*v2[nu]/r2kl)); 
     db2bb=rij/rkl*DDt(gamma,ak,al)
         *(Dt(mu,nu)*DDt(delta,ak,al) + v2[mu]*(DDt(delta,aj,ai)*v1[nu]/r2ij+DDt(delta,al,ak)*v2[nu]/r2kl)); 
     dB2=-CalDZ(ith,k,ib,delta,nu)*B2b - Yk*(db2ba+db2bb);
   }
   DDZ=DA*(B1+B2) + A*(dB1 +dB2); //DA*(B1+B2) + A*(dB1+ dB2)
   

   return(DDZ);
}

void SetHessianIR(long NMolecules,double ***R,double **H) {
   
   long k,gamma,delta,alpha,beta,j,mu,nu,i,ai,aj;
   double Roa[3],roa;
   
   for(k=0;k<NMolecules;k++) {
   //for(k=0;k<1;k++) {
      for(j=0;j<np[0];j++) {
	 ai=*(*(*(ato+0)+0)+j); 
	 aj=*(*(*(ato+0)+1)+j); 
         roa=GetR(*(*(R+k)+aj),*(*(R+k)+ai),Roa);
         for(i=0;i<3;i++) Roa[i]=Roa[i]/roa;
         for(alpha=0;alpha<2;alpha++) {
            for(beta=0;beta<2;beta++) {
	       gamma=(alpha==0 ? ai: aj);
	       delta=(beta ==0 ? ai: aj);
               for(mu=0;mu<3;mu++) {
                  for(nu=0;nu<3;nu++) {
                //     roa=GetR(*(*(R+k)+gamma),*(*(R+k)+delta),Roa);
                     *(*(H +k*3*NAtom +gamma*3 +mu+1) +k*3*NAtom+delta*3+nu+1) +=
                        (DDPhiIR(roa,coe,j)-DPhiIR(roa,coe,j)/roa)*
                         Roa[mu]*Roa[nu]*
			 DDt(gamma,ai,aj)*DDt(delta,ai,aj)

                       + DPhiIR(roa,coe,j)*Dt(mu,nu)*
		         DDt(gamma,ai,aj)*DDt(delta,ai,aj)/roa;
                  }
               }
            }
         }
      }
   }
}
void SetHessianLJ14(long NMolecules,double ***R,double **H) {
   
   long k,gamma,delta,alpha,beta,j,mu,nu,i,ai,aj;
   double Roa[3],roa;
   
   for(k=0;k<NMolecules;k++) {
      for(j=0;j<np[3];j++) {
	 ai=*(*(*(ato+3)+0)+j); 
	 aj=*(*(*(ato+3)+1)+j); 
         roa=GetR(*(*(R+k)+ai),*(*(R+k)+aj),Roa);
         for(i=0;i<3;i++) Roa[i]=Roa[i]/roa;
         for(alpha=0;alpha<2;alpha++) {
            for(beta=0;beta<2;beta++) {
	       gamma=(alpha==0 ? ai: aj);
	       delta=(beta ==0 ? ai: aj);
               for(mu=0;mu<3;mu++) {
                  for(nu=0;nu<3;nu++) {
                     *(*(H +k*3*NAtom +gamma*3 +mu+1) +k*3*NAtom+delta*3+nu+1) +=
                        (DDPhiLJ14(roa,coe,j)-DPhiLJ14(roa,coe,j)/roa)*
                         Roa[mu]*Roa[nu]*
                        (Dt(gamma,ai)-Dt(gamma,aj))*
                        (Dt(delta,ai)-Dt(delta,aj)) 
                       + DPhiLJ14(roa,coe,j)*Dt(mu,nu)*
                        (Dt(gamma,ai)-Dt(gamma,aj))*
                        (Dt(delta,ai)-Dt(delta,aj))/roa;
                  }
               }
            }
         }
      }
   }
}


double potLJSR(long NMolecules,double ***R){
   
   long k,gamma,mu,l,delta,nu,i,j,lig;
   long ai,aj,ak,al;
   double rji,Rji[3],rkl,Rkl[3],rkj,rjk,Rkj[3],Rjk[3],r;
   double Rjit[3],Rklt[3],Rjin[3],Rkln[3];
   double sij,eij,tmp,dot,the,the0,thedot,cc,c[5],nji,nkl;
   tmp=0.0;
   for(k=0;k<NMolecules;k++) {
      lig=k+1;
      for(j=0;j<np[2];j++) {
	  ai=*(*(*(ato+2)+0)+j); //0,1,2,3 ;IR,IA,ID,LJ14
	  aj=*(*(*(ato+2)+1)+j);
	  ak=*(*(*(ato+2)+2)+j);
	  al=*(*(*(ato+2)+3)+j);
	  //the0= *(*(*(coe+1)+0)+j);
	  //cc  = *(*(*(coe+1)+1)+j);
          //the=acos(CalY(R,k,j)); 
	 // the0 *= M_PI/180.0;
	  //printf("the=%8.2f ,the0=%8.2f\n",the,the0);

	  for(i=0;i<4;i++)c[i+1]=*(*(*(coe+2)+i)+j);
	  the=acos(CalT(k,j));
	  eij=0.0;
	  for(i=1;i<=4;i++)eij+= c[i]*pow(cos(the),i*1.0);

	//  eij=0.5*(c[1]*(1.0+cos(the))+c[2]*(1.0-cos(2.0*the))+c[3]*(1.0+cos(3.0*the))+c[4]*(1.0+cos(4.0*the)));
	  printf("LIG=%5ld, the=%12.2e \n",lig,CalT(k,j));
         //riajb=GetR(*(*(R+k)+ai),*(*(R+k)+aj),Riajb);
	 //eij=PhiIR(riajb,coe,j);
	 //eij=PhiLJ14(riajb,coe,j);
//	 eij=0.5*cc*(the-the0)*(the-the0);
//	 printf("%5ld ,bond %3ld,ai=%3ld aj=%3ld r=%8.4lf eij=%8.4lf\n",k,j,ai,aj,riajb,eij);
        // for(i=0;i<3;i++) Riajb[i]=Riajb[i]/riajb;
	// for(gamma=0;gamma<NAtom;gamma++){
	//     for(delta=0;delta<NAtom;delta++){
        //	riajb=GetR(*(*(R+k)+gamma),*(*(R+l)+delta),Riajb);
	//	if(l==k && ecl[gamma][delta]==1)continue;
	// 	if(l==k && delta<gamma)continue;
	// 	if(riajb>1.2)continue;
	// 	sij=sqrt(*(sig+gamma)* *(sig+delta));
	// 	eij=sqrt(*(eps+gamma)* *(eps+delta));
	        //tmp += PhiLJ(riajb,sij,eij);
	        tmp += eij;
	 //}
	 //}
      }
   }
   return tmp;
}

void SetHessianLJNE(long NMolecules,double ***R,double **H) {
   long k,gamma,mu,l,delta,nu,i;
   double riajb,Riajb[3];
   double sij,eij;
   for(k=0;k<NMolecules;k++) {
   for(gamma=0;gamma<NAtom;gamma++){
          for(l=k;l<NMolecules;l++) {
	  for(delta=0;delta<NAtom;delta++){
              riajb=GetR(*(*(R+k)+gamma),*(*(R+l)+delta),Riajb);
	      if(l==k && ecl[gamma][delta]==1)continue;
	      if(l==k && delta<=gamma)continue;// non-equal LJ term
	      if(riajb>1.2)continue;
	      sij=sqrt(*(sig+gamma)* *(sig+delta));
	      eij=sqrt(*(eps+gamma)* *(eps+delta));

         	for(i=0;i<3;i++)Riajb[i]=Riajb[i]/riajb;
         	for(mu=0;mu<3;mu++) {
             	for(nu=0;nu<3;nu++) {
                *(*(H+k*3*NAtom+gamma*3+mu+1)+l*3*NAtom+delta*3+nu+1) -=
                   (DDPhiLJ(riajb,sij,eij)-DPhiLJ(riajb,sij,eij)/riajb)*
                   Riajb[mu]*Riajb[nu]+
                   DPhiLJ(riajb,sij,eij)*Dt(mu,nu)/riajb;
             	}
         	}
      	  }
	  }
   }
   }
}

void SetHessianLJE(long NMolecules,double ***R,double **H) {
   long k,gamma,mu,j,delta,nu,i,l;
   double riajb,Riajb[3],sij,eij;
   for(k=0;k<NMolecules;k++) {
   for(gamma=0;gamma<NAtom;gamma++){
      for(l=0;l<NMolecules;l++) {
      for(delta=0;delta<NAtom;delta++){
          riajb=GetR(*(*(R+k)+gamma),*(*(R+l)+delta),Riajb);

	  if(l==k && ecl[gamma][delta]==1)continue;
	  if(l==k && delta==gamma)continue;
	  if(riajb>1.2)continue;
	  sij=sqrt(*(sig+gamma)* *(sig+delta));
	  eij=sqrt(*(eps+gamma)* *(eps+delta));
          for(i=0;i<3;i++)Riajb[i]=Riajb[i]/riajb;

          for(mu=0;mu<3;mu++) {
             for(nu=0;nu<3;nu++) {
              *(*(H+k*3*NAtom+gamma*3+mu+1)+k*3*NAtom+gamma*3+nu+1) +=
                   (DDPhiLJ(riajb,sij,eij)-DPhiLJ(riajb,sij,eij)/riajb)*
                   Riajb[mu]*Riajb[nu]
                  +DPhiLJ(riajb,sij,eij)*Dt(mu,nu)/riajb;
             }
          }
     }
     }
   }
   }
}

void SetHessianCE(long NMolecules,double ***R,double **H) {

   long k,gamma,mu,j,delta,beta,nu,i,l;
   double riajb,Riajb[3],cg,cd;

   for(k=0;k<NMolecules;k++) {
       for(l=0;l<NMolecules;l++) {
     	for(gamma=0;gamma<NAtom;gamma++) {
       	for(delta=0;delta<NAtom;delta++) {
	     cg=Charges[gamma];cd=Charges[delta];
             riajb=GetR(*(*(R+k)+gamma),*(*(R+l)+delta),Riajb);
	      if(l==k && ecl[gamma][delta]==1)continue;
	      if(l==k && delta<=gamma)continue;// non-equal LJ term
             for(i=0;i<3;i++)Riajb[i]=Riajb[i]/riajb;
             for(mu=0;mu<3;mu++){
               for(nu=0;nu<3;nu++) {
                 *(*(H+k*3*NAtom+gamma*3+mu+1)+k*3*NAtom+gamma*3+nu+1) +=
                   (DDPhiC(cg,cd,riajb)-DPhiC(cg,cd,riajb)/riajb)*
                   Riajb[mu]*Riajb[nu]+
                   DPhiC(cg,cd,riajb)*Dt(mu,nu)/riajb;
               }
             }
       	}
     	}
       }
   }
}

void SetHessianCNE(long NMolecules, double ***R,double **H) {
   
   long k,gamma,mu,l,delta,nu,i;
   double riajb,Riajb[3],cg,cd;
   
   for(k=0;k<NMolecules;k++) {
     for(l=k;l<NMolecules;l++) {
       for(gamma=0;gamma<NAtom;gamma++) {
         for(delta=0;delta<NAtom;delta++) {
	   cg=Charges[gamma];cd=Charges[delta];  
           riajb=GetR(*(*(R+k)+gamma),*(*(R+l)+delta),Riajb);
	      if(l==k && ecl[gamma][delta]==1)continue;
	      if(l==k && delta<=gamma)continue;// non-equal LJ term
           for(i=0;i<3;i++)Riajb[i]=Riajb[i]/riajb;
           for(mu=0;mu<3;mu++) {
             for(nu=0;nu<3;nu++) {
               *(*(H+k*3*NAtom+gamma*3+mu+1)+l*3*NAtom+delta*3+nu+1) -=
                 (DDPhiC(cg,cd,riajb)-DPhiC(cg,cd,riajb)/riajb)*
                 Riajb[mu]*Riajb[nu]+
                 DPhiC(cg,cd,riajb)*Dt(mu,nu)/riajb;
             }
           }
         }
       }
     }
   }
}

double GetR(double *Ri,double *Rj,double *Rr){
   
   long i;
   double rij;
   
   rij=0.0;
   for(i=0;i<3;i++) {
      *(Rr+i)=*(Ri+i)-*(Rj+i);
      if(*(Rr+i)>HBL)   *(Rr+i)=*(Rr+i)-BL;
      if(*(Rr+i)<NHBL)  *(Rr+i)=*(Rr+i)+BL;
      rij=rij+*(Rr+i) * *(Rr+i);
   }
   return(sqrt(rij));
}

double Dt(long i,long j) {
   if(i==j) return(1.0);
   return(0.0);
}
double DDt(long a,long i,long j) {
   double tmp;
   tmp=Dt(a,i)-Dt(a,j);
   return(tmp);
}

void SetZero(double **H,long N) {
   long i,j;
   
   for(i=0;i<3*NAtom*N+1;i++)
      for(j=0;j<3*NAtom*N+1;j++)
         *(*(H+i)+j)=0.0;
}


void GetConf(FILE *Conf,long NMolecules,double ***R) {
   long i,alpha,j;
   
   for(i=0;i<NMolecules;i++){
       for(alpha=0;alpha<NAtom;alpha++){
         for(j=0;j<3;j++){
	     fscanf(Conf,"%lf",*(*(R+i)+alpha)+j);
	     //printf("%8.5f",*(*(*(R+i)+alpha)+j));
       }
       //printf("\n");	 
       }
   }
}
void GetCoeff(FILE *Coeff,long ***ato,double ***coe,long *np) {
   int s;
   long i,alpha,j,k;
   char name[20];
   for(i=0;i<4;i++){
       for(j=0;j<4;j++){
	   for(k=0;k<18;k++){
	   *(*(*(ato+i)+j)+k)=0;
	   *(*(*(coe+i)+j)+k)=0.0;
	   }
       }
   }
   for(i=0;i<4;i++){
	 fscanf(Coeff,"%ld%s",&np[i],name);
	 //printf("%ld %s\n",np[i],name);
         for(j=0;j<np[i];j++){
	     if(i==0)fscanf(Coeff,"%ld%ld%lf%lf",*(*(ato+i)+0)+j,*(*(ato+i)+1)+j,*(*(coe+i)+0)+j,*(*(coe+i)+1)+j);
	     if(i==1)fscanf(Coeff,"%ld%ld%ld%lf%lf",*(*(ato+i)+0)+j,*(*(ato+i)+1)+j,*(*(ato+i)+2)+j,*(*(coe+i)+0)+j,*(*(coe+i)+1)+j);
	     if(i==2)fscanf(Coeff,"%ld%ld%ld%ld%lf%lf%lf%lf",*(*(ato+i)+0)+j,*(*(ato+i)+1)+j,*(*(ato+i)+2)+j,*(*(ato+i)+3)+j,*(*(coe+i)+0)+j,*(*(coe+i)+1)+j,*(*(coe+i)+2)+j,*(*(coe+i)+3)+j);
	     if(i==3)fscanf(Coeff,"%ld%ld%lf%lf",*(*(ato+i)+0)+j,*(*(ato+i)+1)+j,*(*(coe+i)+0)+j,*(*(coe+i)+1)+j);
	     //printf("%8.5f",*(*(*(R+i)+alpha)+j));
	   //  printf("%3ld%3ld%3ld%3ld%15.5e%15.5e%15.5e%15.5e\n",*(*(*(ato+i)+0)+j),*(*(*(ato+i)+1)+j),*(*(*(ato+i)+2)+j),*(*(*(ato+i)+3)+j),*(*(*(coe+i)+0)+j),*(*(*(coe+i)+1)+j),*(*(*(coe+i)+2)+j),*(*(*(coe+i)+3)+j));
       }
       }
}

void GetParam(double *Charges,double *M,double *sig,double *eps,long **ecl) {
   long i,j,d;
   for(i=0;i<NAtom;i++){
   //getAdouble("Mass of H? ",M+i);
   scanf("%lf%lf%lf%lf",Charges+i,M+i,sig+i,eps+i);
   //printf("%lf\n",*(M+i));
   *(M+i)=1.0/sqrt(*(M+i));
   }
   for(i=0;i<NAtom;i++){
     for(j=0;j<NAtom;j++)scanf("%ld",*(ecl+i)+j);
     //for(j=0;j<NAtom;j++)printf("%2ld",*(*(ecl+i)+j));
     //printf("\n");
   }
}
