#include<stdio.h>
#include<stdlib.h>
#include<math.h>
int main(){
    //int NS=5184,ncf=1;
    int NS=7200,ncf=1;
    int i,j,k,count=0,NN=600,PN=2000,NL=NS*ncf;
    double x,Ndos[NN],Pdos[PN],ev,dev=2.0;
    double unit=1.0;
    FILE *fp;
    fp=fopen("eg_1.dat","r");
    for(i=0;i<PN;i++)Pdos[i]=0.0;
    for(i=0;i<NN;i++)Ndos[i]=0.0;
    for(i=0;i<NL;i++){
	fscanf(fp,"%lf",&ev);
	//if(ev<1e-2 && ev>-1e-2)continue;
	if(ev>0.0){
	    k=ev/dev;
	    Pdos[k] += 1.0;
	//    count ++;
	//    printf("%8d, P %8d\n",i,k);
	}
	if(ev<0.0){
	    ev=(-ev);
	    k=ev/dev;
	    Ndos[k] += 1.0;
	//    count ++;
	//    printf("%8d, N %8d\n",i,k);
	}
    }
    for(i=1;i<PN;i++)Pdos[i] /= (1.0*dev*NL);
    for(i=1;i<NN;i++)Ndos[i] /= (1.0*dev*NL);
    for(i=NN-1;i>=1;i--){
	x=-i*dev;
	printf("%10.2f%12.4e\n",x,Ndos[i]);
    }
    for(i=1;i<PN;i++){
	x=i*dev;
	printf("%10.2f%12.4e\n",x,Pdos[i]);
    }
    return 0;
}
