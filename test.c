#include<stdio.h>
#include<stdlib.h>
int main(){
    int i,j,d;
    FILE *fp;
    fp=fopen("excels.dat","r");
    for(i=0;i<12;i++){
	for(j=0;j<12;j++){
	    fscanf(fp,"%d",&d);
	    printf("%2d",d);
	}
	printf("\n");
    }
    return 0;
}

