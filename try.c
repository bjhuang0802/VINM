#include<stdio.h>
#include<stdlib.h>
int main(){
    int i;
    for(i=0;i<10;i++)printf("%d %d\n",i,i==1?99: 0);
    return 0;
} 
