#include <stdio.h>
#include <stdlib.h>

void *AllocateMemory1(long,long);
void **AllocateMemory2(long,long,long);
void ***AllocateMemory3(long,long,long,long);
void ****AllocateMemory4(long,long,long,long,long);
void *****AllocateMemory5(long,long,long,long,long,long);
void Free2(void **,long);
void Free3(void ***,long,long);
void Free4(void ****,long,long,long);
void Free5(void *****,long,long,long,long);
void getAlong(char *,long *);
void getAdouble(char *,double *);
FILE *getAfile(char *,char *,char *);
