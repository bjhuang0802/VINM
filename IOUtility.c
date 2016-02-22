#include <IOUtility.h>

void getAlong(char *str,long *p) {
    printf("%s",str);
    scanf("%ld",p);
}

void getAdouble(char *str,double *p) {
    printf("%s",str);
    scanf("%lf",p);
}
FILE *getAfile(char *str,char *name,char *state) {
    FILE *fp;
    
    printf("%s",str);
    scanf("%s",name);
    fp=fopen(name,state);
    if(fp==NULL) {
        printf("I can't open %s.\n",name);
        exit(0);
    }
    return(fp);
}

void *AllocateMemory1(long Unit,long N) {
    void *p;
    
    p=(void *)malloc(Unit*N);
    if(p==NULL) {
        printf("I can't allocate memory!\n");
        exit(0);
    }
    return(p);
}

void **AllocateMemory2(long Unit,long N0,long N1) {
    long i;
    void **p;
    
    p=(void **)malloc(sizeof(void *)*N0);
    if(p==NULL) {
        printf("I can't allocate memory!\n");
        exit(0);
    }
    for(i=0;i<N0;i++)
        *(p+i)=(void *)AllocateMemory1(Unit,N1);
    return(p);
}

void ***AllocateMemory3(long Unit,long N0,long N1,long N2) {
    long i;
    void ***p;
    
    p=(void ***)malloc(sizeof(void **)*N0);
    if(p==NULL) {
        printf("I can't allocate memory!\n");
        exit(0);
    }
    for(i=0;i<N0;i++)
        *(p+i)=(void **)AllocateMemory2(Unit,N1,N2);
    return(p);
}

void ****AllocateMemory4(long Unit,long N0,long N1,long N2,long N3) {
    long i;
    void ****p;
    
    p=(void ****)malloc(sizeof(void ***)*N0);
    if(p==NULL) {
        printf("I can't allocate memory!\n");
        exit(0);
    }
    for(i=0;i<N0;i++)
        *(p+i)=(void ***)AllocateMemory3(Unit,N1,N2,N3);
    return(p);
}

void *****AllocateMemory5(long Unit,long N0,long N1,long N2,long N3,long N4) {
    long i;
    void *****p;
    
    p=(void *****)malloc(sizeof(void ****)*N0);
    if(p==NULL) {
        printf("I can't allocate memory!\n");
        exit(0);
    }
    for(i=0;i<N0;i++)
        *(p+i)=(void ****)AllocateMemory4(Unit,N1,N2,N3,N4);
    return(p);
}

void Free2(void **p,long N0) {
    long i;
    
    for(i=0;i<N0;i++)
        free(*(p+i));
    free(p);
}

void Free3(void ***p,long N0,long N1) {
    long i;
    
    for(i=0;i<N0;i++)
        Free2(*(p+i),N1);
    free(p);
}

void Free4(void ****p,long N0,long N1,long N2) {
    long i;
    
    for(i=0;i<N0;i++)
        Free3(*(p+i),N1,N2);
    free(p);
}

void Free5(void *****p,long N0,long N1,long N2,long N3) {
    long i;
    
    for(i=0;i<N0;i++)
        Free4(*(p+i),N1,N2,N3);
    free(p);
}
