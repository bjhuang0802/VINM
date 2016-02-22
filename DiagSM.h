static double dsqrarg;
#define DSQR(a) ((dsqrarg=(a)) == 0.0 ? 0.0 : dsqrarg*dsqrarg)
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
void tqli(double *, double *, long , double **);
double pythag(double, double);
void tred2(double **, long, double *, double *);
