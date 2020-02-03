#include <stdlib.h>
#include <math.h>

extern int comparedim;

extern double globaldiscr;

double fmax(double a,double b);

//int cmpkeyk(double **pt1, double **pt2);
int cmpkeyk(const void *pt1_, const void *pt2_);

int intpoints(double **pointset, int dim, int npoints, double *base);

double poly_discr(double **pointset, int dim, int npoints);
