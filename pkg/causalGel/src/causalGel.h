#ifndef R_CAUSALGEL_H
#define R_CAUSALGEL_H

#include <R_ext/RS.h>


void F77_SUB(cvfct)(double *p, double *y, int *n,
		    double *h, int *kern, int *nh, 
		    double *cv);
		 
void F77_SUB(llr)(double *p1, double *p2, double *y2, int *n1,
		       int *n2, double *h, int *kern,
		       int *info, double *y2h);

void F77_SUB(gridcv)(double *p, double *y, int *n,
		     int *kern, double *lh, double *uh,
		     int *nh, int *info, double *h, double *cv);

void F77_SUB(brentcv)(double *p1, double *y1, int *n,
		      int *kern, double *tol, int *maxit, double *h,
		      double *cv, int *info,
		      double *minh, double *mincv);

void F77_SUB(findnn)(double *x1, double *x2, double *y2,
		     double *tol, int *n1, int *n2, int *k,
		     int *minn, double *rk, int *nn,
		     int *ind, double *y2h);

#endif		     
