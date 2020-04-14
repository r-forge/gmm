#ifndef R_SEIR_H
#define R_SEIR_H

#include <R_ext/RS.h>


void F77_SUB(solveseir)(int *n, double *c, 
			double *sigma, double *gamma, 
			double *rzero, double *y0,
			double *x, double *y, double *h);
		 
void F77_SUB(sysseir)(double *x, double *y, double *par,
		   double *dy);

#endif		     
