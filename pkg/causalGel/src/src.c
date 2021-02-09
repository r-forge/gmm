#include <stdlib.h>
#include <R_ext/Rdynload.h>
#include "causalGel.h"

static const R_FortranMethodDef fortranMethods[] = {
  {"cvfct", (DL_FUNC) &F77_SUB(cvfct),7},
  {"llr", (DL_FUNC) &F77_SUB(llr),9},
  {"gridcv", (DL_FUNC) &F77_SUB(gridcv),10},
  {"brentcv", (DL_FUNC) &F77_SUB(brentcv),11},
  {"findnn", (DL_FUNC) &F77_SUB(findnn),12},  
  {NULL, NULL, 0}
};

void R_init_causalGel(DllInfo *dll)
   {
     R_registerRoutines(dll,
			NULL, NULL, 
			fortranMethods, NULL);
     R_useDynamicSymbols(dll, FALSE);
   }

