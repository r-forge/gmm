#include <stdlib.h>
#include <R_ext/Rdynload.h>
#include "SEIR.h"

static const R_FortranMethodDef fortranMethods[] = {
  {"sysseir", (DL_FUNC) &F77_SUB(sysseir), 9},
  {"solveseir", (DL_FUNC) &F77_SUB(solveseir), 15},
  {NULL, NULL, 0}
};

void R_init_SEIR(DllInfo *dll)
   {
     R_registerRoutines(dll,
			NULL, NULL, 
			fortranMethods, NULL);
     R_useDynamicSymbols(dll, FALSE);
   }

