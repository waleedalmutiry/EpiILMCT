// RegisteringDynamic Symbols and defining infinity values

#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* .Fortran calls */
extern void F77_NAME(datasimulation)(int *, int *, int *, int *, double *, double *, double *, int *,
         double *, double *, int *, double *, double *, double *, double *, double *,
         double *, double *, double *, double *, double *);
extern void F77_NAME(datasimulationsinr)(int *, int *, int *, int *, double *, double *, double *, int *,
             double *, double *, int *, double *, double *, double *, double *, double *,
             double *, double *, double *, double *, double *, double *, double *, double *);
extern void F77_NAME(loglikcontilmsinr)(int *, int *, int *, int *, int *, double *, double *, double *,
            double *, double *, double *, double *, double *, double *, double *, double *,
            double *, double *, double *, double *, double *, double *);
extern void F77_NAME(loglikcontilm)(int *, int *, int *, int *, int *, double *, double *, double *,
            double *, double *, double *, double *, double *, double *, double *, double *,
            double *, double *, double *);
extern void F77_NAME(mcmcsir)(int *, int *, int *, int *, int *, int *, int *, double *,
              double *, double *, int *, int *, int *, int *, int *, int *,
              int *, double *, double *, double *, double *, double *, double *, double *,
              double *, double *, double *, double *, double *, double *, double *, double *,
              double *, double *, double *, double *, double *, double *, double *, double *,
              double *, double *, double *, double *, double *, double *, double *, double *,
              double *, double *, double *, double *, double *, double *, double *, double *,
              double *);
extern void F77_NAME(mcmcsinr)(int *, int *, int *, int *, int *, int *, int *, double *,
               double *, double *, int *, int *, int *, int *, int *, int *,
               int *, int *, double *, double *, double *, double *, double *, double *,
               double *, double *, double *, double *, double *, double *, double *, double *,
               double *, double *, double *, double *, double *, double *, double *, double *,
               double *, double *, double *, double *, double *, double *, double *, double *,
               double *, double *, double *, double *, double *, double *, double *, double *,
               double *, double *, double *, double *, double *, double *, double *, double *,
               double *);


void F77_SUB(infinity_value)(double *infval){
    *infval = INFINITY;
}

static const R_FortranMethodDef FortranEntries[] = {
    {"datasimulation",      (DL_FUNC) &F77_NAME(datasimulation),        21},
    {"datasimulationsinr",  (DL_FUNC) &F77_NAME(datasimulationsinr),    24},
    {"loglikcontilmsinr",   (DL_FUNC) &F77_NAME(loglikcontilmsinr),     22},
    {"loglikcontilm",       (DL_FUNC) &F77_NAME(loglikcontilm),         19},
    {"mcmcsir",             (DL_FUNC) &F77_NAME(mcmcsir),               57},
    {"mcmcsinr",            (DL_FUNC) &F77_NAME(mcmcsinr),              65},
    {NULL, NULL, 0}
};

void R_init_EpiILMCT(DllInfo *info){
    R_registerRoutines(info, NULL, NULL, FortranEntries, NULL);
    R_useDynamicSymbols(info, FALSE);
}

