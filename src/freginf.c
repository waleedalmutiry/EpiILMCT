// RegisteringDynamic Symbols and defining infinity values

#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>
#include <R_ext/RS.h>

/* .Fortran calls */

extern void F77_NAME(datasimulation_f)(int *n, int *anum, int *num, int *observednum, double *observedepi,
                                     double *tmax, int *temp, double *suspar, int *nsuspar, double *powersus,
                                     double *transpar, int *ntranspar, double *powertrans, double *kernelpar, double *spark,
                                     double *delta1, double *delta2, double *suscov, double *transcov, double *cc, double *d3,
                                     double *epidat);

extern void F77_NAME(datasimulationsinr_f)(int *n, int *anum, int *num, int *observednum,
                                         double *observedepi, double *tmax, int *temp, double *suspar, int *nsuspar, double *powersus,
                                         double *transpar, int *ntranspar, double *powertrans, double *kernelpar,
                                         double *spark, double *gamma, double *deltain1, double *deltain2, double *deltanr1,
                                         double *deltanr2, double *suscov, double *transcov, double *cc, double *d3, double *epidat);


extern void F77_NAME(loglikcontilmsinr_f)(int *n, int *ninfected, int *num, int *nsuspar, int *ntranspar,
                                        double *cc, double *d3, double *epidat, double *suscov, double *transcov, double *suspar,
                                        double *powersus, double *transpar, double *powertrans, double *kernelpar, double *spark,
                                        double *gamma, double *deltain1, double *deltain2, double *deltanr1, double *deltanr2,
                                        double *likk);

extern void F77_NAME(loglikcontilm_f)(int *n, int *ninfected, int *num, int *nsuspar, int *ntranspar,
                                    double *cc, double *d333, double *epidat, double *suscov, double *transcov, double *suspar,
                                    double *powersus, double *transpar, double *powertrans, double *kernelpar,
                                    double *spark, double *deltain1, double *deltain2, double *likk);


extern void F77_NAME(mcmcsinr_f)(int *n, int *nsim, int *ni, int *temp, int *num, int *anum2, int *nsuspar,
    int *ntranspar, double *net, double *dis, double *epidat, int *blockupdate, int *priordistsuspar,
    int *priordisttranspar, int *priordistkernelparpar, int *priordistpowersus, int *priordistpowertrans,
    int *priordistsparkpar, int *priordistgammapar, double *suspar, double *suscov, double *powersus,
    double *transpar, double *transcov, double *powertrans, double *kernelpar, double *spark,
    double *gamma, double *deltain, double *deltanr, double *kernelparproposalvar,
    double *sparkproposalvar, double *gammaproposalvar, double *susproposalvar,
    double *powersusproposalvar, double *transproposalvar,
    double *powertransproposalvar, double *infperiodproposalin, double *infperiodproposalnr,
    double *priorpar1sus, double *priorpar2sus, double *priorpar1powersus, double *priorpar2powersus,
    double *priorpar1trans, double *priorpar2trans, double *priorpar1powertrans,
    double *priorpar2powertrans, double *kernelparprior, double *sparkprior, double *gammaprior,
    double *deltain2prior, double *deltanr2prior, double *susparop, double *powersusparop,
    double *transparop, double *powertransparop, double *kernelparop, double *sparkop, double *gammaop,
    double *deltain2op, double *deltanr2op, double *epidatmctim, double *epidatmcrem, double *loglik);

extern void F77_NAME(mcmcsir_f)(int *n, int *nsim, int *ni, int *num, int *anum2, int *temp, int *nsuspar,
   int *ntranspar, double *net, double *dis, double *epidat, int *blockupdate,
   int *priordistsuspar, int *priordisttranspar, int *priordistkernelparpar, int *priordistsparkpar,
   int *priordistpowersus, int *priordistpowertrans, double *suspar, double *suscov,
   double *powersus, double *transpar, double *transcov, double *powertrans, double *kernelpar,
   double *spark, double *delta1, double *kernelparproposalvar, double *sparkproposalvar,
   double *susproposalvar, double *powersusproposalvar, double *transproposalvar,
   double *powertransproposalvar, double *infperiodproposal, double *priorpar1sus,
   double *priorpar2sus, double *priorpar1powersus, double *priorpar2powersus,
   double *priorpar1trans, double *priorpar2trans, double *priorpar1powertrans,
   double *priorpar2powertrans, double *kernelparprior, double *sparkprior,
   double *delta2prior, double *susparop, double *powersusparop, double *transparop,
   double *powertransparop, double *kernelparop, double *sparkop,
   double *delta2op, double *epidatmctim, double *epidatmcrem,
   double *loglik);

void F77_SUB(infinity_value)(double *infval){
    *infval = R_PosInf;
}

static const R_FortranMethodDef FortranEntries[] = {
    {"datasimulation_f",      (DL_FUNC) &F77_NAME(datasimulation_f),        22},
    {"datasimulationsinr_f",  (DL_FUNC) &F77_NAME(datasimulationsinr_f),    25},
    {"loglikcontilmsinr_f",   (DL_FUNC) &F77_NAME(loglikcontilmsinr_f),     22},
    {"loglikcontilm_f",       (DL_FUNC) &F77_NAME(loglikcontilm_f),         19},
    {"mcmcsir_f",             (DL_FUNC) &F77_NAME(mcmcsir_f),               55},
    {"mcmcsinr_f",            (DL_FUNC) &F77_NAME(mcmcsinr_f),              64},
    {NULL, NULL, 0}
};

void R_init_EpiILMCT(DllInfo *dll){
    R_registerRoutines(dll, NULL, NULL, FortranEntries, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
