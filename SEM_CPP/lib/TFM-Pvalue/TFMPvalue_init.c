#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP lazyScore(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP pv2sc(SEXP, SEXP, SEXP, SEXP);
extern SEXP sc2pv(SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
  {"lazyScore", (DL_FUNC) &lazyScore, 5},
  {"pv2sc",     (DL_FUNC) &pv2sc,     4},
  {"sc2pv",     (DL_FUNC) &sc2pv,     4},
  {NULL, NULL, 0}
};

void R_init_TFMPvalue(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}