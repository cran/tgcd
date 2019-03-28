#include <R_ext/RS.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Fortran calls */

extern void F77_NAME(calcam)(void *, void *, void *, void *, void *, void *, void *);

extern void F77_NAME(calcei)(void *, void *, void *);

extern void F77_NAME(qeotor)(void *, void *, void *, void *, void *, void *, void *, 
                             void *, void *, void *, void *);

extern void F77_NAME(savgol_filter)(void *, void *, void *, void *, void *, void *, void *);

extern void F77_NAME(simpeak)(void *, void *, void *, void *, void *, void *, void *, 
                             void *, void *, void *, void *);

extern void F77_NAME(tgcd_drive)(void *, void *, void *, void *, void *, void *, void *,
                                 void *, void *, void *, void *, void *, void *, void *,
                                 void *, void *, void *, void *, void *, void *, void *);

extern void F77_NAME(wrightomega)(void *, void *);






static const R_FortranMethodDef FortranEntries[] = {
    {"calcam",        (DL_FUNC) &F77_NAME(calcam),        7},
    {"calcei",        (DL_FUNC) &F77_NAME(calcei),        3},
    {"qeotor",        (DL_FUNC) &F77_NAME(qeotor),       11},
    {"savgol_filter", (DL_FUNC) &F77_NAME(savgol_filter), 7},
    {"simpeak",       (DL_FUNC) &F77_NAME(simpeak),      11},
    {"tgcd_drive",    (DL_FUNC) &F77_NAME(tgcd_drive),   21},
    {"wrightomega",   (DL_FUNC) &F77_NAME(wrightomega),   2},
    {NULL, NULL, 0}
};

void R_init_tgcd(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, NULL, FortranEntries, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
