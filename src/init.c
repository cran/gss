#include <R_ext/RS.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .C calls */
extern void quad_smolyak(void *, void *, void *, void *);
extern void size_smolyak(void *, void *, void *);

/* .Fortran calls */
extern void F77_NAME(cdennewton)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(cdennewton10)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(cdenrkl)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(copu2newton)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void*);
extern void F77_NAME(coxaux)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(dcrdr)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(dmudr0)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(dnewton)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(dnewton10)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(drkl)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(dsidr0)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(dsms)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(gaussq)(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(hrkl)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(hzdaux1)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(hzdaux101)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(hzdaux2)(void *, void *, void *, void *, void *, void *);
extern void F77_NAME(hzdnewton)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(hzdnewton10)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(llrmaux)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(llrmnewton)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(llrmrkl)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(reg)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(regaux)(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(dmcdc)(void *, void *, void *, void *, void *, void *);

static const R_CMethodDef CEntries[] = {
    {"quad_smolyak", (DL_FUNC) &quad_smolyak, 4},
    {"size_smolyak", (DL_FUNC) &size_smolyak, 3},
    {NULL, NULL, 0}
};

static const R_FortranMethodDef FortranEntries[] = {
    {"cdennewton",   (DL_FUNC) &F77_NAME(cdennewton),   19},
    {"cdennewton10", (DL_FUNC) &F77_NAME(cdennewton10), 15},
    {"cdenrkl",      (DL_FUNC) &F77_NAME(cdenrkl),      20},
    {"copu2newton",  (DL_FUNC) &F77_NAME(copu2newton),  34},
    {"coxaux",       (DL_FUNC) &F77_NAME(coxaux),       16},
    {"dcrdr",        (DL_FUNC) &F77_NAME(dcrdr),        18},
    {"dmudr0",       (DL_FUNC) &F77_NAME(dmudr0),       23},
    {"dnewton",      (DL_FUNC) &F77_NAME(dnewton),      19},
    {"dnewton10",    (DL_FUNC) &F77_NAME(dnewton10),    15},
    {"drkl",         (DL_FUNC) &F77_NAME(drkl),         14},
    {"dsidr0",       (DL_FUNC) &F77_NAME(dsidr0),       20},
    {"dsms",         (DL_FUNC) &F77_NAME(dsms),         12},
    {"gaussq",       (DL_FUNC) &F77_NAME(gaussq),        9},
    {"hrkl",         (DL_FUNC) &F77_NAME(hrkl),         19},
    {"hzdaux1",      (DL_FUNC) &F77_NAME(hzdaux1),      13},
    {"hzdaux101",    (DL_FUNC) &F77_NAME(hzdaux101),    10},
    {"hzdaux2",      (DL_FUNC) &F77_NAME(hzdaux2),       6},
    {"hzdnewton",    (DL_FUNC) &F77_NAME(hzdnewton),    19},
    {"hzdnewton10",  (DL_FUNC) &F77_NAME(hzdnewton10),  17},
    {"llrmaux",      (DL_FUNC) &F77_NAME(llrmaux),      16},
    {"llrmnewton",   (DL_FUNC) &F77_NAME(llrmnewton),   19},
    {"llrmrkl",      (DL_FUNC) &F77_NAME(llrmrkl),      21},
    {"reg",          (DL_FUNC) &F77_NAME(reg),          18},
    {"regaux",       (DL_FUNC) &F77_NAME(regaux),        9},
    {"dmcdc",        (DL_FUNC) &F77_NAME(dmcdc),         6},
    {NULL, NULL, 0}
};

void R_init_gss(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, FortranEntries, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
