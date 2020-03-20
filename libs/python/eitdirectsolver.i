/* File : eitdirectsolver.i */
%module eitdirectsolver
%{
/* Put headers and other declarations here */
#include "solver.h"
#include "problem.h"
#include "solution.h"
#include "observations.h"
extern void    init(const char *meshfilename,const  char *currentfilename);
#define SWIG_FILE_WITH_INIT
extern void    setconds(double* cond, int n);
extern void    solve(double* cond, int n, int patterno);
%}

%include "numpy.i"

%init %{
    import_array();
%}

extern void    init(const char *meshfilename,const  char *currentfilename);
%apply (double* INPLACE_ARRAY1, int DIM1) {(double* cond, int n)}
extern void    setconds(double* cond, int n);
%apply (double* INPLACE_ARRAY1, int DIM1, int patterno) {(double* potentials, int n, int patterno)}
extern void    solve(double* cond, int n, int patterno);