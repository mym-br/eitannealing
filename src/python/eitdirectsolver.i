/* File : eitdirectsolver.i */
%module eitdirectsolver
%{
/* Put headers and other declarations here */
#include "solver.h"
#include "problem.h"
#include "solution.h"
#include "observations.h"
extern int    solve(const char *meshfilename,const  char *currentfilename);
%}

extern int    solve(const char *meshfilename,const  char *currentfilename);