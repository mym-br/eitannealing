/* File : eitdirectsolver.i */
%module eitdirectsolver
%{
/* Put headers and other declarations here */
#include "solver.h"
#include "problem.h"
#include "solution.h"
#include "observations.h"
extern int    solve();
%}

extern int    solve();