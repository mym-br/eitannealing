#ifndef PYEITSOLVER_H
#define PYEITSOLVER_H

#include <string>
#include <vector>
#include <map>

std::map<std::string, int> init(const char* meshfilename, const  char* currentfilename);
std::pair<int, std::vector<double>> solveForwardProblem(std::vector<double> conds);
std::vector<double> solveFullForwardProblem(std::vector<double> conds);

#endif // PYEITSOLVER_H