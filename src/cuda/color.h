#ifndef COLOR_H
#define COLOR_H

#include "settings.h"

#include <list>
#include <map>
#include <memory>

#include <Eigen/Core>
#include "../basematrix.h"

bool graphColoring(Eigen::SparseMatrix<numType, Eigen::ColMajor> * data, std::unique_ptr<int[]> &colors);

#endif /* COLOR_H */