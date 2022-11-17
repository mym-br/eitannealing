#ifndef CONVERSIONS_H
#define CONVERSIONS_H

#include <string>
#include "basematrix.h"
#include "../src/cuda/settings.h"

namespace EITFILECONVERSIONS {
	std::string convertMeshFile(const std::string filename);
	void saveMtx(Eigen::SparseMatrix<double> *m, FILE *f, bool symmetric = true);
}

#endif // CONVERSIONS_H