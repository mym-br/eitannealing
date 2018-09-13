#ifndef CONVERSIONS_H
#define CONVERSIONS_H

#include <string>
#include "basematrix.h"

namespace EITFILECONVERSIONS {
	std::string convertMeshFile(const std::string filename);
	void saveMtx(matrix *m, FILE *f, bool symmetric = true);
}

#endif // CONVERSIONS_H