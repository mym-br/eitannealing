#ifndef COLOR_H
#define COLOR_H

#include "settings.h"

#include <list>
#include <map>

bool graphColoring(int dim, int *indices, int *rl, int rowMax, int * colors);

bool graphColoring(int dim, numType * data, int * colors);

#endif /* COLOR_H */