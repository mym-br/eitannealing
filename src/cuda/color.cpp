#include "color.h"

#include <iostream>

#define NULL 0
#define RED 1
#define BLUE 2
#define GREEN 3
#define YELLOW 4
#define PINK 5

//#define V 5
#define COLORS 5

/* A utility function to check if the current color assignment
is safe for vertex v */
bool isSafe(Eigen::SparseMatrix<numType, Eigen::ColMajor> * data, std::unique_ptr<int[]> &colors, int color, int row) {
	// works only with full matrix. Is there a fast way to do this with a lower triangular matrix?
	for (Eigen::SparseMatrix<numType, Eigen::ColMajor>::InnerIterator it(*data, row); it; ++it) {
		if (colors[it.col()] == color || colors[it.row()] == color) {
			return false;
		}
	}
	return true;
}

/* A recursive utility function to solve m coloring problem */
bool graphColoringUtil(Eigen::SparseMatrix<numType, Eigen::ColMajor> * data, std::unique_ptr<int[]> &colors, int row) {
	int dim = data->cols();
	/* base case: If all vertices are assigned a color then
	return true */
	if (row == dim) {
		return true;
	}

	for (int row = 0; row < dim; row++) {
		/* Consider this vertex v and try different colors */
		for (int c = 1; c <= dim; c++) {
			/* Check if assignment of color c to v is fine*/
			if (isSafe(data, colors, c, row)) {
				colors[row] = c;
				break;
				///* recur to assign colors to rest of the vertices */
				//if (graphColoringUtil(dim, data, colors, row + 1) == true) {
				//	return true;
				//}

				///* If assigning color c doesn't lead to a solution
				//then remove it */
				//colors[row] = 0;
			}
		}

		/* If no color can be assigned to this vertex then return false */
		//return false;
	}
	return true;
}

bool graphColoring(Eigen::SparseMatrix<numType, Eigen::ColMajor> * data, std::unique_ptr<int[]> &colors) {
	int dim = data->cols();
	// Initialize all color values as 0. This initialization is needed
	// correct functioning of isSafe()
	for (int i = 0; i < dim; i++) {
		colors[i] = 0;
	}

	// Call graphColoringUtil() for vertex 0
	if (graphColoringUtil(data, colors, 0) == false) {
		printf("Solution does not exist");
		return false;
	}

	return true;
}