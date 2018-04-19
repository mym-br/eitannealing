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
bool isSafe(int dim, int *indices, int *rl, int rowMax, int * colors, int color, int row) {
	for (int i = 0; i < dim; i++) { // searching column
		// searching for elements that depends on current row - rows that have non-zeros on the same column
		for (int k = 0; k < rl[i]; k++) { // find elements in rest of matrix with same index (same column)
			if (indices[i * rowMax + k] == row && colors[i] == color) { // if indices in rest of matrix equals current row
				return false;
			}
		}
	}
	return true;
}

/* A recursive utility function to solve m coloring problem */
bool graphColoringUtil(int dim, int *indices, int *rl, int rowMax, int * colors, int row) {
	/* base case: If all vertices are assigned a color then
	return true */
	if (row == dim) {
		return true;
	}

	for (int row = 0; row < dim; row++) {
		/* Consider this vertex v and try different colors */
		for (int c = 1; c <= dim; c++) {
			/* Check if assignment of color c to v is fine*/
			if (isSafe(dim, indices, rl, rowMax, colors, c, row)) {
				colors[row] = c;
				break;
				///* recur to assign colors to rest of the vertices */
				//if (graphColoringUtil(dim, indices, rl, rowMax, colors, row + 1) == true) {
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

/* This function solves the m Coloring problem using Backtracking.
It mainly uses graphColoringUtil() to solve the problem. It returns
false if the m colors cannot be assigned, otherwise return true and
prints assignments of colors to all vertices. Please note that there
may be more than one solutions, this function prints one of the
feasible solutions.*/
bool graphColoring(int dim, int *indices, int *rl, int rowMax, int * colors) {
	// Initialize all color values as 0. This initialization is needed
	// correct functioning of isSafe()
	/*int *color = new int[V];
	for (int i = 0; i < V; i++)
	color[i] = 0;*/

	// Call graphColoringUtil() for vertex 0
	if (graphColoringUtil(dim, indices, rl, rowMax, colors, 0) == false) {
		printf("Solution does not exist");
		return false;
	}

	return true;
}

/* A utility function to check if the current color assignment
is safe for vertex v */
bool isSafe(int dim, numType * data, int * colors, int color, int row) {
	for (int i = 0; i < dim; i++) { // searching column
		if (MOD(data[row * dim + i]) > EPS && colors[i] == color) {
			return false;
		}
	}
	return true;
}

/* A recursive utility function to solve m coloring problem */
bool graphColoringUtil(int dim, numType * data, int * colors, int row) {
	/* base case: If all vertices are assigned a color then
	return true */
	if (row == dim) {
		return true;
	}

	for (int row = 0; row < dim; row++) {
		/* Consider this vertex v and try different colors */
		for (int c = 1; c <= dim; c++) {
			/* Check if assignment of color c to v is fine*/
			if (isSafe(dim, data, colors, c, row)) {
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

/* This function solves the m Coloring problem using Backtracking.
It mainly uses graphColoringUtil() to solve the problem. It returns
false if the m colors cannot be assigned, otherwise return true and
prints assignments of colors to all vertices. Please note that there
may be more than one solutions, this function prints one of the
feasible solutions.*/
bool graphColoring(int dim, numType * data, int * colors) {
	// Initialize all color values as 0. This initialization is needed
	// correct functioning of isSafe()
	for (int i = 0; i < dim; i++) {
		colors[i] = 0;
	}

	// Call graphColoringUtil() for vertex 0
	if (graphColoringUtil(dim, data, colors, 0) == false) {
		printf("Solution does not exist");
		return false;
	}

	return true;
}

/* A utility function to check if the current color assignment
is safe for vertex v */
bool isSafe(matrix * data, int * colors, int color, int row) {
	for (int k = 0; k<data->outerSize(); ++k)
		for (matrix::InnerIterator it(*data, k); it; ++it)
		{
			if ((it.row() == row && MOD(it.value()) > EPS && colors[it.col()] == color) || 
				(it.col() == row && MOD(it.value()) > EPS && colors[it.row()] == color) ) {
				return false;
			}
		}
	return true;
}

/* A recursive utility function to solve m coloring problem */
bool graphColoringUtil(matrix * data, int * colors, int row) {
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

bool graphColoring(matrix * data, int * colors) {
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