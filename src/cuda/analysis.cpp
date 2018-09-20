#include "analysis.h"

#include "utils.h"
#include "color.h"
#include <iostream>

namespace analysis {
	int find(std::unique_ptr<int[]> &arr, int n, int val, int offset) {
		for (int i = offset; i < n; i++) {
			if (arr[i] == val) {
				return i;
			}
		}
		return -1;
	}
	int findMax(std::unique_ptr<int[]> &arr, int n, int offset, int span) {
		int max = offset;
		for (int i = offset; i < span + offset && i < n; i++) {
			if (arr[i] > arr[max]) {
				max = i;
			}
		}
		return max;
	}
}

void fixIndicesMapPadding(int size, std::unique_ptr<int[]> &reorderIdx, std::unique_ptr<int[]> &unorderIdx,
	int sizePadded, std::unique_ptr<int[]> &reorderIdxFixed, std::unique_ptr<int[]> &unorderIdxFixed,
	int colorCount, std::unique_ptr<int[]> &colorOff, std::shared_ptr<int[]> colorsOffPadded) {
	for (int i = 0; i < size; i++) { // initiating
		unorderIdx[i] = -1;
	}
	for (int i = 0; i < size; i++) { // initiating
		unorderIdx[reorderIdx[i]] = i;
	}
	for (int i = 0; i < size; i++) { // initiating
		if (unorderIdx[i] < 0) {
			std::cout << "damn, something went wrong!" << std::endl;
		}
	}

	for (int i = 0; i < size; i++) {
		for (int j = 0; j < colorCount; j++) {
			if (unorderIdx[i] >= colorOff[j] && unorderIdx[i] < colorOff[j + 1]) {
				unorderIdxFixed[i] = unorderIdx[i] + (colorsOffPadded[j] - colorOff[j]);
				break;
			}
		}
	}
	for (int i = 0; i < sizePadded; i++) { // initiating
		reorderIdxFixed[i] = -1;
	}
	for (int i = 0; i < size; i++) { // translating
		reorderIdxFixed[unorderIdxFixed[i]] = i;
	}
}

void swap(Eigen::SparseMatrix<numType, Eigen::ColMajor> * arr, Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic> &perm,  int a, int b, std::unique_ptr<int[]> &newIdx) {
	int permCoeffA = perm.indices().coeff(a);
	perm.indices().coeffRef(a) = perm.indices().coeff(b);
	perm.indices().coeffRef(b) = permCoeffA;

	if (newIdx != NULL) {
		int aux_int = newIdx[a];
		newIdx[a] = newIdx[b];
		newIdx[b] = aux_int;
	}
}

std::unique_ptr<int[]> colorSort(Eigen::SparseMatrix<numType, Eigen::ColMajor> * data, int &colorCount, std::unique_ptr<int[]> &newIdx) {
	return colorSort(data, colorCount, newIdx, true);
}

std::unique_ptr<int[]> colorSort(Eigen::SparseMatrix<numType, Eigen::ColMajor> * data, int &colorCount, std::unique_ptr<int[]> &newIdx, bool reorderByLowerTriangular) {
	int n = data->cols();
	for (int i = 0; i < n; i++) {
		newIdx[i] = i;
	}

	std::unique_ptr<int[]> colors(new int[n]);
	bool result = graphColoring(data, colors);
	if (!result) return NULL;
	//LOGV2(colors, n, "Color:", LOGCPU);

	// find max color index (which equals to the number of colors)
	int count = 0;
	for (int i = 0; i < n; i++) {
		if (colors[i] >(colorCount)) {
			(colorCount) = colors[i];
		}
	}

	// count each color's size
	std::unique_ptr<int[]> colorsCount(new int[(colorCount)]);
	for (int i = 0; i < (colorCount); i++) {
		colorsCount[i] = 0;
	}
	for (int i = 0; i < n; i++) {
		colorsCount[colors[i] - 1]++;
	}
	//LOGV2(colorsCount, (*colorCount), "Colors count:", LOGCPU);

	// reorder colors by count
	std::list<int> reorderL(1, 0);
	for (int i = 1; i < (colorCount); i++) {
		int size = reorderL.size();
		for (std::list<int>::iterator it = reorderL.begin(); it != reorderL.end(); it++) {
			if (colorsCount[i] > colorsCount[(*it)]) {
				reorderL.insert(it, i);
				break;
			}
		}
		if (size == reorderL.size()) {
			reorderL.push_back(i);
		}
	}
	std::unique_ptr<int[]> reorderA(new int[(colorCount)]); int idx = 0;
	for (std::list<int>::iterator it = reorderL.begin(); it != reorderL.end(); it++, idx++) {
		reorderA[idx] = (*it) + 1;
	}
	//LOGV2(reorderA, (*colorCount), "Reordered colors:", LOGCPU);

	// populate colors sizes array
	std::unique_ptr<int[]> colorsOff(new int[(colorCount) + 1]);
	Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic> perm(data->cols());
	perm.setIdentity();
	int offset = 0;
	for (int i = 0; i < (colorCount); i++) {
		colorsOff[i] = offset;
		for (int j = 0; j < colorsCount[reorderA[i] - 1]; j++) {
			if (colors[offset + j] != reorderA[i]) {
				int idx = analysis::find(colors, n, reorderA[i], offset + j + 1);
				if (idx > -1) {
					swap(data, perm, offset + j, idx, newIdx);
					int aux = colors[offset + j];
					colors[offset + j] = colors[idx];
					colors[idx] = aux;
				}
			}
		}
		offset += colorsCount[reorderA[i] - 1];
	}
	colorsOff[(colorCount)] = n;
	*data = (*data).selfadjointView<Eigen::Lower>().twistedBy(perm.inverse());
	*data = (*data).triangularView<Eigen::Lower>();
	
	std::unique_ptr<int[]> reorderIdx(new int[n]);
	std::unique_ptr<int[]> reorderRL(new int[n]);
	for (int i = 0; i < n; i++) {
		reorderIdx[i] = -1;
		reorderRL[i] = 0;
	}
	if (reorderByLowerTriangular) { // matrix should be reordered by its lower triangular rows' length
									// calculate lower triangular RL		
		for (int col = 0; col<data->outerSize(); ++col)
			for (Eigen::SparseMatrix<numType, Eigen::ColMajor>::InnerIterator it(*data, col); it; ++it) {
				if (MOD(it.value()) > EPS) reorderRL[it.row()]++;
			}
	}
	else { // matrix should be reordered by full matrix rows' length
		   // calculate full matrix RL
		for (int col = 0; col<data->outerSize(); ++col)
			for (Eigen::SparseMatrix<numType, Eigen::ColMajor>::InnerIterator it(*data, col); it; ++it) {
				if (MOD(it.value()) > EPS) {
					reorderRL[it.row()]++;
					if (it.row() != it.col()) reorderRL[it.col()]++;
				}
			}
	}

	int maxIdx = 0;
	for (int i = 0; i < (colorCount); i++) {
		int offset = colorsOff[i],
			offsetNext = colorsOff[i + 1];
		for (int j = offset; j < offsetNext; j++) {
			maxIdx = analysis::findMax(reorderRL, n, offset, offsetNext - offset);
			if (maxIdx > -1) {
				reorderIdx[j] = maxIdx;
				reorderRL[maxIdx] = -1;
			}
		}
	}
	for (int i = 0; i < n; i++) {
		if (reorderIdx[i] == -1) {
			std::cout << "invalid reordered index" << std::endl;
		}
	}
	
	std::unique_ptr<int[]> swapIdx(new int[n]);
	perm.setIdentity();
	for (int i = 0; i < n; i++) {
		swapIdx[i] = i;
	}
	for (int i = 0; i < n; i++) {
		if (swapIdx[i] != reorderIdx[i]) { // if not desired row...
			int whereIdx = analysis::find(swapIdx, n, reorderIdx[i], 0); // find where the desired row is
			if (whereIdx != -1) { // if index was found...
				int a = i,
					b = whereIdx;
				swap(data, perm, a, b, newIdx);

				int aux = swapIdx[whereIdx];
				swapIdx[whereIdx] = swapIdx[i];
				swapIdx[i] = aux;
			}
		}
	}
	*data = (*data).selfadjointView<Eigen::Lower>().twistedBy(perm.inverse());
	*data = (*data).triangularView<Eigen::Lower>();

	return colorsOff;
}


std::unique_ptr<Eigen::SparseMatrix<numType>> fillPadding(Eigen::SparseMatrix<numType, Eigen::ColMajor> * data, int colorCount, std::unique_ptr<int[]> &colorOff, std::shared_ptr<int[]> colorsOffPadded, int &sizePadding) {
	// Copy original matrix
	int oldSize = data->cols();
	std::unique_ptr<Eigen::SparseMatrix<numType>> paddedMatrix(new Eigen::SparseMatrix<numType>(*data));

	colorsOffPadded[0] = 0;
	for (int i = 1; i < colorCount+1; i++) {
		// Round up color sizes to closest warp multiple
		colorsOffPadded[i] = colorsOffPadded[i - 1] + ceil((double)(colorOff[i] - colorOff[i - 1]) / WARP_SIZE) * WARP_SIZE;
	}
	int newSize = colorsOffPadded[colorCount];

	// Enlarge matrix due to color padding
	paddedMatrix->conservativeResize(newSize, newSize);
	// Fill new matrix (add 1 to extra diagonals)
	for (int i = oldSize; i < newSize; i++) paddedMatrix->coeffRef(i, i) = 1.0;

	// Create permutation vertex for the color padding (insert rows/collumns between colors)
	Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic> perm(newSize);
	//perm.setIdentity();
	int accPadding = 0;
	int lastZeroRowIdx = oldSize;
	for (int k = 0; k < colorCount; k++) {
		for (int i = colorOff[k]; i < colorOff[k + 1]; i++) {
			perm.indices().coeffRef(i) = colorsOffPadded[k] + i - colorOff[k];
			//std::cout << "perm(" << i + 1 << ") = " << paddedColorOff[k] + i - colorOff[k] + 1 << ";" << std::endl;
		}
		for (int i = colorOff[k+1]+ accPadding; i < colorsOffPadded[k + 1]; i++) {
			perm.indices().coeffRef(lastZeroRowIdx++) = i;
			//std::cout << "perm(" << lastZeroRowIdx << ") = " << i+1 << ";" << std::endl;
		}
		accPadding = colorsOffPadded[k + 1] - colorOff[k + 1];
	}
	// Apply permutation to matrix
	*paddedMatrix = paddedMatrix->selfadjointView<Eigen::Lower>().twistedBy(perm);
	*paddedMatrix = paddedMatrix->triangularView<Eigen::Lower>();

	sizePadding = newSize;
	return paddedMatrix;
}