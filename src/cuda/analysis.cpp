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
							int colorCount, std::unique_ptr<int[]> &colorOff) {
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

	int * newOff = new int[colorCount + 1];
	int newSize = 0;
	for (int i = 0; i < colorCount; i++) {
		newOff[i] = newSize;
		newSize += ceil((double)(colorOff[i + 1] - colorOff[i]) / WARP_SIZE) * WARP_SIZE;
	}
	newOff[colorCount] = newSize;
	//LOGV2(newOff, colorCount + 1, "Matrix padding insert position:", LOGCPU);

	if (newSize != sizePadded) {
		LOG("Incorrect padded size!");
	}

	for (int i = 0; i < size; i++) {
		for (int j = 0; j < colorCount; j++) {
			if (unorderIdx[i] >= colorOff[j] && unorderIdx[i] < colorOff[j + 1]) {
				unorderIdxFixed[i] = unorderIdx[i] + (newOff[j] - colorOff[j]);
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

	delete newOff;
}

void swap(matrix * arr, int a, int b, std::unique_ptr<int[]> &newIdx) {
	Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic> perm(arr->cols());
	perm.setIdentity();
	perm.indices().coeffRef(a) = b;
	perm.indices().coeffRef(b) = a;
	*arr = (*arr).selfadjointView<Eigen::Lower>().twistedBy(perm);
	*arr = (*arr).triangularView<Eigen::Lower>();

	if (newIdx != NULL) {
		int aux_int = newIdx[a];
		newIdx[a] = newIdx[b];
		newIdx[b] = aux_int;
	}
}

std::unique_ptr<int[]> colorSort(matrix * data, int &colorCount, std::unique_ptr<int[]> &newIdx) {
	return colorSort(data, colorCount, newIdx, true);
}

std::unique_ptr<int[]> colorSort(matrix * data, int &colorCount, std::unique_ptr<int[]> &newIdx, bool reorderByLowerTriangular) {
	int n = data->cols();
	for (int i = 0; i < n; i++) {
		newIdx[i] = i;
	}

	std::unique_ptr<int[]> colors(new int[n]);
	bool result = graphColoring(data, colors);
	if (!result) return NULL;
	//std::vector<int> teste = std::vector<int>(colors, colors + n);
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
	int offset = 0;
	for (int i = 0; i < (colorCount); i++) {
		colorsOff[i] = offset;
		for (int j = 0; j < colorsCount[reorderA[i] - 1]; j++) {
			if (colors[offset + j] != reorderA[i]) {
				int idx = analysis::find(colors, n, reorderA[i], offset + j + 1);
				if (idx > -1) {
					swap(data, offset + j, idx, newIdx);
					int aux = colors[offset + j];
					colors[offset + j] = colors[idx];
					colors[idx] = aux;
				}
			}
		}
		offset += colorsCount[reorderA[i] - 1];
	}
	colorsOff[(colorCount)] = n;

	std::unique_ptr<int[]> reorderIdx(new int[n]);
	std::unique_ptr<int[]> reorderRL(new int[n]);
	for (int i = 0; i < n; i++) {
		reorderIdx[i] = -1;
		reorderRL[i] = 0;
	}
	if (reorderByLowerTriangular) { // matrix should be reordered by its lower triangular rows' length
									// calculate lower triangular RL		
		for (int col = 0; col<data->outerSize(); ++col)
			for (matrix::InnerIterator it(*data, col); it; ++it) {
				if (MOD(it.value()) > EPS) reorderRL[it.row()]++;
			}
	}
	else { // matrix should be reordered by full matrix rows' length
		   // calculate full matrix RL
		for (int col = 0; col<data->outerSize(); ++col)
			for (matrix::InnerIterator it(*data, col); it; ++it) {
				if (MOD(it.value()) > EPS) {
					reorderRL[it.row()]++;
					if (it.row() != it.col()) reorderRL[it.col()]++;
				}
			}
	}
	std::vector<int> teste = std::vector<int>(reorderRL.get(), reorderRL.get() + n);

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
	for (int i = 0; i < n; i++) {
		swapIdx[i] = i;
	}
	for (int i = 0; i < n; i++) {
		if (swapIdx[i] != reorderIdx[i]) { // if not desired row...
			int whereIdx = analysis::find(swapIdx, n, reorderIdx[i], 0); // find where the desired row is
			if (whereIdx != -1) { // if index was found...
				int a = i,
					b = whereIdx;
				swap(data, a, b, newIdx);

				int aux = swapIdx[whereIdx];
				swapIdx[whereIdx] = swapIdx[i];
				swapIdx[i] = aux;
			}
		}
	}

	return colorsOff;
}


std::unique_ptr<numType[]> fillPadding(matrix * data, int colorCount, std::unique_ptr<int[]> &colorOff, int &sizePadding) {
	int size = data->cols();
	std::unique_ptr<int[]> newColorCount(new int[colorCount]);
	std::unique_ptr<int[]> newOff(new int[colorCount + 1]);
	std::unique_ptr<int[]> insCount(new int[colorCount]);
	int newSize = 0, oldSize = 0;
	for (int i = 0; i < colorCount; i++) {
		newOff[i] = newSize;
		newColorCount[i] = ceil((double)(colorOff[i + 1] - colorOff[i]) / WARP_SIZE) * WARP_SIZE;
		newSize += newColorCount[i];
		oldSize += (colorOff[i + 1] - colorOff[i]);
		insCount[i] = newColorCount[i] - (colorOff[i + 1] - colorOff[i]);
	}
	newOff[colorCount] = newSize;
	//LOGV2(newOff, colorCount + 1, "Matrix padding insert position:", LOGCPU);
	//LOGV2(insCount, colorCount, "Matrix padding count:", LOGCPU);

	std::unique_ptr<numType[]> newData(new numType[newSize * newSize]);
	for (int i = 0; i < newSize; i++) {
		for (int j = 0; j < newSize; j++) {
			newData[i * newSize + j] = 0;
		}
		newData[i * newSize + i] = 1;
	}

	for (int col = 0; col<data->outerSize(); ++col)
		for (matrix::InnerIterator it(*data, col); it; ++it)
		{
			int row = it.row();
			newData[row * newSize + col] = it.value();
			newData[col * newSize + row] = it.value();
		}
	//LOGM2(newData, newSize, newSize, "Copied padded matrix:", LOGCPU);

	// insert rows
	for (int k = colorCount - 1; k >= 0; k--) {
		for (int i = size - 1; i >= 0; i--) {
			if (i >= colorOff[k] && i < colorOff[k + 1]) {
				for (int j = 0; j < newSize; j++) {
					int newRow = (i + newOff[k] - colorOff[k]);
					if (newRow != i) {
						// move value to correct position
						newData[newRow * newSize + j] = newData[i * newSize + j];
						// reset previous position
						newData[i * newSize + j] = 0;
					}
				}
			}
		}
	}
	//LOGM2(newData, newSize, newSize, "Padded reordered matrix (rows):", LOGCPU);
	// insert columns
	for (int k = colorCount - 1; k >= 0; k--) {
		for (int j = size - 1; j >= 0; j--) {
			if (j >= colorOff[k] && j < colorOff[k + 1]) {
				for (int i = 0; i < newSize; i++) {
					int newColumn = (j + newOff[k] - colorOff[k]);
					if (newColumn != j) {
						// move value to correct position
						newData[i * newSize + newColumn] = newData[i * newSize + j];
						// reset previous position
						newData[i * newSize + j] = 0;
					}
				}
			}
		}
	}
	//LOGM2(newData, newSize, newSize, "Padded reordered matrix (cols):", LOGCPU);
	// "fix" main diagonal
	for (int i = 0; i < newSize; i++) {
		if (MOD(newData[i * newSize + i]) < EPS) {
			newData[i * newSize + i] = 1;
		}
	}
	//LOGM2(newData, newSize, newSize, "Padded reordered matrix:", LOGCPU);

	sizePadding = newSize;
	return newData;
}