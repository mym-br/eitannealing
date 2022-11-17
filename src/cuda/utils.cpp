#include "utils.h"

#include <cstdio>
#include <iostream>
#include <fstream>
#include <sstream>

#ifdef USECUDA
#include "cuda_runtime.h"
#endif

std::ofstream myfile;

void logInit(char * filename) {
#ifdef DEBUG
	myfile.open(filename);
#endif
}

void logClose() {
#ifdef DEBUG
	myfile.close();
#endif
}

void logMessage(char * msg) {
#ifdef DEBUG
	std::ostringstream str;
	str << "\n" << msg << "\n";
	str.flush();
	printf(&str.str()[0]);
	myfile << "\n" << msg << "\n";
#endif
}

void logNumber(double * data, char * msg, int mode) {
#ifdef DEBUG

#ifdef USECUDA
	if (mode == LOGAUTO || mode == LOGCUDA) {
		// copy from CUDA device memory to host memory
		double * data_h = new double[1];
		cudaMemcpy(data_h, data, (size_t) sizeof (double), cudaMemcpyDeviceToHost);
		data = data_h;
	}
#endif
	double val = data[0];
	if (MOD(val) < EPS) {
		val = 0;
	}

	std::cout << "\n" << msg << " : " << val << std::endl;
	myfile << "\n" << msg << " : " << val << "\n";
	myfile.flush();

#ifdef USECUDA
	if (mode == LOGAUTO || mode == LOGCUDA) {
		// clear host memory
		delete data;
	}
#endif
#endif
}

void logVector(double * data, int size, char * msg, int mode) {
#ifdef DEBUG

#ifdef USECUDA
	if (mode == LOGAUTO || mode == LOGCUDA) {
		// copy from CUDA device memory to host memory
		double * data_h = new double[size];
		cudaMemcpy(data_h, data, (size_t)size * sizeof (double), cudaMemcpyDeviceToHost);
		data = data_h;
	}
#endif
	double val;
	if (size < SIZETHRESH) {
		printf("\n");
		printf(msg);
		printf("\n");
		for (int i = 0; i < size; i++) {
			val = data[i];
			if (MOD(val) < EPS) {
				val = 0;
			}
			printf("%.4f\n", val);
		}
	}
	else {
		myfile << "\n" << msg << "\n";
		for (int i = 0; i < size; i++) {
			val = data[i];
			if (MOD(val) < EPS) {
				val = 0;
			}
			myfile << val << "\n";
		}
		myfile.flush();
	}
#ifdef USECUDA
	if (mode == LOGAUTO || mode == LOGCUDA) {
		// clear host memory
		delete data;
	}
#endif
#endif
}

void logVector(int * data, int size, char * msg, int mode) {
#ifdef DEBUG
#ifdef USECUDA
	if (mode == LOGAUTO || mode == LOGCUDA) {
		// copy from CUDA device memory to host memory
		int * data_h = new int[size];
		cudaMemcpy(data_h, data, (size_t)size * sizeof (int), cudaMemcpyDeviceToHost);
		data = data_h;
	}
#endif
	if (size < SIZETHRESH) {
		printf("\n");
		printf(msg);
		printf("\n");
		for (int i = 0; i < size; i++) {
			printf("%d\n", data[i]);
		}
	}
	else {
		myfile << "\n" << msg << "\n";
		for (int i = 0; i < size; i++) {
			myfile << data[i] << "\n";
		}
		myfile.flush();
	}
#ifdef USECUDA
	if (mode == LOGAUTO || mode == LOGCUDA) {
		// clear host memory
		delete data;
	}
#endif
#endif
}

void logMatrix(double * data, int rows, int cols, char * msg, int mode) {
#ifdef DEBUG
#ifdef USECUDA
	if (mode == LOGAUTO || mode == LOGCUDA) {
		// copy from CUDA device memory to host memory
		double * data_h = new double[rows * cols];
		cudaMemcpy(data_h, data, (size_t)rows * cols * sizeof (double), cudaMemcpyDeviceToHost);
		data = data_h;
	}
#endif
	double val = 0;
	if (rows < SIZETHRESH && cols < SIZETHRESH) {
		printf("\n");
		printf(msg);
		printf("\n");
		for (int i = 0; i < rows; i++) {
			for (int j = 0; j < cols; j++) {
				val = data[i * cols + j];
				if (MOD(val) < EPS) {
					val = 0;
				}
				printf("%.4f\t", val);
			}
			printf("\n");
		}
	}
	else {
		myfile << "\n" << msg << "\n";
		for (int i = 0; i < rows; i++) {
			for (int j = 0; j < cols; j++) {
				val = data[i * cols + j];
				if (MOD(val) < EPS) {
					val = 0;
				}
				myfile << val << "\t";
			}
			myfile << "\n";
		}
		myfile.flush();
	}
#ifdef USECUDA
	if (mode == LOGAUTO || mode == LOGCUDA) {
		// clear host memory
		delete data;
	}
#endif
#endif
}

void logMatrix(double * data, int rows, int * rm, char * msg, int mode) {
#ifdef DEBUG
#ifdef USECUDA
	if (mode == LOGAUTO || mode == LOGCUDA) {
		// copy from CUDA device memory to host memory
		int * rm_h = new int[rows];
		cudaMemcpy(rm_h, rm, (size_t)rows * sizeof (int), cudaMemcpyDeviceToHost);
		rm = rm_h;

		int elCount = 0;
		for (int i = 0; i < rows; i++) {
			elCount += rm_h[i];
		}
		// copy from CUDA device memory to host memory
		double * data_h = new double[elCount];
		cudaMemcpy(data_h, data, (size_t)elCount * sizeof (double), cudaMemcpyDeviceToHost);
		data = data_h;
	}
#endif
	double val = 0;
	int offset = 0;
	if (rows < SIZETHRESH) {
		printf("\n");
		printf(msg);
		printf("\n");
		for (int i = 0; i < rows; i++) {
			for (int j = 0; j < rm[i]; j++) {
				val = data[offset + j];
				if (MOD(val) < EPS) {
					val = 0;
				}
				printf("%.4f\t", val);
			}
			printf("\n");
			offset += rm[i];
		}
	}
	else {
		myfile << "\n" << msg << "\n";
		for (int i = 0; i < rows; i++) {
			for (int j = 0; j < rm[i]; j++) {
				val = data[offset + j];
				if (MOD(val) < EPS) {
					val = 0;
				}
				myfile << val << "\t";
			}
			myfile << "\n";
			offset += rm[i];
		}
		myfile.flush();
	}
#ifdef USECUDA
	if (mode == LOGAUTO || mode == LOGCUDA) {
		// clear host memory
		delete data;
		delete rm;
	}
#endif
#endif
}

void logMatrix(int * data, int rows, int cols, char * msg, int mode) {
#ifdef DEBUG
#ifdef USECUDA
	if (mode == LOGAUTO || mode == LOGCUDA) {
		// copy from CUDA device memory to host memory
		int * data_h = new int[rows * cols];
		cudaMemcpy(data_h, data, (size_t)rows * cols * sizeof (int), cudaMemcpyDeviceToHost);
		data = data_h;
	}
#endif
	if (rows < SIZETHRESH && cols < SIZETHRESH) {
		printf("\n");
		printf(msg);
		printf("\n");
		for (int i = 0; i < rows; i++) {
			for (int j = 0; j < cols; j++) {
				printf("%d\t", data[i * cols + j]);
			}
			printf("\n");
		}
	}
	else {
		myfile << "\n" << msg << "\n";
		for (int i = 0; i < rows; i++) {
			for (int j = 0; j < cols; j++) {
				myfile << data[i * cols + j] << "\t";
			}
			myfile << "\n";
		}
		myfile.flush();
	}
#ifdef USECUDA
	if (mode == LOGAUTO || mode == LOGCUDA) {
		// clear host memory
		delete data;
	}
#endif
#endif
}

void logMatrix(int * data, int rows, int * rm, char * msg, int mode) {
#ifdef DEBUG
#ifdef USECUDA
	if (mode == LOGAUTO || mode == LOGCUDA) {
		// copy from CUDA device memory to host memory
		int * rm_h = new int[rows];
		cudaMemcpy(rm_h, rm, (size_t)rows * sizeof (int), cudaMemcpyDeviceToHost);
		rm = rm_h;

		int elCount = 0;
		for (int i = 0; i < rows; i++) {
			elCount += rm_h[i];
		}
		// copy from CUDA device memory to host memory
		int * data_h = new int[elCount];
		cudaMemcpy(data_h, data, (size_t)elCount * sizeof (int), cudaMemcpyDeviceToHost);
		data = data_h;
	}
#endif
	int offset = 0;
	if (rows < SIZETHRESH) {
		printf("\n");
		printf(msg);
		printf("\n");
		for (int i = 0; i < rows; i++) {
			for (int j = 0; j < rm[i]; j++) {
				printf("%d\t", data[offset + j]);
			}
			printf("\n");
			offset += rm[i];
		}
	}
	else {
		myfile << "\n" << msg << "\n";
		for (int i = 0; i < rows; i++) {
			for (int j = 0; j < rm[i]; j++) {
				myfile << data[offset + j] << "\t";
			}
			myfile << "\n";
			offset += rm[i];
		}
		myfile.flush();
	}
#ifdef USECUDA
	if (mode == LOGAUTO || mode == LOGCUDA) {
		// clear host memory
		delete data;
		delete rm;
	}
#endif
#endif
}
