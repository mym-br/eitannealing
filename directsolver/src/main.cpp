#include <stdio.h>
#include <iostream>
#include "cudacg.h"
#include "args/args.hxx"
#include "basematrix.h"
#include "conversions.h"
#include "incomplete_cholesky.h"
#include "solver.h"
#include "cuda/solvercuda.h"
#include "cuda/solvercublas.h"
#include <chrono>
#include <fstream>
#include <algorithm>
#include <numeric>

using namespace EITFILECONVERSIONS;

extern "C" {
	#include "mm/mmio.h"
}

#define DEFAULTMAXITS 100

struct raw_matrix {
	raw_matrix() {}
	raw_matrix(int _M, int _N) : M(_M), N(_N) {}
	void addElement(std::tuple<int, int, double> el) { elements.push_back(el); }
	std::vector<std::tuple<int, int, double>>::iterator begin() { return elements.begin(); }
	std::vector<std::tuple<int, int, double>>::iterator end() { return elements.end(); }
	std::vector<std::tuple<int, int, double>>::const_iterator cbegin() { return elements.cbegin(); }
	std::vector<std::tuple<int, int, double>>::const_iterator cend() { return elements.cend(); }
	int M, N;
	std::vector<std::tuple<int, int, double>> elements;
};
typedef std::vector<double> raw_vector;


void saveVals(const char* fname, Eigen::SparseMatrix<numType> &mat, bool symm = false) {
	FILE *f;
	if ((f = fopen(fname, "w")) == NULL) { std::cerr << "Could not open market matrix file to write" << *fname << std::endl; return; }
	EITFILECONVERSIONS::saveMtx(&mat, f, symm);
	fclose(f);
}

/*
* Get File extension from File path or File Name
*/
std::string getFileExtension(std::string filePath)
{
	// Find the last position of '.' in given string
	std::size_t pos = filePath.rfind('.');

	// If last '.' is found
	if (pos != std::string::npos) {
		// return the substring
		return filePath.substr(pos);
	}
	// In case of no extension return empty string
	return "";
}

void saveDenseVectorMtx(const std::string filename, raw_vector &vec) {
	std::ofstream myfile;
	myfile.open(filename);
	myfile << "%%MatrixMarket matrix array real general" << std::endl << vec.size() << " 1" << std::endl;
	for (auto val : vec) myfile << val << std::endl;
	myfile.close();
}

std::tuple<long, long> runEigenCGTest(raw_matrix &Araw, raw_vector &braw, raw_vector &x, double res, int maxit);
std::tuple<long, long> runCudaCGTest(raw_matrix &Araw, raw_vector &braw, raw_vector &x, double res, int maxit, bool isConsolidated);
std::tuple<long, long> runCusparseCublasCGTest(raw_matrix &Araw, raw_vector &braw, raw_vector &x, double res, int maxit);

int main(int argc, char *argv[])
{
	// Parse arguments
	args::ArgumentParser parser("This is a performance test program for the CG implementation of eitannealingtest.", "No comments.");
	args::HelpFlag help(parser, "help", "Display this help menu", { 'h', "help" });
	args::CompletionFlag completion(parser, { "complete" });
	args::ValueFlag<std::string> bfname(parser, "filename", "b vector file", { 'b' });
	args::ValueFlag<std::string> xfname(parser, "filename", "Output x filename prefix", { "output" });
	args::ValueFlag<std::string> resultsfname(parser, "compilation", "Results compilation filename prefix (append)", { "compilation" });
	args::ValueFlag<double> res(parser, "number", "Residual for CG convergence", { "res" });
	args::ValueFlag<int> maxits(parser, "number", "Maximum number of CG iterations", { "maxits" });
	args::Flag convertonly(parser, "flag", "Skip tests (only converts mesh file to mtx)", { "conversiononly" });
	args::Positional<std::string> Afname(parser, "filename", "A matrix file. Supported files: .mtx or .msh (with automatic conversion)");
	try { parser.ParseCLI(argc, argv); }
	catch (args::Completion e) { std::cout << e.what(); return 0; }
	catch (args::Help) { std::cout << parser; return 0; }
	catch (args::ParseError e) { std::cerr << e.what() << std::endl << parser; return 1; }

	// Open matrix file
	std::string fileName = args::get(Afname);
	raw_matrix A;
	std::string ext = getFileExtension(fileName);
	if (ext == ".msh") {
		// convert file
		fileName = convertMeshFile(fileName);
	}
	else if (ext != ".mtx") { std::cerr << "Wrong input file extension"; return 1; }
	
	if (fileName == "") { std::cerr << "Error converting file " << args::get(Afname); return 1; }
	if (convertonly) { std::cout << "Finished conversion part, skipping the remainder of the program"; return 0; }
	{
		FILE *f;
		if ((f = fopen(fileName.c_str(), "r")) == NULL) { std::cerr << "Could not read market matrix file " << fileName << std::endl; return 1; }
		MM_typecode matcode;
		if (mm_read_banner(f, &matcode) != 0) { std::cerr << "Could not process Matrix Market banner" << std::endl; return 1; }

		// Check matrix type
		if (mm_is_complex(matcode)) { std::cerr << "Sorry, this application does not support Market Market type: [" << mm_typecode_to_str(matcode) << "]" << std::endl; return 1; }

		// Get matrix size
		int ret_code, M, N, nz;
		if ((ret_code = mm_read_mtx_crd_size(f, &M, &N, &nz)) != 0 || M != N) { std::cerr << "Incorrect matrix size"; return 1; }
		std::cout << "Detected A matrix with size " << M << "x" << N << " from file " << fileName << std::endl;

		// Read matrix
		int row, col;  Scalar val;
		A = raw_matrix(M, N);
		for (int i = 0; i < nz; i++) {
			if (mm_is_pattern(matcode)) { fscanf(f, "%d %d\n", &row, &col, &val); val = 1.0; }
			else fscanf(f, "%d %d %lg\n", &row, &col, &val);
			A.addElement({ row - 1, col - 1, val });  /* adjust from 1-based to 0-based */
		}
		// Close file
		fclose(f);
	}

	// Create vector
	raw_vector b(A.M);
	if (bfname) { 
		FILE *f;
		MM_typecode matcode;

		// Open vector file
		if ((f = fopen(args::get(bfname).c_str(), "r")) == NULL) { std::cerr << "Could not read file " << args::get(bfname) << std::endl; return 1; }
		if (mm_read_banner(f, &matcode) != 0) { std::cerr << "Could not process Matrix Market banner.\n" << std::endl; return 1; }

		// Check vector type
		if (mm_is_complex(matcode)) { std::cerr << "Sorry, this application does not support Market Market type: [" << mm_typecode_to_str(matcode) << "]" << std::endl; return 1; }

		// Get vector size
		int ret_code, M, N, nz;
		if (
			(mm_is_coordinate(matcode) && (ret_code = mm_read_mtx_crd_size(f, &M, &N, &nz)) != 0) ||
			(mm_is_array(matcode) && (ret_code = mm_read_mtx_array_size(f, &M, &N)) != 0) ||
			N != 1) {
			std::cerr << "Incorrect vector size"; return 1; 
		}
		std::cout << "Detected b vector with size " << M << "x" << N << std::endl;

		// Read vector to Eigen type
		int row, col;  Scalar val;
		if (mm_is_coordinate(matcode)) {
			std::fill(b.begin(), b.end(), 0);
			for (int i = 0; i < nz; i++)
			{
				fscanf(f, "%d %d %lg\n", &row, &col, &val);
				b[row - 1] = val;  /* adjust from 1-based to 0-based */
			}
		}
		else {
			for (int i = 0; i < M; i++) {
				fscanf(f, "%lg\n", &val);
				b[i] = val;
			}
		}
	}
	else {
		// TODO: Generate random rhs
		std::fill(b.begin(), b.end(), 1);
	}

	raw_vector x(A.M);
	auto[analysertime, executiontime] = runEigenCGTest(A, b, x, res ? args::get(res) : -1, maxits ? args::get(maxits) : DEFAULTMAXITS);
	if(xfname) saveDenseVectorMtx(args::get(xfname) + "_serial.mtx", x);
	auto[analysertimeCuda, executiontimeCuda] = runCudaCGTest(A, b, x, res ? args::get(res) : -1, maxits ? args::get(maxits) : DEFAULTMAXITS, false);
	if (xfname) saveDenseVectorMtx(args::get(xfname) + "_cuda.mtx", x);
	auto[analysertimeCCuda, executiontimeCCuda] = runCudaCGTest(A, b, x, res ? args::get(res) : -1, maxits ? args::get(maxits) : DEFAULTMAXITS, true);
	if (xfname) saveDenseVectorMtx(args::get(xfname) + "_ccuda.mtx", x);
	auto[analysertimeCublas, executiontimeCublas] = runCusparseCublasCGTest(A, b, x, res ? args::get(res) : -1, maxits ? args::get(maxits) : DEFAULTMAXITS);
	if (xfname) saveDenseVectorMtx(args::get(xfname) + "_cusparse.mtx", x);

	if(resultsfname) {
		// Append execution times to compilation file args::get(resultsfname)
		std::ofstream outfile(args::get(resultsfname), std::ios_base::app);
		outfile << fileName << "\t" << A.N << "\t" << A.elements.size() << "\t" << analysertime << "\t" << executiontime << "\t" << analysertimeCuda  << "\t" << executiontimeCuda  << "\t" << analysertimeCCuda  << "\t" << executiontimeCCuda << "\t" << analysertimeCublas << "\t" << executiontimeCublas  << std::endl;
	}

	return 0;
}

std::tuple<long, long> runEigenCGTest(raw_matrix &Araw, raw_vector &braw, raw_vector &x, double res, int maxit) {
	// Convert matrix data
	std::vector<Eigen::Triplet<Scalar>> tripletList;
	for (auto el : Araw) tripletList.push_back(Eigen::Triplet<Scalar>(std::get<0>(el), std::get<1>(el), std::get<2>(el)));
	auto t1 = std::chrono::high_resolution_clock::now();
	vectorx b(braw.size());
	for (int i = 0; i < b.size(); i++) b[i] = braw[i];

	// Create sparse matrix
	matrix A(Araw.M, Araw.N);
	A.setFromTriplets(tripletList.begin(), tripletList.end());
	A.makeCompressed();
	// Create preconditioner
	SparseIncompleteLLT L(A);
	auto t2 = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> time_analyser = t2 - t1;
	std::cout << "Serial analyser on A used " << std::chrono::duration_cast<std::chrono::microseconds>(time_analyser).count() << " us." << std::endl;

	std::cout << "Starting CG with res = " << res << " and max iterations = " << maxit << std::endl;
	res = res > 0 ? res * res : -1;
	// Create solver
	t1 = std::chrono::high_resolution_clock::now();
	CG_Solver solver(A, b, L, res);
	// Execute solver iterations
	int totalIts = solver.getIteration();
	double curRes = solver.getResidueSquaredNorm();
	//for (int i = 0; i < 100; i++) {
	while(curRes > res && totalIts < maxit) {
		solver.do_iteration();
		curRes = solver.getResidueSquaredNorm();
		totalIts++;
	}
	t2 = std::chrono::high_resolution_clock::now();

	///************************/
	///* now write out result */
	///************************/
	vectorx xeig = solver.getX();
	std::chrono::duration<double> time_executor = t2 - t1;
	std::cout << "Serial executor on A used " << std::chrono::duration_cast<std::chrono::microseconds>(time_executor).count() << " us. Final residual is " << sqrt(curRes) << " after " <<  totalIts << " iterations." << std::endl;
	for (int i = 0; i < x.size(); i++) x[i] = xeig[i];

	return std::make_tuple(std::chrono::duration_cast<std::chrono::microseconds>(time_analyser).count(), std::chrono::duration_cast<std::chrono::microseconds>(time_executor).count());
}

std::tuple<long, long> runCudaCGTest(raw_matrix &Araw, raw_vector &braw, raw_vector &x, double res, int maxit, bool isConsolidated) {
	// Convert to full matrix
	std::vector<Eigen::Triplet<numType>> tripletList;
	for (auto el : Araw) {
		tripletList.push_back(Eigen::Triplet<numType>(std::get<0>(el), std::get<1>(el), std::get<2>(el)));
		if(std::get<0>(el) != std::get<1>(el)) tripletList.push_back(Eigen::Triplet<numType>(std::get<1>(el), std::get<0>(el), std::get<2>(el)));
	}
	auto t1 = std::chrono::high_resolution_clock::now();
	// Create sparse matrix
	Eigen::SparseMatrix<numType> A(Araw.M, Araw.N);
	A.setFromTriplets(tripletList.begin(), tripletList.end());
	A.makeCompressed();
	std::unique_ptr<MatrixCPJDS> stiffness = std::unique_ptr<MatrixCPJDS>(new MatrixCPJDS);
	std::unique_ptr<MatrixCPJDSManager> mgr = std::unique_ptr<MatrixCPJDSManager>(CGCUDA_Solver::createManager(&A, stiffness.get()));
	// Create vector on GPU
	std::unique_ptr<numType[]> bdata(new numType[braw.size()]);
	for (int i = 0; i < braw.size(); i++) bdata[i] = braw[i];
	Vector *b = CGCUDA_Solver::createCurrentVector(bdata.get() , *mgr, stiffness->matrixData.n, braw.size());
	// Create preconditioner
	numType lINFinityNorm = CGCUDA_Solver::createPreconditioner(*stiffness, stiffness->cpuData.data);
	auto t2 = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> time_analyser = t2 - t1;
	std::cout << "Cuda analyser on A used " << std::chrono::duration_cast<std::chrono::microseconds>(time_analyser).count() << " us." << std::endl;
	//// Output cpjds matrices and vector to mtx files
	//Eigen::SparseMatrix<numType> mPrint = CGCUDA_Solver::getCpjdsStiffness(*stiffness, stiffness->cpuData.data); saveVals("Acpjds.mtx", mPrint, true);
	//Eigen::SparseMatrix<numType> precondPrint = CGCUDA_Solver::getCpjdsStiffness(*stiffness, stiffness->cpuData.precond); saveVals("Lcpjds.mtx", precondPrint);
	//Eigen::SparseMatrix<numType> bPrint(stiffness->matrixData.n, 1); for (int i = 0; i < braw.size(); i++) bPrint.coeffRef(mgr->original2PaddedIdx[i], 0) = braw[i]; saveVals("bcpjds.mtx", bPrint);

	std::cout << "Starting " << (isConsolidated ? "consolidated" : "") << " Cuda CG with res = " << res << " and max iterations = " << maxit << std::endl;
	res = res > 0 ? res * res : -1;
	// Create solver
	t1 = std::chrono::high_resolution_clock::now();
	CGCUDA_Solver solvercuda(stiffness.get(), mgr.get(), b, lINFinityNorm, res, isConsolidated);
	// Execute solver iterations
	int totalIts = solvercuda.getIteration();
	double curRes = solvercuda.getResidueSquaredNorm();
	//for (int i = 0; i < 100; i++) {
	while (curRes > res && totalIts < maxit) {
		solvercuda.do_iteration();
		curRes = solvercuda.getResidueSquaredNorm();
		totalIts++;
	}
	t2 = std::chrono::high_resolution_clock::now();

	///************************/
	///* now write out result */
	///************************/
	try {
		Eigen::Matrix<numType, Eigen::Dynamic, 1> xeig = solvercuda.getX();
		std::chrono::duration<double> time_executor = t2 - t1;
		std::cout << (isConsolidated ? "Consolidated " : "") << "Cuda executor on A used " << std::chrono::duration_cast<std::chrono::microseconds>(time_executor).count() << " us. Final residual is " << sqrt(curRes) << " after " << totalIts << " iterations." << std::endl;
		for (int i = 0; i < x.size(); i++) x[i] = xeig[i];
		return std::make_tuple(std::chrono::duration_cast<std::chrono::microseconds>(time_analyser).count(), std::chrono::duration_cast<std::chrono::microseconds>(time_executor).count());
	}
	catch (const std::exception& e) { 
		for (int i = 0; i < x.size(); i++) x[i] = 0;
		std::cerr << "Failed to process " << (isConsolidated ? "consolidated " : "") << "cuda executor. Message: " << e.what() << std::endl;
		return std::make_tuple(-1,-1);
	}
}

std::tuple<long, long> runCusparseCublasCGTest(raw_matrix &Araw, raw_vector &braw, raw_vector &x, double res, int maxit) {
	std::vector<int> unsorted_mmrow;
	std::vector<int> unsorted_J;
	std::vector<float> unsorted_val;
	{
		int row, col;
		double curVal;
		for (auto el : Araw) {
			row = std::get<0>(el);
			col = std::get<1>(el);
			curVal = std::get<2>(el);
			unsorted_mmrow.push_back(row); unsorted_J.push_back(col); unsorted_val.push_back(curVal);
			if (col != row) { unsorted_mmrow.push_back(col); unsorted_J.push_back(row); unsorted_val.push_back(curVal); }
		}
	}
	int nz = unsorted_val.size();
	std::vector<std::size_t> permvec(nz);
	std::iota(permvec.begin(), permvec.end(), 0);
	std::sort(permvec.begin(), permvec.end(), [&](std::size_t i, std::size_t j) { if (unsorted_mmrow[i] == unsorted_mmrow[j]) return unsorted_J[i] < unsorted_J[j]; return unsorted_mmrow[i] < unsorted_mmrow[j]; });

	std::vector<int> mm_row(nz);
	std::transform(permvec.begin(), permvec.end(), mm_row.begin(), [&](std::size_t i) { return unsorted_mmrow[i]; });
	std::vector<int> J(nz);
	std::transform(permvec.begin(), permvec.end(), J.begin(), [&](std::size_t i) { return unsorted_J[i]; });
	std::vector<float> val(nz);
	std::transform(permvec.begin(), permvec.end(), val.begin(), [&](std::size_t i) { return unsorted_val[i]; });

	std::vector<int> I(Araw.N + 1);
	I[0] = 0;
	int currow = 0;
	int curJ = 0;
	for (auto curmm_row : mm_row) {
		while (currow < curmm_row) {
			I[++currow] = I[currow - 1] + curJ;
			curJ = 0;
		}
		curJ++;
	}
	I[++currow] = I[currow - 1] + curJ;

	std::vector<float> rhs(Araw.N);
	for (int i = 0; i < braw.size(); i++) rhs[i] = braw[i];

	std::vector<float> xcublas(Araw.N);
	auto ans = runCusparseCublasCG(I, J, val, rhs, xcublas, Araw.M, Araw.N, nz, res, maxit);

	for (int i = 0; i < x.size(); i++) x[i] = xcublas[i];

	return ans;
}