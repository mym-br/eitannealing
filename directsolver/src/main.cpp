#include <stdio.h>
#include <iostream>
#include "args/args.hxx"
#include "basematrix.h"
#include "incomplete_cholesky.h"
#include "solver.h"
#include "cuda/solvercuda.h"
#include <chrono>

extern "C" {
	#include "mm/mmio.h"
}

#define DEFAULTMAXITS 100

struct raw_matrix {
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

void runEigenCGTest(raw_matrix &Araw, raw_vector &braw, double res, int maxit);
void runCudaCGTest(raw_matrix &Araw, raw_vector &braw, double res, int maxit);

int main(int argc, char *argv[])
{
	// Parse arguments
	args::ArgumentParser parser("This is a performance test program for the CG implementation of eitannealingtest.", "No comments.");
	args::HelpFlag help(parser, "help", "Display this help menu", { 'h', "help" });
	args::CompletionFlag completion(parser, { "complete" });
	args::ValueFlag<std::string> bfname(parser, "filename", "b vector file", { 'b' });
	args::ValueFlag<double> res(parser, "number", "Residual for CG convergence", { "res" });
	args::ValueFlag<int> maxits(parser, "number", "Maximum number of CG iterations", { "maxits" });
	args::Positional<std::string> Afname(parser, "filename", "A matrix file");
	try { parser.ParseCLI(argc, argv); }
	catch (args::Completion e) { std::cout << e.what(); return 0; }
	catch (args::Help) { std::cout << parser; return 0; }
	catch (args::ParseError e) { std::cerr << e.what() << std::endl << parser; return 1; }

	// Open matrix file
	FILE *f;
	if ((f = fopen(args::get(Afname).c_str(), "r")) == NULL) { std::cerr << "Could not read file " << args::get(Afname) << std::endl; return 1; }
	MM_typecode matcode;
	if (mm_read_banner(f, &matcode) != 0) { std::cerr << "Could not process Matrix Market banner.\n" << std::endl; return 1; }
	
	// Check matrix type
	if (mm_is_complex(matcode)) { std::cerr << "Sorry, this application does not support Market Market type: [" << mm_typecode_to_str(matcode) << "]" << std::endl; return 1; }

	// Get matrix size
	int ret_code, M, N, nz;
	if ((ret_code = mm_read_mtx_crd_size(f, &M, &N, &nz)) != 0 || M != N) { std::cerr << "Incorrect matrix size"; return 1; }
	std::cout << "Detected A matrix with size " << M << "x" << N << std::endl;

	// Read matrix
	int row, col;  Scalar val;
	raw_matrix A(M,N);
	for (int i = 0; i < nz; i++) {
		fscanf(f, "%d %d %lg\n", &row, &col, &val);
		A.addElement({ row - 1, col - 1, val });  /* adjust from 1-based to 0-based */
	}
	// Close file
	fclose(f);

	///************************/
	///* now write out matrix */
	///************************/
	//mm_write_banner(stdout, matcode);
	//mm_write_mtx_crd_size(stdout, M, N, nz);
	//for (int k = 0; k<A.outerSize(); ++k)
	//	for (matrix::InnerIterator it(A, k); it; ++it)
	//		fprintf(stdout, "%d %d %20.19g\n", it.row() + 1, it.col() + 1, it.value());

	// Create vector
	raw_vector b(M);
	if (bfname) { 
		// Open vector file
		if ((f = fopen(args::get(bfname).c_str(), "r")) == NULL) { std::cerr << "Could not read file " << args::get(bfname) << std::endl; return 1; }
		if (mm_read_banner(f, &matcode) != 0) { std::cerr << "Could not process Matrix Market banner.\n" << std::endl; return 1; }

		// Check vector type
		if (mm_is_complex(matcode)) { std::cerr << "Sorry, this application does not support Market Market type: [" << mm_typecode_to_str(matcode) << "]" << std::endl; return 1; }

		// Get vector size
		if ((ret_code = mm_read_mtx_crd_size(f, &M, &N, &nz)) != 0 || N != 1) { std::cerr << "Incorrect vector size"; return 1; }
		std::cout << "Detected b vector with size " << M << "x" << N << std::endl;

		// Read vector to Eigen type
		std::fill(b.begin(), b.end(), 0);
		for (int i = 0; i < nz; i++)
		{
			fscanf(f, "%d %d %lg\n", &row, &col, &val);
			b[row - 1] = val;  /* adjust from 1-based to 0-based */
		}
	}
	else {
		// TODO: Generate random rhs
		std::fill(b.begin(), b.end(), 1);
	}

	///************************/
	///* now write out vector */
	///************************/
	//mm_write_banner(stdout, matcode);
	//mm_write_mtx_crd_size(stdout, M, N, nz);
	//std::cout << b.transpose() << std::endl;

	runEigenCGTest(A, b, res ? args::get(res) : -1, maxits ? args::get(maxits) : DEFAULTMAXITS);
	runCudaCGTest(A, b, res ? args::get(res) : -1, maxits ? args::get(maxits) : DEFAULTMAXITS);

	return 0;
}

void runEigenCGTest(raw_matrix &Araw, raw_vector &braw, double res, int maxit) {
	std::cout << "Creating Eigen sparse matrix and preconditioner" << std::endl;
	
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
	// Create solver
	CG_Solver solver(A, b, L);
	// Execute solver iterations
	t1 = std::chrono::high_resolution_clock::now();
	int totalIts = 3; 
	double curRes = std::numeric_limits<double>::max();
	//for (int i = 0; i < 100; i++) {
	while(curRes > res && totalIts < maxit) {
		solver.do_iteration();
		curRes = sqrt(solver.getResidueSquaredNorm());
		totalIts++;
	}
	t2 = std::chrono::high_resolution_clock::now();

	///************************/
	///* now write out result */
	///************************/
	vectorx x = solver.getX();
	std::chrono::duration<double> time_executor = t2 - t1;
	std::cout << "Serial executor on A used " << std::chrono::duration_cast<std::chrono::microseconds>(time_executor).count() << " us. Final residual is " << curRes << " after " <<  totalIts << " iterations." << std::endl;
	//for (int i = 0; i < M; i++) {
	//	std::cout << i + 1 << "\t" << x[i] << std::endl;
	//}
}

void runCudaCGTest(raw_matrix &Araw, raw_vector &braw, double res, int maxit) {
	// Convert matrix data
	std::vector<Eigen::Triplet<Scalar>> tripletList;
	for (auto el : Araw) tripletList.push_back(Eigen::Triplet<Scalar>(std::get<0>(el), std::get<1>(el), std::get<2>(el)));
	auto t1 = std::chrono::high_resolution_clock::now();
	// Create sparse matrix
	matrix A(Araw.M, Araw.N);
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

	std::cout << "Starting Cuda CG with res = " << res << " and max iterations = " << maxit << std::endl;
	// Create solver
	CGCUDA_Solver solvercuda(stiffness.get(), mgr.get(), b, lINFinityNorm, false);

	// Execute solver iterations
	t1 = std::chrono::high_resolution_clock::now();
	int totalIts = 3;
	double curRes = std::numeric_limits<double>::max();
	//for (int i = 0; i < 100; i++) {
	while (curRes > res && totalIts < maxit) {
		solvercuda.do_iteration();
		curRes = sqrt(solvercuda.getResidueSquaredNorm());
		totalIts++;
	}
	t2 = std::chrono::high_resolution_clock::now();

	///************************/
	///* now write out result */
	///************************/
	Eigen::VectorXf x = solvercuda.getX();
	std::chrono::duration<double> time_executor = t2 - t1;
	std::cout << "Cuda executor on A used " << std::chrono::duration_cast<std::chrono::microseconds>(time_executor).count() << " us. Final residual is " << curRes << " after " << totalIts << " iterations." << std::endl;
}