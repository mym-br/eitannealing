#include <stdio.h>
#include <iostream>
#include "args/args.hxx"
#include "basematrix.h"
#include "incomplete_cholesky.h"
#include "solver.h"
#include <chrono>

extern "C" {
	#include "mm/mmio.h"
}

#define DEFAULTMAXITS 100

void runCGTest(matrix &A, vectorx &b, double res, int maxit, SparseIncompleteLLT &L);

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

	// Read matrix to Eigen type
	std::cout << "Creating Eigen sparse matrix" << std::endl;
	int row, col;  Scalar val;
	std::vector<Eigen::Triplet<Scalar>> tripletList;
	for (int i = 0; i < nz; i++) {
		fscanf(f, "%d %d %lg\n", &row, &col, &val);
		tripletList.push_back(Eigen::Triplet<Scalar>(row-1, col-1, val));  /* adjust from 1-based to 0-based */
	}
	auto t1 = std::chrono::high_resolution_clock::now();
	// Create sparse matrix
	matrix A(M, N);
	A.setFromTriplets(tripletList.begin(), tripletList.end());
	A.makeCompressed();
	// Create preconditioner
	SparseIncompleteLLT L(A);
	auto t2 = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> time_analyser = t2 - t1;
	std::cout << "Serial analyser on A used " << std::chrono::duration_cast<std::chrono::microseconds>(time_analyser).count() << " us." << std::endl;
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
	vectorx b(M);
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
		b.setZero();
		for (int i = 0; i < nz; i++)
		{
			fscanf(f, "%d %d %lg\n", &row, &col, &val);
			b[row - 1] = val;  /* adjust from 1-based to 0-based */
		}
	}
	else {
		// TODO: Generate random rhs
		b.setOnes();
	}

	///************************/
	///* now write out vector */
	///************************/
	//mm_write_banner(stdout, matcode);
	//mm_write_mtx_crd_size(stdout, M, N, nz);
	//std::cout << b.transpose() << std::endl;

	runCGTest(A, b, res ? args::get(res) : -1, maxits ? args::get(maxits) : DEFAULTMAXITS, L);
	
	return 0;
}

void runCGTest(matrix &A, vectorx &b, double res, int maxit, SparseIncompleteLLT &L) {
	// TODO: Read res and maxit parameters
	std::cout << "Starting CG with res = " << res << " and max iterations = " << maxit << std::endl;

	// Create solver
	CG_Solver solver(A, b, L);

	// Execute solver iterations
	auto t1 = std::chrono::high_resolution_clock::now();
	int totalIts = 3; 
	double curRes = std::numeric_limits<double>::max();
	//for (int i = 0; i < 100; i++) {
	while(curRes > res && totalIts < maxit) {
		solver.do_iteration();
		curRes = sqrt(solver.getResidueSquaredNorm());
		totalIts++;
	}
	auto t2 = std::chrono::high_resolution_clock::now();

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