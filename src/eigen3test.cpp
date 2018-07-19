#include "solver.h"
#include "problem.h"
#include "solution.h"
#include <QTableView>
#include <QApplication>
#include <QCommandLineParser>
#include "parameters\parametersparser.h"
#include "matrixview.h"
#include <iostream>
#include <Eigen/SparseCholesky>
#include "twodim/problem2D.h"
#include "threedim/problem3D.h"
#include "observations.h"
#include "cuda/HighResClock.h"
#include "cuda/solvercuda.h"
#include "cuda/solvercublas.h"

void saveVals(const char* fname, matrix &mat, bool symm = false) {
	std::ofstream myfile;
	myfile.open(fname, std::ios::binary);
	for (int i = 0; i < mat.rows(); i++) {
		for (int j = 0; j < mat.cols(); j++) {
			double valre = j < i && symm ? mat.coeff(i, j) : mat.coeff(j, i);
			double valim = 0.0;
			myfile.write((char*)&valre, sizeof(double));
			myfile.write((char*)&valim, sizeof(double));
		}
	}
	myfile.close();
}

void saveVals(const char* fname, const Eigen::VectorXd &vec) {
	std::ofstream myfile;
	myfile.open(fname, std::ios::binary); myfile;
	for (int i = 0; i < vec.size(); i++) {
		double valre = vec.coeff(i);
		double valim = 0.0;
		myfile.write((char*)&valre, sizeof(double));
		myfile.write((char*)&valim, sizeof(double));
	}
	myfile.close();
}

int main(int argc, char *argv[])
{
	QApplication app(argc, argv);
	QApplication::setApplicationName("EIT Annealing Test");
	QApplication::setApplicationVersion("0.1");

	// --> Parse command line arguments
	QCommandLineParser parser;
	parser.setApplicationDescription("EIT Annealing Test.");
	EitAnnealingArgs params;
	QString errorMessage;
	switch (parseCommandLine(parser, &params, &errorMessage)) {
	case CommandLineOk:
		break;
	case CommandLineError:
		fputs(qPrintable(errorMessage), stderr);
		fputs("\n\n", stderr);
		fputs(qPrintable(parser.helpText()), stderr);
		return 1;
	case CommandLineVersionRequested:
		printf("%s %s\n", qPrintable(QCoreApplication::applicationName()),
			qPrintable(QCoreApplication::applicationVersion()));
		return 0;
	case CommandLineHelpRequested:
		parser.showHelp();
		Q_UNREACHABLE();
	}

	observations<double> *readings = new observations<double>;

	bool is2dProblem;
	std::string meshfname = params.inputMesh.toStdString();
	std::string currentsfname = params.inputCurrents.toStdString();
	std::string tensionsfname = params.inputTensions.toStdString();
	std::shared_ptr<problem> input = problem::createNewProblem(meshfname.c_str(), is2dProblem);
	input->setGroundNode(params.ground);
	input->initProblem(meshfname.c_str());
	const char *currentsfnamecstr = currentsfname.c_str();
	readings->initObs(&currentsfnamecstr, tensionsfname.c_str(), input->getNodesCount(), input->getGenericElectrodesCount());
	input->buildNodeCoefficients();
	input->prepareSkeletonMatrix();
	input->createCoef2KMatrix();
	
    matrix *m;
	Eigen::VectorXd vcond(input->getNumCoefficients());
	for (int i = 0; i < vcond.rows(); i++) vcond[i] = 0.3815;
	input->assembleProblemMatrix(&vcond[0], &m);
	input->postAssembleProblemMatrix(&m);

    QTableView matrixView;
    matrixView.setModel(makeMatrixTableModel(m->selfadjointView<Eigen::Lower>()));
    matrixView.setWindowTitle("Stiffness");
    matrixView.show();

	Eigen::VectorXd currents;
	QTableView vectorbView;
	vectorbView.setWindowTitle("Currents");

	Eigen::VectorXd x;
	QTableView vectorView;
	vectorView.setWindowTitle("Tensions");

	//saveVals("A.txt", *m, true);

	// Create symmetric matrix with filled upper left and float values
	std::vector<Eigen::Triplet<Scalar>> tripletList;
	m->makeCompressed();
	int col = -1;
	for (int i = 0; i < m->nonZeros(); ++i) {
		if (m->outerIndexPtr()[col + 1] <= i) col++;
		int row = m->innerIndexPtr()[i];
		tripletList.push_back(Eigen::Triplet<Scalar>(row, col, m->valuePtr()[i]));
		if (row != col) tripletList.push_back(Eigen::Triplet<Scalar>(col, row, m->valuePtr()[i]));
	}
	Eigen::SparseMatrix<float, Eigen::ColMajor> msymm(m->rows(), m->cols());
	msymm.setFromTriplets(tripletList.begin(), tripletList.end());
	msymm.makeCompressed();
	Cublas::Matrix *Acublas = Cublas::Matrix::createCublasMatrix(&msymm);

	// Create preconditioner
	SparseIncompleteLLT precond(*m);
	//matrix L = precond.matrixL();
	//saveVals("L.txt", L);

	// Create CUDA preconditioner
	matrix mCpjds = *m;
	std::unique_ptr<MatrixCPJDS> stiffness = std::unique_ptr<MatrixCPJDS>(new MatrixCPJDS);
	std::unique_ptr<MatrixCPJDSManager> mgr = std::unique_ptr<MatrixCPJDSManager>(CGCUDA_Solver::createManager(&mCpjds, stiffness.get(), input->getNodeCoefficients(), input->getNodesCount(), input->getNumCoefficients()));
	numType lINFinityNorm  = CGCUDA_Solver::createPreconditioner(*stiffness, stiffness->cpuData.data);

	// Create Cublas preconditioner
	Cublas::Precond *precondcublas = Cublas::Precond::createPrecond(Acublas);

	std::vector<Eigen::VectorXd> solutions, solutionscublas;
	std::vector<Eigen::VectorXf> solutionscuda;
	for (int patterno = 0; patterno < readings->getCurrentsCount(); patterno++) {
		// Get current vector
		currents = input->getCurrentVector(patterno, readings);
		#ifndef BLOCKGND
		currents[input->getGroundNode()] = 0;
		#endif
		//saveVals(("b" + std::to_string(patterno + 1) + ".txt").c_str(), currents);

		// Create CUDA current vector
		numType *currentsData = new numType[m->cols()];
		for (int i = 0; i < m->cols(); i++) currentsData[i] = currents[i];
		Vector *bVec = CGCUDA_Solver::createCurrentVector(currentsData, *mgr, stiffness->matrixData.n, m->cols());
		//saveVals(("b" + std::to_string(patterno + 1) + "_cuda.txt").c_str(), currentsData, n, 1);
		
		// Current GUI view
		vectorbView.setModel(makeMatrixTableModel(currents));
		vectorbView.show();

		// Solve the direct problem
		HighResClock::time_point ts1 = HighResClock::now();
		CG_Solver solver(*m, currents, precond);
		for (int i = 0; i < 100; i++) solver.do_iteration();
		HighResClock::time_point ts2 = HighResClock::now();
		x = solver.getX();
		//saveVals(("x" + std::to_string(patterno + 1) + ".txt").c_str(), currents);
		
		// CUDA solver for the direct problem
		HighResClock::time_point tc1 = HighResClock::now();
		CGCUDA_Solver solvercuda(stiffness.get(), mgr.get(), bVec, lINFinityNorm);
		for (int i = 0; i < 100; i++) solvercuda.do_iteration();
		//std::vector<numType> xcuda = solvercuda.getX();
		Eigen::VectorXf xcudavec = solvercuda.getX();
		HighResClock::time_point tc2 = HighResClock::now();
		//Eigen::VectorXd xcudavec(m->cols()); for (int i = 0; i <  m->cols(); i++) xcudavec[i] = xcuda[i];
		//saveVals(("x" + std::to_string(patterno + 1) + "_cuda.txt").c_str(), xcudavec);

		// Cublas solver for the direct problem
		HighResClock::time_point tb1 = HighResClock::now();
		Cublas::CG_Solver solvercublas(Acublas, currentsData, precondcublas);
		for (int i = 0; i < 100; i++) solvercublas.doIteration();
		float *xcublas = solvercublas.getX();
		HighResClock::time_point tb2 = HighResClock::now();
		Eigen::VectorXd xcublasvec(m->cols()); for (int i = 0; i <  m->cols(); i++) xcublasvec[i] = xcublas[i];
		//saveVals(("x" + std::to_string(patterno + 1) + "_cublas.txt").c_str(), xcublasvec);

		// Potential GUI view
		vectorView.setModel(makeMatrixTableModel(x.selfadjointView<Eigen::Lower>()));
		vectorView.show();
		
		// Store solutions
		solutions.push_back(x);
		solutionscuda.push_back(xcudavec);
		solutionscublas.push_back(xcublasvec);
		std::cout << "Finished solution " << patterno + 1 << " of " << readings->getCurrentsCount() << ". Times: " << std::chrono::duration_cast<std::chrono::microseconds>(ts2 - ts1).count()  <<
			"us (serial), " << std::chrono::duration_cast<std::chrono::microseconds>(tc2 - tc1).count()  << "us (cjpds-cuda), " << std::chrono::duration_cast<std::chrono::microseconds>(tb2 - tb1).count()  << "us (cublas)." << std::endl;
	}

	// Save solutions to gmsh files
	solution::savePotentials(solutions, params.outputMesh.toStdString().c_str(), input, readings);
	std::string refname(params.outputMesh.toStdString()); std::size_t dotfound = refname.find_last_of("."); refname.replace(dotfound, 1, "_cuda."); solution::savePotentials(solutionscuda, refname.c_str(), input, readings);
	std::string refname2(params.outputMesh.toStdString()); refname2.replace(dotfound, 1, "_cublas."); solution::savePotentials(solutionscublas, refname2.c_str(), input, readings);

	return app.exec();
}