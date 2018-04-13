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
//#include "cuda/matrix-cpjds.h"
//#include "cuda/settings.h"
//#include "cuda/solver-pcg.h"
//#include "cuda/HighResClock.h"
#include "cuda/solvercuda.h"

numType * buildMeshStiffness(matrix &mat, int m, int n) {
	numType * stiffnessData = new numType[n * n];
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < m; j++) {
			double valre = j < i ? mat.coeff(i, j) : mat.coeff(j, i);
			stiffnessData[i * n + j] = valre;
		}
	}
	return  stiffnessData;
}

//void saveVals(const char* fname, numType *mat, int m, int n) {
//	std::ofstream myfile;
//	myfile.open(fname, std::ios::binary);
//	for (int i = 0; i < n; i++) {
//		for (int j = 0; j < m; j++) {
//			double valre = mat[i *n + j];
//			double valim = 0.0;
//			myfile.write((char*)&valre, sizeof(double));
//			myfile.write((char*)&valim, sizeof(double));
//		}
//	}
//	myfile.close();
//}

//void saveVals(const char* fname, std::vector<numType> &mat, int m, int n) {
//	std::ofstream myfile;
//	myfile.open(fname, std::ios::binary);
//	for (int i = 0; i < n; i++) {
//		for (int j = 0; j < m; j++) {
//			double valre = mat[i *n + j];
//			double valim = 0.0;
//			myfile.write((char*)&valre, sizeof(double));
//			myfile.write((char*)&valim, sizeof(double));
//		}
//	}
//	myfile.close();
//
//}
void saveVals(const char* fname, std::vector<numType> vec, int n) {
	std::ofstream myfile;
	myfile.open(fname, std::ios::binary);
	for (int i = 0; i < n; i++) {
		double valre = vec[i];
		double valim = 0.0;
		myfile.write((char*)&valre, sizeof(double));
		myfile.write((char*)&valim, sizeof(double));
	}
	myfile.close();
}

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
	//for (int i = 0; i < vcond.rows(); i++) vcond[i] = 0.3815;
	//for (int i = 0; i < 32; i++) vcond[i] = 6; for (int i = 32; i < vcond.rows(); i++) vcond[i] = 0.29;
	for (int i = 0; i < vcond.rows(); i++) vcond[i] = 1;
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

	int n = (*m).rows();
	//numType * stiffnessData = buildMeshStiffness(*m, n, n);
	//saveVals("A_cuda.txt", stiffnessData, n, n);
	//MatrixCPJDS stiffness;
	//MatrixCPJDSManager mgr(stiffnessData, n);
	//mgr.buidMatrixCPJDS(&stiffness);
	//int size = stiffness.matrixData.n;
	//std::cout << "Original size: " << n << std::endl;
	//std::cout << "CPJDSMatrix size: " << size << std::endl;
	//int repetitions = 100;
	//std::cout << "\n--- Conjugate Gradients ---" << std::endl;
	//std::cout << "...running..." << std::endl;
	//Vector ** currentsVec = new Vector*[32];
	//numType * vecArr = new numType[size];
	//for (int k = 0; k < 32; k++) {
	//	for (int i = 0; i < size; i++) {
	//		vecArr[i] = 0;
	//	}
	//	for (int i = 0; i < n; i++) {
	//		currents = input->getCurrentVector(k, readings);
	//		vecArr[mgr.original2PaddedIdx[i]] = currents[i];
	//	}
	//	currentsVec[k] = new Vector(vecArr, size);
	//}
	//cudaMemcpy(stiffness.preconditionedData, stiffness.cpuData.precond, (size_t)stiffness.matrixData.elCount * sizeof(numType), cudaMemcpyHostToDevice);
	//PCGSolverCPJDS ** solvers;
	//solvers = new PCGSolverCPJDS*[32];
	//for (int i = 0; i < 32; i++) solvers[i] = new PCGSolverCPJDS(&mgr, &stiffness, currentsVec[i]);
	//HighResClock::time_point t11 = HighResClock::now();
	//for (int pad = 0; pad < 32; pad++) {
	//	solvers[pad]->init();
	//}
	//HighResClock::time_point t9 = HighResClock::now();
	//for (int pad = 0; pad < 32; pad++) {
	//	for (int i = 0; i < repetitions; i++) {
	//		solvers[pad]->doIteration();
	//	}
	//}
	//cudaDeviceSynchronize();
	//HighResClock::time_point t10 = HighResClock::now();
	//HighResClock::time_point t12 = HighResClock::now();
	//auto duration910 = std::chrono::duration_cast<std::chrono::microseconds>(t10 - t9).count();
	//std::cout << "...done in " << (((float)duration910) / repetitions) << " us!" << std::endl;
	//auto duration1112 = std::chrono::duration_cast<std::chrono::microseconds>(t12 - t11).count();
	//std::cout << "...done (+init) in " << (((float)duration1112) / repetitions) << " us!" << std::endl;

	std::vector<Eigen::VectorXd> solutions, solutioncuda;
	for (int patterno = 0; patterno < readings->getCurrentsCount(); patterno++) {
		//currents = (*m).selfadjointView<Eigen::Lower>() * input->getCurrents()[patterno];
		currents = input->getCurrentVector(patterno, readings);

		numType *currentsData = new numType[n];
		for (int i = 0; i < n; i++) currentsData[i] = currents[i];
		//saveVals(("b" + std::to_string(patterno + 1) + "_cuda.txt").c_str(), currentsData, n, 1);

		//#ifndef BLOCKGND
		//currents[input->getGroundNode()] = 0;
		//#endif

		//std::ofstream myfile;
		//myfile.open("currents" + std::to_string(patterno+1) + ".txt");
		//for (int i = 0; i < currents.size(); i++) myfile << currents.coeff(i) << std::endl;
		//myfile.close();
		//saveVals(("b" + std::to_string(patterno + 1) + ".txt").c_str(), currents);
		
		vectorbView.setModel(makeMatrixTableModel(currents));
		vectorbView.show();

		SparseIncompleteLLT precond(*m);
		matrix L = precond.matrixL();
		//saveVals("L.txt", L);
		CG_Solver solver(*m, currents, precond);
		for (int i = 0; i < 100; i++) solver.do_iteration();
		x = solver.getX();
		//#ifndef BLOCKGND
		//// Correct potentials
		//double avg = 0;
		//for (int i = input->getNodesCount() - input->getGenericElectrodesCount(); i < input->getNodesCount(); i++) avg += x[i];
		//avg /= input->getGenericElectrodesCount();
		//for (int i = 0; i < input->getNodesCount(); i++) x[i] -= avg;
		//#endif
		numType * stiffnessData = buildMeshStiffness(*m, n, n); // FIXME: more efficiently copy data, also do it only once
		CGCUDA_Solver solvercuda(stiffnessData, currentsData, input->getNodeCoefficients(), input->getNodesCount(), input->getNumCoefficients(), NULL, n);
		for (int i = 0; i < 100; i++) solvercuda.doIteration();
		std::vector<numType> xcuda = solvercuda.getX();
		Eigen::VectorXd xcudavec(n);
		for (int i = 0; i < n; i++) xcudavec[i] = xcuda[i];
		//saveVals(("x" + std::to_string(patterno + 1) + "_cuda.txt").c_str(), xcudavec);

		vectorView.setModel(makeMatrixTableModel(x.selfadjointView<Eigen::Lower>()));
		vectorView.show();

		//myfile.open("tensions" + std::to_string(patterno+1) + ".txt");
		//for (int i = 0; i < x.size(); i++) myfile << x.coeff(i) << std::endl;
		//myfile.close();
		//saveVals(("x" + std::to_string(patterno + 1) + ".txt").c_str(), currents);

		solutions.push_back(x);
		solutioncuda.push_back(xcudavec);
		std::cout << "Finished solution " << patterno + 1 << " of " << readings->getCurrentsCount() << std::endl;
	}


	solution::savePotentials(solutions, params.outputMesh.toStdString().c_str(), input, readings);

	std::string refname(params.outputMesh.toStdString());
	std::size_t dotfound = refname.find_last_of(".");
	refname.replace(dotfound, 1, "_cuda.");
	solution::savePotentials(solutioncuda, refname.c_str(), input, readings);

	return app.exec();
}