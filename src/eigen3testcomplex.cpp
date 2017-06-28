#include "solver.h"
#include "solvercomplex.h"
#include "solutioncomplex.h"
#include "problem.h"
//#include "solution.h"
#include <QTableView>
#include <QApplication>
#include <QCommandLineParser>
#include "parameters\parametersparser.h"
#include "matrixview.h"
#include <iostream>
#include <Eigen/SparseCholesky>
//#include "incomplete_choleskycomplex.h"
#include "random.h"
#include "twodim/problem2D.h"
#include "threedim/problem3D.h"
#include "observations.h"

void saveVals(const char* fname, matrixcomplex &mat, bool symm = false) {
	std::ofstream myfile;
	myfile.open(fname, std::ios::binary);
	for (int i = 0; i < mat.rows(); i++) {
		for (int j = 0; j < mat.cols(); j++) {
			double valre = j < i && symm ? mat.coeff(i, j).real() : mat.coeff(j, i).real();
			double valim = j < i && symm ? mat.coeff(i, j).imag() : mat.coeff(j, i).imag();
			//double valre = mat.coeff(j, i).real();
			//double valim = mat.coeff(j, i).imag();
			myfile.write((char*)&valre, sizeof(double));
			myfile.write((char*)&valim, sizeof(double));
		}
	}
	myfile.close();
}

void saveVals(const char* fname, const Eigen::VectorXcd &vec) {
	std::ofstream myfile;
	myfile.open(fname, std::ios::binary); myfile;
	for (int i = 0; i < vec.size(); i++) {
		double valre = vec.coeff(i).real();
		double valim = vec.coeff(i).imag();
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

	observations<std::complex<double>> *readingsComplex = new observations<std::complex<double>>;

	bool is2dProblem;
	std::string meshfname = params.inputMesh.toStdString();
	std::string currentsfname = params.inputCurrents.toStdString();
	std::string tensionsfname = params.inputTensions.toStdString();
	std::shared_ptr<problem> input;// = problem::createNewProblem(meshfname.c_str(), is2dProblem);
	is2dProblem = problem::isProblem2D(meshfname.c_str());
	if (is2dProblem) input = std::shared_ptr<problem>(new problem2D(meshfname.c_str()));
	else input = std::shared_ptr<problem>(new problem3D(meshfname.c_str()));
	input->setGroundNode(params.ground);
	input->initProblem(meshfname.c_str());
	//input->initObs(currentsfname.c_str(), tensionsfname.c_str());
	readingsComplex->initObs(currentsfname.c_str(), tensionsfname.c_str(), input->getNodesCount());
	input->buildNodeCoefficients();
	input->prepareSkeletonMatrix();
	input->createCoef2KMatrix();
	input->setCapacitance(80E-12);
	input->setCurrentFreq(275000);

	matrixcomplex *m;

	unsigned long long seed = 4939495;
	init_genrand64(seed);
	std::complex<double> *sol = solutioncomplex::getNewRandomSolution(input);
	matrixcomplex *m1 = solutioncomplex::getNewStiffness(sol, &m, input);

	saveVals("A.txt", *m, true);
	saveVals("A_H.txt", (matrixcomplex)(*m1).selfadjointView<Eigen::Lower>(), false);
	//saveVals("A_H.txt", *m1, false);

    QTableView matrixView;
	//matrixView.setModel(makeMatrixTableModel(m->selfadjointView<Eigen::Lower>()));
    matrixView.setWindowTitle("Stiffness");
    matrixView.show();

	Eigen::VectorXcd currents;
	QTableView vectorbView;
	vectorbView.setWindowTitle("Currents");

	Eigen::VectorXcd x;
	QTableView vectorView;
	vectorView.setWindowTitle("Tensions");

	SparseIncompleteLLTComplex precond(*m1);
	matrixcomplex L = precond.matrixL();
	saveVals("L.txt", L);

	std::vector<Eigen::VectorXcd> solutions;
	for (int patterno = 0; patterno < readingsComplex->getCurrentsCount(); patterno++) {
		currents = input->getConjugatedCurrentVector(patterno, m, readingsComplex);

		saveVals(("b" + std::to_string(patterno + 1) + ".txt").c_str(), readingsComplex->getCurrents()[patterno]);
		saveVals(("b_H" + std::to_string(patterno + 1) + ".txt").c_str(), currents);

		//vectorbView.setModel(makeMatrixTableModel(currents));
		vectorbView.show();

		CG_SolverComplex solver(*m1, currents, precond);
		for (int i = 0; i < 100 && solver.getResidueSquaredNorm() != 0.0; i++)
			solver.do_iteration();
		
		x = solver.getX();
		saveVals(("x" + std::to_string(patterno + 1) + ".txt").c_str(), x);
		#ifndef BLOCKGND
		// Correct potentials
		std::complex<double> avg = 0;
		for (int i = input->getNodesCount() - input->getGenericElectrodesCount(); i < input->getNodesCount(); i++) avg += x[i];
		avg /= input->getGenericElectrodesCount();
		for (int i = 0; i < input->getNodesCount(); i++) x[i] -= avg;
		#endif

		//vectorView.setModel(makeMatrixTableModel(x.selfadjointView<Eigen::Lower>()));
		vectorView.show();

		//std::ofstream myfile("tensions" + std::to_string(patterno + 1) + ".txt", std::ios::binary);
		//for (int i = 0; i < x.size(); i++) {
		//	double valre = x.coeff(i).real(); double valim = x.coeff(i).imag();
		//	myfile.write((char*)&valre, sizeof(double)); myfile.write((char*)&valim, sizeof(double));
		//}
		//myfile.close();

		solutions.push_back(x);
		std::cout << "Finished solution " << patterno + 1 << " of " << readingsComplex->getCurrentsCount() << std::endl;
	}

	solutioncomplex::savePotentials(solutions, params.outputMesh.toStdString().c_str(), input, readingsComplex);

	return app.exec();
}