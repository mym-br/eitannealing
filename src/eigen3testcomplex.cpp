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
	std::string currentsinfname = params.inputCurrents.toStdString();
	std::string currentsoutfname = params.inputCurrentsOut.toStdString();
	std::string tensionsfname = params.inputTensions.toStdString();
	std::shared_ptr<problem> input = problem::createNewProblem(meshfname.c_str(), is2dProblem);
	input->setGroundNode(params.ground);
	input->initProblem(meshfname.c_str());
	const char **currentspair = new const char*[2]; currentspair[0] = currentsinfname.c_str(); currentspair[1] = currentsoutfname.c_str();
	readingsComplex->initObs(currentspair, tensionsfname.c_str(), input->getNodesCount(), input->getGenericElectrodesCount());
	input->buildNodeCoefficients();
	input->prepareSkeletonMatrix();
	input->createCoef2KMatrix();
	input->setCapacitance(80E-12);
	input->setCurrentFreq(275000);

	matrixcomplex *m;

	unsigned long long seed = 4939495;
	init_genrand64(seed);
	std::vector<std::complex<double>> electrodesCoeffs;
	electrodesCoeffs.push_back(std::complex<double>(8063.9, 5.82505e-008));
	electrodesCoeffs.push_back(std::complex<double>(7684.01, 2.64247e-008));
	electrodesCoeffs.push_back(std::complex<double>(6543.7, 7.99271e-008));
	electrodesCoeffs.push_back(std::complex<double>(9286.28, 2.3217e-008));
	electrodesCoeffs.push_back(std::complex<double>(4122.22, 6.75443e-008));
	electrodesCoeffs.push_back(std::complex<double>(9368.93, 4.48715e-008));
	electrodesCoeffs.push_back(std::complex<double>(8894.94, 5.52865e-008));
	electrodesCoeffs.push_back(std::complex<double>(1448.95, 7.93369e-008));
	electrodesCoeffs.push_back(std::complex<double>(7445.41, 2.7454e-008));
	electrodesCoeffs.push_back(std::complex<double>(9393.67, 4.66246e-009));
	electrodesCoeffs.push_back(std::complex<double>(9007.99, 4.6255e-008));
	electrodesCoeffs.push_back(std::complex<double>(7695.65, 5.22615e-008));
	electrodesCoeffs.push_back(std::complex<double>(8644.59, 3.46799e-008));
	electrodesCoeffs.push_back(std::complex<double>(8428.48, 8.74252e-008));
	electrodesCoeffs.push_back(std::complex<double>(9367.94, 1.45588e-008));
	electrodesCoeffs.push_back(std::complex<double>(9303.94, 3.81534e-008));
	electrodesCoeffs.push_back(std::complex<double>(7815.82, 9.14307e-008));
	electrodesCoeffs.push_back(std::complex<double>(9832.13, 7.96253e-008));
	electrodesCoeffs.push_back(std::complex<double>(5466.56, 9.52926e-008));
	electrodesCoeffs.push_back(std::complex<double>(8004.29, 9.38392e-008));
	electrodesCoeffs.push_back(std::complex<double>(9994.72, 5.12151e-008));
	electrodesCoeffs.push_back(std::complex<double>(9164.13, 7.21611e-009));
	electrodesCoeffs.push_back(std::complex<double>(8979.89, 5.32692e-009));
	electrodesCoeffs.push_back(std::complex<double>(7348.07, 9.7583e-008));
	electrodesCoeffs.push_back(std::complex<double>(8796.92, 4.86889e-008));
	electrodesCoeffs.push_back(std::complex<double>(8767.28, 7.96476e-008));
	electrodesCoeffs.push_back(std::complex<double>(9015.4, 4.61737e-008));
	electrodesCoeffs.push_back(std::complex<double>(9242.83, 8.51366e-008));
	electrodesCoeffs.push_back(std::complex<double>(7116.48, 7.5862e-008));
	electrodesCoeffs.push_back(std::complex<double>(1114.67, 7.72724e-008));
	electrodesCoeffs.push_back(std::complex<double>(7003.08, 4.15017e-008));
	electrodesCoeffs.push_back(std::complex<double>(9617.3, 5.9162e-008));
	std::complex<double> *sol = solutioncomplex::getNewRandomSolution(input, electrodesCoeffs);
	matrixcomplex *m1 = solutioncomplex::getNewStiffness(sol, &m, input);

	matrixcomplex Afull = (*m) + ((matrixcomplex)((matrixcomplex)((*m).selfadjointView<Eigen::Lower>())).triangularView<Eigen::StrictlyUpper>()).conjugate();
	matrixcomplex Atfull = Afull.conjugate();
	saveVals("A.txt", *m, true);
	//matrixcomplex m_h = (matrixcomplex)(*m1).selfadjointView<Eigen::Lower>();
	saveVals("A_H.txt", (matrixcomplex)(*m1).selfadjointView<Eigen::Lower>(), false);
	//saveVals("A_H.txt", *m1, false);

    QTableView matrixView;
	//matrixView.setModel(makeMatrixTableModel(m->selfadjointView<Eigen::Lower>()));
    matrixView.setWindowTitle("Stiffness");
    matrixView.show();

	Eigen::VectorXcd currents;
	QTableView vectorbView;
	vectorbView.setWindowTitle("Currents");

	Eigen::VectorXcd x;// , x2;
	QTableView vectorView;
	vectorView.setWindowTitle("Tensions");

	SparseIncompleteLLTComplex precond(*m1);
	//SparseIncompleteLLTComplex2 precond2(*m);
	matrixcomplex L = precond.matrixL();
	//matrixcomplex L2 = precond2.matrixL();
	saveVals("L.txt", L);
	//saveVals("L2.txt", L2);

	std::vector<Eigen::VectorXcd> solutions, solutions2;
	for (int patterno = 0; patterno < readingsComplex->getCurrentsCount(); patterno++) {
		currents = input->getConjugatedCurrentVector(patterno, m, readingsComplex);

		saveVals(("b" + std::to_string(patterno + 1) + ".txt").c_str(), readingsComplex->getCurrents()[patterno]);
		saveVals(("b_H" + std::to_string(patterno + 1) + ".txt").c_str(), currents);

		//vectorbView.setModel(makeMatrixTableModel(currents));
		vectorbView.show();

		CG_SolverComplex solver(*m1, currents, precond);
		int i = 0;
		for (; i < 1500 && solver.getResidueSquaredNorm() > 1e-19; i++)
			solver.do_iteration();
		std::cout << "Num its = " << solver.getIteration() << ". Rnorm = " << solver.getResidueSquaredNorm() << std::endl;

		//CG_SolverComplex solver2(Afull, Atfull, currents, precond2);
		//for (int i = 0; i < 900; i++)
		//	solver2.do_iteration(Atfull);

		x = solver.getX();
		//x2 = solver2.getX();
		saveVals(("x" + std::to_string(patterno + 1) + ".txt").c_str(), x);
		//saveVals(("xi" + std::to_string(patterno + 1) + ".txt").c_str(), x2);
		#ifndef BLOCKGND
		// Correct potentials
		std::complex<double> avg = 0;
		for (int i = input->getNodesCount() - input->getGenericElectrodesCount(); i < input->getNodesCount(); i++) avg += x[i];
		avg /= input->getGenericElectrodesCount();
		for (int i = 0; i < input->getNodesCount(); i++) x[i] -= avg;
		#endif

		//vectorView.setModel(makeMatrixTableModel(x.selfadjointView<Eigen::Lower>()));
		vectorView.show();

		solutions.push_back(x);
		//solutions2.push_back(x2);
		std::cout << "Finished solution " << patterno + 1 << " of " << readingsComplex->getCurrentsCount() << std::endl;
	}

	solutioncomplex::savePotentials(solutions, params.outputMesh.toStdString().c_str(), input, readingsComplex);
	//solutioncomplex::savePotentials(solutions2, (params.outputMesh.toStdString() + "2").c_str(), input, readingsComplex);

	return app.exec();
}