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

	bool is2dProblem;
	std::string meshfname = params.inputMesh.toStdString();
	std::string currentsfname = params.inputCurrents.toStdString();
	std::string tensionsfname = params.inputTensions.toStdString();
	std::shared_ptr<problem> input = problem::createNewProblem(meshfname.c_str(), is2dProblem);
	input->setGroundNode(params.ground);
	input->initProblem(meshfname.c_str());
	input->initObsComplex(currentsfname.c_str(), tensionsfname.c_str());
	input->buildNodeCoefficients();
	input->prepareSkeletonMatrix();
	input->createCoef2KMatrix();

	matrixcomplex *m;
	//Eigen::VectorXd vcond(input->getNumCoefficients()), vperm(input->getNumCoefficients());
	//for (int i = 0; i < input->getNumCoefficients(); i++) {
	//	vcond[i] = 0.3815; vperm[i] = 0.00000000070922044418976;
	//}
	//input->assembleProblemMatrix(&vcond[0], &vperm[0], &m);

	//Eigen::VectorXcd vcond(input->getNumCoefficients());// , vperm(input->getNumCoefficients());
	//for (int i = 0; i < input->getNumCoefficients(); i++) {
	//	vcond[i] = std::complex<double>(0.3815, 0.00000000070922044418976);
	//}
	//input->assembleProblemMatrix(&vcond[0], &m);

	unsigned long long seed = 4939495;
	init_genrand64(seed);
	std::complex<double> *sol = solutioncomplex::getNewRandomSolution(input);
	matrixcomplex *m1 = solutioncomplex::getNewStiffness(sol, &m, input);
	//#ifndef BLOCKGND
	//for (int i = 0; i < input->getGroundNode(); i++) *(&m1->coeffRef(input->getGroundNode(), i)) = std::complex<double>(0, 0);
	//*(&m1->coeffRef(input->getGroundNode(), input->getGroundNode())) = std::complex<double>(1, 0);
	//#endif

    QTableView matrixView;
	//matrixView.setModel(makeMatrixTableModel(m->selfadjointView<Eigen::Lower>()));
    matrixView.setWindowTitle("Stiffness");
    matrixView.show();

	//matrixcomplex *m1 = new matrixcomplex(m->rows(), m->cols());
	//*m1 = (*m).conjugate().selfadjointView<Eigen::Lower>() * (matrixcomplex)(*m).selfadjointView<Eigen::Lower>();
	SparseIncompleteLLTComplex precond(*m1);

	Eigen::VectorXcd currents;
	QTableView vectorbView;
	vectorbView.setWindowTitle("Currents");

	Eigen::VectorXcd x;
	QTableView vectorView;
	vectorView.setWindowTitle("Tensions");

	saveVals("A.txt", *m, true);
	saveVals("A_H.txt", *m1, true);
	matrixcomplex L = precond.matrixL();
	saveVals("L.txt", L);
	//SparseIncompleteLLTComplex precondalt(*m); matrixcomplex Lalt = precondalt.matrixL(); saveVals("Lalt.txt", Lalt);


	std::vector<Eigen::VectorXcd> solutions;
	int a = input->getCurrentsComplexCount();
	for (int patterno = 0; patterno < input->getCurrentsComplexCount(); patterno++) {
		//currents = (*m).conjugate().selfadjointView<Eigen::Lower>() * input->getCurrents()[patterno];
		currents = *input->getCurrentsComplex(m, patterno);
		//#ifndef BLOCKGND
		//currents[input->getGroundNode()] = 0;
		//#endif

		//saveVals(("b" + std::to_string(patterno + 1) + ".txt").c_str(), input->getCurrents()[patterno]);
		saveVals(("b_H" + std::to_string(patterno + 1) + ".txt").c_str(), currents);

		//vectorbView.setModel(makeMatrixTableModel(currents));
		vectorbView.show();

		CG_SolverComplex solver(*m1, currents, precond);
		for (int i = 0; i < 100 && solver.getResidueSquaredNorm() != 0.0; i++)
			solver.do_iteration();
		
		x = solver.getX();
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
		std::cout << "Finished solution " << patterno + 1 << " of " << input->getCurrentsComplexCount() << std::endl;
	}

	solutioncomplex::savePotentials(solutions, params.outputMesh.toStdString().c_str(), input);

	return app.exec();
}