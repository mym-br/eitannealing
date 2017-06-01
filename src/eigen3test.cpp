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
//#include "incomplete_choleskycomplex.h"

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
	input->initObs(currentsfname.c_str(), tensionsfname.c_str());
	input->buildNodeCoefficients();
	input->prepareSkeletonMatrix();
	input->createCoef2KMatrix();
	
    matrix *m;
	Eigen::VectorXd vcond(input->getNumCoefficients());
	for (int i = 0; i < vcond.rows(); i++) vcond[i] = 0.3815;
	input->assembleProblemMatrix(&vcond[0], &m);

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

	matrix *m1 = new matrix(m->rows(), m->cols());
	//*m1 = (*m).selfadjointView<Eigen::Lower>() * (matrix)(*m).selfadjointView<Eigen::Lower>();
	*m1 = *m;
	#ifndef BLOCKGND
	for (int i = 0; i < input->getGroundNode(); i++) *(&m1->coeffRef(input->getGroundNode(), i)) = 0;
	*(&m1->coeffRef(input->getGroundNode(), input->getGroundNode())) = 1.0/32.0;
	#endif

	saveVals("A.txt", *m, true);
	saveVals("A_H.txt", *m1, true);

	std::vector<Eigen::VectorXd> solutions;
	for (int patterno = 0; patterno < input->getCurrentsCount(); patterno++) {
		//currents = (*m).selfadjointView<Eigen::Lower>() * input->getCurrents()[patterno];
		currents = input->getCurrents()[patterno];
		#ifndef BLOCKGND
		currents[input->getGroundNode()] = 0;
		#endif

		std::ofstream myfile;
		myfile.open("currents" + std::to_string(patterno+1) + ".txt");
		for (int i = 0; i < currents.size(); i++) myfile << currents.coeff(i) << std::endl;
		myfile.close();
		
		vectorbView.setModel(makeMatrixTableModel(currents));
		vectorbView.show();

		SparseIncompleteLLT precond(*m1);
		matrix L = precond.matrixL();
		saveVals("L.txt", L);
		CG_Solver solver(*m1, currents, precond);
		for (int i = 0; i < 100; i++) solver.do_iteration();
		
		x = solver.getX();
		#ifndef BLOCKGND
		// Correct potentials
		double avg = 0;
		for (int i = input->getNodesCount() - input->getGenericElectrodesCount(); i < input->getNodesCount(); i++) avg += x[i];
		avg /= input->getGenericElectrodesCount();
		for (int i = 0; i < input->getNodesCount(); i++) x[i] -= avg;
		#endif

		vectorView.setModel(makeMatrixTableModel(x.selfadjointView<Eigen::Lower>()));
		vectorView.show();

		myfile.open("tensions" + std::to_string(patterno+1) + ".txt");
		for (int i = 0; i < x.size(); i++) myfile << x.coeff(i) << std::endl;
		myfile.close();

		solutions.push_back(x);
		std::cout << "Finished solution " << patterno + 1 << " of " << input->getCurrentsCount() << std::endl;
	}

	solution::savePotentials(solutions, params.outputMesh.toStdString().c_str(), input);

	return app.exec();
}