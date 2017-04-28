#include "solver.h"
#include "problem.h"
#include "solution.h"
#include <QTableView>
#include <QApplication>
#include <QCommandLineParser>
#include "parameters\parametersparser.h"
#include "matrixview.h"
#include <iostream>

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
	

    matrix *m1;
	Eigen::VectorXd v(input->getNumCoefficients());
	for (int i = 0; i < v.rows(); i++) v[i] = 1;
	input->assembleProblemMatrix(&v[0], &m1);  
    
    QTableView matrixView;
    matrixView.setModel(makeMatrixTableModel(m1->selfadjointView<Eigen::Lower>()));
    matrixView.setWindowTitle("Stiffness");
    matrixView.show();

	std::ofstream myfile;
	myfile.open("stiffness.txt");
	for (int i = 0; i < m1->rows(); i++) {
		for (int j = 0; j < m1->cols(); j++)
			if (j < i) myfile << m1->coeff(i, j) << " ";
			else myfile << m1->coeff(j, i) << " ";
		myfile << std::endl;
	}
	myfile.close();

	Eigen::VectorXd currents;
	QTableView vectorbView;
	vectorbView.setWindowTitle("Currents");

	Eigen::VectorXd x;
	QTableView vectorView;
	vectorView.setWindowTitle("Tensions");

	std::vector<Eigen::VectorXd> solutions;
	for (int patterno = 0; patterno < input->getCurrentsCount(); patterno++) {
		currents = input->getCurrents()[patterno];

		myfile;
		myfile.open("currents" + std::to_string(patterno+1) + ".txt");
		for (int i = 0; i < currents.size(); i++) myfile << currents.coeff(i) << std::endl;
		myfile.close();
		
		vectorbView.setModel(makeMatrixTableModel(currents));	
		vectorbView.show();

		SparseIncompleteLLT precond(*m1);
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

		myfile;
		myfile.open("tensions" + std::to_string(patterno+1) + ".txt");
		for (int i = 0; i < x.size(); i++) myfile << x.coeff(i) << std::endl;
		myfile.close();

		solutions.push_back(x);
		std::cout << "Finished solution " << patterno + 1 << " of " << input->getCurrentsCount() << std::endl;
	}

	solution::savePotentials(solutions, params.outputMesh.toStdString().c_str(), input);

	return app.exec();
}