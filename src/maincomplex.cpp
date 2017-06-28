/*
 * main.cpp
 *
 *  Created on: Jun 25, 2010
 *      Author: thiago
 * 
 *   FIXME: This file now is a mess. Its content needs refactoring!
 */

#include <QApplication>
#include <QCoreApplication>
#include <QTableView>
#include <QListView>
#include <QThread>
#include <QAction>
#include <thread>
#include <memory>
#include <ctime>
#include <iostream>
#include <QCommandLineParser>
//#include "problemdescription.h"
#include "graphics.h"
#include "solvercomplex.h"
//#include "nodecoefficients.h"
#include "solutioncomplex.h"
#include "problem.h"
#include "twodim/problem2D.h"
#include "threedim/problem3D.h"
//#include "solution_lb.h"
//#include "observations.h"
#include "random.h"
//#include "sparseincompletelq.h"
#include "gradientnormregularisationcomplex.h"
#include "gmsh\gmshgraphics.h"
#include "parameters\parametersparser.h"
#define _USE_MATH_DEFINES
#include <math.h>

solutionView *viewre, *viewim, *viewabs, *viewang;

#include <Eigen/QR>

matrix *stiffness;
//Eigen::VectorXd *sqDiag;
matrix *stiffness0;

float *currentSolution;

float param;

std::shared_ptr<problem> input;
observations<std::complex<double>> *readings;

unsigned long seed;

void workProc() 
{

	/*float *sol = new float[65];
	int i;
	sol[0] = 1.0;
	for(i=1;i<65;i++) {
		sol[i] = 2.0;
	}

	sol[28] = 1.0;
	sol[29] = 1.0;
	sol[36] = 1.0;
	sol[37] = 1.1;

	solution newsol(sol);

	for(i=0;i<numNodes*31;i++) {
		newsol.improve();
		std::cout << i << " - " << newsol.getDEstimate() << std::endl;
	}*/

/*
	float *sol = new float[65];
	int i;
	sol[0] = 1.0;
	for(i=1;i<65;i++) {
		sol[i] = 1.0;
	}

	obs::initObsProblem();
	matrix *big = obs::buildObsProblemMatrix(sol);

	SparseIncompleteLLT pbig(*big);
	matrix *small;
	assembleProblemMatrix(sol, &small);
	SparseIncompleteLLT psmall(*small);

	Eigen::VectorXd currentBig(obs::numNodes-1);
	Eigen::VectorXd currentSmall(numNodes-1);

	currentBig.fill(0);
	currentSmall.fill(0);

	/*
	for(i=0;i<7;i++) {
		currentBig[i] = -1;
		currentSmall[i] = -1;
	}

	for(i=0;i<8;i++) {
		currentBig[i+15] = 1;
		currentSmall[i+15] = 1;
	}*/




	// Simulated annealing
	double *solre = new double[input->getNumCoefficients()];
	double *solim = new double[input->getNumCoefficients()];
	std::unique_ptr<solutioncomplex> current, next;
	float kt = 0.05f;
	
	int totalit;
	int acceptit;
	shuffleData sdata;
	shufflercomplex sh(input);
	current.reset(new solutioncomplex(input, readings));
	
	std::cout.flush();
	
	int iterations;
	int solutions;
	double e;
	double r;
	double sqe;
	iterations = 0;
	int no_avance_count = 0;
	double prevE = 10000000000.0;
	while(kt > 0.00000000005 && no_avance_count < 3) {
		e = sqe = r = 0;
		totalit = acceptit = 0;
		solutions = 0;
		iterations = 0;		
		while(totalit<15000 && acceptit < 3000) {
                  
			next.reset(current->shuffle(&sdata, sh));
			//next.reset(new solutioncomplex(input));
			bool decision;
			decision = current->compareWith(*next, kt, 1-param);
			if(decision) {
				iterations += current->getTotalIt();
				solutions++;
				sh.addShufflerFeedback(sdata, true);
				current = std::move(next);
				acceptit++;
			} else {
				iterations += next->getTotalIt();
				solutions++;
				sh.addShufflerFeedback(sdata, false);
			}
			e += current->getDEstimate();
			r += current->getRegularisationValue();
			sqe += current->getDEstimate()*current->getDEstimate();

			totalit++;
            if(totalit % 100 == 0) {
				//std::cout << current->getDEstimate() << ":" << current->getRegularisationValue() << std::endl;
				for (int kk = 0; kk < input->getNumCoefficients(); kk++) {
					solre[kk] = current->getSolution()[kk].real();
					solim[kk] = current->getSolution()[kk].imag();
				}
				viewre->setCurrentSolution(solre);
				viewim->setCurrentSolution(solim);
			}
		}
		double eav = e/solutions;
		double rav = r/solutions;
		double sige = sqrt(sqe/solutions - eav*eav);
		//solution probe(current->getSolution());
		//probe.saturate();
		std::cout << kt << ":" << totalit << ":" << eav << ":" << sige << ":" << rav << ":" << ((float)iterations) / (readings->getNObs()*solutions) << ":" << seed << std::endl;
		//std::cout << "last:" << current->getDEstimate() << " real:" << probe.getDEstimate() <<  std::endl;
		/*for(int it=0;it<numcoefficients;it++) {
		    std::cout << it << ":" << current->getSolution()[it] << std::endl;
		}*/
		
		/*solution_lb probe(current->getSolution());
                probe.saturate();
                std::cout << "last (max):" << current->getDMax() << "last (min):" << current->getDMin() << " LB:" << probe.getDEstimate() <<  std::endl;
                *//*for(int kk=0;kk<9000;kk++) {
                  std::cout << (kk/32) << " Dest:" << current->getDEstimate() << std::endl;
                  current->improve();
                }*/
                  
                
                
                
		kt *= 0.9f;
		double variation = fabs(prevE-current->getDEstimate())/prevE;
		//std::cout << "totalit:" << iterations << std::endl;
		//std::cout << "variation: " << variation << std::endl;
		if((fabs(prevE-current->getDEstimate())/prevE) < 2.0e-15)
		  no_avance_count++;
		else
		  no_avance_count = 0;		
		prevE = current->getDEstimate();  
		
	}
	//probe.saturate();
	//std::cout << "real:" << probe.getDEstimate() << " iterations: " << iterations << std::endl;

#ifdef GGGGGGGGGGGG

	boost::mt11213b rng(std::clock());
	
	// Prepare a current vector, flowing from left to right
	Eigen::VectorXd current(numNodes-1);
	int i;
	Eigen::VectorXd tension(numNodes-1);
	for(i=0;i<numNodes-1;i++) tension[i] = i%10 - 5;
	current = *stiffness * tension;
	/*// Left electrodes (but the top left)
	for(i=0;i<7;i++) current[i] = -1;
	// Bottom
	for(;i<7+8;i++) current[i] = 0;
	// Right
	for(;i<7+8+8;i++) current[i] = 1;
	// Top and the rest
	for(;i<numNodes-1;i++) current[i] = 0;*/
	
	// Now solve for tensions


	Eigen::VectorXd error1, error2;
	Eigen::VectorXd error;
	//Eigen::VectorXd perror;
	SparseIncompleteLLT precond(*stiffness);
	CG_Solver solver(*stiffness, current, precond);
	//CG_PrecondSolver psolver(*stiffness, current);
/*
	matrix jacobi = solver.buildJacobiMatrx();
	Eigen::QR<Eigen::MatrixXd> jacobi_qr(jacobi.toDense());
	Eigen::VectorXd e1(jacobi.rows());
	e1.fill(0); e1[0] = 1.0;
	Eigen::VectorXd w(jacobi_qr.matrixR().transpose().solveTriangular(e1));
	std::cout << jacobi << std::endl;
	std::cout << "W:\n";
	std::cout << w << std::endl;
	std::cout << "norm:" << w.squaredNorm();
	std::cout << "Using jacobi:    " << current.squaredNorm()*w.squaredNorm() << std::endl;*/
	Eigen::VectorXd preTension(tension);
	solver.getPrecond().halfMultInPlace(preTension);
	std::cout << "Using stiffness: " << preTension.squaredNorm() << std::endl;
	//std::cout << "Using stiffness: " << sqrt(tension.dot(preTension)) << std::endl;
	std::cout << "Using stiffness (l2): " << tension.squaredNorm() << std::endl;



	for(i=0;i<numNodes;i++) {
		error = preTension - solver.getY();
		//perror = tension - psolver.getX();
			
		solver.do_iteration();
		//psolver.do_iteration();
		//if((solver.getIteration() % 30)==0)
			//solver.setrefresh();


		std::cout << solver.getIteration() << ":" << sqrt(error.squaredNorm()) << ":" <<  sqrt(solver.getErrorl2Estimate()) << std::endl;
	}
	//matrix jacobi = solver.buildJacobiMatrx();




/*
	CG_Solver solver0(*stiffness0, current, tension);

	for(i=0;i<numNodes+30;i++) {
		solver0.do_iteration();
	}

	CG_Solver solver1(*stiffness0, current, tension);
	CG_PrecondSolver psolver1(*stiffness0, current, tension);
	Eigen::VectorXd oldtension(tension);
	tension = solver0.getX();

	for(i=0;i<numNodes+30;i++) {
		error1 = tension - solver1.getX();
		error2 = tension - psolver1.getX();

		solver1.do_iteration();
		psolver1.do_iteration();

		if((psolver1.getIteration() % 30)==0)
			psolver1.setrefresh();

		std::cout << solver1.getIteration() << ":" << error1.lpNorm<2>() << ":" << solver1.getErrorl2Estimate() << std::endl;
		//std::cout << solver.getIteration() << ":" << norm.rows() << "," << norm.cols() << std::endl;
		//std::cout << solver.getIteration() << ":" << solver.getResidueSquaredNorm() << std::endl;
	}

	double totalerror = 0;
	for(i=0;i<31;i++) {
		totalerror += (solver1.getX()[i] - oldtension[i])*(solver1.getX()[i] - oldtension[i]);
	}
	totalerror = sqrt(totalerror);
	std::cout << "Total error:" << totalerror << std::endl;*/
#endif //#ifdef GGGGGGGGGGGG
	QApplication::quit();
}

void setSol(float *sol);

#ifdef WIN32
#include <windows.h>
#else
#include <time.h>
#endif

double get_time()
{
#ifdef WIN32
	SYSTEMTIME time;
	GetSystemTime(&time);
	return (double) (time.wSecond * 1000) + time.wMilliseconds;
#else
    struct timespec t;
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &t);
    return ((double)t.tv_sec + (double)t.tv_nsec*1e-9)*1000;
#endif
}

unsigned long getSeed() {
	// Windows specific
#ifdef _WIN32
	return GetTickCount();
#else // _WIN32
	// Linux specific
	struct timeval tv1;
	gettimeofday(&tv1, (struct timezone*)0);
	return (unsigned long)(tv1.tv_sec*1.E3 + tv1.tv_usec);
#endif
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

	if (!params.isSeedSpecified()) seed = getSeed();
	else seed = params.getSeed();
	seed = 4939495;
	init_genrand64(seed);
	param = params.peParam;
	bool is2dProblem;
	std::string meshfname = params.inputMesh.toStdString();
	std::string currentsfname = params.inputCurrents.toStdString();
	std::string tensionsfname = params.inputTensions.toStdString();
	//input = problem::createNewProblem(meshfname.c_str(), is2dProblem);
	is2dProblem = problem::isProblem2D(meshfname.c_str());
	if (is2dProblem) input = std::shared_ptr<problem>(new problem2D(meshfname.c_str()));
	else input = std::shared_ptr<problem>(new problem3D(meshfname.c_str()));
	input->setGroundNode(params.ground);
	// TODO: read parameter from commanline
	input->setCurrentFreq(275000);
	input->initProblem(meshfname.c_str());
	readings = new observations<std::complex<double>>;
	readings->initObs(currentsfname.c_str(), tensionsfname.c_str(), input->getNodesCount(), input->getGenericElectrodesCount());
	input->buildNodeCoefficients();
	input->prepareSkeletonMatrix();
	input->createCoef2KMatrix();

	gradientNormRegularisationComplex::initInstance(input);

	qRegisterMetaType<QModelIndex>("QModelIndex");
	qRegisterMetaType<QModelIndex>("QVector<int>");

	viewre = new solutionView(input->getNumCoefficients());
	QTableView list;
	list.setModel(viewre);
	list.setWindowTitle("Sol Real");
	QAction *copyDataAction = new QAction("Copy", &list);
	TableViewCopyDataPopupMenu::getInstance()->connect(copyDataAction, SIGNAL(triggered()), SLOT(actionFired()));
	list.addAction(copyDataAction);
	list.setContextMenuPolicy(Qt::ActionsContextMenu);
	list.show();

	viewim = new solutionView(input->getNumCoefficients());
	QTableView listim;
	listim.setModel(viewim);
	listim.setWindowTitle("Sol Imag");
	QAction *copyDataActionim = new QAction("Copy", &listim);
	TableViewCopyDataPopupMenu::getInstance()->connect(copyDataActionim, SIGNAL(triggered()), SLOT(actionFired()));
	listim.addAction(copyDataActionim);
	listim.setContextMenuPolicy(Qt::ActionsContextMenu);
	listim.show();

	double w = 2 * M_PI * input->getCurrentFreq();
	viewportcomplex graphics(600, 600, "Reverse Problem Real", std::dynamic_pointer_cast<problem2D>(input), mincond, maxcond);
	viewportcomplex graphicsim(600, 600, "Reverse Problem Imaginary", std::dynamic_pointer_cast<problem2D>(input), w*minperm, w*maxperm);
	//gmshviewport graphics_gmsh("eitannealingtest", params.outputMesh.toStdString().c_str(), params.gmeshAddress.toStdString().c_str(), input);
	if (!params.gmeshAddress.isEmpty()) {
		//graphics_gmsh.connect(view, SIGNAL(dataChanged(QModelIndex, QModelIndex)), SLOT(solution_updated(QModelIndex, QModelIndex)));
	}
	if (is2dProblem) {
		graphics.show();
		graphics.connect(viewre, SIGNAL(dataChanged(QModelIndex, QModelIndex)), SLOT(solution_updated(QModelIndex, QModelIndex)));
		graphicsim.show();
		graphicsim.connect(viewim, SIGNAL(dataChanged(QModelIndex, QModelIndex)), SLOT(solution_updated(QModelIndex, QModelIndex)));
	}

   double *sol = new double[input->getNumCoefficients()];
   for (int i = 0; i<input->getNumCoefficients(); i++) sol[i] = 1.0;
   viewre->setCurrentSolution(sol);
   viewim->setCurrentSolution(sol);
   delete[] sol;
   std::thread worker(workProc);
   
   int retval =  app.exec();
   worker.join();

   delete viewre;
   delete viewim;
   return 0;
 }
 

