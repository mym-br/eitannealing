/*
 * main.cpp
 *
 *  Created on: Jun 25, 2010
 *      Author: thiago
 */

#include <QApplication>
#include <QCoreApplication>
#include <QTableView>
#include <QListView>
#include <QThread>
#include <QAction>
#include <Eigen/Array>
#include <boost/thread.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <ctime>
#include "problemdescription.h"
#include "graphics.h"
#include "solver.h"
#include "nodecoefficients.h"
#include "solution.h"
#include "observations.h"
#include "random.h"

#include "init_obs_problem.h"


solutionView *view;

#include <Eigen/QR>

matrix *stiffness;
//Eigen::VectorXd *sqDiag;
matrix *stiffness0;

float *currentSolution;

float param;

bool e2test = false;

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
	std::auto_ptr<solution> current, next;
	float kt =  0.005;
	int totalit;
	int acceptit;
	shuffleData sdata;
	shuffler sh;
	current.reset(new solution());
	
	int iterations;
	int solutions;
	double e;
	double sqe;
	iterations = 0;
	int no_avance_count = 0;
	double prevE = 10000000000.0;
	while(kt > 0.0005 && no_avance_count < 3) {
		e = sqe = 0;
		totalit = acceptit = 0;
		solutions = 0;
		//iterations = 0;		
		while(totalit<12000 && acceptit < 2000) {
			next.reset(current->shuffle(&sdata, sh));
			bool decision;
			if(e2test)
			  decision = current->compareWithMaxE2(*next, kt, param);
			else
			  decision = current->compareWith(*next, kt, param);
			if(decision) {
			//if(current->compareWithMinIt(*next, kt, 13)) {
			//if(current->compareWithMaxE2(*next, kt, param)) {
				iterations += current->getTotalIt();
				solutions++;
				sh.addShufflerFeedback(sdata, true);
				current = next;
				acceptit++;
			} else {
				iterations += next->getTotalIt();
				solutions++;
				sh.addShufflerFeedback(sdata, false);
			}
			e += current->getDEstimate();
			sqe += current->getDEstimate()*current->getDEstimate();

			totalit++;
			if(totalit % 25 == 0) {
				view->setCurrentSolution(current->getSolution());				
			}
		}
		double eav = e/solutions;
		double sige = sqrt(sqe/solutions - eav*eav);
		//solution probe(current->getSolution());
		//probe.saturate();
		std::cout << kt << ":" << totalit << ":" << eav << ":" << sige << ":" << ((float)iterations)/(nobs*solutions) << std::endl;
		//std::cout << "last:" << current->getDEstimate() << " real:" << probe.getDEstimate() <<  std::endl;
		/*for(int it=0;it<numcoefficients;it++) {
		    std::cout << it << ":" << current->getSolution()[it] << std::endl;
		}*/
		
		
		kt *= 0.95;
		double variation = fabs(prevE-current->getDEstimate())/prevE;
		std::cout << "totalit:" << iterations << std::endl;
		//std::cout << "variation: " << variation << std::endl;
		if((fabs(prevE-current->getDEstimate())/prevE) < 2.0e-15)
		  no_avance_count++;
		else
		  no_avance_count = 0;		
		prevE = current->getDEstimate();
	}
	
	solution probe(current->getSolution());
	probe.saturate();
	std::cout << "real:" << probe.getDEstimate() << " iterations: " << iterations << std::endl;
		

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

}

void setSol(float *sol);

#include <time.h>

int main(int argc, char *argv[])
 {
  //  struct timespec time;
  //  clock_gettime(CLOCK_REALTIME, &time);
  // init_genrand64(time.tv_nsec);
   if(argc > 4)
     param = atof(argv[4]);
   else
     param = 0.875f;
   if(argc > 5)
     e2test = true;
   QApplication app(argc, argv);
     initProblem(argv[1]);

	 initObs(argv[2], argv[3]);
	 buildNodeCoefficients();
/*
	 float *solution = new float[numcoefficients];
	 for(int i=0;i<numcoefficients;i++) solution[i] = 1/4.76;


	 //solution[32] = 1.5;
	 //solution[node2coefficient[288]] = 1.5;
	 //solution[node2coefficient[358]] = 1.0;
	 //solution[node2coefficient[346]] = 1.0;	 
	 
     /*
	 
     assembleProblemMatrix(solution, &stiffness0);
	 
	 SparseIncompleteLLT precond(*stiffness0);

	 CG_Solver solver(*stiffness0, currents[29], precond);

	 for(int i=0;i<900;i++)
		 solver.do_iteration();

	 std::cout << "currents:\n";
	 std::cout << currents[29].end(31) << std::endl;
	 
	 std::cout << "tensions:\n";
	 std::cout << solver.getX().end(31) << std::endl;*/
     

     //int i;
     //float *coefficients = new float[65];
     //// Chessboard
     //for(i=0;i<65;i++) coefficients[i] = 1.0;
     //coefficients[28] = 1;
     //coefficients[29] = 1;
     //coefficients[36] = 1;
     //coefficients[37] = 1;
     //     assembleProblemMatrix(coefficients, &stiffness);
    // coefficients[1] = 2;
     //stiffness0 = assembleProblemMatrix(coefficients);

     //QTableView matrixView;
     
     //float sol[65];
     //for(int i=0;i<65;i++) sol[i] = 1.0;
     //matrix *stiffness = obs::buildObsProblemMatrix(sol);
     //matrixView.setModel(new matrixViewModel(*stiffness));
     //matrixView.setWindowTitle("Stiffness");
     //matrixView.show();
     
     //QTableView matrixView;
     //matrixView.setModel(new matrixViewModel(*stiffness0));
     //matrixView.setWindowTitle("Stiffness");
     //matrixView.show();

     float *sol = new float[numcoefficients];
     for(int i=0;i<numcoefficients;i++) sol[i] = 1.0;
     matrix *Kii, *Kcc;
     matrix2 *Kic;
     assembleProblemMatrix_lb(sol, &Kii, &Kic, &Kcc, 32);
     
	 SparseIncompleteLLT precond(*Kii);
	 Eigen::VectorXd J(currents[0].end(31));
     Eigen::VectorXd Phi(tensions[0].end(31));
     LB_Solver solver(Kii, Kic, Kcc, J, Phi, precond);
     
     for(int i=0;i<100;i++) {
       std::cout << solver.getIteration() << ":" <<  solver.getErrorl2Estimate() << std::endl;
       solver.do_iteration();       
     }
     
     QTableView matrixView1, matrixView2, matrixView3;
     matrixView1.setModel(new matrixViewModel(*Kii));
     matrixView1.setWindowTitle("Kii");
     matrixView1.show();
     matrixView2.setModel(new matrix2ViewModel(*Kic));
     matrixView2.setWindowTitle("Kic");
     matrixView2.show();
     matrixView3.setModel(new matrixViewModel(*Kcc));
     matrixView3.setWindowTitle("Kcc");
     matrixView3.show();
     
     
     qRegisterMetaType<QModelIndex>("QModelIndex");
     
     view = new solutionView(numcoefficients);
     QTableView list;
     list.setModel(view);
     QAction *copyDataAction = new QAction("Copy", &list);
     TableViewCopyDataPopupMenu::getInstance()->connect(copyDataAction, SIGNAL(triggered()), SLOT(actionFired()));
     list.addAction(copyDataAction);
     list.setContextMenuPolicy(Qt::ActionsContextMenu);
     list.show();
     viewport graphics(600, 600, argc>3?argv[4]:"Reverse Problem");
     graphics.show();
     graphics.connect(view, SIGNAL(dataChanged(QModelIndex,QModelIndex)), SLOT(solution_updated(QModelIndex,QModelIndex)));
     
          
     //boost::thread worker(workProc);

     int retval =  app.exec();
     //worker.join();
     return 0;
 }
 

