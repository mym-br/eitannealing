/*
 * main.cpp
 *
 *  Created on: Jun 25, 2010
 *      Author: thiago
 */

#include <QApplication>
#include <QTableView>
#include <QThread>
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
#include <sstream>

#include "init_obs_problem.h"


viewport *view;

#include <Eigen/QR>

matrix *stiffness;
//Eigen::VectorXd *sqDiag;
matrix *stiffness0;

float *currentSolution;

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
	float kt = 9;
	int totalit;
	int acceptit;
	shuffleData sdata;
	shuffler sh;
	
	current.reset(new solution());
	std::cout.flush();
	
	int iterations;
	int solutions;
	double e;
	double sqe;
	int i;
	int filecount = 0;
	std::stringstream filename;
	float solAv[65];
	float solVa[65];
	
	while(kt > 0.0000001) {
		for(i=0;i<65;i++) solAv[i] =0;
		for(i=0;i<65;i++) solVa[i] =0;
		e = sqe = 0;
		totalit = acceptit = 0;
		iterations = solutions = 0;
		while(totalit<12288 && acceptit < 4096) {
			next.reset(current->shuffle(&sdata, sh));
			if(current->compareWith(*next, kt, 0.01)) {
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
			
			for(i=0;i<65;i++) {
			  float v = current->getSolution()[i];
			  solAv[i] += v;
			  solVa[i] += v*v;
			}			
			if(totalit % 1024 == 0) {
			  
			      for(i=0;i<65;i++) {
				solAv[i] /= 1024;
				solVa[i] /= 1024;
				solVa[i] -= solAv[i]*solAv[i];
				if(solVa[i]>0) {
				  solVa[i] = std::sqrt(solVa[i])/0.28867513;			
				} else solVa[i]=0;
			      }
			      view->setCurrentSolution(solAv, solVa);
      			      for(i=0;i<65;i++) {
				solAv[i] = 0;
				solVa[i] = 0;
			      }
			      
			      filename.str("");
			      filename << "img/sol" << std::setw(5) << std::setfill('0') << filecount++ << ".png";
			      view->saveImage(filename.str());
			}
		}
		double eav = e/solutions;
		double sige = sqrt(sqe/solutions - eav*eav);
		std::cout << kt << ":" << totalit << ":" << eav << ":" << sige << ":" << ((float)iterations)/(nobs*solutions) << std::endl;
		kt *= 0.95;
	}

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


 int main(int argc, char *argv[])
 {
     QApplication app(argc, argv);
     initProblem();
     buildNodeCoefficients();
     initObs();

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
     //obs::initObsProblem();
     //float sol[65];
     //for(int i=0;i<65;i++) sol[i] = 1.0;
     //matrix *stiffness = obs::buildObsProblemMatrix(sol);
     //matrixView.setModel(new matrixViewModel(*stiffness));
     //matrixView.setWindowTitle("Stiffness");
     //matrixView.show();

     view =new viewport(400, 400, "Reverse Problem");
     view->show();
     drawElements(*view);
     
	 boost::thread worker(workProc);

     return app.exec();
    // worker.join();
 }
 
 