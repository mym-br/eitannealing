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
	float kt = 400;
	int totalit;
	int acceptit;
	shuffleData sdata;
	shuffler sh;
	current.reset(new solution());
	int iterations;
	int solutions;
	double e;
	double sqe;
	while(kt > 0.000001) {
		e = sqe = 0;
		totalit = acceptit = 0;
		iterations = solutions = 0;
		while(totalit<15000 && acceptit < 3000) {
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
			if(totalit % 25 == 0) {
				view->setCurrentSolution(current->getSolution());
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

void setSol(float *sol);

 int main(int argc, char *argv[])
 {
     QApplication app(argc, argv);
     initProblem(argv[1]);
	 initObs(argv[2], argv[3]);
	 buildNodeCoefficients();

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


     view =new viewport(600, 600, argc>3?argv[4]:"Reverse Problem");
	 view->setCurrentSolution(solution);
     view->show();
          
     boost::thread worker(workProc);

     return app.exec();
     worker.join();
 }
 

 
void setSol(float *sol)
{
	int i;
	sol[++i] = 0.118859;
	sol[++i] = 0.091111399;
	sol[++i] = 0.100251;
	sol[++i] = 0.106117;
	sol[++i] = 0.160247;
	sol[++i] = 0.105798;
	sol[++i] = 0.106278;
	sol[++i] = 0.104063;
	sol[++i] = 0.088288397;
	sol[++i] = 0.0978375;
	sol[++i] = 0.0876754;
	sol[++i] = 0.079351403;
	sol[++i] = 0.078733303;
	sol[++i] = 0.084508903;
	sol[++i] = 0.075319901;
	sol[++i] = 0.078695297;
	sol[++i] = 0.079169698;
	sol[++i] = 0.083236098;
	sol[++i] = 0.0828951;
	sol[++i] = 0.092743203;
	sol[++i] = 0.101078;
	sol[++i] = 0.100138;
	sol[++i] = 0.102309;
	sol[++i] = 0.112525;
	sol[++i] = 0.125682;
	sol[++i] = 0.087949798;
	sol[++i] = 0.121403;
	sol[++i] = 0.11835;
	sol[++i] = 0.120835;
	sol[++i] = 0.095488697;
	sol[++i] = 0.132771;
	sol[++i] = 0.118369;
	sol[++i] = 0.35421553;
	sol[++i] = 0.054601301;
	sol[++i] = 0.078092746;
	sol[++i] = 0.39992017;
	sol[++i] = 0.3999936;
	sol[++i] = 0.39985937;
	sol[++i] = 0.3993668;
	sol[++i] = 0.39980483;
	sol[++i] = 0.39742151;
	sol[++i] = 0.39999518;
	sol[++i] = 0.39999586;
	sol[++i] = 0.39995813;
	sol[++i] = 0.39999288;
	sol[++i] = 0.39998108;
	sol[++i] = 0.0010174469;
	sol[++i] = 0.39970458;
	sol[++i] = 0.39991745;
	sol[++i] = 0.20288259;
	sol[++i] = 0.39999634;
	sol[++i] = 0.39996758;
	sol[++i] = 0.39999619;
	sol[++i] = 0.39863312;
	sol[++i] = 0.39999369;
	sol[++i] = 0.39997411;
	sol[++i] = 0.39999986;
	sol[++i] = 0.39994782;
	sol[++i] = 0.39992842;
	sol[++i] = 0.39999962;
	sol[++i] = 0.39998999;
	sol[++i] = 0.39997974;
	sol[++i] = 0.39995131;
	sol[++i] = 0.3999697;
	sol[++i] = 0.39997223;
	sol[++i] = 0.39982533;
	sol[++i] = 0.39996183;
	sol[++i] = 0.39997625;
	sol[++i] = 0.39998993;
	sol[++i] = 0.39990193;
	sol[++i] = 0.31986132;
	sol[++i] = 0.34098965;
	sol[++i] = 0.399921;
	sol[++i] = 0.39998892;
	sol[++i] = 0.39998636;
	sol[++i] = 0.39988568;
	sol[++i] = 0.39999491;
	sol[++i] = 0.0011092508;
	sol[++i] = 0.39999205;
	sol[++i] = 0.39999107;
	sol[++i] = 0.39999676;
	sol[++i] = 0.39999563;
	sol[++i] = 0.39999962;
	sol[++i] = 0.39999881;
	sol[++i] = 0.3999902;
	sol[++i] = 0.39978328;
	sol[++i] = 0.39997184;
	sol[++i] = 0.39998195;
	sol[++i] = 0.39965218;
	sol[++i] = 0.0010428318;
	sol[++i] = 0.19049102;
	sol[++i] = 0.0011390392;
	sol[++i] = 0.39996099;
	sol[++i] = 0.39996141;
	sol[++i] = 0.39999294;
	sol[++i] = 0.39953852;
	sol[++i] = 0.39999306;
	sol[++i] = 0.39996853;
	sol[++i] = 0.39999583;
	sol[++i] = 0.39994928;
	sol[++i] = 0.0010730034;
	sol[++i] = 0.39996752;
	sol[++i] = 0.3999494;
	sol[++i] = 0.39854196;
	sol[++i] = 0.39999577;
	sol[++i] = 0.39997965;
	sol[++i] = 0.39979729;
	sol[++i] = 0.39998087;
	sol[++i] = 0.39980671;
	sol[++i] = 0.0020898296;
	sol[++i] = 0.39998797;
	sol[++i] = 0.39981347;
	sol[++i] = 0.39996654;
	sol[++i] = 0.39997408;
	sol[++i] = 0.39999658;
	sol[++i] = 0.39997998;
	sol[++i] = 0.39986727;
	sol[++i] = 0.39996868;
	sol[++i] = 0.39999387;
	sol[++i] = 0.39990991;
	sol[++i] = 0.39999226;
	sol[++i] = 0.39988443;
	sol[++i] = 0.39997697;
	sol[++i] = 0.39855841;
	sol[++i] = 0.39983752;
	sol[++i] = 0.39998797;
	sol[++i] = 0.39999831;
	sol[++i] = 0.39992708;
	sol[++i] = 0.39998728;
	sol[++i] = 0.39996022;
	sol[++i] = 0.39998612;
	sol[++i] = 0.39998898;
	sol[++i] = 0.39999723;
	sol[++i] = 0.39998898;
	sol[++i] = 0.39999962;
	sol[++i] = 0.39969119;
	sol[++i] = 0.39999124;
	sol[++i] = 0.39998105;
	sol[++i] = 0.39964896;
	sol[++i] = 0.39999911;
	sol[++i] = 0.39999616;
	sol[++i] = 0.39998847;
	sol[++i] = 0.39999321;
	sol[++i] = 0.39999932;
	sol[++i] = 0.39997113;
	sol[++i] = 0.39973491;
	sol[++i] = 0.39999938;
	sol[++i] = 0.39997441;
	sol[++i] = 0.39844882;
	sol[++i] = 0.39994127;
	sol[++i] = 0.39984819;
	sol[++i] = 0.39998332;
	sol[++i] = 0.39996198;
	sol[++i] = 0.39999422;
	sol[++i] = 0.39995658;
	sol[++i] = 0.39999804;
	sol[++i] = 0.39999989;
	sol[++i] = 0.39979181;
	sol[++i] = 0.39983839;
	sol[++i] = 0.39990836;
	sol[++i] = 0.39998832;
	sol[++i] = 0.39997208;
	sol[++i] = 0.39979103;
	sol[++i] = 0.39992902;
	sol[++i] = 0.0010252512;
	sol[++i] = 0.39994299;
	sol[++i] = 0.39999968;
	sol[++i] = 0.3999728;
	sol[++i] = 0.39999697;
}