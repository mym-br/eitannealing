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
#include "problemdescription.h"
#include "graphics.h"
#include "solver.h"
#include "nodecoefficients.h"
#include "solution.h"
#include "solution_lb.h"
#include "observations.h"
#include "random.h"
//#include "sparseincompletelq.h"
#include "gradientnormregularisation.h"
#include "intcoef.h"


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

	// Simulated annealing
	std::unique_ptr<solution_lb> current, next;
	float kt = 0.001;

	int totalit;
	int acceptit;
	shuffleData sdata;
	shuffler sh;
	current.reset(new solution_lb);

	std::cout.flush();

	int iterations;
	int solutions;
	double e;
	double r;
	double sqe;
	iterations = 0;
	int no_avance_count = 0;
	double prevE = 10000000000.0;
	//while(kt >0.00515377 && no_avance_count < 3) {
	while(kt >0.00000000515377 && no_avance_count < 3) {
		e = sqe = r = 0;
		totalit = acceptit = 0;
		solutions = 0;
		iterations = 0;
		while(totalit<15000 && acceptit < 3000) {

			next.reset(current->shuffle(&sdata, sh));
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
				view->setCurrentSolution(current->getSolution());
			}
		}

		double eav = e/solutions;
		double rav = r/solutions;
		double sige = sqrt(sqe/solutions - eav*eav);
		std::cout << kt << ":" << totalit << ":" << eav << ":" << sige << ":" << rav << ":" << ((float)iterations)/(nobs*solutions) << std::endl;

		kt *= 0.9;
		double variation = fabs(prevE-current->getDEstimate())/prevE;
		//std::cout << "totalit:" << iterations << std::endl;
		//std::cout << "variation: " << variation << std::endl;
		if((fabs(prevE-current->getDEstimate())/prevE) < 2.0e-15)
		  no_avance_count++;
		else
		  no_avance_count = 0;
		prevE = current->getDEstimate();
	}
	std::cout << "**********************\n";
	std::unique_ptr<solution_lb> current_lb, next_lb;
	current_lb.reset(new solution_lb(current->getSolution()));
	kt = 4.17456e-06;
	iterations = 0;
	no_avance_count = 0;
	prevE = 10000000000.0;
	while( no_avance_count < 3) {
		e = sqe = r = 0;
		totalit = acceptit = 0;
		solutions = 0;
		iterations = 0;
		while(totalit<15000 && acceptit < 3000) {

			next_lb.reset(current_lb->shuffle(&sdata, sh));
			bool decision;
			decision = current_lb->compareWith(*next_lb, kt, 1-param);
			if(decision) {
				iterations += current_lb->getTotalIt();
				solutions++;
				sh.addShufflerFeedback(sdata, true);
				current_lb = std::move(next_lb);
				acceptit++;
			} else {
				iterations += next_lb->getTotalIt();
				solutions++;
				sh.addShufflerFeedback(sdata, false);
			}
			e += current_lb->getDEstimate();
			r += current_lb->getRegularisationValue();
			sqe += current_lb->getDEstimate()*current_lb->getDEstimate();

			totalit++;
												if(totalit % 100 == 0) {
				//std::cout << current->getDEstimate() << ":" << current->getRegularisationValue() << std::endl;
				view->setCurrentSolution(current_lb->getSolution());
			}
		}


		double eav = e/solutions;
		double rav = r/solutions;
		double sige = sqrt(sqe/solutions - eav*eav);
		std::cout << kt << ":" << totalit << ":" << eav << ":" << sige << ":" << rav << ":" << ((float)iterations)/(nobs*solutions) << std::endl;

		kt *= 0.9;
		double variation = fabs(prevE-current_lb->getDEstimate())/prevE;
		//std::cout << "totalit:" << iterations << std::endl;
		//std::cout << "variation: " << variation << std::endl;
		if((fabs(prevE-current_lb->getDEstimate())/prevE) < 2.0e-15)
			no_avance_count++;
		else
			no_avance_count = 0;
		prevE = current_lb->getDEstimate();
	}
}

void setSol(float *sol);

/*#include <time.h>

double get_time()
{
    struct timespec t;
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &t);
    return ((double)t.tv_sec + (double)t.tv_nsec*1e-9)*1000;
}*/

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
	 prepareSkeletonMatrix();
         createCoef2KMatrix();
	 gradientNormRegularisation::initInstance();
	 intCoef::initInstance();

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

     double *sol = new double[numcoefficients];
     for(int i=0;i<numcoefficients;i++) sol[i]=1.0;
     view->setCurrentSolution(sol);


     std::thread worker(workProc);

     int retval =  app.exec();
     worker.join();
     return retval;
 }
