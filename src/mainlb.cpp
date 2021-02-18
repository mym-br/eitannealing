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
#include "solution_lb.h"
#include "solution_lb_complex.h"
#include "problem.h"
#include "twodim/problem2D.h"
#include "threedim/problem3D.h"
#include "observations.h"
#include "random.h"
//#include "sparseincompletelq.h"
#include "gradientnormregularisation.h"
#include "gradientnormregularisationcomplex.h"
#include "gmsh/gmshgraphics.h"
#include "parameters/parametersparser.h"
#define _USE_MATH_DEFINES
#include <math.h>
#include <QStyledItemDelegate>

solutionView *viewre, *viewim, *viewabs, *viewang;

#include <Eigen/QR>

matrix *stiffness;
//Eigen::VectorXd *sqDiag;a
matrix *stiffness0;

float *currentSolution;

float param;

std::shared_ptr<problem> input;
observations<std::complex<double>> *readingsComplex;
observations<double> *readingsScalar;
bool isComplexProblem;
unsigned long seed;
float kt;

void workProc()
{





	// Simulated annealing
	double *solre = new double[input->getNumCoefficients()];
	double *solim = new double[input->getNumCoefficients()];
	std::unique_ptr<solution_lb_complex> currentComplex, nextComplex;
	std::unique_ptr<solution_lb> currentScalar, nextScalar;

	int totalit;
	int acceptit;
	shuffleData sdata;
	std::unique_ptr<shuffler> sh;

	shuffleData sdata_i;
	std::unique_ptr<shuffler> sh_i;
	bool shuffleR = true;
	//sh = isComplexProblem ? std::unique_ptr<shuffler>(new shuffler(input, readingsComplex)) : new std::unique_ptr<shuffler>(shuffler(input, readingsScalar));
	if (isComplexProblem) {
		sh.reset(new shuffler(input, readingsScalar));
		sh_i.reset(new shuffler(input, readingsScalar));
		std::vector<std::complex<double>> electrodesCoeffs;
		currentComplex.reset(new solution_lb_complex(input, *readingsComplex));
	}
	else {
		std::vector<double> electrodesCoeffs;
		//for (int j = 0; j < 32; j++) electrodesCoeffs.push_back(0.002);
		sh.reset(new shuffler(input, readingsScalar));
		currentScalar.reset(new solution_lb(input, *readingsScalar));
	}

	std::cout.flush();

	int iterations;
	int solutions;
	double e;
	double r;
	double v;
	double sqe;
	iterations = 0;
	int no_avance_count = 0;
	double prevE = 10000000000.0;
	while (kt > 0.00000000005 && no_avance_count < 3) {
		e = sqe = r = 0; v = 0;
		totalit = acceptit = 0;
		solutions = 0;
		iterations = 0;
		while (totalit<15000 && acceptit < 3000) {

			isComplexProblem ? nextComplex.reset(currentComplex->shuffle(&sdata, *sh, &sdata_i, *sh_i, &shuffleR)) : nextScalar.reset(currentScalar->shuffle(&sdata, *sh));
			//next.reset(current->shuffle(&sdata, sh));
			bool decision;
			decision = isComplexProblem ? currentComplex->compareWith(*nextComplex, kt, 1 - param) : currentScalar->compareWith(*nextScalar, kt, 1 - param);
			int curits = isComplexProblem ? currentComplex->getTotalIt() : currentScalar->getTotalIt();
			if (decision) {
				iterations += isComplexProblem ? currentComplex->getTotalIt() : currentScalar->getTotalIt();
				solutions++;
				if(shuffleR) sh->addShufflerFeedback(sdata, true);
				else sh_i->addShufflerFeedback(sdata_i, true);
				if (isComplexProblem) currentComplex = std::move(nextComplex); else currentScalar = std::move(nextScalar);
				acceptit++;
			}
			else {
				iterations += isComplexProblem ? nextComplex->getTotalIt() : nextScalar->getTotalIt();
				solutions++;
				if(shuffleR) sh->addShufflerFeedback(sdata, false);
				else sh_i->addShufflerFeedback(sdata_i, false);
			}
			if (isComplexProblem) {
				e += currentComplex->getDEstimate();
				r += currentComplex->getRegularisationValue();
				sqe += currentComplex->getDEstimate()*currentComplex->getDEstimate();
			}
			else {
				e += currentScalar->getDEstimate();
				r += currentScalar->getRegularisationValue();
				sqe += currentScalar->getDEstimate()*currentScalar->getDEstimate();
			}

			totalit++;
			if (totalit % 100 == 0) {
				//std::cout << current->getDEstimate() << ":" << current->getRegularisationValue() << std::endl;
				//std::cout << totalit << ":" << acceptit << ":" << (isComplexProblem ? currentComplex->getDEstimate() : currentScalar->getDEstimate()) << ":" << (isComplexProblem ? currentComplex->getRegularisationValue() : currentScalar->getRegularisationValue()) << std::endl;
				double w = 2 * M_PI * input->getCurrentFreq();
				if (isComplexProblem) {
					for (int kk = 0; kk < input->getNumCoefficients(); kk++) {
						solre[kk] = currentComplex->getSolution_R()[kk];
						solim[kk] = currentComplex->getSolution_I()[kk] / w;
					}
					viewre->setCurrentSolution(solre);
					viewim->setCurrentSolution(solim);
				}
				else viewre->setCurrentSolution(currentScalar->getSolution());
			}
		}
		double eav = e / solutions;
		double rav = r / solutions;
		double vav = v / solutions;
		double sige = sqrt(sqe / solutions - eav*eav);
		//solution probe(current->getSolution());
		//probe.saturate();
		int nObs = isComplexProblem ? readingsComplex->getNObs() : readingsScalar->getNObs();
		std::cout << kt << ":" << totalit << ":" << eav << ":" << sige << ":" << rav << ":" << vav << ":" << ((float)iterations) / (nObs*solutions) << ":" << seed << std::endl;
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
		double variation = isComplexProblem ? fabs(prevE - currentComplex->getDEstimate()) / prevE : fabs(prevE - currentScalar->getDEstimate()) / prevE;
		//std::cout << "totalit:" << iterations << std::endl;
		//std::cout << "variation: " << variation << std::endl;
		//if ((fabs(prevE - current->getDEstimate()) / prevE) < 2.0e-15)
		if ((variation) < 2.0e-15)
			no_avance_count++;
		else
			no_avance_count = 0;
		prevE = isComplexProblem ? currentComplex->getDEstimate() : currentScalar->getDEstimate();

	}
	//probe.saturate();
	//std::cout << "real:" << probe.getDEstimate() << " iterations: " << iterations << std::endl;

#ifdef GGGGGGGGGGGG

	boost::mt11213b rng(std::clock());

	// Prepare a current vector, flowing from left to right
	Eigen::VectorXd current(numNodes - 1);
	int i;
	Eigen::VectorXd tension(numNodes - 1);
	for (i = 0; i<numNodes - 1; i++) tension[i] = i % 10 - 5;
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



	for (i = 0; i<numNodes; i++) {
		error = preTension - solver.getY();
		//perror = tension - psolver.getX();

		solver.do_iteration();
		//psolver.do_iteration();
		//if((solver.getIteration() % 30)==0)
		//solver.setrefresh();


		std::cout << solver.getIteration() << ":" << sqrt(error.squaredNorm()) << ":" << sqrt(solver.getErrorl2Estimate()) << std::endl;
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
#endif

double get_time()
{
#ifdef WIN32
	SYSTEMTIME time;
	GetSystemTime(&time);
	return (double)(time.wSecond * 1000) + time.wMilliseconds;
#else
	struct timespec t;
	clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &t);
	return ((double)t.tv_sec + (double)t.tv_nsec*1e-9) * 1000;
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
	init_genrand64(seed);
	param = params.peParam;
	bool is2dProblem;
	std::string meshfname = params.inputMesh.toStdString();
	std::string currentsinfname = params.inputCurrents.toStdString();
	std::string currentsoutfname = params.inputCurrentsOut.toStdString();
	std::string tensionsfname = params.inputTensions.toStdString();
	input = problem::createNewProblem(meshfname.c_str(), &is2dProblem);
	input->setGroundNode(params.ground);
	kt = params.kt;
	input->electrodevar = params.electrodevar;
	input->regularizationFactor = params.regularizationFactor;
	isComplexProblem = !currentsoutfname.empty();
	if (isComplexProblem) {
		// TODO: read parameters from commanline
		input->setCapacitance(80E-12);
		input->setCurrentFreq(275000);
	}
	input->initProblem(meshfname.c_str());
	if (params.calibrationMode) {
		input->setCalibrationMode(params.calibrationMode == 2);
	}
	const char **currentspair;

	if (isComplexProblem) {
		readingsComplex = new observations<std::complex<double>>;
		currentspair = new const char*[2]; currentspair[0] = currentsinfname.c_str(); currentspair[1] = currentsoutfname.c_str();
		readingsComplex->initObs(currentspair, tensionsfname.c_str(), input->getNodesCount(), input->getGenericElectrodesCount());
	}
	else {
		readingsScalar = new observations<double>;
		const char *currentsfnamecstr = currentsinfname.c_str();
		readingsScalar->initObs(&currentsfnamecstr, tensionsfname.c_str(), input->getNodesCount(), input->getGenericElectrodesCount());
	}

	input->buildNodeCoefficients();
	input->prepareSkeletonMatrix();
	input->createCoef2KMatrix();

	gradientNormRegularisation::initInstance(input);
	gradientNormRegularisationComplex::initInstance(input);
	gradientNormRegularisationComplex::initCalibrationInstance(input);

	qRegisterMetaType<QModelIndex>("QModelIndex");
	//qRegisterMetaType<QModelIndex>("QVector<int>");

	double minvalre, maxvalre, minvalim, maxvalim;
	minvalre = input->getCalibrationMode() != 0 ? mincondint : mincond; maxvalre = input->getCalibrationMode() != 0 ? maxcondint : maxcond;
	minvalim = input->getCalibrationMode() != 0 ? minpermint : minperm; maxvalim = input->getCalibrationMode() != 0 ? maxpermint : maxperm;

	viewport graphics(600, 600, "Reverse Problem Real", std::dynamic_pointer_cast<problem2D>(input), minvalre, maxvalre);
	viewport graphicsim(600, 600, "Reverse Problem Imaginary", std::dynamic_pointer_cast<problem2D>(input), minvalim, maxvalim);
	// Proccess mesh file name
	std::string outputMeshRe(params.outputMesh.toStdString()), outputMeshIm(params.outputMesh.toStdString());
	std::size_t dotfound = params.outputMesh.toStdString().find_last_of(".");
	outputMeshRe.replace(dotfound, 1, "_re."); outputMeshIm.replace(dotfound, 1, "_im.");
	// TODO: Proccess gmesh second address
	gmshviewport graphics_gmshre("eitannealingtest", outputMeshRe.c_str(), "Condutivity", params.gmeshAddress.toStdString().c_str(), input);
	gmshviewport graphics_gmshim("eitannealingtest", outputMeshIm.c_str(), "Permittivity", params.gmeshAddress.toStdString().c_str(), input);
	if (is2dProblem) {
		graphics.show();
		if (isComplexProblem) {
			graphicsim.show();
			graphicsim.move(graphics.pos() + QPoint(graphics.width() + 1, 0));
		}
	}

	viewre = new solutionView(input->getNumCoefficients());
	QTableView list;
	list.setModel(viewre);
	class ListViewDelegateRe : public QStyledItemDelegate {
	protected: QString displayText(const QVariant &value, const QLocale &locale) const { return locale.toString(value.toDouble(), 'f', 4); }
	} *redelegate = new ListViewDelegateRe;
	list.setItemDelegate(redelegate);
	list.setWindowTitle("Sol Real");
	QAction *copyDataAction = new QAction("Copy", &list);
	TableViewCopyDataPopupMenu::getInstance()->connect(copyDataAction, SIGNAL(triggered()), SLOT(actionFired()));
	list.addAction(copyDataAction);
	list.setContextMenuPolicy(Qt::ActionsContextMenu);
	list.show();
	list.resize(QSize(graphics.width(), graphics.height()));
	list.move(graphics.pos() + QPoint(0, graphics.height() + QApplication::style()->pixelMetric(QStyle::PM_TitleBarHeight)+8));

	QTableView listim;
	if (isComplexProblem) {
		viewim = new solutionView(input->getNumCoefficients());
		listim.setModel(viewim);
		class ListViewDelegateIm : public QStyledItemDelegate {
		protected: QString displayText(const QVariant &value, const QLocale &locale) const { return locale.toString(value.toDouble(), 'e', 4); }
		} *imdelegate = new ListViewDelegateIm;
		listim.setItemDelegate(imdelegate);
		listim.setWindowTitle("Sol Imag");
		QAction *copyDataActionim = new QAction("Copy", &listim);
		TableViewCopyDataPopupMenu::getInstance()->connect(copyDataActionim, SIGNAL(triggered()), SLOT(actionFired()));
		listim.addAction(copyDataActionim);
		listim.setContextMenuPolicy(Qt::ActionsContextMenu);
		listim.show();
		listim.resize(QSize(graphics.width(), graphics.height()));
		listim.move(graphicsim.pos() + QPoint(0, graphics.height() + QApplication::style()->pixelMetric(QStyle::PM_TitleBarHeight) + 8));
	}

	if (!params.gmeshAddress.isEmpty()) {
		graphics_gmshre.connect(viewre, SIGNAL(dataChanged(QModelIndex, QModelIndex)), SLOT(solution_updated(QModelIndex, QModelIndex)));
		if (isComplexProblem) graphics_gmshim.connect(viewim, SIGNAL(dataChanged(QModelIndex, QModelIndex)), SLOT(solution_updated(QModelIndex, QModelIndex)));
	}
	if (is2dProblem) {
		graphics.connect(viewre, SIGNAL(dataChanged(QModelIndex, QModelIndex)), SLOT(solution_updated(QModelIndex, QModelIndex)));
		if (isComplexProblem)  graphicsim.connect(viewim, SIGNAL(dataChanged(QModelIndex, QModelIndex)), SLOT(solution_updated(QModelIndex, QModelIndex)));
	}

	double *sol = new double[input->getNumCoefficients()];
	for (int i = 0; i<input->getNumCoefficients(); i++) sol[i] = 1.0;
	viewre->setCurrentSolution(sol);
	if (isComplexProblem) viewim->setCurrentSolution(sol);
	delete[] sol;
	std::thread worker(workProc);

	int retval = app.exec();
	worker.join();

	delete viewre;
	if (isComplexProblem) delete viewim;
	delete currentspair;
	return 0;
}
