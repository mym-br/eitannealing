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
#include "solution.h"
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
	std::unique_ptr<solutionbase<std::complex<double>>> currentComplex, nextComplex;
	std::unique_ptr<solutionbase<double>> currentScalar, nextScalar;

	int totalit;
	int acceptit;
	shuffleData sdata;
	std::unique_ptr<shuffler> sh;
	//sh = isComplexProblem ? std::unique_ptr<shuffler>(new shuffler(input, readingsComplex)) : new std::unique_ptr<shuffler>(shuffler(input, readingsScalar));
	if (isComplexProblem) {
		sh.reset(new shuffler(input, readingsComplex));
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
		if (input->getCalibrationMode()) currentComplex.reset(new solutioncomplexcalibration(input, readingsComplex, electrodesCoeffs));
		else currentComplex.reset(new solutioncomplex(input, readingsComplex, electrodesCoeffs));
	}
	else {
		std::vector<double> electrodesCoeffs;
		//for (int j = 0; j < 32; j++) electrodesCoeffs.push_back(0.002);
		sh.reset(new shuffler(input, readingsScalar));
		currentScalar.reset(new solution(input, readingsScalar, electrodesCoeffs));
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

			isComplexProblem ? nextComplex.reset(currentComplex->shuffle(&sdata, *sh)) : nextScalar.reset(currentScalar->shuffle(&sdata, *sh));
			//next.reset(current->shuffle(&sdata, sh));
			bool decision;
			decision = isComplexProblem ? currentComplex->compareWith(*nextComplex, kt, 1 - param) : currentScalar->compareWith(*nextScalar, kt, 1 - param);
			int curits = isComplexProblem ? currentComplex->getTotalIt() : currentScalar->getTotalIt();
			if (decision) {
				iterations += isComplexProblem ? currentComplex->getTotalIt() : currentScalar->getTotalIt();
				solutions++;
				sh->addShufflerFeedback(sdata, true);
				if (isComplexProblem) currentComplex = std::move(nextComplex); else currentScalar = std::move(nextScalar);
				acceptit++;
			}
			else {
				iterations += isComplexProblem ? nextComplex->getTotalIt() : nextScalar->getTotalIt();
				solutions++;
				sh->addShufflerFeedback(sdata, false);
			}
			if (isComplexProblem) {
				e += currentComplex->getDEstimate();
				r += currentComplex->getRegularisationValue();
				sqe += currentComplex->getDEstimate()*currentComplex->getDEstimate();
			}
			else {
				e += currentScalar->getDEstimate();
				r += currentScalar->getRegularisationValue() - currentScalar->getElectrodeVariance();
				v += currentScalar->getElectrodeVariance();
				sqe += currentScalar->getDEstimate()*currentScalar->getDEstimate();
			}

			totalit++;
			if (totalit % 100 == 0) {
				//std::cout << current->getDEstimate() << ":" << current->getRegularisationValue() << std::endl;
				//std::cout << totalit << ":" << acceptit << ":" << (isComplexProblem ? currentComplex->getDEstimate() : currentScalar->getDEstimate()) << ":" << (isComplexProblem ? currentComplex->getRegularisationValue() : currentScalar->getRegularisationValue()) << std::endl;
				double w = 2 * M_PI * input->getCurrentFreq();
				if (isComplexProblem) {
					for (int kk = 0; kk < input->getNumCoefficients(); kk++) {
						solre[kk] = currentComplex->getSolution()[kk].real();
						solim[kk] = currentComplex->getSolution()[kk].imag() / w;
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

		int nObs = isComplexProblem ? readingsComplex->getNObs() : readingsScalar->getNObs();
		std::cout << kt << ":" << totalit << ":" << eav << ":" << sige << ":" << rav << ":" << vav << ":" << ((float)iterations) / (nObs*solutions) << ":" << seed << std::endl;

		kt *= 0.9f;
		double variation = isComplexProblem ? fabs(prevE - currentComplex->getDEstimate()) / prevE : fabs(prevE - currentScalar->getDEstimate()) / prevE;
		if ((variation) < 2.0e-15)
			no_avance_count++;
		else
			no_avance_count = 0;
		prevE = isComplexProblem ? currentComplex->getDEstimate() : currentScalar->getDEstimate();
	}
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
		readingsComplex->initObs(currentspair, tensionsfname.c_str(), input->getNodesCount(), input->getGenericElectrodesCount(), input->getGroundNode());
	}
	else {
		readingsScalar = new observations<double>;
		const char *currentsfnamecstr = currentsinfname.c_str();
		readingsScalar->initObs(&currentsfnamecstr, tensionsfname.c_str(), input->getNodesCount(), input->getGenericElectrodesCount(), input->getGroundNode());
	}

	input->buildNodeCoefficients();
	input->prepareSkeletonMatrix();
	input->createCoef2KMatrix();

	gradientNormRegularisation::initInstance(input);
	gradientNormRegularisationComplex::initInstance(input);
	gradientNormRegularisationComplex::initCalibrationInstance(input);
	qRegisterMetaType<QModelIndex>("QModelIndex");

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
	delete[] sol;
	std::thread worker(workProc);

	int retval = app.exec();
	worker.join();

	delete viewre;
	if (isComplexProblem) delete viewim;
	delete currentspair;
	return 0;
}
