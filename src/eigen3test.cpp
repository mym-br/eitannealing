#include <iostream>
#include <vector>
#include "basematrix.h"
#include "solver.h"
#include "incomplete_cholesky.h"
#include "incomplete_LDLT.h"
#include "assembleMatrix.h"
#include "problemdescription.h"
#include "gradientnormregularisation.h"
#include <sys/time.h>
#include <QTableView>
#include <QApplication>
#include "matrixview.h"
typedef Eigen::Triplet<double> T;


inline double gettime ()
{
    timeval tv;
    gettimeofday (&tv, NULL);
    return double (tv.tv_sec) + 0.000001 * tv.tv_usec;
}

int main(int argc, char *argv[])
{
    initProblem(argv[1]);
    buildNodeCoefficients();
    prepareSkeletonMatrix();
    createCoef2KMatrix();
    matrix *m1, *m2;
    double startTime, endTime;
    Eigen::VectorXd v(numcoefficients);
    for(int i=0; i<v.rows(); i++)
      v[i] = 1.0+0.0001*(i%11);
    assembleProblemMatrix(&v[0], &m1);
    
    QApplication app(argc, argv);
    
    QTableView matrixView;
    matrixView.setModel(makeMatrixTableModel(m1->selfadjointView<Eigen::Lower>()));
    matrixView.setWindowTitle("Stiffness");
    matrixView.show();
     
        
    return app.exec();
}