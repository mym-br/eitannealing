#include <ctime>
#include <iostream>
#include "../src/lbmatrixbuilder.h"
#include "../src/solver_lb.h"
#include "../src/solver_lb_complex2.h"
#include "../src/intcoef.h"

typedef LBMatrixBuilder<eigen_double_engine::scalar, eigen_double_engine::symMatrix, eigen_double_engine::matrix> realMatrxiBuilder;
typedef LBMatrixBuilder<eigen_complexdouble_engine::scalar, eigen_complexdouble_engine::symMatrix, eigen_complexdouble_engine::matrix> complexMatrxiBuilder;

int main(int argc, char *argv[])
{
    bool is2d;
    if(argc < 2) {
        std::cerr << "need mesh filename, currents and tensions\n";
        return 1;
    }
    std::cout << "Parsing mesh file..." << std::flush;
    std::shared_ptr<problem> input(problem::createNewProblem(argv[1], &is2d));
    input->initProblem(argv[1]);
    input->buildNodeCoefficients();
    std::cout << "\n\t" << input->getNumCoefficients() << " coefficients detected\n";
    std::cout << "Done.\n" << std::flush;

    std::cout << "Preparing dummy solution (real)..." << std::flush;
    std::vector<double> coefficients;
    for(int i=0;i<input->getNumCoefficients();i++) coefficients.push_back(mincond  + i*(maxcond-mincond)/input->getNumCoefficients());
    std::cout << "Done.\n" << std::flush;

    std::cout << "Assembling matrices(old)... " << std::flush;
    std::unique_ptr<matrix>  Aii, Aic, Acc;
    {
        matrix *_Aii, *_Aic, *_Acc;
        assembleProblemMatrix_lb(coefficients.data(), &_Aii, &_Aic, &_Acc, *input);
        Aii.reset(_Aii);
        Aic.reset(_Aic);
        Acc.reset(_Acc);
    }
    std::cout << "Done.\n" << std::flush;

    std::cout << "Preparing new matrix builder..." << std::flush;
    realMatrxiBuilder builder(*input);
    std::cout << "Done\n" << std::flush;

    std::cout << "Assembling matrices(new)... " << std::flush;
    std::unique_ptr<eigen_double_engine::symMatrix> nAii(builder.buildAiiMatrix(coefficients));
    std::unique_ptr<eigen_double_engine::symMatrix> nAic(builder.buildAicMatrix(coefficients));
    std::unique_ptr<eigen_double_engine::matrix> nAcc(builder.buildAccMatrix(coefficients));
    std::cout << "Done.\n" << std::flush;

    std::cout << "Aii (new) - Aii (old) residual: " << ((*Aii) - (*nAii)).squaredNorm() << "\n";
    std::cout << "Aic (new) - Aic (old) residual: " << ((*Aic) - (*nAic)).squaredNorm() << "\n";
    std::cout << "Acc (new) - Acc (old) residual: " << ((*Acc) - (*nAcc)).squaredNorm() << "\n";

    std::cout << "Preparing dummy solution (complex)..." << std::flush;
    std::vector<double> icoefficients;
    std::vector<eigen_complexdouble_engine::scalar> ccoefficients;
    for(int i=0;i<input->getNumCoefficients();i++) {
        double im = maxcond  - i*(maxcond-mincond)/input->getNumCoefficients();
        icoefficients.push_back(im);
        ccoefficients.emplace_back(coefficients[i], im);
    }
    std::cout << "Done.\n" << std::flush;

    std::cout << "Assembling imaginary matrices(old)... " << std::flush;
    std::unique_ptr<matrix>  iAii, iAic, iAcc;
    {
        matrix *_iAii, *_iAic, *_iAcc;
        assembleProblemMatrix_lb(icoefficients.data(), &_iAii, &_iAic, &_iAcc, *input);
        iAii.reset(_iAii);
        iAic.reset(_iAic);
        iAcc.reset(_iAcc);
    }
    std::cout << "Done.\n" << std::flush;

    std::cout << "Preparing new complex matrix builder..." << std::flush;
    complexMatrxiBuilder cbuilder(*input);
    std::cout << "Done\n" << std::flush;

    std::cout << "Assembling matrices(new)... " << std::flush;
    std::unique_ptr<eigen_complexdouble_engine::symMatrix> ncAii(cbuilder.buildAiiMatrix(ccoefficients));
    std::unique_ptr<eigen_complexdouble_engine::symMatrix> ncAic(cbuilder.buildAicMatrix(ccoefficients));
    std::unique_ptr<eigen_complexdouble_engine::matrix> ncAcc(cbuilder.buildAccMatrix(ccoefficients));
    std::cout << "Done.\n" << std::flush;

    std::cout << "Re(Aii) (new) - Re(Aii) (old) residual: " << ((*Aii) - (ncAii->real())).squaredNorm() << "\n";
    std::cout << "Im(Aii) (new) - Im(Aii) (old) residual: " << ((*iAii) - (ncAii->imag())).squaredNorm() << "\n";
    std::cout << "Re(Aic) (new) - Re(Aic) (old) residual: " << ((*Aic) - (ncAic->real())).squaredNorm() << "\n";
    std::cout << "Im(Aic) (new) - Im(Aic) (old) residual: " << ((*iAic) - (ncAic->imag())).squaredNorm() << "\n";
    std::cout << "Re(Acc) (new) - Re(Acc) (old) residual: " << ((*Acc) - (ncAcc->real())).squaredNorm() << "\n";
    std::cout << "Im(Acc) (new) - Im(Acc) (old) residual: " << ((*iAcc) - (ncAcc->imag())).squaredNorm() << "\n";



    return 0;
}
