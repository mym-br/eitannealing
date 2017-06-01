// complex constructor example
#include <iostream>     // std::cout
#include <complex>      // std::complex
#include <Eigen/Core>
#include "basematrix.h"

int main()
{
	//Eigen::Matrix3Xcd compA(3, 3);
	matrixcomplex compA(3, 3);
	std::vector<Eigen::Triplet<std::complex<double>>> coefficients;
	coefficients.push_back(Eigen::Triplet<std::complex<double>>(0, 0, std::complex<double>(1.0, 2.0)));
	coefficients.push_back(Eigen::Triplet<std::complex<double>>(0, 1, std::complex<double>(2.0, 3.0)));
	coefficients.push_back(Eigen::Triplet<std::complex<double>>(0, 2, std::complex<double>(3.0, 4.0)));
	coefficients.push_back(Eigen::Triplet<std::complex<double>>(1, 0, std::complex<double>(4.0, 5.0)));
	coefficients.push_back(Eigen::Triplet<std::complex<double>>(1, 1, std::complex<double>(5.0, 6.0)));
	coefficients.push_back(Eigen::Triplet<std::complex<double>>(1, 2, std::complex<double>(6.0, 7.0)));
	coefficients.push_back(Eigen::Triplet<std::complex<double>>(2, 0, std::complex<double>(7.0, 8.0)));
	coefficients.push_back(Eigen::Triplet<std::complex<double>>(2, 1, std::complex<double>(8.0, 9.0)));
	coefficients.push_back(Eigen::Triplet<std::complex<double>>(2, 2, std::complex<double>(9.0, 11.0)));
	compA.setFromTriplets(coefficients.begin(), coefficients.end());
	compA.makeCompressed();

	matrixcomplex compAconj(3, 3);
	compAconj = compA.conjugate();

	matrixcomplex compAherm(3, 3);
	compAherm = compAconj*compA;
	//compA(0, 0) = std::complex<double>(1.0,2.0);
	//compA(0, 1) = std::complex<double>(2.0, 3.0);
	//compA(0, 2) = std::complex<double>(3.0, 4.0);
	//compA(1, 0) = std::complex<double>(4.0, 5.0);
	//compA(1, 1) = std::complex<double>(5.0, 6.0);
	//compA(1, 2) = std::complex<double>(6.0, 7.0);
	//compA(2, 0) = std::complex<double>(7.0, 8.0);
	//compA(2, 1) = std::complex<double>(8.0, 9.0);
	//compA(2, 2) = std::complex<double>(9.0, 11.0);

	//Eigen::VectorXcd compb(3);
	//compb(0) = std::complex<double>(6.0, 9.0);
	//compb(1) = std::complex<double>(15.0, 18.0);
	//compb(2) = std::complex<double>(24.0, 28.0);

	//Eigen::VectorXcd compteste(3);
	//compteste(0) = std::complex<double>(1.0, 1.0);
	//compteste(1) = std::complex<double>(2.0, 2.0);
	//compteste(2) = std::complex<double>(3.0, 3.0);
	//
	std::cout << "A" << '\n';
	std::cout << compA << '\n';
	//std::cout << "b" << '\n';
	//std::cout << compb << '\n';
	//std::cout << "teste" << '\n';
	//std::cout << compteste << '\n';
	//std::cout << "b + 2*teste" << '\n';
	//std::cout << compb + 2*compteste << '\n';
	std::cout << "conj(A)*A" << '\n';
	std::cout << compAherm << '\n';

	//std::complex<double> first(2.0, 2.0);
	//std::complex<double> second(first);
	//std::complex<long double> third(second);

	//std::cout << third << '\n';

	//std::cout << first + second << '\n';

	//std::cout << first * second << '\n';

	return 0;
}