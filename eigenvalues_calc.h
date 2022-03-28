#ifndef EIGENVALUES_CALC
#define EIGENVALUES_CALC

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include <array>
#include <complex>
#include <iostream>

template <size_t DIM>
void eigenvalues_calc(const std::array< std::array<double, DIM>, DIM >& monodromy_mat) {
	Eigen::Matrix<double, DIM, DIM> eig_mat;
	for (size_t i = 0; i < DIM; ++i) {
		for (size_t j = 0; j < DIM; ++j) {
			eig_mat(i, j) = monodromy_mat[i][j];
		}
	}
	Eigen::Matrix<std::complex<double>, DIM, 1> eigenvals( eig_mat.eigenvalues() );
	std::cout << eigenvals << "\n";
	std::cout << "abs values: " << std::abs(eigenvals(0)) << " " << std::abs(eigenvals(1)) << "\n";

}

#endif
