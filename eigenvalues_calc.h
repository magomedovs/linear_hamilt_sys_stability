#ifndef EIGENVALUES_CALC
#define EIGENVALUES_CALC

#include "isapprox.h"

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include <array>
#include <complex>
#include <iostream>

template <size_t DIM>
std::array< std::complex<double>, DIM > eigenvalues_calc(const std::array< std::array<double, DIM>, DIM >& monodromy_mat) {
	Eigen::Matrix< double, DIM, DIM > eig_mat;
	for (size_t i = 0; i < DIM; ++i) {
		for (size_t j = 0; j < DIM; ++j) {
			eig_mat(i, j) = monodromy_mat[i][j];
		}
	}
	// Eigenvalues calculation
	Eigen::Matrix< std::complex<double>, DIM, 1 > eigenvals( eig_mat.eigenvalues() );
//	std::cout << eigenvals << "\n";
//	std::cout << "abs values: " << std::abs(eigenvals(0)) << " " << std::abs(eigenvals(1)) << "\n";

	std::array< std::complex<double>, DIM > eigenvalsarr;
	for (size_t i = 0; i < DIM; ++i) {
		eigenvalsarr[i] = eigenvals(i);
	}
	return eigenvalsarr;
}

template <size_t DIM>
void stability_investigation(const std::array< std::array<double, DIM>, DIM >& monodromy_mat) {
	Eigen::Matrix< double, DIM, DIM > eig_mat;
	for (size_t i = 0; i < DIM; ++i) {
		for (size_t j = 0; j < DIM; ++j) {
			eig_mat(i, j) = monodromy_mat[i][j];
		}
	}	
	// Eigenvalues calculation
	Eigen::Matrix< std::complex<double>, DIM, 1 > eigenvals( eig_mat.eigenvalues() );
	// Eigenvectors calculation
	Eigen::EigenSolver< Eigen::Matrix< double, DIM, DIM > > es(eig_mat);
	Eigen::Matrix< std::complex<double>, DIM, DIM > eigenvectors( es.eigenvectors() );
//	std::cout << eigenvectors << "\n";
//	std::cout << eigenvectors.col(0) << "\n";

	for (size_t i = 0; i < DIM; ++i) {
		if ( !isapprox( std::abs(eigenvals(i)), 1.0 ) ) {
			// unstable case
			// ...
			std::cout << "unstable solution" << "\n";
			return ;
		}
	}

	std::cout << "probably stable solution" << "\n";
	// investigate eigenvectors and make conclusion about stability !

}

#endif
