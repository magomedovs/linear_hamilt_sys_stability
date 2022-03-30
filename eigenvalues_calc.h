#ifndef EIGENVALUES_CALC
#define EIGENVALUES_CALC

#include "isapprox.h"

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include <array>
#include <complex>
#include <iostream>
#include <iomanip>

template <size_t DIM>
std::array< std::complex<double>, DIM > eigenvalues_calc(const std::array< std::array<double, DIM>, DIM >& monodromy_mat) {
	Eigen::Matrix< double, DIM, DIM > eig_mat;
	for (size_t i = 0; i < DIM; ++i) {
		for (size_t j = 0; j < DIM; ++j) {
			eig_mat(i, j) = monodromy_mat[i][j];
		}
	}
	/* Eigenvalues calculation */
	Eigen::Matrix< std::complex<double>, DIM, 1 > eigenvals( eig_mat.eigenvalues() );
//	std::cout << eigenvals << "\n";
//	std::cout << "abs values: " << std::abs(eigenvals(0)) << " " << std::abs(eigenvals(1)) << "\n";

	std::array< std::complex<double>, DIM > eigenvalsarr;
	for (size_t i = 0; i < DIM; ++i) {
		eigenvalsarr[i] = eigenvals(i);
	}
	return eigenvalsarr;
}

enum class Solution_type {
	stable,
	unstable,
	double_minus_one,
	double_plus_one
};

template <size_t DIM>
void printEig(const std::array< std::complex<double>, DIM >& eigenvals) {
	for (const auto& eig : eigenvals) {
		std::cout << eig << " "; 
	}
}

template <size_t DIM>
Solution_type stability_investigation(const std::array< std::complex<double>, DIM >& eigenvals) {
//	std::cout << std::fixed << std::setprecision(8);
	for (size_t i = 0; i < DIM; ++i) {
		if ( !isapprox( std::abs(eigenvals[i]), 1.0 ) ) {
			/* unstable case */
//			printEig(eigenvals);
//			std::cout << "unstable" << "\n";
			return Solution_type::unstable;
		}
	}

	/* There are only complex conjugate roots lying on unit circle and double multiple roots +1 or -1. */
	/* 2D case */
	if ( !isapprox(eigenvals[0].imag(), 0.) ) {	// complex conjugate roots
//		printEig(eigenvals);
//		std::cout << "stable" << "\n";
		return Solution_type::stable;
	} else if ( isapprox(eigenvals[0].real(), -1.0) ) {	// double -1
//		printEig(eigenvals);
//		std::cout << "double_minus_one" << "\n";
		return Solution_type::double_minus_one;
	} else {	// double +1
//		printEig(eigenvals);
//		std::cout << "double_plus_one" << "\n";
		return Solution_type::double_plus_one;	
	}
	/* For multiple roots one needs to investigate eigenvectors to make conclusion about stability (do they form the basis of the vector space). */
}

#endif

	/* Eigenvectors calculation */
//	Eigen::EigenSolver< Eigen::Matrix< double, DIM, DIM > > es(eig_mat);
//	Eigen::Matrix< std::complex<double>, DIM, DIM > eigenvectors( es.eigenvectors() );

