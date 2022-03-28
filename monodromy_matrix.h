#ifndef MONODROMY_MATRIX
#define MONODROMY_MATRIX

#include "integrate_system.h"

#include <array>

template <typename ODE_obj_T, size_t DIM>
std::array<std::array<double, DIM>, DIM> Monodromy_matrix(const ODE_obj_T& ode_system_obj, 
		Integrator_Type integrator_type=Integrator_Type::dopri5, const double dt=0.0001) {

	// creating (DIM * DIM) matrix and initializing it as identity matrix
	std::array<std::array<double, DIM>, DIM> monodr_mat;
	for (size_t i = 0; i < DIM; ++i) {
		for (size_t j = 0; j < DIM; ++j) {
			monodr_mat[i][j] = (i == j) ? 1  : 0;
		}
	}
	// calculating the monodromy matrix
	for (std::array<double, DIM>& sol_row : monodr_mat) {
		IntegrateSystem(ode_system_obj, sol_row, 0., ode_system_obj.T, dt, integrator_type);
	}

	// actually the obtained matrix is the transposed monodromy matrix
	return monodr_mat;
}

#endif

