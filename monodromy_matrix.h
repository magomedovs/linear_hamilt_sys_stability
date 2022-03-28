#ifndef MONODROMY_MATRIX
#define MONODROMY_MATRIX

#include "ode_system_classes.h"
#include "integrate_system.h"

#include <array>

template <typename ODE_obj_T>
std::array<state_type, DIM> Monodromy_matrix(const ODE_obj_T& ode_system_obj, 
		Integrator_Type integrator_type=Integrator_Type::dopri5, const double dt=0.0001) {

	// creating (DIM * DIM) matrix and initializing it as identity matrix
	std::array<state_type, DIM> monodr_mat;
	for (size_t i = 0; i < DIM; ++i) {
		for (size_t j = 0; j < DIM; ++j) {
			monodr_mat[i][j] = (i == j) ? 1  : 0;
		}
	}
	// calculating the monodromy matrix
	for (state_type& sol_row : monodr_mat) {
		IntegrateSystem(ode_system_obj, sol_row, 0., ode_system_obj.T, dt, integrator_type);
	}

	// actually the obtained matrix is the transposed monodromy matrix
	return monodr_mat;
}

#endif

