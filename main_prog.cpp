#include "eigenvalues_calc.h"
#include "linspace.h"
#include "jacobi_am_from_boost.h"
#include "ode_system_classes.h"
#include "integrate_system.h"
#include "monodromy_matrix.h"

#include <iostream>
#include <fstream>
#include <array>
#include <complex>

int main()
{
/*	test for monodromy matrix */
//	const Oscillation_system obj(0.5, 0.65, 0.);
	const Rotation_system obj(0.5, 0.85, 1.5);
	const size_t DIM(obj.DIM);
	typedef std::array<double, DIM> state_type;

//	std::array<state_type, DIM> monodr_mat = Monodromy_matrix<Oscillation_system, DIM>(obj);
	std::array<state_type, DIM> monodr_mat = Monodromy_matrix<Rotation_system, DIM>(obj);

//	std::array< std::complex<double>, DIM > eigvalsarr( eigenvalues_calc(monodr_mat) );
//	std::cout << eigvalsarr[0] << " " << eigvalsarr[1] << "\n";

	stability_investigation(monodr_mat);

/*	test for ode integration */
//	Oscillation_system Rotation_system
//	typedef std::array<double, Base_system::DIM> state_type;
//	const Rotation_system obj(0.5, 0.85, 1.5);
//	const Oscillation_system obj(0.5, 0.65, -0.5);
//	state_type x;
//	x[0] = 1.;
//	x[1] = 0.;
//	const double t1 = 0;
//	const double t2 = 1. * obj.T;
//	const double dt = 0.0001;
//
//	IntegrateSystem(obj, x, t1, t2, dt, Integrator_Type::dopri5, true);
//
//	std::cout << x[0] << " " << x[1] << "\n";

	return 0;
}
