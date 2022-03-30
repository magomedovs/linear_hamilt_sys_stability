#include "eigenvalues_calc.h"
#include "linspace.h"
#include "jacobi_am_from_boost.h"
#include "ode_system_classes.h"
#include "integrate_system.h"
#include "monodromy_matrix.h"

#include <iostream>
#include <fstream>
#include <iomanip>
#include <array>
#include <complex>
#include <string>

template <typename ODE_obj_T, size_t DIM>
void CalculateAndWriteToFile(std::ofstream& stream, const ODE_obj_T& obj) {
	typedef std::array<double, DIM> state_type;
    std::array<state_type, DIM> monodr_mat = Monodromy_matrix<ODE_obj_T, DIM>(obj, Integrator_Type::fehlberg78);
	std::array< std::complex<double>, DIM > eigvalsarr( eigenvalues_calc(monodr_mat) );
	Solution_type stability_status = stability_investigation(eigvalsarr);
	stream << std::fixed << std::setprecision(8);
	stream << obj.param.alpha << " " 
		<< obj.param.beta << " "
		<< obj.param.h << " ";
	for (const auto& eig : eigvalsarr) {
		stream << eig << " ";
	}
	stream << Solution_type_to_string(stability_status) << "\n";
}


int main()
{
/*	test for monodromy matrix */
	const Oscillation_system osc_obj1(0.5, 0.65, 0.);
	const Oscillation_system osc_obj2(0.5, 0.65, -0.5);
	const Rotation_system rot_obj(0.5, 0.65, 1.5);

	std::ofstream output("output_stability_data.txt");
	CalculateAndWriteToFile<Oscillation_system, osc_obj1.DIM>(output, osc_obj1);
	CalculateAndWriteToFile<Oscillation_system, osc_obj2.DIM>(output, osc_obj2);
	CalculateAndWriteToFile<Rotation_system, rot_obj.DIM>(output, rot_obj);
	output.close();

	return 0;
}

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


