#include "linspace.h"
#include "jacobi_am_from_boost.h"
#include "ode_system_classes.h"
#include "integrate_system.h"

#include <iostream>
#include <fstream>
#include <array>

/*
class Harmonic_oscillator {
    const double coef;
public:
    Harmonic_oscillator(double new_coef) : coef(new_coef) {}

    void operator() (const state_type& x, state_type& dxdt, const double t) {
        dxdt[0] = x[1];
        dxdt[1] = -x[0]; // -coef * std::sin(x[0]) + 0.1 * std::sin(t); // -x[0] - coef * x[1];
    }

};
*/

int main()
{
//	Oscillation_system Rotation_system
//	const Rotation_system< state_type > obj(0.5, 0.65, 1.5);
	const Oscillation_system< state_type > obj(0.5, 0.65, -0.5);
	state_type x;
	x[0] = 0.02;
	x[1] = 0.03;
	const double t1 = 0;
	const double t2 = 10 * obj.T;
	const double dt = 0.0001;

	IntegrateSystem(obj, x, t1, t2, dt, Integrator_Type::dopri5, true);

	std::cout << x[0] << " " << x[1] << "\n";

/*	test for ode class */
//	const size_t steps_number = 10000;
//	const double interval_length = 2;
//	const double step = interval_length / steps_number;
//	std::ofstream output("output.txt");
//
//	std::array<double, steps_number> h_arr;
//
//	Oscillation_system< std::array<double, 2> > osc_obj(3, 4, 0.3);
//
//	const double tol = std::pow(10, -4);
//	linspace(-1+tol, 1-tol, h_arr);
//	for (double h : h_arr) {
//		output << h << " " << osc_obj.k_func(h) << " " << osc_obj.omega_func(osc_obj.k_func(h)) << "\n";
//	}

//	const double k = 0.9;
//	Rotation_system rot_obj(k);
//	std::array<double, steps_number> u_arr;
//	linspace(0, 5 * boost::math::ellint_1(k), u_arr);
//	for (double u : u_arr) {
//		output << u << " " << rot_obj.q1_func(u) << "\n";
//	}

//	output.close();
//
//	std::cout << osc_obj.T << "\n";
//	std::cout << osc_obj.param.alpha << " " << osc_obj.param.beta << " " << osc_obj.param.h << "\n";
//
//	Rotation_system< std::array<double, 2> > rot_obj(5, 6, 13);
//	std::cout << rot_obj.T << "\n";
//	std::cout << rot_obj.param.alpha << " " << rot_obj.param.beta << " " << rot_obj.param.h << "\n";

/*	test for integrator	*/
//	state_type x;
//	x[0] = 2.;
//	x[1] = 0.;
//
//	const Harmonic_oscillator ho(1.15);
//	const double t1 = 0;
//	const double t2 = 2. * M_PI;
//	const double dt = 0.0001;

//	IntegrateSystem(ho, x, t1, t2, dt, Integrator_Type::dopri5, true);
//
//	std::cout << x[0] << " " << x[1] << "\n";

	return 0;
}
