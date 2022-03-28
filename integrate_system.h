#ifndef INTEGRATE_SYSTEM
#define INTEGRATE_SYSTEM

#include "save_solution_function.h"

#include <boost/numeric/odeint.hpp>

#include <vector>
#include <array>
#include <cmath>

template <typename state_type>
class Push_back_state_time {
public:
    std::vector<state_type>& m_states;
    std::vector<double>& m_times;

    Push_back_state_time(std::vector<state_type>& states, std::vector<double>& times)
        : m_states(states), m_times(times) {}

    void operator() (const state_type& x, double t) {
        m_states.push_back(x);
        m_times.push_back(t);
    }
};

// type of RK integrator
enum class Integrator_Type {
	dopri5,
	cash_karp54,
	rosenbrock4,
	fehlberg78
};

// write saving to file instructions !!!
template <typename ODE_obj_T, typename state_type>
void IntegrateSystem(const ODE_obj_T& ode_system_obj, 
		state_type& x, const double t_begin, const double t_end, const double dt,
		Integrator_Type integrator_type=Integrator_Type::dopri5, bool save_to_file_flag=false,
		const double abs_er_tol=1.0e-13, const double rel_er_tol=1.0e-12) {

	if (save_to_file_flag) {

		std::vector<state_type> x_vec;
		std::vector<double> t_vec;

		using namespace boost::numeric::odeint;

		if (integrator_type == Integrator_Type::dopri5) {
			typedef runge_kutta_dopri5< state_type > error_stepper_type;
			integrate_adaptive( make_controlled< error_stepper_type >(abs_er_tol, rel_er_tol),
                ode_system_obj, x, t_begin, t_end, dt, Push_back_state_time< state_type >(x_vec, t_vec) );
		} else if (integrator_type == Integrator_Type::cash_karp54) {
			typedef runge_kutta_cash_karp54< state_type > error_stepper_type;
			integrate_adaptive( make_controlled< error_stepper_type >(abs_er_tol, rel_er_tol),
                ode_system_obj, x, t_begin, t_end, dt, Push_back_state_time< state_type >(x_vec, t_vec) );
		} else if (integrator_type == Integrator_Type::rosenbrock4) {
//			typedef rosenbrock4< state_type > error_stepper_type;
//			integrate_adaptive( make_controlled< error_stepper_type >(abs_er_tol, rel_er_tol),
//                ode_system_obj, x, t_begin, t_end, dt, Push_back_state_time< state_type >(x_vec, t_vec) );
		} else { // Integrator_Type::fehlberg78
			typedef runge_kutta_fehlberg78< state_type > error_stepper_type;
			integrate_adaptive( make_controlled< error_stepper_type >(abs_er_tol, rel_er_tol),
                ode_system_obj, x, t_begin, t_end, dt, Push_back_state_time< state_type >(x_vec, t_vec) );
		}
//		integrate_adaptive( make_controlled< error_stepper_type >(abs_er_tol, rel_er_tol), 
//				ode_system_obj, x, t_begin, t_end, dt, Push_back_state_time< state_type >(x_vec, t_vec) );
	
		SaveSolutionIntoFile(x_vec, t_vec);

	} else {
	
		using namespace boost::numeric::odeint;
 
       	if (integrator_type == Integrator_Type::dopri5) {
			typedef runge_kutta_dopri5< state_type > error_stepper_type;
			integrate_adaptive( make_controlled< error_stepper_type >(abs_er_tol, rel_er_tol),
                ode_system_obj, x, t_begin, t_end, dt );
		} else if (integrator_type == Integrator_Type::cash_karp54) {
			typedef runge_kutta_cash_karp54< state_type > error_stepper_type;
			integrate_adaptive( make_controlled< error_stepper_type >(abs_er_tol, rel_er_tol),
                ode_system_obj, x, t_begin, t_end, dt );
		} else if (integrator_type == Integrator_Type::rosenbrock4) {
//			typedef rosenbrock4< state_type > error_stepper_type;
//			integrate_adaptive( make_controlled< error_stepper_type >(abs_er_tol, rel_er_tol),
//                ode_system_obj, x, t_begin, t_end, dt );
		} else { // Integrator_Type::fehlberg78
			typedef runge_kutta_fehlberg78< state_type > error_stepper_type;
			integrate_adaptive( make_controlled< error_stepper_type >(abs_er_tol, rel_er_tol),
                ode_system_obj, x, t_begin, t_end, dt );
		}
//		integrate_adaptive( make_controlled< error_stepper_type >(abs_er_tol, rel_er_tol), 
//				ode_system_obj, x, t_begin, t_end, dt );

	}
	
}

#endif

//  runge_kutta_cash_karp54
//  runge_kutta_dopri5
//  runge_kutta_fehlberg78

