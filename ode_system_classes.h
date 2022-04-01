#ifndef ODE_SYSTEM_CLASSES
#define ODE_SYSTEM_CLASSES

// Split this file into several files !!!

#include "jacobi_am_from_boost.h"

#include <boost/math/special_functions/jacobi_elliptic.hpp>
#include <boost/math/special_functions/ellint_1.hpp>

#include <cmath>
#include <vector>

/* rhs functional coefficients for ode system */
double F_20_func(double q1, double p1, double alpha, double beta) {
    return 1./2. * ( (alpha + (beta - alpha) * std::pow(std::sin(q1), 2)) * std::pow(p1, 2) + std::cos(q1) );
}
double F_11_func(double q1, double p1, double alpha, double beta) {
    return (beta - alpha) * p1 * std::sin(q1) * std::cos(q1);
}
double F_02_func(double q1, double p1, double alpha, double beta) {
    return 1./2. * ( beta + (alpha - beta) * std::pow(std::sin(q1), 2) );
}

/* Functions for calculation of the left and right constraints of the parameter beta depending on the parameter alpha */
double BetaLeftConstrCalc(double alpha) {
	return alpha / (alpha + 1);
}
double BetaRightConstrCalc(double alpha) {
	return alpha / std::abs(alpha - 1);
}

/* Parameters of the problem */
struct Params {
    double alpha;	// > 0 and != 1
    double beta;	// > 0
    double h;	// > -1

	double beta_left_constraint;
	double beta_right_constraint;
    Params() {}
    Params(double alpha_new, double beta_new, double h_new) : alpha(alpha_new), beta(beta_new), h(h_new), 
	beta_left_constraint( BetaLeftConstrCalc(alpha) ), beta_right_constraint( BetaRightConstrCalc(alpha) ) {}
};

/* Base class for oscillations and rotations systems */
class Base_system {
public:
    Base_system(double period) : T(period) {}
    Base_system(double period, double alpha_new, double beta_new, double h_new) : T(period), param(alpha_new, beta_new, h_new) {}

	virtual double k_func(double h) const = 0;
	virtual double omega_func(double k) const = 0;
	virtual double u_func(double psi) const = 0;
	virtual double q1_func(double u) const = 0;
	virtual double p1_func(double u) const = 0;
	
	double f_11_func(double psi) const {
        double q1 = q1_func(u_func(psi));
        double p1 = p1_func(u_func(psi));
        return F_11_func(q1, p1, param.alpha, param.beta);
    }
    double f_02_func(double psi) const {
        double q1 = q1_func(u_func(psi));
        double p1 = p1_func(u_func(psi));
        return F_02_func(q1, p1, param.alpha, param.beta);
    }
    double f_20_func(double psi) const {
       double q1 = q1_func(u_func(psi));
       double p1 = p1_func(u_func(psi));
       return F_20_func(q1, p1, param.alpha, param.beta);
    }	

	static const size_t DIM = 2;
	const double T; // period
    const Params param; // parameters
};

/* Oscillations. k = k1 */
class Oscillation_system: public Base_system {
    const double k;
    const double omega;
public:
    Oscillation_system() : Base_system(M_PI), k(0), omega(0) {}
    Oscillation_system(double k_new) : Base_system(M_PI), k(k_new), omega(omega_func(k)) {}
    Oscillation_system(double alpha_new, double beta_new, double h_new) : Base_system(M_PI, alpha_new, beta_new, h_new),
        k(k_func(h_new)), omega(omega_func(k)) {
    }
    double k_func(double h) const {   // -1 < h < 1
        return std::sqrt((h + 1.)/2.);
    }
    double omega_func(double k) const {
        return M_PI / (2. * boost::math::ellint_1(k));
    }

    double u_func(double psi) const {
        return 2. * boost::math::ellint_1(k) / M_PI * psi;
    }
    double q1_func(double u) const {
        return 2. * std::asin(k * boost::math::jacobi_sn(k, u));
    }
    double p1_func(double u) const {
        return 2. * k * boost::math::jacobi_cn(k, u);
    }

	/* x[0] := q_2, x[1] := p_2 */
    void operator() (const std::array<double, DIM> &x , std::array<double, DIM> &dxdt , const double t) {
		dxdt[0] = 1./omega * (this->f_11_func(t) * x[0] + 2. * this->f_02_func(t) * x[1]);
		dxdt[1] = -1./omega * (2. * this->f_20_func(t) * x[0] + this->f_11_func(t) * x[1]);
    }

//  void getOmega() const {
//      std::cout << omega << "\n";
//  }

};

/* Rotations. k = k2 */
class Rotation_system: public Base_system {
    const double k;
    const double omega;
public:
    Rotation_system() : Base_system(2. * M_PI), k(0), omega(0) {}
    Rotation_system(double k_new) : Base_system(2. * M_PI), k(k_new), omega(omega_func(k)) {}
    Rotation_system(double alpha_new, double beta_new, double h_new) : Base_system(2. * M_PI, alpha_new, beta_new, h_new),
        k(k_func(h_new)), omega(omega_func(k)) {
    }
    double k_func(double h) const { // 1 < h
        return std::sqrt(2./(h + 1.));
    }
    double omega_func(double k) const {
        return M_PI / (k * boost::math::ellint_1(k));
    }

    double u_func(double psi) const {
        return boost::math::ellint_1(k) / M_PI * psi;
    }
    double q1_func(double u) const {
        return 2. * jacobi_am(k, u);
    }
    double p1_func(double u) const {
        return 2. / k * boost::math::jacobi_dn(k, u);
    }

	/* x[0] := q_2, x[1] := p_2 */
    void operator() (const std::array<double, DIM> &x , std::array<double, DIM> &dxdt , const double t) {
		dxdt[0] = 1./omega * (this->f_11_func(t) * x[0] + 2. * this->f_02_func(t) * x[1]);
		dxdt[1] = -1./omega * (2. * this->f_20_func(t) * x[0] + this->f_11_func(t) * x[1]);
    }

//  void GetK() const {
//      std::cout << k << "\n";
//  }

};

#endif
