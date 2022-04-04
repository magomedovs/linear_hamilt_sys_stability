#include "jacobi_am_from_boost.h"
#include "ode_system_classes.h"

#include <boost/math/special_functions/jacobi_elliptic.hpp>
#include <boost/math/special_functions/ellint_1.hpp>

#include <cmath>
#include <array>

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

/* struct Params */
Params::Params() {}
Params::Params(double alpha_new, double beta_new, double h_new) : alpha(alpha_new), beta(beta_new), h(h_new),
    beta_left_constraint( BetaLeftConstrCalc(alpha) ), beta_right_constraint( BetaRightConstrCalc(alpha) ) {}

/* class Base_system */
Base_system::Base_system(double period) : T(period) {}
Base_system::Base_system(double period, double alpha_new, double beta_new, double h_new) : T(period), param(alpha_new, beta_new, h_new) {}

double Base_system::f_11_func(double psi) const {
	double q1 = q1_func(u_func(psi));
	double p1 = p1_func(u_func(psi));
	return F_11_func(q1, p1, param.alpha, param.beta);
}
double Base_system::f_02_func(double psi) const {
	double q1 = q1_func(u_func(psi));
	double p1 = p1_func(u_func(psi));
	return F_02_func(q1, p1, param.alpha, param.beta);
}
double Base_system::f_20_func(double psi) const {
	double q1 = q1_func(u_func(psi));
	double p1 = p1_func(u_func(psi));
	return F_20_func(q1, p1, param.alpha, param.beta);
}

/* class Oscillation_system */
Oscillation_system::Oscillation_system() : Base_system(M_PI), k(0), omega(0) {}
Oscillation_system::Oscillation_system(double k_new) : Base_system(M_PI), k(k_new), omega(omega_func(k)) {}
Oscillation_system::Oscillation_system(double alpha_new, double beta_new, double h_new) : Base_system(M_PI, alpha_new, beta_new, h_new), 
	k(k_func(h_new)), omega(omega_func(k)) {}

double Oscillation_system::k_func(double h) const {   // -1 < h < 1
	return std::sqrt((h + 1.)/2.);
}
double Oscillation_system::omega_func(double k) const {
    return M_PI / (2. * boost::math::ellint_1(k));
}

double Oscillation_system::u_func(double psi) const {
    return 2. * boost::math::ellint_1(k) / M_PI * psi;
}
double Oscillation_system::q1_func(double u) const {
    return 2. * std::asin(k * boost::math::jacobi_sn(k, u));
}
double Oscillation_system::p1_func(double u) const {
    return 2. * k * boost::math::jacobi_cn(k, u);
}

/* x[0] := q_2, x[1] := p_2 */
void Oscillation_system::operator() (const std::array<double, DIM> &x , std::array<double, DIM> &dxdt , const double t) {
    dxdt[0] = 1./omega * (f_11_func(t) * x[0] + 2. * f_02_func(t) * x[1]);
    dxdt[1] = -1./omega * (2. * f_20_func(t) * x[0] + f_11_func(t) * x[1]);
}

/* class Rotation_system */
Rotation_system::Rotation_system() : Base_system(2. * M_PI), k(0), omega(0) {}
Rotation_system::Rotation_system(double k_new) : Base_system(2. * M_PI), k(k_new), omega(omega_func(k)) {}
Rotation_system::Rotation_system(double alpha_new, double beta_new, double h_new) : Base_system(2. * M_PI, alpha_new, beta_new, h_new),
    k(k_func(h_new)), omega(omega_func(k)) {}

double Rotation_system::k_func(double h) const { // 1 < h
    return std::sqrt(2./(h + 1.));
}
double Rotation_system::omega_func(double k) const {
    return M_PI / (k * boost::math::ellint_1(k));
}

double Rotation_system::u_func(double psi) const {
    return boost::math::ellint_1(k) / M_PI * psi;
}
double Rotation_system::q1_func(double u) const {
    return 2. * jacobi_am(k, u);
}
double Rotation_system::p1_func(double u) const {
    return 2. / k * boost::math::jacobi_dn(k, u);
}

/* x[0] := q_2, x[1] := p_2 */
void Rotation_system::operator() (const std::array<double, DIM> &x , std::array<double, DIM> &dxdt , const double t) {
    dxdt[0] = 1./omega * (f_11_func(t) * x[0] + 2. * f_02_func(t) * x[1]);
    dxdt[1] = -1./omega * (2. * f_20_func(t) * x[0] + f_11_func(t) * x[1]);
}

