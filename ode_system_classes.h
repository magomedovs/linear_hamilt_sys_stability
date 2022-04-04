#ifndef ODE_SYSTEM_CLASSES
#define ODE_SYSTEM_CLASSES

#include <array>

/* rhs functional coefficients for ode system */
double F_20_func(double q1, double p1, double alpha, double beta); 
double F_11_func(double q1, double p1, double alpha, double beta); 
double F_02_func(double q1, double p1, double alpha, double beta);

/* Functions for calculation of the left and right constraints of the parameter beta depending on the parameter alpha */
double BetaLeftConstrCalc(double alpha);
double BetaRightConstrCalc(double alpha);

/* Parameters of the problem */
struct Params {
    double alpha;	// > 0 and != 1
    double beta;	// > 0
    double h;	// > -1

	double beta_left_constraint;
	double beta_right_constraint;
    Params();
    Params(double alpha_new, double beta_new, double h_new);
};

/* Base class for oscillations and rotations systems */
class Base_system {
public:
    Base_system(double period);
    Base_system(double period, double alpha_new, double beta_new, double h_new);

	virtual double k_func(double h) const = 0;
	virtual double omega_func(double k) const = 0;
	virtual double u_func(double psi) const = 0;
	virtual double q1_func(double u) const = 0;
	virtual double p1_func(double u) const = 0;
	
	double f_11_func(double psi) const;
	double f_02_func(double psi) const;
    double f_20_func(double psi) const;	

	static const size_t DIM = 2;
	const double T; // period
    const Params param; // parameters
};

/* Oscillations. k = k1 */
class Oscillation_system: public Base_system {
    const double k;
    const double omega;
public:
    Oscillation_system();
    Oscillation_system(double k_new);
    Oscillation_system(double alpha_new, double beta_new, double h_new);

    double k_func(double h) const;	// -1 < h < 1
    double omega_func(double k) const;

    double u_func(double psi) const;
    double q1_func(double u) const;
    double p1_func(double u) const;

	/* x[0] := q_2, x[1] := p_2 */
    void operator() (const std::array<double, DIM> &x , std::array<double, DIM> &dxdt , const double t);
};

/* Rotations. k = k2 */
class Rotation_system: public Base_system {
    const double k;
    const double omega;
public:
    Rotation_system();
    Rotation_system(double k_new);
    Rotation_system(double alpha_new, double beta_new, double h_new);

    double k_func(double h) const;	// 1 < h
    double omega_func(double k) const;

    double u_func(double psi) const;
    double q1_func(double u) const;
    double p1_func(double u) const;

	/* x[0] := q_2, x[1] := p_2 */
    void operator() (const std::array<double, DIM> &x , std::array<double, DIM> &dxdt , const double t);
};

#endif
