#ifndef CALCULATE
#define CALCULATE

#include "eigenvalues_calc.h"
#include "linspace.h"
#include "ode_system_classes.h"
#include "monodromy_matrix.h"

#include <iostream>
#include <fstream>
#include <iomanip>
#include <array>
#include <vector>
#include <complex>
#include <string>

const size_t PARAM_W = 14;
const size_t EIG_W = 34;
const size_t STABIL_STAT_W = 19;

template <typename ODE_obj_T, size_t DIM>
void CalculateForPoint(const ODE_obj_T& obj) {
    typedef std::array<double, DIM> state_type;
    std::array<state_type, DIM> monodr_mat = Monodromy_matrix<ODE_obj_T, DIM>(obj, Integrator_Type::fehlberg78);
    std::array< std::complex<double>, DIM > eigvalsarr( eigenvalues_calc(monodr_mat) );
    Solution_type stability_status = stability_investigation(eigvalsarr);

	std::cout << "alpha = " << obj.param.alpha << "\n"
		<< "beta = " << obj.param.beta << "\n"
		<< "h = " << obj.param.h << "\n";

	std::cout << "eigenvalues : ";
	for (const auto& eig : eigvalsarr) {
		std::cout << eig << " ";
    }
	std::cout << "\n";
	std::cout << "status : " << Solution_type_to_string(stability_status) << "\n";
}

//template <size_t DATA_SIZE>
//void WriteToFile(const std::array<double, DATA_SIZE>& alpha_vec, const std::array<double, DATA_SIZE>& beta_vec,
//        const std::array<double, DATA_SIZE>& h_vec,
//        const std::array< std::complex<double>, DATA_SIZE >& eig1_vec, const std::array< std::complex<double>, DATA_SIZE >& eig2_vec,
//        const std::array<std::string, DATA_SIZE>& stab_stat) {

template <typename T>
void WriteToFile(const std::vector<T>& alpha_vec, const std::vector<T>& beta_vec,
        const std::vector<T>& h_vec,
        const std::vector< std::complex<T> >& eig1_vec, const std::vector< std::complex<T> >& eig2_vec,
        const std::vector<std::string>& stab_stat) {

	const size_t DATA_SIZE = alpha_vec.size();

	std::string filename = "output_stability_data_alpha" + std::to_string(alpha_vec[0]) + ".txt";
    std::ofstream output(filename);
    output << std::setw(PARAM_W) << "alpha" << " "
        << std::setw(PARAM_W) << "beta" << " "
        << std::setw(PARAM_W) << "h" << " "
        << std::setw(EIG_W) << "eigenvalue_1" << " "
        << std::setw(EIG_W) << "eigenvalue_2" << " "
        << std::setw(STABIL_STAT_W) << "stability_status" << "\n";

    output << std::fixed << std::setprecision(8);
    for (size_t i = 0; i < DATA_SIZE; ++i) {
        output << std::setw(PARAM_W) << alpha_vec[i] << " "
            << std::setw(PARAM_W) << beta_vec[i] << " "
            << std::setw(PARAM_W) << h_vec[i] << " "
            << std::setw(EIG_W) << eig1_vec[i] << " "
            << std::setw(EIG_W) << eig2_vec[i] << " "
            << std::setw(STABIL_STAT_W) << stab_stat[i] << "\n";
    }

    output.close();
}

//template <typename ODE_obj_T, size_t DIM, size_t DATA_SIZE>
//void Calculate(const ODE_obj_T& obj,
//        std::array< std::complex<double>, DATA_SIZE >& eig1_vec, std::array< std::complex<double>, DATA_SIZE >& eig2_vec,
//        std::array<std::string, DATA_SIZE>& stab_stat, size_t i) {

template <typename ODE_obj_T, size_t DIM>
void Calculate(const ODE_obj_T& obj,
        std::vector< std::complex<double> >& eig1_vec, std::vector< std::complex<double> >& eig2_vec,
        std::vector<std::string>& stab_stat, size_t i) {

	typedef std::array<double, DIM> state_type;
    std::array<state_type, DIM> monodr_mat = Monodromy_matrix<ODE_obj_T, DIM>(obj, Integrator_Type::fehlberg78);
    std::array< std::complex<double>, DIM > eigvalsarr( eigenvalues_calc(monodr_mat) );
    Solution_type stability_status = stability_investigation(eigvalsarr);
    eig1_vec[i] = eigvalsarr[0];
    eig2_vec[i] = eigvalsarr[1];
    stab_stat[i] = Solution_type_to_string(stability_status);
}

template <size_t BETA_P_NUM, size_t H_OSC_P_NUM, size_t H_ROT_P_NUM>
void CalculateAndWriteToFile(const double alpha) { 

    const size_t DATA_SIZE = BETA_P_NUM * (H_OSC_P_NUM + H_ROT_P_NUM);

//	std::array<double, BETA_P_NUM> beta_span;
    std::vector<double> beta_span(BETA_P_NUM);
    linspace(BetaLeftConstrCalc(alpha) + 0.01, BetaRightConstrCalc(alpha) - 0.01, beta_span);

    const double SHIFT = 1.0e-3;

//	std::array<double, H_OSC_P_NUM> h_osc_span;
    std::vector<double> h_osc_span(H_OSC_P_NUM);
    linspace(-1. + SHIFT, 1. - SHIFT, h_osc_span);
 
//	std::array<double, H_ROT_P_NUM> h_rot_span;
	std::vector<double> h_rot_span(H_ROT_P_NUM);
    linspace(1. + SHIFT, 7., h_rot_span);

//	std::array<double, DATA_SIZE> alpha_vec, beta_vec, h_vec;
//	std::array< std::complex<double>, DATA_SIZE > eig1_vec, eig2_vec;
//	std::array<std::string, DATA_SIZE> stab_stat;
    std::vector<double> alpha_vec(DATA_SIZE), beta_vec(DATA_SIZE), h_vec(DATA_SIZE);
    std::vector< std::complex<double> > eig1_vec(DATA_SIZE), eig2_vec(DATA_SIZE);
    std::vector<std::string> stab_stat(DATA_SIZE);

	size_t i = 0;
    for (const double beta : beta_span) {
        for (const double h : h_osc_span) {
            alpha_vec[i] = alpha;
            beta_vec[i] = beta;
            h_vec[i] = h;
            const Oscillation_system osc_obj(alpha, beta, h);
//          Calculate<Oscillation_system, Base_system::DIM, DATA_SIZE>(osc_obj, eig1_vec, eig2_vec, stab_stat, i);
            Calculate<Oscillation_system, Base_system::DIM>(osc_obj, eig1_vec, eig2_vec, stab_stat, i);
            ++i;
        }
        for (const double h : h_rot_span) {
            alpha_vec[i] = alpha;
            beta_vec[i] = beta;
            h_vec[i] = h;
            const Rotation_system rot_obj(alpha, beta, h);
//          Calculate<Rotation_system, Base_system::DIM, DATA_SIZE>(rot_obj, eig1_vec, eig2_vec, stab_stat, i);
            Calculate<Rotation_system, Base_system::DIM>(rot_obj, eig1_vec, eig2_vec, stab_stat, i);
            ++i;
        }
    }

    WriteToFile(alpha_vec, beta_vec, h_vec, eig1_vec, eig2_vec, stab_stat);
}

#endif
