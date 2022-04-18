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
#include <thread>
#include <utility>

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

//template <typename T>
void CalcThread(double alpha, std::vector<double>& beta_vec,
        std::vector<double>& h_vec,
        std::vector< std::complex<double> >& eig1_vec, std::vector< std::complex<double> >& eig2_vec,
        std::vector<std::string>& stab_stat, 
		const std::vector<double>& beta_span, const std::vector<double>& h_osc_span, const std::vector<double>& h_rot_span,
		std::pair<size_t, size_t> interval, size_t i) {

	for (size_t beta_ind = interval.first; beta_ind < interval.second; ++beta_ind) {
        for (const double h : h_osc_span) {
            beta_vec[i] = beta_span[beta_ind];
            h_vec[i] = h;
            const Oscillation_system osc_obj(alpha, beta_span[beta_ind], h);
            Calculate<Oscillation_system, Base_system::DIM>(osc_obj, eig1_vec, eig2_vec, stab_stat, i);
            ++i;
        }
        for (const double h : h_rot_span) {
            beta_vec[i] = beta_span[beta_ind];
            h_vec[i] = h;
            const Rotation_system rot_obj(alpha, beta_span[beta_ind], h);
            Calculate<Rotation_system, Base_system::DIM>(rot_obj, eig1_vec, eig2_vec, stab_stat, i);
            ++i;
        }
    }

}

template <size_t BETA_P_NUM, size_t H_OSC_P_NUM, size_t H_ROT_P_NUM>
void CalculateAndWriteToFile(const double alpha) { 

    const size_t DATA_SIZE = BETA_P_NUM * (H_OSC_P_NUM + H_ROT_P_NUM);

    std::vector<double> beta_span(BETA_P_NUM);
    linspace(BetaLeftConstrCalc(alpha) + 0.01, BetaRightConstrCalc(alpha) - 0.01, beta_span);

    const double SHIFT = 1.0e-3;

    std::vector<double> h_osc_span(H_OSC_P_NUM);
    linspace(-1. + SHIFT, 1. - SHIFT, h_osc_span);
 
	std::vector<double> h_rot_span(H_ROT_P_NUM);
    linspace(1. + SHIFT, 7., h_rot_span);

    std::vector<double> alpha_vec(DATA_SIZE, alpha), beta_vec(DATA_SIZE), h_vec(DATA_SIZE);
    std::vector< std::complex<double> > eig1_vec(DATA_SIZE), eig2_vec(DATA_SIZE);
    std::vector<std::string> stab_stat(DATA_SIZE);

	const size_t NUM_OF_THREADS = 2;
	std::vector< std::thread > threads;//(NUM_OF_THREADS);

	std::vector< std::pair<size_t, size_t> > intervals(NUM_OF_THREADS);
	for (size_t i = 0; i < NUM_OF_THREADS; ++i) {
		size_t from_ind = i * static_cast<size_t>(BETA_P_NUM / NUM_OF_THREADS);
		size_t to_ind = (i + 1) * static_cast<size_t>(BETA_P_NUM / NUM_OF_THREADS);
		
		if (i == NUM_OF_THREADS - 1) {
			to_ind = BETA_P_NUM;
		}

		intervals[i] = std::make_pair(from_ind, to_ind);
	}
	
	for (size_t i = 0; i < NUM_OF_THREADS; ++i) {
		threads[i] = std::thread(CalcThread, 
				alpha, std::ref(beta_vec), std::ref(h_vec), 
				std::ref(eig1_vec), std::ref(eig2_vec), std::ref(stab_stat),
				std::ref(beta_span), std::ref(h_osc_span), std::ref(h_rot_span),
				intervals[i], intervals[i].first * (H_OSC_P_NUM + H_ROT_P_NUM));
	}

	for (size_t i = 0; i < NUM_OF_THREADS; ++i) {
		threads[i].join();
	}

    WriteToFile(alpha_vec, beta_vec, h_vec, eig1_vec, eig2_vec, stab_stat);
}

#endif
