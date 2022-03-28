#include "linspace.h"
#include "ode_system_classes.h"

#include <iostream>
#include <fstream>
#include <vector>
#include <string>

const size_t N = 1000;

template <typename Ode_System_t>
void save_fij(const Ode_System_t& obj, double t0, double t1) {
	std::vector<double> x_vec(N);
	linspace(t0, t1, x_vec);
	{
	std::string str = "f_20";
	std::string filename = "test_output_" + str + ".txt";
	std::ofstream output(filename);
	for (double& num : x_vec) {
		output << num << " " << obj.f_20_func(num) << "\n";
	}
	output.close();
	}
	{
	std::string str = "f_11";
	std::string filename = "test_output_" + str + ".txt";
	std::ofstream output(filename);
	for (double& num : x_vec) {
		output << num << " " << obj.f_11_func(num) << "\n";
	}
	output.close();
	}
	{
	std::string str = "f_02";
	std::string filename = "test_output_" + str + ".txt";
	std::ofstream output(filename);
	for (double& num : x_vec) {
		output << num << " " << obj.f_02_func(num) << "\n";
	}
	output.close();
	}

}

int main()
{
//	Rotation_system Oscillation_system	std::array<double, 2>
	Rotation_system obj(0.5, 0.7, 2.7);

	save_fij(obj, 0, 2. * obj.T);

//	std::cout << obj.f_20_func(2.) << "\n";

	return 0;
}

//void save_fij(double (*f)(double), double t0, double t1, std::string str) {
//	std::vector<double> x_vec(N);
//	linspace(t0, t1, x_vec);
//
//	std::string filename = "test_output_" + str + ".txt";
//	std::ofstream output(filename);
//	for (double& num : x_vec) {
//		output << num << " " << (*f)(num) << "\n";
//	}
//	output.close();
//}


//	save_fij(osc_obj.f_20_func, 0, 2. * osc_obj.T, "f20");
