#include "calculate.h"

#include <vector>
#include <chrono>
#include <iostream>
#include <string>

int main()
{
//	const double alpha = 0.3;
	
	std::vector<double> alpha_vec{0.3, 0.33, 0.37,
		0.41, 0.45, 0.47, 0.49,
		0.51, 0.53, 0.55, 0.57,
		0.61, 0.67, 
		0.7, 0.75};

	const size_t BETA_P_NUM = 150;
	const size_t H_OSC_P_NUM = 67; //41;
	const size_t H_ROT_P_NUM = 201; //101;

	std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

	const size_t NUM_OF_THREADS = 8;

	for (const double alpha : alpha_vec) {
		CalculateAndWriteToFile<BETA_P_NUM, H_OSC_P_NUM, H_ROT_P_NUM>(alpha, NUM_OF_THREADS);
	}

//	const Rotation_system obj(0.5, 0.7, 2.5);
//	CalculateForPoint<Rotation_system, obj.DIM>(obj);

	std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
	
	int s_dur = std::chrono::duration_cast<std::chrono::seconds> (end - begin).count();
	int ms_dur = std::chrono::duration_cast<std::chrono::milliseconds> (end - begin).count();

	std::stringstream time_difference;
	time_difference << "Time difference = " 
		<< s_dur / 60 << " [min] : " 
		<< s_dur % 60 << " [s] : " 
		<< ms_dur % 1000 << " [ms]";
	std::cout << time_difference.str() << "\n";
	
	return 0;
}


