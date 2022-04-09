#include "calculate.h"

#include <chrono>
#include <iostream>

int main()
{
	const double alpha = 0.3;
	
	const size_t BETA_P_NUM = 120;
	const size_t H_OSC_P_NUM = 41;
	const size_t H_ROT_P_NUM = 101;

	std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

	CalculateAndWriteToFile<BETA_P_NUM, H_OSC_P_NUM, H_ROT_P_NUM>(alpha);

//	const Rotation_system obj(0.5, 0.85, 1.5);
//	CalculateForPoint<Rotation_system, obj.DIM>(obj);

	std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
	
	std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::seconds> (end - begin).count() << "[s]" << std::endl;
	std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::minutes> (end - begin).count() << "[min]" << std::endl;
	
	return 0;
}


