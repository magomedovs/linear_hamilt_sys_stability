#include "calculate.h"

#include <iostream>

int main()
{
	const double alpha = 0.33;
	
	const size_t BETA_P_NUM = 120;
	const size_t H_OSC_P_NUM = 41;
	const size_t H_ROT_P_NUM = 101;

	CalculateAndWriteToFile<BETA_P_NUM, H_OSC_P_NUM, H_ROT_P_NUM>(alpha);

//	const Rotation_system obj(0.5, 0.85, 1.5);
//	CalculateForPoint<Rotation_system, obj.DIM>(obj);

	return 0;
}


