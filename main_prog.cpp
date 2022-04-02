#include "calculate.h"

#include <iostream>

int main()
{
	const double alpha = 0.33;
	
	const size_t BETA_P_NUM = 120;
	const size_t H_OSC_P_NUM = 41;
	const size_t H_ROT_P_NUM = 101;

	CalculateAndWriteToFile<BETA_P_NUM, H_OSC_P_NUM, H_ROT_P_NUM>(alpha);
		
	return 0;
}


