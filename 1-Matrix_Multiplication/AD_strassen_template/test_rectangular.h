#ifndef __TEST_Rectangular__
#include "stdlib.h"

double test_rectangular(void (*f)(float **,
	                  float const *const *const,
	                  float const *const *const,
	                  size_t, size_t, size_t), 
	                  float **C, float** A, float **B, size_t Row_A, size_t Col_A, size_t Col_B);

#endif // __TEST_Rectangular__