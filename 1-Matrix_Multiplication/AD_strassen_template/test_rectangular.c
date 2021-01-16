#include "stdlib.h"
#include "time.h"

double test_rectangular(void (*f)(float **,
	                  float const *const *const,
	                  float const *const *const,
	                  size_t, size_t, size_t), 
	                  float **C, float** A, float **B, size_t Row_A, size_t Col_A, size_t Col_B)
{
    struct timespec requestStart, requestEnd;
    double accum;
    size_t rep = 1;

    clock_gettime(CLOCK_REALTIME, &requestStart);
    for (size_t i = 0; i < rep; i++) {
        f(C, (float const *const *const)A,
        (float const *const *const)B, Row_A, Col_A, Col_B);
    }

    clock_gettime(CLOCK_REALTIME, &requestEnd);

    accum = (requestEnd.tv_sec - requestStart.tv_sec) +
            (requestEnd.tv_nsec - requestStart.tv_nsec) / 1E9;
    
    return accum / rep;
}                      