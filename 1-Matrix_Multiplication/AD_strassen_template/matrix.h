#ifndef __MATRIX__
#include <stdlib.h>

void naive_matrix_multiplication(float **C, float const *const *const A,
                                float const *const *const B,
                                const size_t n);

void rectangular_naive_matrix_multiplication(float **C, float const *const *const A,
                                        float const *const *const B,
                                        size_t Row_A, size_t Col_A, size_t Col_B);

int same_matrix(float const *const *const A, float const *const *const B,
                const size_t rows, const size_t cols);

float **allocate_matrix(const size_t rows, const size_t cols);
void deallocate_matrix(float **A, const size_t rows);


float **allocate_random_matrix(const size_t rows, const size_t cols);
#endif //__MATRIX__
