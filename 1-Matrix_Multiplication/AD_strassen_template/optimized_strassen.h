#ifndef __OPTIMIZED__
#include "stdlib.h"

void optimized_strassen_aux(float **C, float const *const *const A,
                  float const *const *const B, 
                  const size_t C_f_row, const size_t C_f_col,
                  const size_t A_f_row, const size_t A_f_col,
                  const size_t B_f_row, const size_t B_f_col,
                  const size_t n);

void optimized_strassen_matrix_multiplication(float **C, float const *const *const A,
                                    float const *const *const B, size_t n);                                  
#endif //__OPTIMIZED__