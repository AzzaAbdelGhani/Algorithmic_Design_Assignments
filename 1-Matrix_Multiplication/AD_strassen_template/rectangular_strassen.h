#ifndef __Rectangular__
#include "stdlib.h"

void rectangular_strassen_aux(float **C, float const *const *const A,
                            float const *const *const B,
                            const size_t C_f_row, const size_t C_f_col,
                            const size_t A_f_row, const size_t A_f_col,
                            const size_t B_f_row, const size_t B_f_col,
                            const size_t Row_A, const size_t Col_A, const size_t Col_B);

void rectangular_strassen_matrix_multiplication(float **C, float const *const *const A,
                                                float const *const *const B, size_t Row_A, 
                                                size_t Col_A, size_t Col_B);


#endif // Rectangular