#ifndef __GENERALIZED__

void general_sub_matrix_blocks(float **C, float const *const *const A,
                  float const *const *const B, 
                  const size_t C_f_row, const size_t C_f_col,
                  const size_t A_f_row, const size_t A_f_col,
                  const size_t B_f_row, const size_t B_f_col,
                  const size_t rows, const size_t cols);
                  
void general_sum_matrix_blocks(float **C, float const *const *const A,
                  float const *const *const B, 
                  const size_t C_f_row, const size_t C_f_col,
                  const size_t A_f_row, const size_t A_f_col,
                  const size_t B_f_row, const size_t B_f_col,
                  const size_t rows, const size_t cols);

void general_naive_aux(float **C, float const *const *const A,
                  float const *const *const B, 
                  const size_t C_f_row, const size_t C_f_col,
                  const size_t A_f_row, const size_t A_f_col,
                  const size_t B_f_row, const size_t B_f_col,
                  const size_t Row_A, const size_t Col_A,
                  const size_t Col_B);

void repairing(float **C, float const *const *const A,
                    float const *const *const B,
                    const size_t C_f_row, const size_t C_f_col,
                    const size_t A_f_row, const size_t A_f_col,
                    const size_t B_f_row, const size_t B_f_col,
                    const size_t n, const size_t m);

void general_strassen_aux(float **C, float const *const *const A,
                  float const *const *const B, 
                  const size_t C_f_row, const size_t C_f_col,
                  const size_t A_f_row, const size_t A_f_col,
                  const size_t B_f_row, const size_t B_f_col,
                  const size_t n);

void general_strassen_matrix_multiplication(float **C, float const *const *const A,
                                    float const *const *const B, size_t n);                                  
#endif //__STRASSEN__