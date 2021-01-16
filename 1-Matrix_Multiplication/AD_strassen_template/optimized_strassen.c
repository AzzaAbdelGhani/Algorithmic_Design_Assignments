#include "matrix.h"
#include "optimized_strassen.h"
#include "generalized_strassen.h"

void optimized_strassen_aux(float **C, float const *const *const A,
                  float const *const *const B, 
                  const size_t C_f_row, const size_t C_f_col,
                  const size_t A_f_row, const size_t A_f_col,
                  const size_t B_f_row, const size_t B_f_col,
                  const size_t n) 
{
    if (n <= 128) 
    {
        general_naive_aux(C,A,B,
                  C_f_row,C_f_col,
                  A_f_row,A_f_col,
                  B_f_row,B_f_col,
                  n, n, n);
        return;
    }
    size_t n2; // this is the size of the blocks
    if (n % 2 == 0)
    {
        n2 = n/2;
    }
    else
    {
        n2 = (n-1)/2;
    }
    
    

    
    float **S0 = allocate_matrix(n2,n2);
    float **P0 = allocate_matrix(n2,n2);
    //S1 = B12 - B22
    general_sub_matrix_blocks(S0, B, B, 0, 0, B_f_row, B_f_col+n2,
                      B_f_row + n2, B_f_col + n2, n2, n2);

    // P1 = A11 * S1
    optimized_strassen_aux(P0, A, (const float* const const*) S0,
                 0, 0, A_f_row, A_f_col, 0, 0, n2);

    deallocate_matrix(S0, n2);

    float **S1 = allocate_matrix(n2,n2);
    // S2 = A11 + A12
    general_sum_matrix_blocks(S1, A, A,
                      0, 0,
                      A_f_row, A_f_col,
                      A_f_row, A_f_col + n2,
                      n2, n2);


    // P2 = S2 x B22
    optimized_strassen_aux(C, (const float* const *const) S1, B,
                 C_f_row, C_f_col + n2,
                 0, 0,
                 B_f_row + n2, B_f_col + n2,
                 n2);
    deallocate_matrix(S1, n2);
    
    float **S2 = allocate_matrix(n2, n2);
    // S3 = A21 + A22
    general_sum_matrix_blocks(S2, A, A,
                      0, 0,
                      A_f_row + n2, A_f_col,
                      A_f_row + n2, A_f_col + n2,
                      n2, n2);
    float **P2 = allocate_matrix(n2, n2);
    // P3 = S3 x B11
    optimized_strassen_aux(P2, (const float* const *const) S2, B,
                 0, 0,
                 0, 0,
                 B_f_row, B_f_col,
                 n2);
    deallocate_matrix(S2, n2);

    float **S3 = allocate_matrix(n2, n2);
    // S4 = B21 - B11
    general_sub_matrix_blocks(S3, B, B,
                      0, 0,
                      B_f_row + n2, B_f_col,
                      B_f_row, B_f_col,
                      n2, n2);

    // P4 = A22 x S4
    optimized_strassen_aux(C, A, (const float* const *const) S3,
                 C_f_row + n2, C_f_col,
                 A_f_row + n2, A_f_col + n2,
                 0, 0,
                 n2);
    deallocate_matrix(S3, n2);

    float **S4 = allocate_matrix(n2, n2);
    // S5 = A11 + A22
    general_sum_matrix_blocks(S4, A, A,
                      0, 0,
                      A_f_row, A_f_col,
                      A_f_row + n2, A_f_col + n2,
                      n2, n2);
    float **S5 = allocate_matrix(n2, n2);
    // S6 = B11 + B22
    general_sum_matrix_blocks(S5, B, B,
                      0, 0,
                      B_f_row, B_f_col,
                      B_f_row + n2, B_f_col + n2,
                      n2, n2);

    // P5 = S5 x S6
    optimized_strassen_aux(C, (const float* const *const) S4,
                 (const float* const *const) S5,
                 C_f_row + n2, C_f_col + n2,
                 0, 0,
                 0, 0,
                 n2);
    deallocate_matrix(S4, n2);
    deallocate_matrix(S5, n2);

    float **S6 = allocate_matrix(n2, n2);
    // S7 = A12 - A22
    general_sub_matrix_blocks(S6, A, A,
                      0, 0,
                      A_f_row, A_f_col + n2,
                      A_f_row + n2, A_f_col + n2,
                      n2, n2);
    float **S7 = allocate_matrix(n2, n2);
    // S8 = B21 + B22
    general_sum_matrix_blocks(S7, B, B,
                      0, 0,
                      B_f_row + n2, B_f_col,
                      B_f_row + n2, B_f_col + n2,
                      n2, n2);

    // P6 = S7 x S8
    optimized_strassen_aux(C, (const float* const *const) S6,
                 (const float* const *const) S7,
                 C_f_row, C_f_col,
                 0, 0,
                 0, 0,
                 n2);
    deallocate_matrix(S6, n2);
    deallocate_matrix(S7, n2);

    float **S8 = allocate_matrix(n2, n2);
    // S9 = A11 - A21
    general_sub_matrix_blocks(S8, A, A,
                      0, 0,
                      A_f_row, A_f_col,
                      A_f_row + n2, A_f_col,
                      n2, n2);
    float **S9 = allocate_matrix(n2, n2);
    // S10 = B11 + B12
    general_sum_matrix_blocks(S9, B, B,
                      0, 0,
                      B_f_row, B_f_col,
                      B_f_row, B_f_col + n2,
                      n2, n2);
    float **P6 = allocate_matrix(n2, n2);
    // P7 = S9 x S10
    optimized_strassen_aux(P6, (const float* const *const) S8,
                 (const float* const *const) S9,
                 0, 0,
                 0, 0,
                 0, 0,
                 n2);
    deallocate_matrix(S8, n2);
    deallocate_matrix(S9, n2);

    // C11 = P5 + P4 - P2 + P6
    general_sum_matrix_blocks(C, (const float* const *const) C,
                      (const float* const *const) C,
                      C_f_row, C_f_col,
                      C_f_row, C_f_col,
                      C_f_row + n2, C_f_col + n2,
                      n2, n2);
    general_sum_matrix_blocks(C, (const float* const *const) C,
                      (const float* const *const) C,
                      C_f_row, C_f_col,
                      C_f_row, C_f_col,
                      C_f_row + n2, C_f_col,
                      n2, n2);
    general_sub_matrix_blocks(C, (const float* const *const) C,
                      (const float* const *const) C,
                      C_f_row, C_f_col,
                      C_f_row, C_f_col,
                      C_f_row, C_f_col + n2,
                      n2, n2);

    // C12 = P1 + P2
    general_sum_matrix_blocks(C, (const float* const *const) C,
                      (const float* const *const) P0,
                      C_f_row, C_f_col+n2,
                      C_f_row, C_f_col + n2,
                      0, 0,
                      n2, n2);

    // C21 = P3 + P4
    general_sum_matrix_blocks(C, (const float* const *const) C,
                      (const float* const *const) P2,
                      C_f_row+n2, C_f_col,
                      C_f_row+n2, C_f_col,
                      0, 0,
                      n2, n2);

    // C22 = P5 + P1 - P3 - P7
    general_sum_matrix_blocks(C, (const float* const *const) C,
                      (const float* const *const) P0,
                      C_f_row+n2, C_f_col+n2,
                      C_f_row+n2, C_f_col+n2,
                      0, 0,
                      n2, n2);
    general_sub_matrix_blocks(C, (const float* const *const) C,
                      (const float* const *const) P2,
                      C_f_row+n2, C_f_col+n2,
                      C_f_row+n2, C_f_col+n2,
                      0, 0,
                      n2, n2);
    general_sub_matrix_blocks(C, (const float* const *const) C,
                      (const float* const *const) P6,
                      C_f_row+n2, C_f_col+n2,
                      C_f_row+n2, C_f_col+n2,
                      0, 0,
                      n2, n2);
    deallocate_matrix(P0, n2);
    deallocate_matrix(P2, n2);
    deallocate_matrix(P6, n2);

    // repair if n is odd 
    if(n%2 != 0){
        repairing(C,  A, B,
                    C_f_row, C_f_col,
                    A_f_row, A_f_col,
                    B_f_row, B_f_col,
                    n,1);
    } 
    
}


void optimized_strassen_matrix_multiplication(float **C, float const *const *const A,
                                    float const *const *const B, size_t n) 
{
 
    // IMPLEMENT THE Generalized STRASSEN'S ALGORITHM HERE
    optimized_strassen_aux(C, A, B,
                  0, 0,
                  0, 0,
                  0, 0,
                  n);

}

