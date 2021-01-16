#include "matrix.h"
#include "generalized_strassen.h"

void rectangular_strassen_aux(float **C, float const *const *const A,
                            float const *const *const B,
                            const size_t C_f_row, const size_t C_f_col,
                            const size_t A_f_row, const size_t A_f_col,
                            const size_t B_f_row, const size_t B_f_col,
                            const size_t Row_A, const size_t Col_A, const size_t Col_B)
{
    if (Row_A == Col_A && Col_A == Col_B){
        general_strassen_aux(C, A, B,
                        C_f_row, C_f_col, 
                        A_f_row, A_f_col,
                        B_f_row, B_f_col,
                        Col_B);
        return; 
    }
    
    //Block size is the minimum among the 3 dimensions 
    size_t min_dim = (Row_A < Col_A)?((Row_A < Col_B)? Row_A : Col_B):((Col_A < Col_B)? Col_A : Col_B); 

    if (min_dim <= 128){
        general_naive_aux(C, A, B,
                    C_f_row, C_f_col,
                    A_f_row, A_f_col, 
                    B_f_row, B_f_col,
                    Row_A, Col_A, Col_B);
        return;
    }
    size_t m = Row_A-min_dim;
    size_t k = Col_A-min_dim;
    size_t n = Col_B-min_dim;
    
    general_strassen_aux(C, A, B,
                    C_f_row, C_f_col, 
                    A_f_row, A_f_col,
                    B_f_row, B_f_col,
                    min_dim); 

    // C11 = A11 * B11 if min_dim == Col_A
    // C11 = A11 * B11 + A12 * B21 if min_dim == Row_A or Col_B
    if ((min_dim == Row_A && Row_A != Col_A) || (min_dim == Col_B && Col_B != Col_A)){
        float **A12A21 = allocate_matrix(min_dim,min_dim);
        rectangular_strassen_aux(A12A21, A, B,
                            0, 0,
                            A_f_row, A_f_col + min_dim,
                            B_f_row + min_dim, B_f_col,
                            min_dim,k,min_dim);
        general_sum_matrix_blocks(C, (const float* const *const)C, (const float* const *const)A12A21,
                            C_f_row, C_f_col,
                            C_f_row, C_f_col,
                            0, 0,
                            min_dim, min_dim);
        deallocate_matrix(A12A21, min_dim);
        if(Row_A == Col_B){
            return;
        }
    }

    if((min_dim == Row_A && Row_A != Col_B) || (min_dim == Col_A && Col_A != Col_B)){
        //C12 = A11 * B12 if min_dim == Col_A or Row_A 
        //C12 = A11 * B12 + A12 * B22 if min_dim == Row_A  & Row_A  != Col_A
       
        rectangular_strassen_aux(C, A, B,
                            C_f_row, C_f_col + min_dim, 
                            A_f_row, A_f_col,
                            B_f_row, B_f_col + min_dim,
                            min_dim, min_dim, n); 
        if(Row_A == Col_A){
            return;
        }

        if (min_dim == Row_A){            
            float **A12B22 = allocate_matrix(min_dim,n);
            rectangular_strassen_aux(A12B22, A, B,
                                0, 0, 
                                A_f_row, A_f_col + min_dim,
                                B_f_row + min_dim, B_f_col + min_dim,
                                min_dim, k, n); 
            general_sum_matrix_blocks(C, (const float* const *const)C, (const float* const *const)A12B22,
                                    C_f_row, C_f_col + min_dim,
                                    C_f_row, C_f_col + min_dim,
                                    0, 0,
                                    min_dim, n);
            deallocate_matrix(A12B22,min_dim);
        }                                                                  
    }

    if((min_dim == Col_B && Col_B != Row_A)|| (min_dim == Col_A && Col_A != Row_A)){
        //C21 = A21 * B11 if min_dim == Col_A
        //C21 = A21 * B11 + A22 * B21 if min_dim == Col_B
        rectangular_strassen_aux(C, A, B,
                            C_f_row + min_dim, C_f_col,
                            A_f_row + min_dim, A_f_col,
                            B_f_row, B_f_col,
                            m, min_dim, min_dim);
        if(Col_B == Col_A){
            return;
        }
        if (min_dim == Col_B)
        {
            float **A22B21 = allocate_matrix(m,min_dim);
            rectangular_strassen_aux(A22B21, A, B,
                                0, 0, 
                                A_f_row + min_dim, A_f_col + min_dim,
                                B_f_row + min_dim, B_f_col,
                                m, k, min_dim); 
            general_sum_matrix_blocks(C, (const float* const *const)C, (const float* const *const)A22B21,
                                C_f_row + min_dim, C_f_col,
                                C_f_row + min_dim, C_f_col,
                                0, 0,
                                m, min_dim);
            deallocate_matrix(A22B21,m);
        }
        
    }
    if(min_dim == Col_A && Col_A != Row_A && Col_A != Col_B){
        //C22 == A21 * B12 if min_dim == Col_A
        rectangular_strassen_aux(C, A, B,
                            C_f_row + min_dim, C_f_col + min_dim,
                            A_f_row + min_dim, A_f_col,
                            B_f_row, B_f_col + min_dim,
                            m, min_dim, n);
    }
}
void rectangular_strassen_matrix_multiplication(float **C, float const *const *const A,
                                                float const *const *const B, size_t Row_A, size_t Col_A, size_t Col_B) 
{
    if (Row_A == Col_A && Col_A == Col_B){
        general_strassen_matrix_multiplication(C, A, B, Col_B);
        return;
    }
    rectangular_strassen_aux(C, A, B,
                        0, 0,
                        0, 0,
                        0, 0,
                        Row_A, Col_A, Col_B);
}
