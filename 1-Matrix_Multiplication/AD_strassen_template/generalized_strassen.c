#include "matrix.h"
#include "generalized_strassen.h"

//In this file we have the generalized strassen algorithm --> for all square matrices of any dimension


void general_sub_matrix_blocks(float **C, float const *const *const A,
                  float const *const *const B, 
                  const size_t C_f_row, const size_t C_f_col,
                  const size_t A_f_row, const size_t A_f_col,
                  const size_t B_f_row, const size_t B_f_col,
                  const size_t rows, const size_t cols) 
{
    for (size_t y = 0; y < rows; y++){
        for (size_t x = 0; x < cols; x++){
            C[y+C_f_row][x+C_f_col] = 
                A[y+A_f_row][x+A_f_col] - B[y+B_f_row][x+B_f_col];
        }
    }
}


void general_sum_matrix_blocks(float **C, float const *const *const A,
                  float const *const *const B, 
                  const size_t C_f_row, const size_t C_f_col,
                  const size_t A_f_row, const size_t A_f_col,
                  const size_t B_f_row, const size_t B_f_col,
                  const size_t rows, const size_t cols) 
{
    for (size_t y = 0; y < rows; y++){ 
        for (size_t x = 0; x < cols; x++){
            C[y+C_f_row][x+C_f_col] = 
                A[y+A_f_row][x+A_f_col] + B[y+B_f_row][x+B_f_col];
        }
    }
}


void general_naive_aux(float **C, float const *const *const A,
                  float const *const *const B, 
                  const size_t C_f_row, const size_t C_f_col,
                  const size_t A_f_row, const size_t A_f_col,
                  const size_t B_f_row, const size_t B_f_col,
                  const size_t Row_A, const size_t Col_A,
                  const size_t Col_B) 
{
    for (size_t y = 0; y < Row_A; y++) 
    {
    for (size_t x = 0; x < Col_B; x++)
    {
      float value = 0.0;
      for (size_t z = 0; z < Col_A; z++)
      {
        value += A[y+A_f_row][z+A_f_col] * B[z+B_f_row][x+B_f_col]; 
      }
      C[y+C_f_row][x+C_f_col] = value;
    }
    
  }
}

void repairing(float **C, float const *const *const A,
                    float const *const *const B,
                    const size_t C_f_row, const size_t C_f_col,
                    const size_t A_f_row, const size_t A_f_col,
                    const size_t B_f_row, const size_t B_f_col,
                    const size_t n, const size_t m)
{   
    const size_t n1 = n-m;

    // C11 = C + A12 * B21 
    float **c11 = allocate_matrix(n1, n1);
    general_naive_aux(c11, A, B,
                      0, 0,
                      A_f_row, A_f_col + n1,
                      B_f_row + n1, B_f_col,
                      n1, m, n1);  
    
    general_sum_matrix_blocks(C, (const float* const *const)C, (const float* const *const)c11,
                        C_f_row, C_f_col,
                        C_f_row, C_f_col,
                        0, 0,
                        n1, n1);
    deallocate_matrix(c11, n1);
    
    // C12 = A11 * B12 + A12 * B22 
    float **c12 =  allocate_matrix(n1, m);
    float **c121 = allocate_matrix(n1, m);
    general_naive_aux(c12, A, B,
                    0, 0,
                    A_f_row, A_f_col,
                    B_f_row, B_f_col + n1,
                    n1, n1, m);   
    general_naive_aux(c121, A, B,
                    0, 0,
                    A_f_row, A_f_col + n1,
                    B_f_row + n-1, B_f_col + n1,
                    n1, m, m);  

    general_sum_matrix_blocks(C, (const float* const *const)c12, (const float* const *const)c121,
                        C_f_row, C_f_col + n1,
                        0, 0,
                        0, 0,
                        n1, m);     

    deallocate_matrix(c12, n1);
    deallocate_matrix(c121, n1); 
    
    //C21 = A21 * B11 + A22 * B21 
    float **c21 =  allocate_matrix(m, n1);
    float **c211 = allocate_matrix(m, n1);
    general_naive_aux(c21, A, B,
                    0, 0,
                    A_f_row + n1, A_f_col,
                    B_f_row, B_f_col,
                    m, n1, n1);   
                    
    general_naive_aux(c211, A, B,
                    0, 0,
                    A_f_row + n1, A_f_col + n1,
                    B_f_row + n1, B_f_col,
                    m, m, n1);  
    
    general_sum_matrix_blocks(C, (const float* const *const)c21, (const float* const *const)c211,
                        C_f_row + n1, C_f_col,
                        0, 0,
                        0, 0,
                        m, n1);     
    
    deallocate_matrix(c21, m);
    deallocate_matrix(c211, m);         
    
    //C22 = A21 * B12 + A22 * B22
    float **c22 =  allocate_matrix(m, m);
    float **c221 = allocate_matrix(m, m);
    general_naive_aux(c22, A, B,
                    0, 0,
                    A_f_row + n-1, A_f_col,
                    B_f_row, B_f_col + n1,
                    m, n1, m);
    general_naive_aux(c221, A, B,
                    0, 0,
                    A_f_row + n1, A_f_col + n1,
                    B_f_row + n1, B_f_col + n1,
                    m, m, m);   
    general_sum_matrix_blocks(C, (const float* const *const)c22, (const float* const *const)c221,
                        C_f_row + n1, C_f_col + n1,
                        0, 0,
                        0, 0,
                        m, m);    
    deallocate_matrix(c22, m);
    deallocate_matrix(c221, m);  
                                 
}


void general_strassen_aux(float **C, float const *const *const A,
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
    
    

    float ***S = (float ***)malloc(sizeof(float **) * 10); // it is an array of matrices
    for (size_t i = 0; i < 10; i++)
    {
        S[i] = allocate_matrix(n2,n2);
    }

    float ***P = (float ***)malloc(sizeof(float **) * 7);
    for (size_t i = 0; i < 7; i++)
    {
        P[i] = allocate_matrix(n2,n2);
    }

    //S1 = B12 - B22
    general_sub_matrix_blocks(S[0], B, B, 0, 0, B_f_row, B_f_col+n2,
                      B_f_row + n2, B_f_col + n2, n2, n2);

    // P1 = A11 * S1
    general_strassen_aux(P[0], A, (const float* const const*) S[0],
                 0, 0, A_f_row, A_f_col, 0, 0, n2);

    // S2 = A11 + A12
    general_sum_matrix_blocks(S[1], A, A,
                      0, 0,
                      A_f_row, A_f_col,
                      A_f_row, A_f_col + n2,
                      n2, n2);


    // P2 = S2 x B22
    general_strassen_aux(P[1], (const float* const *const) S[1], B,
                 0, 0,
                 0, 0,
                 B_f_row + n2, B_f_col + n2,
                 n2);

    // S3 = A21 + A22
    general_sum_matrix_blocks(S[2], A, A,
                      0, 0,
                      A_f_row + n2, A_f_col,
                      A_f_row + n2, A_f_col + n2,
                      n2, n2);

    // P3 = S3 x B11
    general_strassen_aux(P[2], (const float* const *const) S[2], B,
                 0, 0,
                 0, 0,
                 B_f_row, B_f_col,
                 n2);

    // S4 = B21 - B11
    general_sub_matrix_blocks(S[3], B, B,
                      0, 0,
                      B_f_row + n2, B_f_col,
                      B_f_row, B_f_col,
                      n2, n2);

    // P4 = A22 x S4
    general_strassen_aux(P[3], A, (const float* const *const) S[3],
                 0, 0,
                 A_f_row + n2, A_f_col + n2,
                 0, 0,
                 n2);

    // S5 = A11 + A22
    general_sum_matrix_blocks(S[4], A, A,
                      0, 0,
                      A_f_row, A_f_col,
                      A_f_row + n2, A_f_col + n2,
                      n2, n2);

    // S6 = B11 + B22
    general_sum_matrix_blocks(S[5], B, B,
                      0, 0,
                      B_f_row, B_f_col,
                      B_f_row + n2, B_f_col + n2,
                      n2, n2);

    // P5 = S5 x S6
    general_strassen_aux(P[4], (const float* const *const) S[4],
                 (const float* const *const) S[5],
                 0, 0,
                 0, 0,
                 0, 0,
                 n2);

    // S7 = A12 - A22
    general_sub_matrix_blocks(S[6], A, A,
                      0, 0,
                      A_f_row, A_f_col + n2,
                      A_f_row + n2, A_f_col + n2,
                      n2, n2);

    // S8 = B21 + B22
    general_sum_matrix_blocks(S[7], B, B,
                      0, 0,
                      B_f_row + n2, B_f_col,
                      B_f_row + n2, B_f_col + n2,
                      n2, n2);

    // P6 = S7 x S8
    general_strassen_aux(P[5], (const float* const *const) S[6],
                 (const float* const *const) S[7],
                 0, 0,
                 0, 0,
                 0, 0,
                 n2);

    // S9 = A11 - A21
    general_sub_matrix_blocks(S[8], A, A,
                      0, 0,
                      A_f_row, A_f_col,
                      A_f_row + n2, A_f_col,
                      n2, n2);

    // S10 = B11 + B12
    general_sum_matrix_blocks(S[9], B, B,
                      0, 0,
                      B_f_row, B_f_col,
                      B_f_row, B_f_col + n2,
                      n2, n2);

    // P7 = S9 x S10
    general_strassen_aux(P[6], (const float* const *const) S[8],
                 (const float* const *const) S[9],
                 0, 0,
                 0, 0,
                 0, 0,
                 n2);

    // C11 = P5 + P4 - P2 + P6
    general_sum_matrix_blocks(C, (const float* const *const) P[4],
                      (const float* const *const) P[3],
                      C_f_row, C_f_col,
                      0, 0,
                      0, 0,
                      n2, n2);
    general_sub_matrix_blocks(C, (const float* const *const) C,
                      (const float* const *const) P[1],
                      C_f_row, C_f_col,
                      C_f_row, C_f_col,
                      0, 0,
                      n2, n2);
    general_sum_matrix_blocks(C, (const float* const *const) C,
                      (const float* const *const) P[5],
                      C_f_row, C_f_col,
                      C_f_row, C_f_col,
                      0, 0,
                      n2, n2);

    // C12 = P1 + P2
    general_sum_matrix_blocks(C, (const float* const *const) P[0],
                      (const float* const *const) P[1],
                      C_f_row, C_f_col+n2,
                      0, 0,
                      0, 0,
                      n2, n2);

    // C21 = P3 + P4
    general_sum_matrix_blocks(C, (const float* const *const) P[2],
                      (const float* const *const) P[3],
                      C_f_row+n2, C_f_col,
                      0, 0,
                      0, 0,
                      n2, n2);

    // C22 = P5 + P1 - P3 - P7
    general_sum_matrix_blocks(C, (const float* const *const) P[4],
                      (const float* const *const) P[0],
                      C_f_row+n2, C_f_col+n2,
                      0, 0,
                      0, 0,
                      n2, n2);
    general_sub_matrix_blocks(C, (const float* const *const) C,
                      (const float* const *const) P[2],
                      C_f_row+n2, C_f_col+n2,
                      C_f_row+n2, C_f_col+n2,
                      0, 0,
                      n2, n2);
    general_sub_matrix_blocks(C, (const float* const *const) C,
                      (const float* const *const) P[6],
                      C_f_row+n2, C_f_col+n2,
                      C_f_row+n2, C_f_col+n2,
                      0, 0,
                      n2, n2);

    for (size_t i = 0; i < 10; i++) {
      deallocate_matrix(S[i], n2);
    }
    free(S);
    
    for (size_t i = 0; i < 7; i++) {
      deallocate_matrix(P[i], n2);
    }
    free(P); 
    // repair if n is odd 
    if(n%2 != 0){
        repairing(C,  A, B,
                    C_f_row, C_f_col,
                    A_f_row, A_f_col,
                    B_f_row, B_f_col,
                    n,1);
    } 
    
}


void general_strassen_matrix_multiplication(float **C, float const *const *const A,
                                    float const *const *const B, size_t n) 
{
 
    // IMPLEMENT THE Generalized STRASSEN'S ALGORITHM HERE
    general_strassen_aux(C, A, B,
                  0, 0,
                  0, 0,
                  0, 0,
                  n);

}

