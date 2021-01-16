#include <stdio.h>
#include <stdlib.h>

#include "matrix.h"

// A is a pointer to a pointer to a float (it is a matrix)
// the 3 const mean the raw pointer and the column pointer and the contnets must not change
void naive_matrix_multiplication(float **C, float const *const *const A,
                                float const *const *const B,
                                const size_t n) 
{
  // IMPLEMENT THE NAIVE ALGORITHM HERE
  for (size_t y = 0; y < n; y++)
  {
    for (size_t x = 0; x < n; x++)
    {
      float value = 0.0;
      for (size_t z = 0; z < n; z++)
      {
        value += A[y][z] * B[z][x]; 
      }
      C[y][x] = value;
    }
    
  }
  
}

// This function is the naive algorithm for recatngualr matrices 
void rectangular_naive_matrix_multiplication(float **C, float const *const *const A,
                                        float const *const *const B,
                                        size_t Row_A, size_t Col_A, size_t Col_B) 
{
  for (size_t y = 0; y < Row_A; y++){
    for (size_t x = 0; x < Col_B; x++){
      float value = 0.0;
      for (size_t z = 0; z < Col_A; z++){
        value += A[y][z]*B[z][x];
      }
      C[y][x] = value;
    }
  }
}

int same_matrix(float const *const *const A, float const *const *const B,
                const size_t rows, const size_t cols) {
  for (size_t i = 0; i < rows; i++) {
    for (size_t j = 0; j < cols; j++) {
      if (A[i][j] != B[i][j]) {
        return 0;
      }
    }
  }

  return 1;
}

float **allocate_matrix(const size_t rows, const size_t cols) {
  float **M = (float **)malloc(sizeof(float *) * rows);

  for (size_t i = 0; i < rows; i++) {
    M[i] = (float *)malloc(sizeof(float) * cols);
  }

  return M;
}

void deallocate_matrix(float **A, const size_t rows) {
  for (size_t i = 0; i < rows; i++) {
    free(A[i]);
  }

  free(A);
}

float **allocate_random_matrix(const size_t rows, const size_t cols) {
  
  float **A = allocate_matrix(rows, cols);
  
  srand(10);
  for (size_t i = 0; i < rows; i++) {
    for (size_t j = 0; j < cols; j++) {
      A[i][j] = (rand() - RAND_MAX / 2) % 5;
    }
  }

  return A;
}

