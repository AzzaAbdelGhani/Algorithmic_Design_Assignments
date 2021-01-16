#include <stdio.h>
#include "test.h"
#include "test_rectangular.h"
#include "matrix.h"
#include "strassen.h"
#include "generalized_strassen.h"
#include "rectangular_strassen.h"
#include "optimized_strassen.h"

int main(int argc, char *argv[]) {
 
    //First, compare Naive algorithm and Strassen Algo. 

    size_t n = 4000;
 
    float **A1 = allocate_random_matrix(n, n);
    float **B1 = allocate_random_matrix(n, n);
    float **C0 = allocate_matrix(n, n);
    float **C1 = allocate_matrix(n, n);


    printf("n\tStrassen's Alg.\tNaive Alg.\tSame result\n");
    for (size_t j = 1; j <= n; j *= 2) {

      printf("%ld\t", j);
      fflush(stdout);

      printf("%lf\t", test(strassen_matrix_multiplication, C1, A1, B1, j));
      fflush(stdout);
      printf("%lf\t", test(naive_matrix_multiplication, C0, A1, B1, j));
      fflush(stdout);
    
      printf("%d\n", same_matrix((float const *const *const)C0,
                                (float const *const *const)C1, j, j));
    }

    deallocate_matrix(A1, n);
    deallocate_matrix(B1, n);
    deallocate_matrix(C0, n);
    deallocate_matrix(C1, n);


  /////////////////////////////////////////////////////////////////////////////////
  // Secondly, compare generalized strassen (which is for all squared matrices) and Naive alg.
    n = 1800;
    float **A2 = allocate_random_matrix(n, n);
    float **B2 = allocate_random_matrix(n, n);
    float **C2 = allocate_matrix(n, n);
    float **C3 = allocate_matrix(n, n);
    printf("n\tGen_Stras Alg.\tnaive Alg.\tSame result\n");
    for (size_t j = 1; j <= n; j += 90) {

        printf("%ld\t", j);
        fflush(stdout);
        printf("%lf\t", test(general_strassen_matrix_multiplication, C3, A2, B2, j));
        fflush(stdout);
        printf("%lf\t", test(naive_matrix_multiplication, C2, A2, B2, j));
        fflush(stdout);
      
        printf("%d\n", same_matrix((float const *const *const)C2,
                                  (float const *const *const)C3, j, j));
    }
    deallocate_matrix(A2, n);
    deallocate_matrix(B2, n);
    deallocate_matrix(C2, n);
    deallocate_matrix(C3, n);


    ////////////////////////////////////////////////////////////////////////////
    // Thirdly, Compare rectangular strassen alg. and rectangular Naive alg.
    
    size_t m = 5000;
    size_t k = 5000;

    float **A3 = allocate_random_matrix(m, k);
    float **B3 = allocate_random_matrix(k, n);
    float **C4 = allocate_matrix(m, n);
    float **C5 = allocate_matrix(m, n);

    printf("m\tk\tn\tRect_Strassen's \tRect_Naive \tSame result\n");
    for (size_t i = 234; i <= m; i += 1500){
      for (size_t l = 785; l <= k; l += 1500){
        for (size_t j = 574; j <= n; j += 1500) {
          printf("%ld\t%ld\t%ld\t", i, l, j);
          fflush(stdout);
          printf("%lf\t", test_rectangular(rectangular_strassen_matrix_multiplication, C4, A3, B3, i, l, j));
          fflush(stdout);
          printf("%lf\t", test_rectangular(rectangular_naive_matrix_multiplication, C5, A3, B3, i, l, j));
          fflush(stdout);
          
          printf("%d\n", same_matrix((float const *const *const)C4,
                                      (float const *const *const)C5, i, j));
        }
      }
    }

    deallocate_matrix(A3, m);
    deallocate_matrix(B3, k);
    deallocate_matrix(C4, m);
    deallocate_matrix(C5, m);


    /////////////////////////////////////////////////////////////////////////////////
  // Fourthly, compare generalized strassen (which is for all squared matrices) and optimized strassen alg.
    n = 5000;
    float **A4 = allocate_random_matrix(n, n);
    float **B4 = allocate_random_matrix(n, n);
    float **C6 = allocate_matrix(n, n);
    float **C7 = allocate_matrix(n, n);
    printf("n\tGen_Stras \tOptimized stras.\tSame result\n");
    for (size_t j = 600; j <= n; j += 543) {

        printf("%ld\t", j);
        fflush(stdout);
        printf("%lf\t", test(general_strassen_matrix_multiplication, C7, A4, B4, j));
        fflush(stdout);
        printf("%lf\t", test(optimized_strassen_matrix_multiplication, C6, A4, B4, j));
        fflush(stdout);
      
        printf("%d\n", same_matrix((float const *const *const)C6,
                                  (float const *const *const)C7, j, j));
    }
    deallocate_matrix(A4, n);
    deallocate_matrix(B4, n);
    deallocate_matrix(C6, n);
    deallocate_matrix(C7, n);
    return 0;
}
