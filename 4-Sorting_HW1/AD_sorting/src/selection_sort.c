#include "selection_sort.h"
#include "swap.h"

#define ADDR(A, index, key_size) (A+ (index) * (key_size))

void selection_sort(void *A, const unsigned int n, 
                    const size_t elem_size, 
                    total_order leq)
{
    // Implement selection sort here
    unsigned int max_j;
    for (unsigned int i = n-1; i > 0; i--)
    {
        max_j = 1;
        for (unsigned int j = 2; j <= i; j++)
        {
            if (!leq(ADDR(A,j,elem_size), ADDR(A,max_j,elem_size)))
            {
                max_j = j;
            }
            
        }

        swap(ADDR(A,max_j,elem_size), ADDR(A,i,elem_size), elem_size);
    }
    
}