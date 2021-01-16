#include "bubble_sort.h"
#include "swap.h"

#define ADDR(A, index, key_size) (A+ (index) * (key_size))

void bubble_sort(void *A, const unsigned int n, 
                 const size_t elem_size, 
                 total_order leq)
{
    // Implement bubble sort here
    for (unsigned int i = n-1; i > 0; i--)
    {
        for (unsigned int j = 1; j < i ; j++)
        {
            if (!leq(ADDR(A,j,elem_size), ADDR(A,j+1,elem_size)))
            {
                swap(ADDR(A,j,elem_size), ADDR(A,j+1,elem_size), elem_size);
            }
            
        }
        
    }
    
}