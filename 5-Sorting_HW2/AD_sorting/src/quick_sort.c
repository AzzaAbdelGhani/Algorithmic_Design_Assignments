#include "quick_sort.h"
#include "swap.h"

#define ADDR(A, index, key_size) (A+ (index) * (key_size))

unsigned int PARTITION(void *A, unsigned int l, unsigned int r,
                        unsigned int p, const size_t elem_size,
                        total_order leq)
{
    swap(ADDR(A,l,elem_size), ADDR(A,p,elem_size), elem_size);
    p = l ;
    l++;

    while (l<=r)
    {
        if (!leq(ADDR(A,l,elem_size), ADDR(A,p,elem_size)))
        {
            swap(ADDR(A,l,elem_size), ADDR(A,r,elem_size), elem_size);
            r--;
        }
        else
        {
            l++;
        }  
    }

    swap(ADDR(A,p,elem_size), ADDR(A,r,elem_size), elem_size);
    return r; //return index of the pivot
}

void quick_sort_aux(void *A, unsigned int l, unsigned int r,
                    const size_t elem_size, total_order leq)
{
    while (l<r)
    {
        unsigned int P = PARTITION(A,l,r,l,elem_size,leq);
        quick_sort_aux(A,l,P-1,elem_size,leq);
        l = P+1;
    }
    
}

void quick_sort(void *A, const unsigned int n, 
                const size_t elem_size, 
                total_order leq)
{
    // Implement quick sort here
    quick_sort_aux(A, 1, n, elem_size, leq);
}