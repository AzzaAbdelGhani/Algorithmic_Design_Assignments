#include "select.h"
#include "quick_sort.h"
#include "swap.h"

#define ADDR(A, index, key_size) (A+ (index) * (key_size))

pair_type three_way_PARTITION(void *A, unsigned int l, unsigned int r,
                              unsigned int p, const size_t elem_size,
                              total_order leq)
{
    swap(ADDR(A,l,elem_size), ADDR(A,p,elem_size), elem_size);
    p = l ;
    l++;
    pair_type val; 
    int count = 0;

    while (l<=r)
    {
        if (leq(ADDR(A,l,elem_size), ADDR(A,p,elem_size)) && leq(ADDR(A,p,elem_size), ADDR(A,l,elem_size)))
        {
            //if A[l] = A[p]
            p = l;
            l++;
            count++;
        }
        else if (leq(ADDR(A,p,elem_size), ADDR(A,l,elem_size)))
        {
            //if A[l] > A[p]
            swap(ADDR(A,l,elem_size), ADDR(A,r,elem_size), elem_size);
            r--;
        }
        else
        {
            //if A[l] < A[p]
            swap(ADDR(A,l,elem_size), ADDR(A,p-count,elem_size), elem_size);
            p = l;
            l++;
        }
    }

    swap(ADDR(A,p,elem_size), ADDR(A,r,elem_size), elem_size);
    val.first = r - count;
    val.second = r;
    return val; 

}

unsigned int select_pivot(void *A, const unsigned int n, 
                          const size_t elem_size, total_order leq)
{
    // This function must be implemented
    if(n <= 10)
    {
        quick_sort(A,n,elem_size,leq);
        return (n/2);
    }

    unsigned int chunks = n/5;
    unsigned int c_l, c_r; 
    for (unsigned int c = 0; c < chunks; c++)
    {
        c_l = c*5;
        //c_r = 5 + c*5;
		quick_sort(ADDR(A, c_l, elem_size), 5, elem_size, leq);
        swap(ADDR(A, c_l+2, elem_size),ADDR(A, c, elem_size),elem_size); 
    }
    
	return select_index(A, chunks, chunks/2, elem_size, leq);
}

unsigned int select_index(void *A, const unsigned int n, 
                    const unsigned int i, const size_t elem_size, 
                    total_order leq)
{
    if(n <= 10)
    {
		quick_sort(A, n, elem_size, leq);
		return i;		
	}
	
	unsigned int j = select_pivot(A, n, elem_size, leq);
	pair_type k = three_way_PARTITION(A, 0, n-1, j, elem_size, leq);
	
	if(i < k.first)
    {
		return select_index(A, k.first-1, i, elem_size, leq); 
	}
	else if (i > k.second)
    {
        return select_index(ADDR(A, k.second, elem_size), n - k.second -1 , i, elem_size, leq);
    }
    else
    {
		return i;
	}

    
}



void quick_sort_select_aux(void *A, size_t l, size_t r,
                           const size_t elem_size, 
                           total_order leq)
{
    while (l < r)
    {
        unsigned int pivot_idx = l + select_pivot(ADDR(A, l, elem_size), r - l, elem_size, leq);
        pair_type k = three_way_PARTITION(A, l, r - 1, pivot_idx, elem_size, leq);
        quick_sort_select_aux(A, l, k.first, elem_size, leq);
        l = k.second + 1;
    }
}

void quick_sort_select(void *A, const unsigned int n, 
                       const size_t elem_size, 
                       total_order leq)
{
   // This function must be implemented 
   quick_sort_select_aux(A, 0, n, elem_size, leq);
}
