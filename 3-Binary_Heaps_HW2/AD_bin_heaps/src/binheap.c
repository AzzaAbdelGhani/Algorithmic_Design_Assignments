#include <../include/binheap.h>
#include <string.h>
#include <stdio.h>

#define PARENT(node) ((node-1)/2)
#define LEFT_CHILD(node) (2*(node)+1)
#define RIGHT_CHILD(node) (2*(node+1))
#define VALID_NODE(H, node) ((H)->num_of_elem>(node))
#define ADDR(H, node) ((H)->A + (node)*(H)->key_size)
#define INDEX_OF(H, addr) (((addr)-((H)->A))/(H)->key_size)

//here define the ADDR of key_pos and rev_pos 
#define ADDR_KEY(H, idx) (&(H->key_pos[idx]))
#define ADDR_REV(H, idx) (&(H->rev_pos[idx]))

int is_heap_empty(const binheap_type *H)
{
    return H->num_of_elem == 0;
}

const void *min_value(const binheap_type *H)
{
    if (is_heap_empty(H))
    {
        return NULL;
    }

    // the minimum is stored in the root a.k.a A[0]
    // key of the root a.k.a A[key_pos[0]]
    return ADDR(H, *(H->key_pos));
}

//keep rev_pos updated
void update_rev_pos(binheap_type *H)
{
    unsigned int j;
    for (int  i = 0; i < H->num_of_elem ; i++)
    {
        j = H->key_pos[i];
        H->rev_pos[j] = i; 
    }
}
void swap_keys(binheap_type *H, unsigned int pos_a, unsigned int pos_b)
{
    unsigned int *idx_a = ADDR_KEY(H, pos_a);
    unsigned int *idx_b = ADDR_KEY(H, pos_b);
    unsigned int *idx_tmp = malloc(sizeof(unsigned int));

    memcpy(idx_tmp, idx_a,   sizeof(unsigned int));
    memcpy(idx_a,   idx_b,   sizeof(unsigned int));
    memcpy(idx_b,   idx_tmp, sizeof(unsigned int));

    free(idx_tmp);

    update_rev_pos(H);
    
    
}

void heapify(binheap_type *H, unsigned int node)
{
    unsigned int dst_node = node, child;
    unsigned int dst_node_addr, node_addr, child_addr;

    do
    {
        node = dst_node;
        node_addr = *(ADDR_KEY(H,node));
        dst_node_addr = *(ADDR_KEY(H,dst_node));
        child = RIGHT_CHILD(node);
        child_addr = *(ADDR_KEY(H,child));

        if (VALID_NODE(H, child) && 
            H->leq(ADDR(H, child_addr), ADDR(H, dst_node_addr)))
        {
            dst_node = child;
            dst_node_addr = *(ADDR_KEY(H,dst_node));
        }
        
        child = LEFT_CHILD(node);
        child_addr = *(ADDR_KEY(H,child));

        if (VALID_NODE(H, child) && 
            H->leq(ADDR(H, child_addr), ADDR(H, dst_node_addr)))
        {
            dst_node = child;
        }

        if (dst_node != node) 
        {
            swap_keys(H, dst_node, node);
        }
    } while (dst_node != node);
    
}

const void *extract_min(binheap_type *H)
{
    if (is_heap_empty(H))
    {
        return NULL;
    }
    
    swap_keys(H, 0, H->num_of_elem-1);
    unsigned int min_idx = *(ADDR_KEY(H, (H->num_of_elem)-1));

    if (H->num_of_elem > 1)
    {
        H->num_of_elem--;
        heapify(H, 0);
    }
    else
    {
        H->num_of_elem--;
    }

    return ADDR(H, min_idx); 
}

const void *find_the_max(void *A, 
			             const unsigned int num_of_elem,
                         const size_t key_size, 
			             total_order_type leq)
{
    if (num_of_elem == 0)
    {
        return NULL;
    }
    
    const void *max_value = A;
    for (const void *addr = A+key_size; 
	     addr != A+num_of_elem*key_size; 
	    addr+=key_size)
    {
        if (!leq(addr, max_value))
        {
            max_value = addr;
        }
    }

    return max_value;
}

binheap_type *build_heap(void *A, 
                         const unsigned int num_of_elem,
                         const unsigned int max_size,  
                         const size_t key_size, 
                         total_order_type leq)
{
    binheap_type *H = (binheap_type *)malloc(sizeof(binheap_type));
    
    H->A = A;
    H->num_of_elem = num_of_elem;
    H->max_size = max_size;
    H->key_size = key_size;
    H->leq = leq;
    H->max_order_value = malloc(key_size);
    H->key_pos = (unsigned int *)malloc(sizeof(unsigned int)*max_size);
    H->rev_pos = (unsigned int *)malloc(sizeof(unsigned int)*max_size);

    for (unsigned int i = 0; i < num_of_elem; i++)
    {
        H->key_pos[i] = i;
    }
    
    if (num_of_elem == 0)
    {
        return H;
    }

    const void *value = find_the_max(A, num_of_elem,
                                     key_size, leq);

    memcpy(H->max_order_value, value, key_size);
    
    for (unsigned int i = num_of_elem/2; i > 0; i--)
    {
        heapify(H, i);
    }

    heapify(H, 0);
    update_rev_pos(H);

    return H;
}

void delete_heap(binheap_type *H)
{
    free(H->max_order_value);
    free(H->key_pos);
    free(H->rev_pos);
    free(H);
}

const void *decrease_key(binheap_type *H, void *node, const void *value)
{
    unsigned int node_idx = INDEX_OF(H, node);
    unsigned int *node_addr = ADDR_REV(H,node_idx);
    unsigned int  pos_idx = *(node_addr);

    if (!VALID_NODE(H, node_idx) || !(H->leq(value, node)))
    {
        return NULL;
    }

    memcpy(node, value, H->key_size);
    
    if (H->num_of_elem > 1)
    {
        unsigned int parent_idx = PARENT(pos_idx);
        unsigned int *parent_addr = ADDR_KEY(H, parent_idx);
        void *parent = ADDR(H, *(parent_addr));
    
    
        while (pos_idx != 0 && !H->leq(parent, node))
        {
            swap_keys(H, pos_idx, parent_idx);
            
            node = ADDR(H,*(parent_addr));
            pos_idx = parent_idx;
            parent_idx = PARENT(pos_idx);

            if (pos_idx != 0)
            {
                parent_addr = ADDR_KEY(H,parent_idx);
                parent = ADDR(H, *(parent_addr));
            }
            
        }
        update_rev_pos(H);
    }
    
    return node;
}

const void *insert_value(binheap_type *H, const void *value)
{
    if (H->max_size == H->num_of_elem)
    {
        return NULL;
    }
    
    if (H->num_of_elem == 0 || !H->leq(value, H->max_order_value))
    {
        memcpy(H->max_order_value, value, H->key_size);
    }

    void *new_node_addr = ADDR(H, H->num_of_elem);
    unsigned int *new_node_pos_key = ADDR_KEY(H, H->num_of_elem);
    unsigned int idx_pos_key = INDEX_OF(H,new_node_addr); 
    memcpy(new_node_addr, H->max_order_value, H->key_size);
    memcpy(new_node_pos_key, &idx_pos_key, sizeof(unsigned int));

    H->rev_pos[H->num_of_elem] = H->num_of_elem;
    H->num_of_elem++;
    
    return decrease_key(H, new_node_addr, value);
}

void print_heap(const binheap_type *H, 
                void (*key_printer)(const void *value))
{
    unsigned int next_level_node = 1;

    for (unsigned int  node = 0; node < H->num_of_elem; node++)
    {
        if (node == next_level_node)
        {
            printf("\n");
            next_level_node = LEFT_CHILD(node);
        }
        else
        {
            printf("\t");
        }
        unsigned int *key_pos = ADDR_KEY(H,node);
        key_printer(ADDR(H, *(key_pos)));
    }
    printf("\n");
}
