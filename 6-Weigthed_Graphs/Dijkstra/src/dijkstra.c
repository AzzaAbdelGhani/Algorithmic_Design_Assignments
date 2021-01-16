#include<stdlib.h>
#include<stdio.h>
#include<string.h>

#include "../include/dijkstra.h"

void INIT_SSSP(Graph* G)
{
    for(size_t i; i<G->size; i++)
    {
        G->V[i].idx = i;
        G->V[i].dist = INFTY; 
        G->V[i].pred = NULL;  
    }
}

void RELAX(Node* first, Node* second, unsigned int weight)
{
    if((first->dist + weight) < second->dist)
    {
        second->dist = first->dist + weight;
        second->pred = first; 
    }
}

int leq_dist(const void* x, const void* y)
{
    Node* s = (Node*)x;
    Node* d = (Node*)y;
    return leq_int((void*)&(s->dist), (void*)&(d->dist));
}


// Here is the implementation of Array-based version of Dijkstra's Algorithm 
void DIJKSTRA_array(Graph* G, unsigned int src_idx)
{
    INIT_SSSP(G);
    Node* src = get_node(G,src_idx); 
    src->dist = 0; 

    Queue q = build_queue(G->V, G->size);
    Queue* Q = &q;
    while (!(is_empty_queue(Q)))
    {
        Node* M = extract_min_queue(Q);
        int num = get_num_adjacents(G, M);
        Node** adj = get_adjacents(G, M, num);
        for(size_t i=0; i<num; i++)
        {
            RELAX(M, adj[i], get_weight(G, M, adj[i]));
        }
        
        free(adj);
    }
    delete_queue(q);
    
}


// Here is the implementation of Heap-based version of Dijkstra's Algorithm 
void DIJKSTRA_heap(Graph* G, unsigned int src_idx)
{
    INIT_SSSP(G);
    Node* src = get_node(G,src_idx); 
    src->dist = 0; 

    binheap_type* H = build_heap(G->V, G->size, G->size, sizeof(Node), leq_dist);

    while (!(is_heap_empty(H)))
    {
        Node* M = (Node*)extract_min(H);
        int num = get_num_adjacents(G, M);
        Node** adj = get_adjacents(G, M, num);
        for(size_t i=0; i<num; i++)
        {
            RELAX(M, adj[i], get_weight(G, M, adj[i]));
        }
        heapify(H, 0);
        free(adj);
    }
    delete_heap(H);
    
}

