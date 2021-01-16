#include "../include/queue.h"
#include "stdlib.h"


Queue build_queue(Node* nodes, unsigned int s)
{
    Queue Q; 
    Node** adj_mat = (Node**)malloc(sizeof(Node*)*s);
    for (size_t i = 0; i < s; i++)
    {
        adj_mat[i] = &(nodes[i]);
    }
    
    Q.Array = adj_mat;
    Q.size = s;
    return Q;

} 

void DeQueue(Queue* Q, int min_idx)
{
    Node *n = Q->Array[min_idx];
    Q->Array[min_idx] = Q->Array[Q->size-1];
    Q->Array[Q->size-1] = n;
    Q->size--;
    return;
}

Node* extract_min_queue(Queue *Q)
{
    int min_idx = 0 ;
    unsigned int min_value = Q->Array[0]->dist;
    for (size_t i = 0; i < Q->size; i++)
    {
        if (Q->Array[i]->dist < min_value)
        {
            min_value = Q->Array[i]->dist;
            min_idx = i;
        } 
    }
    Node* m = Q->Array[min_idx];
    DeQueue(Q, min_idx);
    return m ;
    
}

int is_empty_queue(Queue *Q)
{
    return(Q->size == 0); 
}

void delete_queue(Queue Q)
{
    free(Q.Array);
}