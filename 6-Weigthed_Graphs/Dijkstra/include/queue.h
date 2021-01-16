#ifndef __QUEUE__

#include "graph.h"

typedef struct Queue
{
    Node ** Array;
    unsigned int size;

}Queue;

Queue build_queue(Node* nodes, unsigned int s); 

Node* extract_min_queue(Queue *Q);

int is_empty_queue(Queue* Q);

void delete_queue(Queue Q);

void DeQueue(Queue* Q, int min_idx);

#endif