#ifndef __DIJKSTRA__


#include "queue.h"
#include "binheap.h"


void INIT_SSSP(Graph* G);

void RELAX(Node* first, Node* second, unsigned int weight);

void DIJKSTRA_array(Graph* G, unsigned int src_idx);

void DIJKSTRA_heap(Graph* G, unsigned int src_idx);

#endif