#include "stdlib.h"
#include "stdio.h"
#include "../include/graph.h"

//This function to create new graph given the nodes , the weights matrix , and the size of the graph 
Graph build_graph(Node* nodes, unsigned int** weights, unsigned int s)
{
    Graph new_graph ; 
    new_graph.V = nodes;
    new_graph.W = weights;
    new_graph.size = s; 
    return new_graph; 
}

// This function is to retrieve a node by its index 
Node* get_node(Graph* G, unsigned int node_idx)
{
    for(size_t i = 0; i<G->size; i++)
    {
        if (G->V[i].idx == node_idx)
        {
            return &(G->V)[i];
        }
        
    }
    return NULL; 
}

// This function returns the weight of the edge between node a and node b
int get_weight(Graph* G, Node* a, Node* b)
{
    return G->W[a->idx][b->idx];
}

int get_num_adjacents(Graph* G, Node* N)
{
    int count = 0; 
    for (size_t i = 0; i < G->size; i++)
    {
        if (N->idx == G->V[i].idx) continue;
        count+=(get_weight(G,N,&(G->V)[i])) < INFTY;
        
    }
    return count;

}

// This functions returns the adjacent nodes for node N 
Node** get_adjacents(Graph* G, Node* N, int size)
{
    Node** tmp = (Node**)malloc(sizeof(Node*)*size);
    int index = 0;
    for (size_t i = 0; i < G->size; i++)
    {
        if (get_weight(G, N, &(G->V)[i])<INFTY && N->idx != G->V[i].idx)
        {
            tmp[index++] = &(G->V)[i]; 
        }
        
    }
    return tmp;
    
}