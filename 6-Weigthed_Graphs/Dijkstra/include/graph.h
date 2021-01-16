#ifndef __GRAPH__

#define INFTY 999999

// Node structure 
typedef struct Node
{
    unsigned int idx; // node's index 
    unsigned int dist; // node's distance from the source
    struct Node *pred; // node's predecessor
    
}Node;

// Graph representation as an adjacency matrix of nodes
typedef struct Graph
{
    Node* V;           // vertices
    unsigned int** W;  // adjacency matrix 
    unsigned int size; // totla number of nodes in Graph 
}Graph;


Graph build_graph(Node* nodes, unsigned int** weights, unsigned int s);

Node* get_node(Graph* G, unsigned int node_idx);

int get_weight(Graph* G, Node* a, Node* b);

int get_num_adjacents(Graph* G, Node* N);

Node** get_adjacents(Graph* G, Node* N, int size);

#endif