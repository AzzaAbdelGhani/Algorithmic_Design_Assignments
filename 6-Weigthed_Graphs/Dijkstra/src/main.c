#include "stdlib.h"
#include "stdio.h"
#include "time.h"

#include "../include/dijkstra.h"

int main()
{
        
    struct timespec start, end; 
    unsigned int s = 30000;
    unsigned int** adj_mat = (unsigned int **)malloc(sizeof(unsigned int *)*s);  
    
    for(size_t i=0; i<s; i++)
    {
        adj_mat[i] = (unsigned int *)malloc(sizeof(unsigned int)*s);
    }

    for (size_t i = 0; i < s; i++)
    {
        for(size_t j=0; j< s; j++)
        {
          if((i+j)<2500)
            adj_mat[i][j]= i+j; 
          if((i+j)>300)
            adj_mat[j][i] = i+j+5;
        }
    }
    
    
    Node*nodes_1 = (Node*)malloc(sizeof(Node)*s);
    Node*nodes_2 = (Node*)malloc(sizeof(Node)*s);

    Graph g1 = build_graph(nodes_1, adj_mat, s);
    Graph g2 = build_graph(nodes_2, adj_mat, s);


    printf("\n\n          Time Performance for each version of Dijkstra's Algorithm             \n\n");
    
    printf("size\tArray\t\tHeap\n");
    for(size_t i=0; i<10; i++)
    {
        g1.size = s/(1<<(9-i));  
        g2.size = s/(1<<(9-i));

        printf("%d", g1.size);
        
        clock_gettime(CLOCK_REALTIME, &start);
        DIJKSTRA_array(&g1, 0);
        clock_gettime(CLOCK_REALTIME, &end);

        printf("\t%lf", (end.tv_sec-start.tv_sec) +
                        (end.tv_nsec-start.tv_nsec)/1E9 );

        clock_gettime(CLOCK_REALTIME, &start);
        DIJKSTRA_heap(&g2, 0);
        clock_gettime(CLOCK_REALTIME, &end);

        printf("\t%lf\n", (end.tv_sec-start.tv_sec) +
                        (end.tv_nsec-start.tv_nsec)/1E9 );
        
    }
    printf("\n");
    
    free(adj_mat);
    free(g1.V);
    free(g2.V);
    
    return 0;
}