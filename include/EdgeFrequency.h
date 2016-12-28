#ifndef EdgeFrequency_H
#define EdgeFrequency_H
#include "Graph.h"

#pragma offload_attribute (push, target(mic)) 

class EdgeFrequency  
{  
    int *array;  
    int u, v;  
  
public:  
    void init(int max_node_label, int max_edge_label)  
    {  
        u = max_node_label + 1;  
        v = u * (max_edge_label + 1);  
        array = new int[u * v];  
    }  
  
    int& operator()(int x, int a, int y) { return array[x * v + a * u + y]; }  
  
    int operator()(int x, int a, int y) const { return array[x * v + a * u + y]; }  
}; 
#pragma offload_attribute (pop)

#endif