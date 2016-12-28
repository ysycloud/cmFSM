#ifndef Traveler_H
#define Traveler_H

#include "Graph.h"

/* traverse a graph and find its minimal DFS code */ 
#pragma offload_attribute (push, target(mic)) 
class Traveler  
{  
    const vector<const Edge *> &s;  
    const Graph &g;
    bool is_min;  
    vector<int> g2s;  
    vector<vector<bool> > f;   
    void DFS(vector<int> &v, int c, int next);  
  
public:  
    Traveler(const vector<const Edge *> &_s, const Graph &_g) : s(_s), g(_g), is_min(true) {}   
    bool isMin() const { return is_min; }   
    void travel();
};  
#pragma offload_attribute (pop)



#endif