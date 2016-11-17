#ifndef SubTraveler_H
#define SubTraveler_H
#include <assert.h> 
#include <vector>
#include <set>
#include "Graph.h"
#include "EdgeFrequency.h"

using namespace std;

/* traverse a subgraph and find its children in graph */  
class SubTraveler  
{  
    const vector<const Edge *> &s; 
    const Graph &g;  
    set<Edge> &c;
	const int &min_support; 	
	const EdgeFrequency &EF;  	
    int next;  
    vector<int> rm;  // rightmost  
    vector<int> s2g;  
    vector<int> g2s;  
    vector<vector<bool> > f;   
    void DFS(int p); 
  
public:  
    SubTraveler(const vector<const Edge *> &_s, const Graph &_g, set<Edge> &_c, const int &_min_support, const EdgeFrequency &_EF)  
        : s(_s), g(_g), c(_c), min_support(_min_support), EF(_EF){}    
    void travel(int _next);
};

#endif