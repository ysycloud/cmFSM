#ifndef GRAPH_H
#define GRAPH_H

#pragma offload_attribute (push, target(mic)) 

#include <stdio.h>  
#include <stdlib.h>  
#include <string.h>  
#include <assert.h> 
#include <vector>  
#include <map>  
#include <set>
#include <unistd.h> 
#include <getopt.h>
#include <vector>
#include <time.h>  
#include <sys/time.h>  
#include <omp.h>
#include "offload.h"

using namespace std;

class Graph
{  
	public:
		vector<int> node_label;  
		vector<int> *edge_next; //adjacent node ids of every node
		vector<int> *edge_label;  //adjacent edge labels of every node
		vector<int> gs;  //the id set of original graph in dataSet. Current graph will appear in them.
	
		void removeEdge(int x, int a, int y);  
		bool hasEdge(int x, int a, int y);    
};

struct GraphData  
{  
    vector<int> nodel;  
    vector<bool> nodev;  
  
    vector<int> edgex;  
    vector<int> edgey;  
    vector<int> edgel;  
    vector<bool> edgev;  
};

struct Edge  
{  
    int ix, iy;  
    int x, a, y;  
  
    Edge(int _ix, int _iy, int _x, int _a, int _y) : ix(_ix), iy(_iy), x(_x), a(_a), y(_y) {}  
  
    bool operator<(const Edge &e) const  
    {  
        if (ix > iy)  
        {  
            if (e.ix < e.iy)  
                return true;  
            if (iy < e.iy || (iy == e.iy && a < e.a))  
                return true;  
        }  
        else if (e.ix < e.iy)  
        {  
            if (ix > e.ix)  
                return true;  
            if (ix == e.ix)  
            {  
                if (x < e.x)  
                    return true;  
                if (x == e.x && (a < e.a || (a == e.a && y < e.y)))  
                    return true;  
            }  
        }  
        return false;  
    }  
};  
  
struct GraphCode  
{  
    vector<const Edge *> seq;  
    vector<int> gs;  
}; 

#pragma offload_attribute (pop)

#include "mpi.h"


#endif