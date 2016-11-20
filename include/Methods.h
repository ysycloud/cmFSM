#ifndef METHODS_H
#define METHODS_H
#include <stdio.h>  
#include <stdlib.h>  
#include <string.h>  
#include <vector>
#include <set>
#include <map>  
#include "Global.h"
#include "Traveler.h"
#include "SubTraveler.h"
#include <omp.h>

using namespace std;

void split_data_equality(int size, int n, int my_rank, int* begin, int* end, int* local_n);
void split_data_increment(int size, int n, int my_rank, int* begin, int* end, int* local_n);
void split_data_single(int size, int n, int my_rank, int* begin, int* end, int* local_n);
void split_data_circle(int size, int n, int my_rank, int* index, int* local_n);
void subgraph_mining(GraphCode &gc, int next);
void one_edge_expansion(GraphCode &gc, int next, vector<GraphCode> &child_gcs, vector<int> &nexts);
void freqGraphMining(GraphCode &gc, int next);

#endif