#ifndef METHODS_H
#define METHODS_H
#include <stdio.h>  
#include <stdlib.h>  
#include <string.h>  
#include <vector>
#include <set>
#include <map> 
#include <time.h>  
#include <sys/time.h>  
#include "Global.h"
#include "Traveler.h"
#include "SubTraveler.h"
#include <omp.h>

using namespace std;


#define GET_TIME(now) { \
   struct timeval t; \
   gettimeofday(&t, NULL); \
   now = t.tv_sec + t.tv_usec/1000000.0; \
}

void split_data_equality(int size, int n, int my_rank, int* begin, int* end, int* local_n);
void split_data_increment(int size, int n, int my_rank, int* begin, int* end, int* local_n);
void split_data_single(int size, int n, int my_rank, int* begin, int* end, int* local_n);
void split_data_circle(int size, int n, int my_rank, int* index, int* local_n);
void pretreatment(int my_rank, int thread_num, const vector<GraphData *> &v_gd, int *freq_node_label, int *freq_edge_label, int *rank_to_node_label, int *rank_to_edge_label, int &max_node_label, int &max_edge_label);
void subgraph_mining(GraphCode &gc, int next);
void one_edge_expansion(GraphCode &gc, int next, vector<GraphCode> &child_gcs, vector<int> &nexts);
void freqGraphMining(GraphCode &gc, int next);
void paraFreqGraphMining(GraphCode &gc, int next, int thread_num);
void singleEdgeGraphMining(const Edge &e, int thread_num);

#endif