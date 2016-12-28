#ifndef METHODS_H
#define METHODS_H

#include "Global.h"
#include "Traveler.h"
#include "SubTraveler.h"


#pragma offload_attribute (push, target(mic)) 

#define GET_TIME(now) { \
   struct timeval t; \
   gettimeofday(&t, NULL); \
   now = t.tv_sec + t.tv_usec/1000000.0; \
}

void split_data_equality(int size, int n, int my_rank, int* begin, int* end, int* local_n);
void split_data_increment(int size, int n, int my_rank, int* begin, int* end, int* local_n);
void split_data_single(int size, int n, int my_rank, int* begin, int* end, int* local_n);
void split_data_circle(int size, int n, int my_rank, int* index, int* local_n);
void pretreatment(int my_rank, int thread_num, const vector<GraphData *> &v_gd,  /* input paras */
				int *freq_node_label, int *freq_edge_label,  /* output paras */
				int *rank_to_node_label, int *rank_to_edge_label,  /* output paras */
				int &max_node_label, int &max_edge_label  /* output paras */
				);
/* next means the number of vertices in current gc */
void subgraph_mining(GraphCode &gc, int next);
void one_edge_expansion(GraphCode &gc, int next, vector<GraphCode> &child_gcs, vector<int> &nexts);
void freqGraphMining(GraphCode &gc, int next);
void paraFreqGraphMining(GraphCode &gc, int next, int thread_num);
void singleEdgeGraphMining(const Edge &e, vector<Edge> &single_edge_graph, int thread_num, int begin, int end);

#pragma offload_attribute (pop)

#endif