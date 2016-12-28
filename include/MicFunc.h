#ifndef MICFUNC
#define MICFUNC

#include "Methods.h"

//__attribute__((target(mic))) void funcheck(int i);


void cmsingleEdgeGraphMining(int my_rank, const Edge &e, vector<Edge> &single_edge_graph, int thread_num, int begin, int end, int mic_thread);
void freqGraphMiningfrom2edgesOnCPU(int my_rank, vector<GraphCode> two_edges_child_gcs, vector<int> two_edges_nexts, int thread_num);
void split_data_interval(int n, int num, int order, int* index, int& local_n);
void freqGraphMiningfrom2edgesOnMIC(int my_rank, vector<GraphCode> two_edges_child_gcs, vector<int> two_edges_nexts, int mic_thread, vector<Graph *> &MIC_S);
void one_edge_expansion_mic(GraphCode &gc, int next, vector<GraphCode> &child_gcs, vector<int> &nexts, vector<Graph *> &MIC_S);

#pragma offload_attribute (push, target(mic)) 

void constructGraphSetOnMIC(int mic_thread, const vector<GraphData *> &v_gd,  /* input paras */
				int *freq_node_label, int *freq_edge_label,  /* input paras */
				int &max_node_label, int &max_edge_label  /* output paras */);		
void write_resultsOnMIC(char *output);

#pragma offload_attribute (pop)

#endif 