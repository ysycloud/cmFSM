#ifndef MICFUNC
#define MICFUNC

#include "Methods.h"

__attribute__((target(mic))) void funcheck(int i);

void cmsingleEdgeGraphMining(int my_rank, const Edge &e, vector<Edge> &single_edge_graph, int thread_num, int begin, int end, int mic_thread);
void freqGraphMiningfrom2edgesOnCPU(int my_rank, vector<GraphCode> two_edges_child_gcs, vector<int> two_edges_nexts, int thread_num);
void freqGraphMiningfrom2edgesOnMIC(int my_rank, vector<GraphCode> two_edges_child_gcs, vector<int> two_edges_nexts, int mic_thread, vector<Graph *> &MIC_S);
void one_edge_expansion_mic(GraphCode &gc, int next, vector<GraphCode> &child_gcs, vector<int> &nexts, vector<Graph *> &MIC_S);

#endif 