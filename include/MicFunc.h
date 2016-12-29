#ifndef MICFUNC
#define MICFUNC

#include "Methods.h"
#include "IO.h"

//__attribute__((target(mic))) void funcheck(int i);



void cmsingleEdgeGraphMining(const Edge &e, vector<Edge> &single_edge_graph, int thread_num, int begin, int end, int mic_thread);
void cmsingleEdgeGraphMining_simulationOnCPU( const Edge &e, vector<Edge> &single_edge_graph, int thread_num, int begin, int end, int mic_thread);
void split_data_interval(int n, int num, int order, int* index, int& local_n);
void freqGraphMiningfrom2edgesOnCPU(vector<GraphCode> &two_edges_child_gcs, vector<int> &two_edges_nexts, int thread_num);
void freqGraphMiningfrom2edgesOnMIC_simulation(vector<GraphCode> &two_edges_child_gcs, vector<int> &two_edges_nexts, int mic_thread, vector<Graph *> &MIC_S);
void one_edge_expansion_mic_simulation(GraphCode &gc, int next, vector<GraphCode> &child_gcs, vector<int> &nexts, vector<Graph *> &MIC_S);
void createDataForOffload(vector<GraphCode> &two_edges_child_gcs_mic, vector<int> &two_edges_nexts_mic,   /* input paras */
				int *ix, int *iy, int *x, int *a, int *y,  /* split edge parameters(n+1)  */ 
				int *gs, /* combination of gss split by -1 */
				int *nexts /* nexts(n) */);


#pragma offload_attribute (push, target(mic)) 
void constructGraphSetOnMIC(int mic_thread, const vector<GraphData *> &v_gd,  /* input paras */
				int *freq_node_label, int *freq_edge_label  /* input paras */);		
void prepareDataOnMIC(char *input, int mic_thread, float min_Support_Rate);
void reconstructDataforMiningOnMIC(int *ix, int *iy, int *x, int *a, int *y,  /* input split edge parameters(offload_len+1)  */ 
				int *gs, int *nexts, /* input combination of gss split by -1(gs_len) and nexts (offload_len) */ 
				int offload_len, int gs_len, /* length of input parameters  */
				vector<GraphCode> &two_edges_child_gcs_mic, vector<int> &two_edges_nexts_mic /* output  */);
void freqGraphMiningfrom2edgesOnMIC(vector<GraphCode> &two_edges_child_gcs, vector<int> &two_edges_nexts, int mic_thread);
void getResultsSizeOnMIC(int &size);
void write_resultsOnMIC(char *output);

#pragma offload_attribute (pop)

#endif 