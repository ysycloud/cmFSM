#include "Global.h"

const int LABEL_MAX = 1000;
const int FILE_NAME_MAX = 100;
int rank_to_node_label[1001];
int rank_to_edge_label[1001];
int min_support;  
int nr_graph; 
struct Graph *GS;     
EdgeFrequency EF;  
vector<Graph *> S;  // graph mining result  