#include "Global.h"

const int LABEL_MAX = 1000;
const int FILE_NAME_MAX = 100;
int min_support;  
int nr_graph; 
struct Graph *GS;     
EdgeFrequency EF;  
vector<Graph *> S;  // graph mining result  