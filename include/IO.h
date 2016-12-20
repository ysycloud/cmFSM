#ifndef IO_H
#define IO_H
#include <stdio.h>  
#include <stdlib.h> 
#include <string.h> 
#include <unistd.h> 
#include <getopt.h>
#include <assert.h>
#include "Global.h" 
#include "mpi.h"

#define ERRM "cmFSM error:"

using namespace std;
 
static const char *USAGE =
"\n"
"Usage:"
"  cmFSM [options]\n"
"\n"
"  general options:\n"
"    -s --support: The minimal value of support rates\n"
"    -d --division: the division strategy among processes[default:2] \n"
"    	0: equality; 1: single; 2: increment; 3: circle; 4:dynamic. \n"
"    -t --thread: the number of threads in per process_num[default:1].\n"
"  input/output options: \n"
"    -i --input: input file of graph set information. \n"
"    -o --output: the output file of frequent subgraph results. \n";
	
//static char * const UNSET = "unset";

void Usage();
void parse_params(int argc, char **argv, int my_rank, char *input, char *output, float &min_Support_Rate, int &division_way, int &thread_num);
void load_data(char *input, vector<GraphData *> &v_gd, int *freq_node_label, int *freq_edge_label);
void write_results(char *output, int my_rank, int *rank_to_node_label, int *rank_to_edge_label);

#endif