#ifndef IO_H
#define IO_H

#include "Global.h" 

#define ERRM "paraGSpan error:"

static const char *USAGE =
"\n"
"Usage:"
"  paraGSpan [options]\n"
"\n"
"  general options:\n"
"    -s --support: The minimal value of support rates\n"
"    -d --division: the division strategy among processes[default:4] \n"
"    	0: equality; 1: single; 2: increment; 3: circle; 4: dynamic. \n"
"    -t --thread: the number of threads in per process_num[default:1].\n"
"    -m --micthread: the number of threads in per mic[default:0 not use mic].\n"
"    -b --bindmic: the number of mics bind to one process[default:1].\n"
"  input/output options: \n"
"    -i --input: input file of graph set information. \n"
"    -o --output: the output file of frequent subgraph results. \n";
	
//static char * const UNSET = "unset";

void Usage();
void parse_params(int argc, char **argv, int my_rank,  /* input paras */
				char *input, char *output, /* output paras */
				float &min_Support_Rate, int &division_way, 
				int &thread_num, int &mic_thread, int &mic_num /* output paras */
				);
				
#pragma offload_attribute (push, target(mic))
void load_data(char *input, /* input paras */
			vector<GraphData *> &v_gd, /* output paras */
			int *freq_node_label, int *freq_edge_label  /* output paras */
			);
#pragma offload_attribute (pop)

void write_results(char *output, int my_rank);



#endif