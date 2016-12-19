#ifndef IO_H
#define IO_H
#include <stdio.h>  
#include <stdlib.h> 
#include <unistd.h> 
#include <getopt.h>
#include <assert.h>
#include "Global.h" 
#include "mpi.h"

using namespace std;
 
const char *USAGE =
"\n"
"Usage:"
"  paraGSpan [options]\n"
"\n"
"  general options:\n"
"    -s --support: The minimal value of support rates\n"
"    -d --division: the division strategy among processes[default:2] \n"
"    	0: equality; 1: single; 2: increment; 3: circle; 4:dynamic. \n"
"    -t --thread: the number of threads in per process_num[default:1].\n"
"  input/output options: \n"
"    -i --input: input file of graph set information. \n"
"    -o --output: the output file of frequent subgraph results. \n";

void Usage()
{
	fprintf(stderr, "%s\n", USAGE);
}

void parse_params(int argc, char **argv, int my_rank, char *input, char *output, float &min_Support_Rate, int &division_way, int &thread_num);
void load_data();
void write_results();

#endif