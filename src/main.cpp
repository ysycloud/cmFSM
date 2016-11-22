#define _CRT_SECURE_NO_WARNINGS 1  
  
#include <stdio.h>  
#include <stdlib.h>  
#include <string.h>  
#include <assert.h> 
#include <vector>  
#include <map>  
#include <set>
#include <unistd.h> 
#include <getopt.h>
#include "Global.h"
#include "Methods.h"
#include "mpi.h"
   
using namespace std;  

#define ERRM "gSpan error:"

const int LABEL_MAX = 1000;
const char *USAGE =
"\n"
"Usage:"
"  paraGSpan [options]\n"
"\n"
"  general options:\n"
"    -s --support: The minimal value of support rates\n"
"    -d --division: the division strategy among processes[default:2] \n"
"    	0: equality; 1: single; 2: increment; 3: circle. \n"
"    -t --thread: the number of threads in per process_num[default:1].\n"
"  input/output options: \n"
"    -i --input: input file of graph set information. \n"
"    -o --output: the output file of frequent subgraph results. \n";

void Usage()
{
	fprintf(stderr, "%s\n", USAGE);
}
  
  
int main(int argc, char **argv)  
{  

	int	my_rank;   /* My process rank           */
    int	p;         /* The number of processes   */
    int tag = 0;
	int parameternum;
	double start,finish,duration;
	
	/* Let the system do what it needs to start up MPI */
    MPI_Init(&argc, &argv);

    /* Get my process rank */
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    /* Find out how many processes are being used */
    MPI_Comm_size(MPI_COMM_WORLD, &p);
	
	/* check parameter*/
	if(my_rank == 0)
	{
		parameternum = argc;
		if(parameternum == 1)
			Usage();		
	}
	MPI_Bcast(&parameternum, 1, MPI_INT, 0, MPI_COMM_WORLD);
	if(parameternum == 1)
	{
		MPI_Finalize();
		exit(0);
	}
	
    GET_TIME(start);
  
    /* parse command line options */
	// Unset flags (value -1).
	float min_Support_Rate = -1;
	int division_way = -1;
	int thread_num = -1;
    // Unset options (value 'UNSET').
	char * const UNSET = "unset";
    char * input = UNSET;
	char * output = UNSET;
	
	if (argc == 1) 
	{
		Usage();
		exit(0);
    }
	
	int c;
	while (1) {
		int option_index = 0;
		static struct option long_options[] = {
			{"input",           required_argument,        0, 'i'},
			{"output",          required_argument,        0, 'o'},
			{"support",         required_argument,        0, 's'},
			{"division",        required_argument,        0, 'd'},
			{"thread",        	required_argument,        0, 't'},
			{0, 0, 0, 0}
		};

		c = getopt_long(argc, argv, "i:o:s:d:t:",
            long_options, &option_index);
	
		if(c==-1)	break;
		
		switch (c) {
		
		case 0:
			// A flag was set. //
			break;

		case 'i':
			if (input == UNSET) 
			{
				input = optarg;
			}
			else 
			{
				if(my_rank==0)
				{
					fprintf(stderr, "%s --input set more than once\n", ERRM);
					Usage();
				}
				MPI_Finalize();	
				exit(0);
			}
			break;
		
		case 'o':
			if (output == UNSET) 
			{
				output = optarg;
			}
			else 
			{
				if(my_rank==0)
				{
					fprintf(stderr, "%s --output set more than once\n", ERRM);
					Usage();
				}
				MPI_Finalize();
				exit(0);
			}
			break;
		
		case 's':
			if (min_Support_Rate < 0) {
				min_Support_Rate = atof(optarg);
				if (min_Support_Rate < 0 || min_Support_Rate >1) {
					if(my_rank==0)
					{
						fprintf(stderr, "%s --support must be a float value between 0-1\n", ERRM);
						Usage();
					}
					MPI_Finalize();
					exit(0);
				}
			}
			else {
				if(my_rank==0)
				{
					fprintf(stderr,"%s --support set more than once\n", ERRM);
					Usage();
				}
				MPI_Finalize();
				exit(0);
			}
			break;
			
		case 'd':
			if (division_way < 0) {
				division_way = atof(optarg);
				if (division_way < 0 || division_way >3) {
					if(my_rank==0)
					{
						fprintf(stderr, "%d --division must be a integer value among 0,1,2,3\n", ERRM);
						Usage();
					}
					MPI_Finalize();
					exit(0);
				}
			}
			else {
				if(my_rank==0)
				{
					fprintf(stderr,"%s --division set more than once\n", ERRM);
					Usage();
				}
				MPI_Finalize();
				exit(0);
			}
			break;
			
		case 't':
			if (thread_num < 0) {
				thread_num = atof(optarg);
				if (thread_num < 0) {
					if(my_rank==0)
					{
						fprintf(stderr, "%d --thread number must be a positive integer value\n", ERRM);
						Usage();
					}
					MPI_Finalize();
					exit(0);
				}
			}
			else {
				if(my_rank==0)
				{
					fprintf(stderr,"%s --thread number set more than once\n", ERRM);
					Usage();
				}
				MPI_Finalize();
				exit(0);
			}
			break;
			
		default:
			// Cannot parse. //
			if(my_rank==0)
				Usage();
			MPI_Finalize();
			exit(0);
		}		
	}
  
    /* read graph data */  
    FILE *fp;
	if((fp=fopen(input,"r"))==NULL)
	{
		if(my_rank==0)
			fprintf(stderr, "[ param error : -i ] can not open input '%s' file\n", input);
		MPI_Finalize();
		exit(0);
	}
	
	if(min_Support_Rate==-1)
	{
		if(my_rank==0)
			fprintf(stderr, "[ param error : -s ] please input the minimal support rates!\n");
		MPI_Finalize();
		exit(0);
	}
	
	if(division_way==-1)
		division_way = 2;
	
	if(thread_num==-1)
		thread_num = 1;
	
	if(output == UNSET)
	{
		if(my_rank==0)
			fprintf(stderr, "[ param error : -o ] please input the output file path!\n");
		MPI_Finalize();
		exit(0);
	}
	
	
	if(my_rank==0)
	{
		GET_TIME(finish);
		duration = finish-start;	
		printf("loading file time: %f seconds\n", duration);
		GET_TIME(start);		
	}
		
    bool occ_node_label[LABEL_MAX + 1], occ_edge_label[LABEL_MAX + 1];  
    int freq_node_label[LABEL_MAX + 1], freq_edge_label[LABEL_MAX + 1];  
    memset(freq_node_label, 0, sizeof(freq_node_label));  
    memset(freq_edge_label, 0, sizeof(freq_edge_label));  
    GraphData *gd = NULL;  
    vector<GraphData *> v_gd;  
    while (1)  
    {  
        static char dummy[10];
		/* there is no more graphs(t)、vertices(v) and edges(e)*/
        if (fscanf(fp, "%s", dummy) <= 0)  
        {  
            if (gd)  
            {  
                v_gd.push_back(gd);  
                for (int i = 0; i <= LABEL_MAX; ++i)  
                {  
                    if (occ_node_label[i])  
                        freq_node_label[i]++;  
                    if (occ_edge_label[i])  
                        freq_edge_label[i]++;  
                }  
            }  
            break;  
        }  
		
		/* get new graph*/
        if (*dummy == 't')  
        {  
            int id;  
            fscanf(fp, "%s%d", dummy, &id);  
			// deal with the last graph 
            if (gd)  
            {  
				// add the last graph to v_gd
                v_gd.push_back(gd);
				//update the frequency of every vertices and edges according to their occurrences in last graph
                for (int i = 0; i <= LABEL_MAX; ++i)  
                {  
                    if (occ_node_label[i])  
                        freq_node_label[i]++;  
                    if (occ_edge_label[i])  
                        freq_edge_label[i]++;  
                }  
            }  
            if (id < 0) break;  
            assert(id == (int)v_gd.size());  
			
			//generate new graphdata and reset the occurrences array to deal with the next graph
            gd = new GraphData;  
            memset(occ_node_label, 0, sizeof(occ_node_label));  
            memset(occ_edge_label, 0, sizeof(occ_edge_label));  
        }  
		
		/*deal with the new vertice*/
        else if (*dummy == 'v')  
        {  
            int id, label;  
            fscanf(fp, "%d%d", &id, &label);  
            assert(id == (int)gd->nodel.size() && label <= LABEL_MAX);  
            gd->nodel.push_back(label);  
            gd->nodev.push_back(true);  
            occ_node_label[label] = true;  
        }  
		/*deal with the new edge*/
        else if (*dummy == 'e')  
        {  
            int x, y, label;  
            fscanf(fp, "%d%d%d", &x, &y, &label);  
            assert(x < (int)gd->nodel.size() && y < (int)gd->nodel.size() && label <= LABEL_MAX);  
            gd->edgex.push_back(x);  
            gd->edgey.push_back(y);  
            gd->edgel.push_back(label);  
            gd->edgev.push_back(true);  
            occ_edge_label[label] = true;  
        }  
        else  
            assert(0);  
    }  
    fclose(fp);  
  
    min_support = (int) (v_gd.size() * min_Support_Rate);
	//min_support	= v_gd.size() * 0.1;
    if (min_support < 1)  
        min_support = 1;
	
	if(my_rank==0)
		printf("%d graphs with minSup = %d\n", (int)v_gd.size(), min_support);  
  
    /* sort labels of vertices and edges in GS by their frequency */  
    int rank_node_label[LABEL_MAX + 1], rank_edge_label[LABEL_MAX + 1];  
    for (int i = 0; i <= LABEL_MAX; ++i)  
    {  
        rank_node_label[i] = i;  
        rank_edge_label[i] = i;  
    }  
    for (int i = LABEL_MAX; i > 0; --i)  
    {  
        for (int j = 0; j < i; ++j)  
        {  
            int tmp;  
            if (freq_node_label[rank_node_label[j]] < freq_node_label[rank_node_label[j + 1]])  
            {  
                tmp = rank_node_label[j];  
                rank_node_label[j] = rank_node_label[j + 1];  
                rank_node_label[j + 1] = tmp;  
            }  
            if (freq_edge_label[rank_edge_label[j]] < freq_edge_label[rank_edge_label[j + 1]])  
            {  
                tmp = rank_edge_label[j];  
                rank_edge_label[j] = rank_edge_label[j + 1];  
                rank_edge_label[j + 1] = tmp;  
            }  
        }  
    }  
  
    /* remove infrequent vertices and edges */  
    /* ralabel the remaining vertices and edges in descending frequency */  
    int max_node_label = 0, max_edge_label = 0;  
    for (int i = 0; i <= LABEL_MAX; ++i)  
    {  
        if (freq_node_label[rank_node_label[i]] >= min_support)  
            max_node_label = i;  
        if (freq_edge_label[rank_edge_label[i]] >= min_support)  
            max_edge_label = i;  
    }

//	if(my_rank==0)	
//		printf("remaining max vertex, edge label: %d, %d\n", max_node_label, max_edge_label);  
  
    int rank_to_node_label[LABEL_MAX + 1], rank_to_edge_label[LABEL_MAX + 1];  
    memcpy(rank_to_node_label, rank_node_label, sizeof(rank_to_node_label));  
    for (int i = 0; i <= LABEL_MAX; ++i)  
        rank_node_label[rank_to_node_label[i]] = i;  
    memcpy(rank_to_edge_label, rank_edge_label, sizeof(rank_to_edge_label));  
    for (int i = 0; i <= LABEL_MAX; ++i)  
        rank_edge_label[rank_to_edge_label[i]] = i;  
  
    for (size_t i = 0; i < v_gd.size(); ++i)  
    {  
        GraphData &gd = *v_gd[i];  
        for (size_t j = 0; j < gd.nodel.size(); ++j)  
        {  
            if (freq_node_label[gd.nodel[j]] < min_support)  
                gd.nodev[j] = false;  
            else  
                gd.nodel[j] = rank_node_label[gd.nodel[j]];  
        }  
        for (size_t j = 0; j < gd.edgel.size(); ++j)  
        {  
            if (!gd.nodev[gd.edgex[j]] || !gd.nodev[gd.edgey[j]])  
            {  
                gd.edgev[j] = false;  
                continue;  
            }  
            if (freq_edge_label[gd.edgel[j]] < min_support)  
                gd.edgev[j] = false;  
            else  
                gd.edgel[j] = rank_edge_label[gd.edgel[j]];  
        }  
  
        /* re-map vertex index */  
        map<int, int> m;  
        int cur = 0;  
        for (size_t j = 0; j < gd.nodel.size(); ++j)  
        {  
            if (!gd.nodev[j]) continue;  
            m[j] = cur++;  
        }  
        for (size_t j = 0; j < gd.edgel.size(); ++j)  
        {  
            if (!gd.edgev[j]) continue;  
            gd.edgex[j] = m[gd.edgex[j]];  
            gd.edgey[j] = m[gd.edgey[j]];  
        }  
    }  
  
    /* build graph set */  
    nr_graph = (int)v_gd.size();  
    GS = new Graph[nr_graph];  
    for (int i = 0; i < nr_graph; ++i)  
    {  
        Graph &g = GS[i];  
        GraphData &gd = *v_gd[i];  
        for (size_t j = 0; j < gd.nodel.size(); ++j)  
            if (gd.nodev[j])  
                g.node_label.push_back(gd.nodel[j]);  
        g.edge_next = new vector<int>[g.node_label.size()];  
        g.edge_label = new vector<int>[g.node_label.size()];  
        for (size_t j = 0; j < gd.edgel.size(); ++j)  
        {  
            if (!gd.edgev[j]) continue;  
            g.edge_label[gd.edgex[j]].push_back(gd.edgel[j]);  
            g.edge_label[gd.edgey[j]].push_back(gd.edgel[j]);  
            g.edge_next[gd.edgex[j]].push_back(gd.edgey[j]);  
            g.edge_next[gd.edgey[j]].push_back(gd.edgex[j]);  
        }  
    }  
  
    /* find all subgraphs with only one vertex for reference(only root process record once) */
	if(my_rank==0)
	{
		for (int i = 0; i <= max_node_label; ++i)  
		{  
			Graph *g = new Graph;  
			g->node_label.push_back(i);  
			for (int j = 0; j < nr_graph; ++j)  
			{  
				const Graph &G = GS[j];  
				for (size_t k = 0; k < G.node_label.size(); ++k)  
				{  
					if (G.node_label[k] == i)  
					{  
						g->gs.push_back(j);  
						break;  
					}  
				}  
			}  
			S.push_back(g);  
		}
		if(my_rank==0)	
			printf("single_vertex_graph_num: %d\n", S.size());  		
	}

  
    /* enumerate all frequent 1-edge graphs in GS */  
    EF.init(max_node_label, max_edge_label);
	
	#pragma omp parallel for num_threads(thread_num) schedule(dynamic)	
    for (int x = 0; x <= max_node_label; ++x)  
    {  
        for (int a = 0; a <= max_edge_label; ++a)  
        {  
            for (int y = x; y <= max_node_label; ++y)  
            {  
                int count = 0;  
                for (int i = 0; i < nr_graph; ++i)  
                    if (GS[i].hasEdge(x, a, y))  
                        count++;  
                EF(x, a, y) = count;  
                EF(y, a, x) = count;  
            }  
        }  
    } 
	
	if(my_rank==0)
	{
		GET_TIME(finish);
		duration = finish-start;	
		printf("prework for mining frequenct subgraph spend: %f seconds\n", duration);
		GET_TIME(start);
	}

	// use frequenct 1-edge graph to submining(edges are sorted, cannot parallel)
	vector<Edge> single_edge_graph;
	for (int x = 0; x <= max_node_label; ++x)
	{  
		for (int a = 0; a <= max_edge_label; ++a)  
		{  
			for (int y = x; y <= max_node_label; ++y)  
			{  
				if (EF(x, a, y) < min_support)  
					continue;  

				Edge e(0, 1, x, a, y);				
				single_edge_graph.push_back(e);
  
			}  
		}  
	} 
	if(my_rank==0)
	{
		printf("single_edge_graph_num: %d\n", single_edge_graph.size());
		switch (division_way)
		{
			case 0:
				printf("the division strategy among processes : equality\n");
				break;
			case 1:
				printf("the division strategy among processes : single\n");
				break;
			case 2:
				printf("the division strategy among processes : increment\n");
				break;
			case 3:
				printf("the division strategy among processes : circle\n");
				break;
			default:
				break;
		}
	}
	
	// using different division strategy to mine frequent subgraph among processes
	if(division_way>=0&&division_way<=2)
	{
		int begin,end,local_n;
		//divide the data
		if(division_way==0)
			split_data_equality(single_edge_graph.size(), p, my_rank, &begin, &end, &local_n);	
		else if(division_way==1)
			split_data_single(single_edge_graph.size(), p, my_rank, &begin, &end, &local_n);	
		else if(division_way==2)
			split_data_increment(single_edge_graph.size(), p, my_rank, &begin, &end, &local_n);	
		for(int i=begin; i<end ; i++)
		{
			//make sure won't over the vector limit when processes are too much
			if( i<single_edge_graph.size() )
			{
				GraphCode gc;
				Edge e = single_edge_graph[i];
				gc.seq.push_back(&e);
				for (int j = 0; j < nr_graph; ++j)  
					if (GS[j].hasEdge(e.x, e.a, e.y))  
						gc.gs.push_back(j);  
				
				//begin to mining frequent subgraph
				if(thread_num == 1)
					freqGraphMining(gc, 2);
				else
					paraFreqGraphMining(gc, 2, thread_num);
				
				// GS <- GS - e 
				for (int j = 0; j < nr_graph; j++)  
					GS[j].removeEdge(e.x, e.a, e.y); 
			}
		}
	}
	else
	{
		int* index = new int[single_edge_graph.size()/p+1];
		int local_n;
		//divide the data circle
		split_data_circle(single_edge_graph.size(), p, my_rank, index, &local_n);	
		for(int i=0; i<local_n ; i++)
		{
			//make sure won't over the vector limit when processes are too much
			if( index[i]<single_edge_graph.size() )
			{
				GraphCode gc;
				Edge e = single_edge_graph[index[i]];
				gc.seq.push_back(&e);
				for (int j = 0; j < nr_graph; ++j)  
					if (GS[j].hasEdge(e.x, e.a, e.y))  
						gc.gs.push_back(j);  
				
				//begin to mining frequent subgraph
				if(thread_num == 1)
					freqGraphMining(gc, 2);
				else
					paraFreqGraphMining(gc, 2, thread_num);
				
				// GS <- GS - e 
				for (int j = 0; j < nr_graph; j++)  
					GS[j].removeEdge(e.x, e.a, e.y); 
			}
		}
		delete index;
	}
	
	MPI_Barrier(MPI_COMM_WORLD);
	
	if(my_rank==0)
	{
		GET_TIME(finish);
		duration = finish-start;	
		printf("Mining frequenct subgraph spend: %f seconds\n", duration);
		GET_TIME(start);
	}
	
	int local_fgraph_number = (int)S.size();
	int gloal_fgraph_number;
	
	MPI_Reduce(&local_fgraph_number, &gloal_fgraph_number, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  
    /* output mining result */  
	if(my_rank==0)
		printf("Found %d frequent subgraphs\n", gloal_fgraph_number);  
  
	char Res[128];
	sprintf(Res,"%s_%d.txt",output,my_rank);
    fp = fopen(Res, "w");  
    assert(fp);  
    for (int i = 0; i < (int)S.size(); ++i)  
    {  
        const Graph &g = *S[i];  
		fprintf(fp, "graph number: %d\n", i);  
        fprintf(fp, "frequent number: %d\n", (int)g.gs.size());  
		fprintf(fp, "oringal dataset graph list where it appeared:\n");  
        for (size_t j = 0; j < g.gs.size(); ++j)  
            fprintf(fp, " %d", g.gs[j]); 
        fprintf(fp, "\n"); 
		fprintf(fp, "node list:\n"); 		
        for (size_t j = 0; j < g.node_label.size(); ++j)  
            fprintf(fp, "v %d %d\n", (int)j, rank_to_node_label[g.node_label[j]]);  
        if (g.node_label.size() < 2)  
        {  
            fprintf(fp, "\n");  
            continue;  
        } 
		
		fprintf(fp, "edge list:\n");	
        for (int j = 0; j < (int)g.node_label.size(); ++j)  
        {  
            for (size_t k = 0; k < g.edge_label[j].size(); ++k)  
                if (j < g.edge_next[j][k])  
                    fprintf(fp, "e %d %d %d\n", j, g.edge_next[j][k], rank_to_edge_label[g.edge_label[j][k]]);  
        }  
        fprintf(fp, "\n");  
    }  
    fclose(fp);  
  
	if(my_rank==0)
	{
		GET_TIME(finish);
		duration = finish-start;	
		printf("output mining results time: %f seconds\n", duration);
	}
  
    MPI_Finalize();
	return 0;  
}  