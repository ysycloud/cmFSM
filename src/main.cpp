#define _CRT_SECURE_NO_WARNINGS 1  

//#pragma offload_attribute (push, target(mic))
#include "IO.h"
#include "MicFunc.h"
#include "Supervisor.h"
//#pragma offload_attribute (pop)
   
using namespace std;  
 
int main(int argc, char **argv)  
{  

	int	my_rank;   /* My process rank           */
    int	p;         /* The number of processes   */
    int tag = 0;
	double start,finish,duration;
	
	/* Let the system do what it needs to start up MPI */
    MPI_Init(&argc, &argv);

    /* Get my process rank */
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    /* Find out how many processes are being used */
    MPI_Comm_size(MPI_COMM_WORLD, &p); 
    
	float min_Support_Rate;
	int division_way;
	int thread_num;
	int mic_thread;
	int mic_num;
    char input[FILE_NAME_MAX];
	char output[FILE_NAME_MAX];
	
	int *sig0,*sig1,*sig2;
	
	/* parse command line options */
	parse_params(argc, argv, my_rank, input, output, min_Support_Rate, division_way, thread_num, mic_thread, mic_num);
	
	GET_TIME(start);
	
	if(mic_thread !=0 )
	{
		for(int i=0; i<mic_num; i++)
		{
			if(i==0)
			{
				#pragma offload target(mic:0) \
					in(input:length(FILE_NAME_MAX)alloc_if(1) free_if(1)) \
					in(min_Support_Rate) in(mic_thread) signal(sig0)												
						prepareDataOnMIC(input, mic_thread, min_Support_Rate);
			}
			else if(i==1)
			{
				#pragma offload target(mic:1) \
					in(input:length(FILE_NAME_MAX)alloc_if(1) free_if(1)) \
					in(min_Support_Rate) in(mic_thread) signal(sig1)												
						prepareDataOnMIC(input, mic_thread, min_Support_Rate);
			}
			else
			{
				#pragma offload target(mic:2) \
					in(input:length(FILE_NAME_MAX)alloc_if(1) free_if(1)) \
					in(min_Support_Rate) in(mic_thread) signal(sig2)												
						prepareDataOnMIC(input, mic_thread, min_Support_Rate);
			}	
		}
	}

//	if(my_rank==0)
//		printf("input:%s\toutput:%s\t%f\t%d\t%d\n",input,output,min_Support_Rate,division_way,thread_num);
		
			    
    int freq_node_label[LABEL_MAX + 1], freq_edge_label[LABEL_MAX + 1];
	//we cannot memset the array in the function, 
	//because the array is transfered into the function as a copying pointer,
	//the sizeof will not work anymore
	memset(freq_node_label, 0, sizeof(freq_node_label));  
    memset(freq_edge_label, 0, sizeof(freq_edge_label));	
    vector<GraphData *> v_gd;
	/* load data from file */ 	
	load_data(input, v_gd, freq_node_label, freq_edge_label);
	
	if(my_rank==0)
	{
		GET_TIME(finish);
		duration = finish-start;	
		printf("loading file time: %f seconds\n", duration);
		GET_TIME(start);		
	}
  
    min_support = (int) (v_gd.size() * min_Support_Rate);
	//min_support = v_gd.size() * 0.1;
    if (min_support < 1)  
        min_support = 1;
	
	if(my_rank==0)
		printf("%d graphs with minSup = %d\n", (int)v_gd.size(), min_support);

	int max_node_label = 0, max_edge_label = 0;
	/******* pretreatment ******/
	/*******
	1.sort labels of vertices and edges in GS by their frequency
	2.remove infrequent vertices and edges
	3.ralabel the remaining vertices and edges in descending frequency
	4.build graph set
	5.find all subgraphs with only one vertex for reference(only root process record once)
	6.enumerate all frequent 1-edge graphs in GS
	********/
	pretreatment(my_rank, thread_num, v_gd, freq_node_label, freq_edge_label, max_node_label, max_edge_label);
	
	/* Synchronize asynchronous preparation tasks */
	#pragma offload_wait target(mic:0) if(mic_num>=1&&mic_thread!=0) wait(sig0)
	#pragma offload_wait target(mic:1) if(mic_num>=2&&mic_thread!=0) wait(sig1)
	#pragma offload_wait target(mic:2) if(mic_num>=3&&mic_thread!=0) wait(sig2)
  
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
			case 4:
				printf("the division strategy among processes : dynamic\n");
				break;
			default:
				break;
		}
	}
	
	int pre = 0; //record untill which single edge graph did not deleted from last iteration in GS 
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
				Edge e = single_edge_graph[i];
				if(mic_thread==0)
					singleEdgeGraphMining(e, single_edge_graph, thread_num, pre, i);
				else
					cmsingleEdgeGraphMining(e, single_edge_graph, thread_num, pre, i, mic_thread, mic_num);
			}
			pre = i;
		}
	}
	else if(division_way == 3)  //circle split
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
				Edge e = single_edge_graph[index[i]];
				if(mic_thread==0)
					singleEdgeGraphMining(e, single_edge_graph, thread_num, pre, index[i]); 
				else
					cmsingleEdgeGraphMining( e, single_edge_graph, thread_num, pre, index[i], mic_thread, mic_num);				
			}
			pre = index[i];
		}
		delete index;
	}
	else  //dynamic slipt
	{
		if( p==1 ) //only one process, no need to use supervisor
		{
			for(int i=0; i<single_edge_graph.size() ; i++)
			{
				Edge e = single_edge_graph[i];
				if(mic_thread==0)
					singleEdgeGraphMining(e, single_edge_graph, thread_num, pre, i); 
				else
					cmsingleEdgeGraphMining( e, single_edge_graph, thread_num, pre, i, mic_thread, mic_num);
				pre = i;				
			}		
		}
		else  //not only one process, startup supervisor in process0
		{
			if(my_rank==0)
			{
				printf("rank:0 ---> start up supervisor!!!\n");
				supervisor(single_edge_graph.size(),p);
			}			
			else
			{
				int pos = request_edge_id(my_rank); 
				while( pos != -1 )
				{
					Edge e = single_edge_graph[pos];
					if(mic_thread==0)
						singleEdgeGraphMining(e, single_edge_graph, thread_num, pre, pos); 
					else
						cmsingleEdgeGraphMining( e, single_edge_graph, thread_num, pre, pos, mic_thread, mic_num);
					pre = pos;
					pos = request_edge_id(my_rank);
				}
			}
		}
	}
	
	MPI_Barrier(MPI_COMM_WORLD);
	
	if(my_rank==0)
	{
		GET_TIME(finish);
		duration = finish-start;	
		printf("Mining frequenct subgraph spend: %f seconds\n", duration);
		GET_TIME(start);
	}

	//get the number of results on MIC
	int mic_size;
	int local_fgraph_number = (int)S.size();
	
	
	if(mic_thread!=0)
	{
		for(int i=0; i<mic_num; i++)
		{
			if(i==0)
			{
				#pragma offload target(mic:0) out(mic_size)											
					getResultsSizeOnMIC(mic_size);
				local_fgraph_number += mic_size;
			}
			else if(i==1)
			{
				#pragma offload target(mic:1) out(mic_size)											
					getResultsSizeOnMIC(mic_size);
				local_fgraph_number += mic_size;
			}
			else
			{
				#pragma offload target(mic:2) out(mic_size)											
					getResultsSizeOnMIC(mic_size);
				local_fgraph_number += mic_size;
			}	
		}
	}
	
	
	int gloal_fgraph_number;
	//reduce to get the whole results
	MPI_Reduce(&local_fgraph_number, &gloal_fgraph_number, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  
    /* output mining result */  
	if(my_rank==0)
		printf("Found %d frequent subgraphs\n", gloal_fgraph_number);  
  
	/*output mining results in MIC of every process*/
	if(mic_thread!=0)
	{
		for(int i=0; i<mic_num; i++)
		{
			if(i==0)
			{
				#pragma offload target(mic:0) in(output:length(FILE_NAME_MAX)alloc_if(1) free_if(1)) signal(sig0)											
					write_resultsOnMIC(output);
			}
			else if(i==1)
			{
				#pragma offload target(mic:1) in(output:length(FILE_NAME_MAX)alloc_if(1) free_if(1)) signal(sig1)											
					write_resultsOnMIC(output);
			}
			else
			{
				#pragma offload target(mic:2) in(output:length(FILE_NAME_MAX)alloc_if(1) free_if(1)) signal(sig2)											
					write_resultsOnMIC(output);
			}	
		}
	}

	/*output mining results in every process*/ 
	write_results(output, my_rank);  
	
	/* Synchronize asynchronous writing tasks */
	#pragma offload_wait target(mic:0) if(mic_num>=1&&mic_thread!=0) wait(sig0)
	#pragma offload_wait target(mic:1) if(mic_num>=2&&mic_thread!=0) wait(sig1)
	#pragma offload_wait target(mic:2) if(mic_num>=3&&mic_thread!=0) wait(sig2)
  
	if(my_rank==0)
	{
		GET_TIME(finish);
		duration = finish-start;	
		printf("output mining results time: %f seconds\n", duration);
	}
  
    MPI_Finalize();
	return 0;  
}  