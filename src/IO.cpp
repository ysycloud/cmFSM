#include "IO.h"

void Usage()
{
	fprintf(stderr, "%s\n", USAGE);
}
 
void parse_params(int argc, char **argv, int my_rank,  /* input paras */
				char *input, char *output, /* output paras */
				float &min_Support_Rate, int &division_way, 
				int &thread_num, int &mic_thread, int &mic_num/* output paras */
				)
{
	
	/* check parameter*/
	if(my_rank == 0)
	{
		if(argc == 1)
			Usage();		
	}
	
	if(argc == 1)
	{
		MPI_Finalize();
		exit(0);
	}
	
	//initiate parameters
	// Unset flags (value -1).
	min_Support_Rate = -1;
	division_way = -1;
	thread_num = -1;
	mic_thread = -1;
	mic_num = -1;
    // Unset options (value 'UNSET').
    strcpy(input,"unset");
	strcpy(output,"unset");
	
//	if(my_rank==0)
//		printf("input:%s\toutput:%s\t%f\t%d\t%d\n",input,output,min_Support_Rate,division_way,thread_num);
	
	int c;
	while (1) {
		int option_index = 0;
		static struct option long_options[] = {
			{"input",           required_argument,        0, 'i'},
			{"output",          required_argument,        0, 'o'},
			{"support",         required_argument,        0, 's'},
			{"division",        required_argument,        0, 'd'},
			{"thread",        	required_argument,        0, 't'},
			{"micthread",       required_argument,        0, 'm'},
			{"bindmic",       	required_argument,        0, 'b'},
			{0, 0, 0, 0}
		};

		c = getopt_long(argc, argv, "i:o:s:d:t:m:b:",
            long_options, &option_index);
	
		if(c==-1)	break;
		
		switch (c) {
		
		case 0:
			// A flag was set. //
			break;

		case 'i':
			if (strcmp(input, "unset")==0) 
			{
				//input = (char *)malloc(sizeof(char)*(strlen(optarg)+1));
				strcpy(input,optarg);
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
			if (strcmp(output, "unset")==0) 
			{
				//output = (char *)malloc(sizeof(char)*(strlen(optarg)+1));
				strcpy(output,optarg);
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
				if (division_way < 0 || division_way >4) {
					if(my_rank==0)
					{
						fprintf(stderr, "%d --division must be a integer value among 0,1,2,3,4\n", ERRM);
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
		
		case 'm':
			if (mic_thread < 0) {
				mic_thread = atof(optarg);
				if (mic_thread < 0) {
					if(my_rank==0)
					{
						fprintf(stderr, "%d --mic thread number must not be a negative integer value\n", ERRM);
						Usage();
					}
					MPI_Finalize();
					exit(0);
				}
			}
			else {
				if(my_rank==0)
				{
					fprintf(stderr,"%s --mic thread number set more than once\n", ERRM);
					Usage();
				}
				MPI_Finalize();
				exit(0);
			}
			break;

		case 'b':
			if (mic_num < 0) {
				mic_num = atof(optarg);
				if (mic_num <= 0 || mic_num > 3) {
					if(my_rank==0)
					{
						fprintf(stderr, "%d --bindmic must be a integer value among 1,2,3\n", ERRM);
						Usage();
					}
					MPI_Finalize();
					exit(0);
				}
			}
			else {
				if(my_rank==0)
				{
					fprintf(stderr,"%s --bindmic set more than once\n", ERRM);
					Usage();
				}
				MPI_Finalize();
				exit(0);
			}
			break;
			
		default:
			// Cannot parse.
			if(my_rank==0)
				Usage();
			MPI_Finalize();
			exit(0);
		}		
	}
  
    /* read graph data to check whether input file is valid */  
    FILE *fp;
	if((fp=fopen(input,"r"))==NULL)
	{
		if(my_rank==0)
			fprintf(stderr, "[ param error : -i ] can not open input '%s' file\n", input);
		MPI_Finalize();
		exit(0);
	}
	fclose(fp);
	
	 /* check whether min_Support_Rate is valid */
	if(min_Support_Rate==-1)
	{
		if(my_rank==0)
			fprintf(stderr, "[ param error : -s ] please input the minimal support rates!\n");
		MPI_Finalize();
		exit(0);
	}
	
	/* check whether division strategy is valid */
	if(division_way==-1)
		division_way = 4;
	
	/* check whether thread_num is valid */
	if(thread_num==-1)
		thread_num = 1;
	
	/* check whether mic thread num is valid */
	if(mic_thread==-1)
		mic_thread = 0;
	
	/* check whether mic num in one process is valid */
	if(mic_num==-1)
		mic_num = 1;
	
	/* check whether output file is set */
	if(strcmp(output, "unset")==0)
	{
		if(my_rank==0)
			fprintf(stderr, "[ param error : -o ] please input the output file path!\n");
		MPI_Finalize();
		exit(0);
	}
	
	//	if(my_rank==0)
	//		printf("input:%s\toutput:%s\t%f\t%d\t%d\n",input,output,min_Support_Rate,division_way,thread_num);
	
}

void load_data(char *input, /* input paras */
			vector<GraphData *> &v_gd, /* output paras */
			int *freq_node_label, int *freq_edge_label  /* output paras */
			)
{
	
	bool occ_node_label[LABEL_MAX + 1], occ_edge_label[LABEL_MAX + 1];	  
    GraphData *gd = NULL; 
	
	FILE *fp;
	fp=fopen(input,"r");
	
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
}

void write_results(char *output, int my_rank)
{
	FILE *fp;
	char Res[128];
	//merge the output path for every process
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
}
