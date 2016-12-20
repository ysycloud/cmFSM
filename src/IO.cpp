#include "IO.h"

void Usage()
{
	fprintf(stderr, "%s\n", USAGE);
}
 
void parse_params(int argc, char **argv, int my_rank, char *input, char *output, float &min_Support_Rate, int &division_way, int &thread_num)
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
			
		default:
			// Cannot parse.
			if(my_rank==0)
				Usage();
			MPI_Finalize();
			exit(0);
		}		
	}
  
    /* read graph data to check*/  
    FILE *fp;
	if((fp=fopen(input,"r"))==NULL)
	{
		if(my_rank==0)
			fprintf(stderr, "[ param error : -i ] can not open input '%s' file\n", input);
		MPI_Finalize();
		exit(0);
	}
	fclose(fp);
	
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
	
	if(strcmp(output, "unset")==0)
	{
		if(my_rank==0)
			fprintf(stderr, "[ param error : -o ] please input the output file path!\n");
		MPI_Finalize();
		exit(0);
	}
	
}

void load_data()
{
}

void write_results()
{

}
