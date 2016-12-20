#include "Supervisor.h"
 
void supervisor(int size, int p)
{
	MPI_Status  status;
	static int count=0;
	int current_rank;
	while(count < size+p-1)
	{
		if(count<size)  //allocate tasks
		{
			MPI_Recv(&current_rank , 1, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);			
			//printf("%d:%d\n",current_rank,count);
			MPI_Send(&count, 1, MPI_INT, current_rank, 0, MPI_COMM_WORLD);		
		}
		else  //tasks is allocated over,reply -1 to every process
		{
			int end=-1;
			MPI_Recv(&current_rank , 1, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);
			//printf("%d:-1\n",current_rank);
			MPI_Send(&end, 1, MPI_INT, current_rank, 0, MPI_COMM_WORLD);
		}
		count++;
	}
}

int request_edge_id(int rank)
{
	MPI_Status  status;
	int task_id;
	MPI_Send(&rank, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);  //request task
	MPI_Recv(&task_id , 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);  //receive task
	return task_id;
}
