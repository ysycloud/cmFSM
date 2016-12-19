#ifndef Supervisor_H
#define Supervisor_H
#include "mpi.h"

using namespace std;
 
void supervisor(int size, int p);
int request_edge_id(int rank);

#endif