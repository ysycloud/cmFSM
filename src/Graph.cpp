#include "Graph.h"

void Graph::removeEdge(int x, int a, int y)  
{  
	for (size_t i = 0; i < node_label.size(); ++i)  
	{  
		int t;  
		if (node_label[i] == x)  
			t = y;  
		else if (node_label[i] == y)  
			t = x;  
		else  
			continue;  
		for (size_t j = 0; j < edge_next[i].size(); ++j)  
		{  
			if (edge_label[i][j] == a && node_label[edge_next[i][j]] == t)  
			{  
				/* remove edge */  
				edge_label[i][j] = edge_label[i].back();  //use the final element to cover
				edge_label[i].pop_back();  
				edge_next[i][j] = edge_next[i].back();  
				edge_next[i].pop_back();  
				j--;  
			}  
		}  
	}  
}   

bool Graph::hasEdge(int x, int a, int y)  
{  
	for (size_t i = 0; i < node_label.size(); ++i)  
	{  
		int t;  
		if (node_label[i] == x)  
			t = y;  
		else if (node_label[i] == y)  
			t = x;  
		else  
			continue;  
		for (size_t j = 0; j < edge_next[i].size(); ++j)  
			if (edge_label[i][j] == a && node_label[edge_next[i][j]] == t)  
				return true;  
	}  
	return false;  
}  