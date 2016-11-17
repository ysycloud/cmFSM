#include "Methods.h"

void split_data_equality(int size, int n, int my_rank, int* begin, int* end, int* local_n)
{
	*local_n = size / n;  
	int leave = size % n;
	if(my_rank < leave){
		(*local_n)++;   
		*begin = my_rank*(*local_n);
	}else{	
		*begin = my_rank*(*local_n) + leave; 
	}	 
	*end = *begin + *local_n;
}


void split_data_increment(int size, int n, int my_rank, int* begin, int* end, int* local_n)
{
	*begin = my_rank*(my_rank+1)/2;
	if(my_rank!=n-1)
		*end = (my_rank+1)*(my_rank+2)/2;
	else
		*end = size;
	*local_n = my_rank+1;
}

void split_data_single(int size, int n, int my_rank, int* begin, int* end, int* local_n)
{
	*begin = my_rank;
	if(my_rank!=n-1)
		*end = my_rank+1;
	else
		*end = size;
	*local_n = 1;
}

void subgraph_mining(GraphCode &gc, int next)
{  
    /* construct graph from DFS code */  
    Graph *g = new Graph;  
    vector<const Edge *> &s = gc.seq;  
    g->node_label.resize(next);  
    g->edge_label = new vector<int>[next];  
    g->edge_next = new vector<int>[next];  
    for (size_t i = 0; i < s.size(); ++i)  
    {  
        int ix = s[i]->ix, iy = s[i]->iy, a = s[i]->a;  
        g->node_label[ix] = s[i]->x;  
        g->node_label[iy] = s[i]->y;  
        g->edge_label[ix].push_back(a);  
        g->edge_label[iy].push_back(a);  
        g->edge_next[ix].push_back(iy);  
        g->edge_next[iy].push_back(ix);  
    }  
  
    /* minimum DFS code pruning stage (4) */  
    Traveler t(s, *g);  
    t.travel();  
    if (!t.isMin())  
    {  
        delete[] g->edge_label;  
        delete[] g->edge_next;  
        delete g;  
        return;  
    }  
  
    S.push_back(g);   //add a new result(frequent subgraph) 
  
    /* enumerate potential children with 1-edge growth */  
    map<Edge, vector<int> > m;  
    for (size_t i = 0; i < gc.gs.size(); ++i)  
    {  
        set<Edge> c;  
        SubTraveler st(s, GS[gc.gs[i]], c, min_support, EF);  
        st.travel(next);  
        for (set<Edge>::const_iterator j = c.begin(); j != c.end(); ++j)  
            m[*j].push_back(gc.gs[i]);    //get the map(edge:gs) for every child's frequency
    }  
  
    /* mining subgraph children */  
    for (map<Edge, vector<int> >::iterator i = m.begin(); i != m.end(); ++i)  
    {  
        const Edge *e = &i->first;  
        vector<int> &spp = i->second;   //gs for current extension
        if ((int)spp.size() < min_support)  
            continue;  
        GraphCode child_gc;  
        child_gc.gs.swap(spp);  
        child_gc.seq = s;  
        child_gc.seq.push_back(e);  //get the correct extend and carry out the next submining
        if (e->iy == next)  
            subgraph_mining(child_gc, next + 1);  
        else  
            subgraph_mining(child_gc, next);  
    }  
  
    g->gs.swap(gc.gs);  
}  