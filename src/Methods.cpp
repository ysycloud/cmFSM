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

void split_data_circle(int size, int n, int my_rank, int* index, int* local_n)
{
	if(size<=n)
	{
		*local_n = 1;
		index[0] = my_rank;
	}
	else
	{
		int step1 = 2*(n-my_rank)-1;
		int step2 = 2*my_rank+1;
		int count=0;
		index[count] = my_rank;
		while(index[count]<size)
		{
			count++;
			index[count] = index[count-1]+step1;
			if(index[count]>=size)
				break;
			else
			{
				count++;
				index[count] = index[count-1]+step2;
			}				
		}
		*local_n = count;
	}
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
  
//	g->gs = gc.gs;	  //final line’s swap operation will finish this work
	#pragma omp critical
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
        if (e->iy == next)   //vextex add, backward extension, continue to next level
            subgraph_mining(child_gc, next + 1);  
        else  //vextex invariant, forward extension, continue to next level
            subgraph_mining(child_gc, next);   
    }  
  
    g->gs.swap(gc.gs);  
} 

void one_edge_expansion(GraphCode &gc, int next, vector<GraphCode> &child_gcs, vector<int> &nexts)
{  
    /* construct graph from DFS code */  
    Graph *g = new Graph;  
    vector<const Edge *> &s = gc.seq;  
    g->node_label.resize(next);  
    g->edge_label = new vector<int>[next];  
    g->edge_next = new vector<int>[next];  
	
//	printf("current edge conditions:\n");
//	for(size_t i=0;i<s.size();i++)
//		printf("current__size:%d;ix:%d;iy:%d;a:%d\n",s.size(),s[i]->ix,s[i]->iy,s[i]->a);
	
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
	
//	g->gs = gc.gs;	  //final line’s swap operation will finish this work，this will be work because g is a point.
	#pragma omp critical
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
		//must new the edge to make sure the spare will not free after function return
        Edge* e = new Edge(i->first.ix,i->first.iy,i->first.x,i->first.a,i->first.y);
		//const Edge *e = &i->first;  
        vector<int> &spp = i->second;   //gs for current extension
        if ((int)spp.size() < min_support)  
            continue;  
        GraphCode child_gc;  
        child_gc.gs.swap(spp);  
        child_gc.seq = s;  
        child_gc.seq.push_back(e);  //get the correct extend and carry out the next subminingW
        
		child_gcs.push_back(child_gc);
		if (e->iy == next)   //vextex add, backward extension, continue to next level
			nexts.push_back(next+1);
        else  //vextex invariant, forward extension, continue to next level
            nexts.push_back(next);   
    } 
    g->gs.swap(gc.gs);  
}   

void freqGraphMining(GraphCode &gc, int next)
{
	vector<GraphCode> current_global_child_gcs;
	vector<int> current_global_nexts;
	vector<GraphCode> tmp_global_child_gcs;
	vector<int> tmp_global_nexts;
	vector<GraphCode> local_child_gcs;
	vector<int> local_nexts;
	
	//first expansion to get the children set
	one_edge_expansion(gc, next, current_global_child_gcs, current_global_nexts);
	
	int len =  current_global_child_gcs.size();
	
	//submining per level
	while(len != 0)
	{	
		for(size_t i=0; i<len; i++)
		{
			//carry out current child's one edge expansion
			one_edge_expansion(current_global_child_gcs[i], current_global_nexts[i], local_child_gcs, local_nexts);
				
			//add the current child's expansion to local result
			tmp_global_child_gcs.insert(tmp_global_child_gcs.end(), local_child_gcs.begin(), local_child_gcs.end());
			tmp_global_nexts.insert(tmp_global_nexts.end(), local_nexts.begin(), local_nexts.end());
		
			//clear current child's expansion to carry out the next child's expansion
			local_child_gcs.clear();
			local_nexts.clear();
		}
		
		//swap tmp global results and global results
		current_global_child_gcs.swap(tmp_global_child_gcs);
		current_global_nexts.swap(tmp_global_nexts);

		//free the space of edges(because the edge is new to generate)
		//for(size_t i=0; i<tmp_global_child_gcs.size(); i++)
		//	for(size_t j=0; j < tmp_global_child_gcs[j].seq.size(); j++)
		//		delete[] tmp_global_child_gcs[i].seq[j];
		
		//clear tmp global results to carry out the next level mining
		tmp_global_child_gcs.clear();
		tmp_global_nexts.clear();
			
		//get the size of global results to judge whether continue next level mining
		len = current_global_child_gcs.size();
	}	
}


void paraFreqGraphMining(GraphCode &gc, int next, int thread_num)
{
	vector<GraphCode> current_global_child_gcs;
	vector<int> current_global_nexts;
	vector<GraphCode> tmp_global_child_gcs;
	vector<int> tmp_global_nexts;
	
	//first expansion to get the children set
	one_edge_expansion(gc, next, current_global_child_gcs, current_global_nexts);
	
	int len =  current_global_child_gcs.size();
	
	//submining per level
	while(len != 0)
	{
		#pragma omp parallel num_threads(thread_num)
		{		
			
			vector<GraphCode> local_child_gcs;
			vector<int> local_nexts;
			vector<GraphCode> current_local_child_gcs;
			vector<int> current_local_nexts;
			
			/*
			int local_t;	//the data number of each thread must hand
			int begin_t,end_t;
			int threadID = omp_get_thread_num();
			// compute the local size, up boundary and down boundary for every thread in 
			split_data_equality(len, thread_num, threadID, &begin_t, &end_t, &local_t);
			for(size_t i=begin_t; i<end_t; i++)
			{
				//carry out current child's one edge expansion
				one_edge_expansion(current_global_child_gcs[i], current_global_nexts[i], local_child_gcs, local_nexts);
					
				//add the current child's expansion to local result
				current_local_child_gcs.insert(current_local_child_gcs.end(), local_child_gcs.begin(), local_child_gcs.end());
				current_local_nexts.insert(current_local_nexts.end(), local_nexts.begin(), local_nexts.end());
			
				//clear current child's expansion to carry out the next child's expansion
				local_child_gcs.clear();
				local_nexts.clear();
			}
			*/
					
			//#pragma omp for schedule(static,1)
			//#pragma omp for schedule(guided)
			#pragma omp for schedule(dynamic)
			for(size_t i=0; i<len; i++)
			{
				//carry out current child's one edge expansion
				one_edge_expansion(current_global_child_gcs[i], current_global_nexts[i], local_child_gcs, local_nexts);
					
				//add the current child's expansion to local result
				current_local_child_gcs.insert(current_local_child_gcs.end(), local_child_gcs.begin(), local_child_gcs.end());
				current_local_nexts.insert(current_local_nexts.end(), local_nexts.begin(), local_nexts.end());
			
				//clear current child's expansion to carry out the next child's expansion
				local_child_gcs.clear();
				local_nexts.clear();
			}
			
			
			//merge local results into tmp global results  
			#pragma omp critical
			{
				tmp_global_child_gcs.insert(tmp_global_child_gcs.end(),current_local_child_gcs.begin(),current_local_child_gcs.end());
				tmp_global_nexts.insert(tmp_global_nexts.end(),current_local_nexts.begin(),current_local_nexts.end());
			}
		}
		
		//swap tmp global results and global results
		current_global_child_gcs.swap(tmp_global_child_gcs);
		current_global_nexts.swap(tmp_global_nexts);

		//free the space of edges(because the edge is new to generate)
		//for(size_t i=0; i<tmp_global_child_gcs.size(); i++)
		//	for(size_t j=0; j < tmp_global_child_gcs[j].seq.size(); j++)
		//		delete tmp_global_child_gcs[i].seq[j];
		
		//clear tmp global results to carry out the next level mining
		tmp_global_child_gcs.clear();
		tmp_global_nexts.clear();
			
		//get the size of global results to judge whether continue next level mining
		len = current_global_child_gcs.size();
	}	
}