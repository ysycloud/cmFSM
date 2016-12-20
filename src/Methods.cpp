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

void pretreatment(int my_rank, int thread_num, const vector<GraphData *> &v_gd,  /* input paras */
				int *freq_node_label, int *freq_edge_label,  /* output paras */
				int *rank_to_node_label, int *rank_to_edge_label,  /* output paras */
				int &max_node_label, int &max_edge_label  /* output paras */
				)
{
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
    for (int i = 0; i <= LABEL_MAX; ++i)  
    {  
        if (freq_node_label[rank_node_label[i]] >= min_support)  
            max_node_label = i;  
        if (freq_edge_label[rank_edge_label[i]] >= min_support)  
            max_edge_label = i;  
    }

//	if(my_rank==0)	
//		printf("remaining max vertex, edge label: %d, %d\n", max_node_label, max_edge_label);  
  
    memcpy(rank_to_node_label, rank_node_label, sizeof(rank_node_label));  
    for (int i = 0; i <= LABEL_MAX; ++i)  
        rank_node_label[rank_to_node_label[i]] = i;  
    memcpy(rank_to_edge_label, rank_edge_label, sizeof(rank_edge_label));  
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
		//there is no need to delete e in mining stage,  
		//because ‘child_gc.seq = s‘ lead every child will use parent seqs.
		//if delete before program end, must error(delete one pointer more than once)
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
		
		//clear tmp global results to carry out the next level mining
		tmp_global_child_gcs.clear();
		tmp_global_nexts.clear();
			
		//get the size of global results to judge whether continue next level mining
		len = current_global_child_gcs.size();
	}	
}

void singleEdgeGraphMining(const Edge &e, int thread_num)
{
	GraphCode gc;
	gc.seq.push_back(&e);
	for (int j = 0; j < nr_graph; ++j)  
		if (GS[j].hasEdge(e.x, e.a, e.y))  
			gc.gs.push_back(j);  
				
		//begin to mining frequent subgraph
		//subgraph_mining(gc, 2);
		if(thread_num == 1)
			freqGraphMining(gc, 2);
		else
			paraFreqGraphMining(gc, 2, thread_num);
				
		// GS <- GS - e 
		for (int j = 0; j < nr_graph; j++)  
			GS[j].removeEdge(e.x, e.a, e.y); 
}